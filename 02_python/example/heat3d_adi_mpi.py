#!/usr/bin/env python3
"""3D heat conduction ADI example on multiple GPUs with CUDA-aware MPI.

Multi-process generalization of heat3d_adi.py, following the parallel structure
of 01_Fortran/examples/convection_3D:
  - 3D cartesian topology (periodic x/z, walls in y) with one 1D
    sub-communicator per direction                       (mpi_topology.f90)
  - para_range domain decomposition; subdomains carry one ghost layer
    on each side                                          (mpi_subdomain.f90)
  - ghost cell update: boundary planes are packed into GPU buffers and
    exchanged directly with CUDA-aware MPI Isend/Irecv    (ghostcell_update_cuda)
  - one TDMA plan per direction, created with that direction's 1D communicator
    (COMM.py2f()); systems split across ranks are solved inside the library
    through device-resident MPI_Alltoall                  (solve_theta_plan_many_cuda)
  - GPU affinity: device = world_rank % (GPUs on the node) (main.f90)

Usage (e.g. two GPUs, domain split in z):
  mpirun -np 2 python heat3d_adi_mpi.py [T_field_all.dat]
Decomposition: HEAT3D_NP="npx,npy,npz" (default "1,1,nprocs"); grid HEAT3D_N
or per-direction HEAT3D_NX/NY/NZ, steps HEAT3D_TMAX. Plan constraint: the two
system-count axes of every plan must be divisible by (HEAT3D_TX, HEAT3D_TY) = (8, 8).

Benchmark mode (HEAT3D_BENCH=1): measures per-step wall time split into
computation / TDMA solve (incl. its internal all-to-all) / ghost exchange,
excludes the first HEAT3D_WARMUP steps (default 2), and prints one line
"BENCH,np,dims,nx,ny,nz,steps,total,compute,solve,ghost,pool_gb,check"
on rank 0 (times are per-step averages, max over ranks). The gather and
reference comparison are skipped.
"""
import os
import socket
import sys

import numpy as np
import cupy as cp

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
# NOTE: pascal_tdma (which loads libpascal_tdma_capi.so) must be imported before
# mpi4py. If mpi4py initializes MPI first, the NVHPC CUDA Fortran runtime fails
# device detection ("No accelerator device found") due to library load order.
from pascal_tdma import PlanManyCuda
from mpi4py import MPI

world = MPI.COMM_WORLD
myrank, nprocs = world.Get_rank(), world.Get_size()

# GPU affinity must precede any CuPy allocation and plan creation. The device
# index comes from the node-local rank: equivalent to main.f90's
# cudaSetDevice(mod(myrank, nDevices)) under block placement, and still correct
# under any rank placement (e.g. mpirun --map-by node).
_local = world.Split_type(MPI.COMM_TYPE_SHARED)
local_rank = _local.Get_rank()
_local.Free()
ndev = cp.cuda.runtime.getDeviceCount()
cp.cuda.Device(local_rank % ndev).use()

# ============================ Parameters (global_inputpara in global.f90) ============================
n = int(os.environ.get("HEAT3D_N", 16))       # mesh count (nx=ny=nz in PARA_INPUT.inp)
nxg = int(os.environ.get("HEAT3D_NX", n))     # per-direction interior grid counts
nyg = int(os.environ.get("HEAT3D_NY", n))
nzg = int(os.environ.get("HEAT3D_NZ", n))
Tmax = int(os.environ.get("HEAT3D_TMAX", 5))  # number of time steps
TX = int(os.environ.get("HEAT3D_TX", 8))      # plan thread-block sizes (thread_in_x/y_pascal)
TY = int(os.environ.get("HEAT3D_TY", 8))
BENCH = int(os.environ.get("HEAT3D_BENCH", 0))    # 1: per-step timing, gather/compare skipped
WARMUP = int(os.environ.get("HEAT3D_WARMUP", 2))  # warmup steps excluded from timing
dt = 5.0e-3       # dtStart

lx = ly = lz = 1.0
Pr, Ra = 5.0, 2.0e2
theta_cold = -1.0
theta_hot = 2.0 + theta_cold                      # = 1.0
alphaG = 1.0
nu = 1.0 / np.sqrt(Ra / (alphaG * Pr * ly**3 * (theta_hot - theta_cold)))
Ct = nu / Pr                                       # thermal diffusivity

dx = lx / nxg                                      # lx/(nx-1); interior points 1..nxg
dy = ly / (nyg + 1)                                # ly/ny
dz = lz / nzg
coefx, coefy, coefz = 0.5 * Ct / dx, 0.5 * Ct / dy, 0.5 * Ct / dz

# ============================ Topology (mpi_topology.f90) ============================
np_dim = [int(s) for s in os.environ.get("HEAT3D_NP", f"1,1,{nprocs}").split(",")]
if np_dim[0] * np_dim[1] * np_dim[2] != nprocs:
    raise SystemExit(f"product of HEAT3D_NP={np_dim} differs from nprocs({nprocs})")
period = [True, False, True]                       # periodic x/z, walls in y
cart = world.Create_cart(np_dim, periods=period, reorder=False)


def sub_comm(dim):
    """1D sub-communicator of direction dim and its (rank, size, west, east) — MPI_Cart_sub + MPI_Cart_shift."""
    keep = [False, False, False]
    keep[dim] = True
    c = cart.Sub(keep)
    west, east = c.Shift(0, 1)
    return c, c.Get_rank(), c.Get_size(), west, east


comm_x, rx, npx, west_x, east_x = sub_comm(0)
comm_y, ry, npy, west_y, east_y = sub_comm(1)
comm_z, rz, npz, west_z, east_z = sub_comm(2)

# ============================ Decomposition (para_range + mpi_subdomain_make) ============================
def para_range(n1, n2, nprocs_dir, rank_dir):
    base, extra = divmod(n2 - n1 + 1, nprocs_dir)
    ista = rank_dir * base + n1 + min(rank_dir, extra)
    iend = ista + base - 1 + (1 if extra > rank_dir else 0)
    return ista, iend


ista, iend = para_range(1, nxg, npx, rx); nx_sub = iend - ista + 2
jsta, jend = para_range(1, nyg, npy, ry); ny_sub = jend - jsta + 2
ksta, kend = para_range(1, nzg, npz, rz); nz_sub = kend - ksta + 2
nxm, nym, nzm = nx_sub - 1, ny_sub - 1, nz_sub - 1     # interior grid counts of the subdomain

# ---- Coordinates and grid spacings (mpi_subdomain_mesh) ----
x_loc = (cp.arange(nx_sub + 1, dtype=cp.float64) - 1.0 + (ista - 1)) * dx
y_loc = (cp.arange(ny_sub + 1, dtype=cp.float64) + (jsta - 1)) * dy
z_loc = (cp.arange(nz_sub + 1, dtype=cp.float64) - 1.0 + (ksta - 1)) * dz

inv_dmx = cp.full(nx_sub + 1, 1.0 / dx)
inv_dmz = cp.full(nz_sub + 1, 1.0 / dz)
inv_dmy = cp.full(ny_sub + 1, 1.0 / dy)                # half cells only at the physical walls
if ry == 0:
    inv_dmy[0] = 2.0 / dy
if ry == npy - 1:
    inv_dmy[ny_sub] = 2.0 / dy

# ---- Flags decoupling the wall-adjacent rows in y (mpi_subdomain_indices) ----
jmbc = cp.ones(ny_sub + 1, dtype=cp.float64)
jpbc = cp.ones(ny_sub + 1, dtype=cp.float64)
if ry == 0:
    jmbc[1] = 0.0
if ry == npy - 1:
    jpbc[ny_sub - 1] = 0.0

IN = (slice(1, nx_sub), slice(1, ny_sub), slice(1, nz_sub))   # interior slice

# ============================ Ghost cell update (ghostcell_update_cuda) ============================
def _bufs(a, b):
    return (cp.empty((a, b)), cp.empty((a, b)), cp.empty((a, b)), cp.empty((a, b)))


sbx0, sbx1, rbx0, rbx1 = _bufs(ny_sub + 1, nz_sub + 1)
sby0, sby1, rby0, rby1 = _bufs(nx_sub + 1, nz_sub + 1)
sbz0, sbz1, rbz0, rbz1 = _bufs(nx_sub + 1, ny_sub + 1)


def _exchange(comm_d, west, east, sb0, sb1, rb0, rb1, periodic):
    """Local swap for a single periodic rank; otherwise Isend/Irecv on GPU buffers (CUDA-aware)."""
    if comm_d.Get_size() == 1 and periodic:
        rb1[...] = sb0
        rb0[...] = sb1
    else:
        cp.cuda.get_current_stream().synchronize()     # cudaStreamSynchronize() before MPI
        reqs = [comm_d.Isend(sb0, dest=west, tag=111),
                comm_d.Irecv(rb1, source=east, tag=111),
                comm_d.Irecv(rb0, source=west, tag=222),
                comm_d.Isend(sb1, dest=east, tag=222)]
        MPI.Request.Waitall(reqs)


def ghost_update(th):
    """Update the ghost layers of th, sweeping x, y, z in order (matches the corner propagation of the Fortran code)."""
    NULL = MPI.PROC_NULL
    # X
    if west_x != NULL:
        sbx0[...] = th[1, :, :]
    if east_x != NULL:
        sbx1[...] = th[nx_sub - 1, :, :]
    _exchange(comm_x, west_x, east_x, sbx0, sbx1, rbx0, rbx1, period[0])
    if west_x != NULL:
        th[0, :, :] = rbx0
    if east_x != NULL:
        th[nx_sub, :, :] = rbx1
    # Y
    if west_y != NULL:
        sby0[...] = th[:, 1, :]
    if east_y != NULL:
        sby1[...] = th[:, ny_sub - 1, :]
    _exchange(comm_y, west_y, east_y, sby0, sby1, rby0, rby1, period[1])
    if west_y != NULL:
        th[:, 0, :] = rby0
    if east_y != NULL:
        th[:, ny_sub, :] = rby1
    # Z
    if west_z != NULL:
        sbz0[...] = th[:, :, 1]
    if east_z != NULL:
        sbz1[...] = th[:, :, nz_sub - 1]
    _exchange(comm_z, west_z, east_z, sbz0, sbz1, rbz0, rbz1, period[2])
    if west_z != NULL:
        th[:, :, 0] = rbz0
    if east_z != NULL:
        th[:, :, nz_sub] = rbz1


def _now():
    """Wall-clock timestamp taken after completing all queued GPU work (benchmark mode)."""
    cp.cuda.get_current_stream().synchronize()
    return MPI.Wtime()


# ============================ RHS (build_RHS_cuda) ============================
idx_i = inv_dmx[1:nxm + 1].reshape(nxm, 1, 1); idx_ip = inv_dmx[2:nxm + 2].reshape(nxm, 1, 1)
idy_j = inv_dmy[1:nym + 1].reshape(1, nym, 1); idy_jp = inv_dmy[2:nym + 2].reshape(1, nym, 1)
idz_k = inv_dmz[1:nzm + 1].reshape(1, 1, nzm); idz_kp = inv_dmz[2:nzm + 2].reshape(1, 1, nzm)
jep = jpbc[1:nym + 1].reshape(1, nym, 1); jem = jmbc[1:nym + 1].reshape(1, nym, 1)


def build_rhs(th, bc3, bc4):
    """Build the RHS, shape (nxm, nym, nzm): viscous + ebc - eRHS.

    bc3/bc4 are the (0:nx_sub, 0:nz_sub) boundary-value planes of the lower and
    upper y-walls (wall temperature on wall ranks, neighbor ghost otherwise).
    """
    tijk = th[1:nx_sub, 1:ny_sub, 1:nz_sub]
    tip, tim = th[2:nx_sub + 1, 1:ny_sub, 1:nz_sub], th[0:nx_sub - 1, 1:ny_sub, 1:nz_sub]
    tjp, tjm = th[1:nx_sub, 2:ny_sub + 1, 1:nz_sub], th[1:nx_sub, 0:ny_sub - 1, 1:nz_sub]
    tkp, tkm = th[1:nx_sub, 1:ny_sub, 2:nz_sub + 1], th[1:nx_sub, 1:ny_sub, 0:nz_sub - 1]

    dedx1 = (tijk - tim) * idx_i; dedx2 = (tip - tijk) * idx_ip
    dedy3 = (tijk - tjm) * idy_j; dedy4 = (tjp - tijk) * idy_jp
    dedz5 = (tijk - tkm) * idz_k; dedz6 = (tkp - tijk) * idz_kp
    viscous = coefx * (dedx2 - dedx1) + coefy * (dedy4 - dedy3) + coefz * (dedz6 - dedz5)

    b3 = bc3[1:nx_sub, 1:nz_sub][:, cp.newaxis, :]     # (nxm, 1, nzm)
    b4 = bc4[1:nx_sub, 1:nz_sub][:, cp.newaxis, :]
    ebc = (1.0 - jem) * (coefy * idy_j * b3) + (1.0 - jep) * (coefy * idy_jp * b4)

    eAPI, eAMI, eACI = -coefx * idx_ip, -coefx * idx_i, coefx * (idx_ip + idx_i)
    eAPK, eAMK, eACK = -coefz * idz_kp, -coefz * idz_k, coefz * (idz_kp + idz_k)
    eAPJ, eAMJ, eACJ = -coefy * idy_jp * jep, -coefy * idy_j * jem, coefy * (idy_jp + idy_j)
    eRHS = (eAPK * tkp + eACK * tijk + eAMK * tkm
            + eAPJ * tjp + eACJ * tijk + eAMJ * tjm
            + eAPI * tip + eACI * tijk + eAMI * tim)
    return viscous + ebc - eRHS


def fc(shape):
    return cp.empty(shape, dtype=cp.float64, order="F")


def main():
    ref_path = sys.argv[1] if len(sys.argv) > 1 else None

    if myrank == 0:
        print(f"grid={nxg}x{nyg}x{nzg}  np_dim={np_dim}  Tmax={Tmax}  Ct={Ct:.6f}  dt={dt}", flush=True)
    print(f"[rank {myrank}] host={socket.gethostname()} cart=({rx},{ry},{rz}) "
          f"sub=({nxm},{nym},{nzm}) dev={local_rank % ndev}", flush=True)

    # ---- Initial condition (mpi_subdomain_initialization): linear profile + sine perturbation ----
    X = x_loc.reshape(-1, 1, 1); Y = y_loc.reshape(1, -1, 1); Z = z_loc.reshape(1, 1, -1)
    theta = ((theta_cold - theta_hot) / ly * Y + theta_hot
             + cp.sin(4.0 * np.pi / lx * X) * cp.sin(4.0 * np.pi / lz * Z)
             * cp.sin(4.0 * np.pi / ly * Y)) * cp.ones((nx_sub + 1, ny_sub + 1, nz_sub + 1))
    theta = cp.asfortranarray(theta)
    if ry == 0:
        theta[:, 0, :] = theta_hot                     # bottom wall
    if ry == npy - 1:
        theta[:, ny_sub, :] = theta_cold               # top wall
    ghost_update(theta)

    # ---- Boundary-value planes (mpi_subdomain_boundary): wall temperature on wall ranks, neighbor ghost otherwise ----
    thetaBC3 = theta[:, 0, :].copy()
    thetaBC4 = theta[:, ny_sub, :].copy()
    if ry == 0:
        thetaBC3[...] = theta_hot
    if ry == npy - 1:
        thetaBC4[...] = theta_cold

    # ---- One plan per direction, created with that direction's 1D communicator (solve_theta_plan_many_cuda) ----
    plan_z = PlanManyCuda(nxm, nym, nzm, tx=TX, ty=TY, myrank=rz, nprocs=npz, comm=comm_z.py2f())
    plan_y = PlanManyCuda(nxm, nzm, nym, tx=TX, ty=TY, myrank=ry, nprocs=npy, comm=comm_y.py2f())
    plan_x = PlanManyCuda(nym, nzm, nxm, tx=TX, ty=TY, myrank=rx, nprocs=npx, comm=comm_x.py2f())

    # LHS coefficients viewed along the solving (third) axis
    idz_k3, idz_kp3 = inv_dmz[1:nzm + 1].reshape(1, 1, nzm), inv_dmz[2:nzm + 2].reshape(1, 1, nzm)
    idy_j3, idy_jp3 = inv_dmy[1:nym + 1].reshape(1, 1, nym), inv_dmy[2:nym + 2].reshape(1, 1, nym)
    jem3, jep3 = jmbc[1:nym + 1].reshape(1, 1, nym), jpbc[1:nym + 1].reshape(1, 1, nym)
    idx_i3, idx_ip3 = inv_dmx[1:nxm + 1].reshape(1, 1, nxm), inv_dmx[2:nxm + 2].reshape(1, 1, nxm)

    steps_rec = []                                     # per measured step: (total, solve, ghost)
    for step in range(Tmax):
        if BENCH:
            world.Barrier()
            t_step = _now(); t_solve = 0.0; t_ghost = 0.0

        rhs = build_rhs(theta, thetaBC3, thetaBC4)     # (i, j, k)

        # ----- z-direction (cyclic): build_LHSz + solve_cycle + copy_ijk2ijk -----
        am = fc((nxm, nym, nzm)); am[...] = (-coefz * dt) * idz_k3
        ap = fc((nxm, nym, nzm)); ap[...] = (-coefz * dt) * idz_kp3
        ac = fc((nxm, nym, nzm)); ac[...] = coefz * dt * (idz_kp3 + idz_k3) + 1.0
        ad = cp.asfortranarray(rhs * dt)
        if BENCH: t0 = _now()
        plan_z.solve(am, ac, ap, ad, cyclic=True)      # solve(a=lower, b=main, c=upper, d=rhs)
        if BENCH: t_solve += _now() - t0
        rhs = ad                                       # copy_ijk2ijk

        # ----- y-direction (non-cyclic): build_LHSy + solve + transpose_ikj2ijk -----
        ad = cp.asfortranarray(rhs.transpose(0, 2, 1))     # ad(i, k, j)
        am = fc((nxm, nzm, nym)); am[...] = (-coefy * dt) * idy_j3 * jem3
        ap = fc((nxm, nzm, nym)); ap[...] = (-coefy * dt) * idy_jp3 * jep3
        ac = fc((nxm, nzm, nym)); ac[...] = coefy * dt * (idy_jp3 + idy_j3) + 1.0
        if BENCH: t0 = _now()
        plan_y.solve(am, ac, ap, ad, cyclic=False)
        if BENCH: t_solve += _now() - t0
        rhs = ad.transpose(0, 2, 1)                    # transpose_ikj2ijk

        # ----- x-direction (cyclic): build_LHSx + solve_cycle + update_theta -----
        ad = cp.asfortranarray(rhs.transpose(1, 2, 0))     # ad(j, k, i)
        am = fc((nym, nzm, nxm)); am[...] = (-coefx * dt) * idx_i3
        ap = fc((nym, nzm, nxm)); ap[...] = (-coefx * dt) * idx_ip3
        ac = fc((nym, nzm, nxm)); ac[...] = coefx * dt * (idx_ip3 + idx_i3) + 1.0
        if BENCH: t0 = _now()
        plan_x.solve(am, ac, ap, ad, cyclic=True)
        if BENCH: t_solve += _now() - t0
        theta[IN] += ad.transpose(2, 0, 1)             # update_theta

        if BENCH: t0 = _now()
        ghost_update(theta)
        if BENCH:
            t_ghost += _now() - t0
            if step >= WARMUP:
                steps_rec.append((_now() - t_step, t_solve, t_ghost))

    for p in (plan_z, plan_y, plan_x):
        p.destroy()

    if BENCH:
        # Per-step averages per category; compute = total - solve - ghost. Max over ranks.
        rec = np.array(steps_rec)
        loc = np.array([rec[:, 0].mean(),
                        (rec[:, 0] - rec[:, 1] - rec[:, 2]).mean(),
                        rec[:, 1].mean(),
                        rec[:, 2].mean(),
                        cp.get_default_memory_pool().total_bytes() / 2**30])
        gmax = np.empty_like(loc)
        world.Allreduce(loc, gmax, op=MPI.MAX)
        finite = np.array([1.0 if bool(cp.all(cp.isfinite(theta[IN])).get()) else 0.0])
        gfin = np.empty_like(finite)
        world.Allreduce(finite, gfin, op=MPI.MIN)
        if myrank == 0:
            print(f"BENCH,{nprocs},{np_dim[0]}x{np_dim[1]}x{np_dim[2]},{nxg},{nyg},{nzg},"
                  f"{rec.shape[0]},{gmax[0]:.6e},{gmax[1]:.6e},{gmax[2]:.6e},{gmax[3]:.6e},"
                  f"{gmax[4]:.2f},{'OK' if gfin[0] > 0 else 'NAN'}", flush=True)
        return 0 if gfin[0] > 0 else 1

    # ---- Gather interior theta on rank 0 and assemble the global field (post-assembly IO of field_file_write) ----
    sol = cp.asnumpy(theta[IN])
    pieces = world.gather((ista, jsta, ksta, sol), root=0)

    rc = 0
    if myrank == 0:
        glob = np.empty((nxg, nyg, nzg))
        for (i0, j0, k0, arr) in pieces:
            si, sj, sk = arr.shape
            glob[i0 - 1:i0 - 1 + si, j0 - 1:j0 - 1 + sj, k0 - 1:k0 - 1 + sk] = arr

        if ref_path is None or not os.path.exists(ref_path):
            print(f"[WARN] reference file not found: {ref_path} — comparison skipped.")
            print(f"solved theta: min={glob.min():.6f} max={glob.max():.6f} mean={glob.mean():.6f}")
        else:
            vals = []
            with open(ref_path) as f:
                for line in f:
                    p = line.split()
                    if len(p) >= 4:
                        vals.append(float(p[3].replace("D", "E")))
            ref = np.array(vals).reshape((nxg, nyg, nzg), order="F")   # (i, j, k), i fastest
            err = np.max(np.abs(glob - ref))
            rel = err / max(np.max(np.abs(ref)), 1e-300)
            print(f"theta range: solved [{glob.min():.6f}, {glob.max():.6f}]  ref [{ref.min():.6f}, {ref.max():.6f}]")
            print(f"max|python(GPU,np={nprocs}) - fortran(ref)| = {err:.3e}   (relative {rel:.3e})")
            rc = 0 if err < 1e-9 else 1
            print("PASS" if rc == 0 else "FAIL")

    return world.bcast(rc, root=0)


if __name__ == "__main__":
    raise SystemExit(main())
