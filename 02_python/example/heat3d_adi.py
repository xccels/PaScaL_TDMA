#!/usr/bin/env python3
"""3D time-dependent heat conduction solved with a factored Crank-Nicolson ADI scheme.

Python/CuPy port of 01_Fortran/examples/convection_3D (single process).

Problem: theta_t = Ct * laplacian(theta) in a unit cube, periodic in x and z,
Dirichlet walls in y (hot bottom, cold top). Each time step performs batched
TDMA sweeps in the z, y, x directions: z and x use cyclic TDMA (periodic),
y uses regular TDMA (wall rows are decoupled via jmbc/jpbc and the boundary
values are moved to the right-hand side).

All field arrays are device-resident (CuPy). The RHS/LHS builds and the ghost
cell update reproduce the CUDA kernels of solve_theta.f90 with CuPy slice
operations; only the tridiagonal solves call libpascal_tdma_capi.so.

Usage:
  python heat3d_adi.py [T_field_all.dat]
The optional argument is a reference output of the Fortran example to compare
against. Grid size and step count: environment variables HEAT3D_N, HEAT3D_TMAX.
"""
import os
import sys

import numpy as np
import cupy as cp

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from pascal_tdma import PlanManyCuda

# ============================ Parameters (global_inputpara in global.f90) ============================
n = int(os.environ.get("HEAT3D_N", 16))       # input mesh count (nx=ny=nz in PARA_INPUT.inp)
Tmax = int(os.environ.get("HEAT3D_TMAX", 5))  # number of time steps
dt = 5.0e-3       # dtStart

lx = ly = lz = 1.0
Pr, Ra = 5.0, 2.0e2
theta_cold = -1.0
theta_hot = 2.0 + theta_cold                      # = 1.0
alphaG = 1.0
nu = 1.0 / np.sqrt(Ra / (alphaG * Pr * ly**3 * (theta_hot - theta_cold)))
Ct = nu / Pr                                       # thermal diffusivity

# Grid: global nx=ny=nz=n+1, interior points nxm=n; arrays hold one ghost layer
# on each side, (n+2)^3 points with interior indices 1..n.
NG = n + 2                                          # array size per direction (0..n+1)
dx = lx / n                                         # lx/(nx-1)
dy = ly / (n + 1)                                   # ly/ny
dz = lz / n
coefx, coefy, coefz = 0.5 * Ct / dx, 0.5 * Ct / dy, 0.5 * Ct / dz

# Node coordinates (mpi_subdomain_mesh with ista=jsta=ksta=1)
i_idx = cp.arange(NG, dtype=cp.float64)
x_sub = (i_idx - 1.0) * dx                          # x_sub(i) = (i-1)*dx
y_sub = i_idx * dy                                  # y_sub(j) = j*dy
z_sub = (i_idx - 1.0) * dz                          # z_sub(k) = (k-1)*dz

# Inverse grid spacings (dmx, dmz uniform; dmy has half cells dy/2 at the walls)
inv_dmx = cp.full(NG, 1.0 / dx); inv_dmz = cp.full(NG, 1.0 / dz)
inv_dmy = cp.full(NG, 1.0 / dy); inv_dmy[0] = 2.0 / dy; inv_dmy[n + 1] = 2.0 / dy

# Flags decoupling the wall-adjacent rows in the y-direction (mpi_subdomain_indices)
jpbc = cp.ones(NG, dtype=cp.float64); jpbc[n] = 0.0      # jpbc_index(ny_sub-1) = 0
jmbc = cp.ones(NG, dtype=cp.float64); jmbc[1] = 0.0      # jmbc_index(1) = 0

# Wall boundary values (single process: thetaBC3 = hot bottom, thetaBC4 = cold top)
thetaBC3, thetaBC4 = theta_hot, theta_cold

IN = slice(1, n + 1)                                # interior indices 1..n


def wrap_ghost(th):
    """Periodic wrap of the x and z ghost layers (ghostcell_update_cuda); y walls stay fixed."""
    th[0, :, :] = th[n, :, :]; th[n + 1, :, :] = th[1, :, :]    # x periodic: 0 <- n, n+1 <- 1
    th[:, :, 0] = th[:, :, n]; th[:, :, n + 1] = th[:, :, 1]    # z periodic


def build_rhs(th):
    """Build the RHS (build_RHS_cuda), shape (n, n, n) indexed (i, j, k): viscous + ebc - eRHS."""
    tijk = th[1:n+1, 1:n+1, 1:n+1]
    tip, tim = th[2:n+2, 1:n+1, 1:n+1], th[0:n, 1:n+1, 1:n+1]
    tjp, tjm = th[1:n+1, 2:n+2, 1:n+1], th[1:n+1, 0:n, 1:n+1]
    tkp, tkm = th[1:n+1, 1:n+1, 2:n+2], th[1:n+1, 1:n+1, 0:n]

    idx_i = inv_dmx[1:n+1].reshape(n, 1, 1); idx_ip = inv_dmx[2:n+2].reshape(n, 1, 1)
    idy_j = inv_dmy[1:n+1].reshape(1, n, 1); idy_jp = inv_dmy[2:n+2].reshape(1, n, 1)
    idz_k = inv_dmz[1:n+1].reshape(1, 1, n); idz_kp = inv_dmz[2:n+2].reshape(1, 1, n)
    jep = jpbc[1:n+1].reshape(1, n, 1); jem = jmbc[1:n+1].reshape(1, n, 1)

    dedx1 = (tijk - tim) * idx_i; dedx2 = (tip - tijk) * idx_ip
    dedy3 = (tijk - tjm) * idy_j; dedy4 = (tjp - tijk) * idy_jp
    dedz5 = (tijk - tkm) * idz_k; dedz6 = (tkp - tijk) * idz_kp
    viscous = coefx * (dedx2 - dedx1) + coefy * (dedy4 - dedy3) + coefz * (dedz6 - dedz5)

    ebc = (1.0 - jem) * (coefy * idy_j * thetaBC3) + (1.0 - jep) * (coefy * idy_jp * thetaBC4)

    eAPI, eAMI, eACI = -coefx * idx_ip, -coefx * idx_i, coefx * (idx_ip + idx_i)
    eAPK, eAMK, eACK = -coefz * idz_kp, -coefz * idz_k, coefz * (idz_kp + idz_k)
    eAPJ, eAMJ, eACJ = -coefy * idy_jp * jep, -coefy * idy_j * jem, coefy * (idy_jp + idy_j)
    eRHS = (eAPK * tkp + eACK * tijk + eAMK * tkm
            + eAPJ * tjp + eACJ * tijk + eAMJ * tjm
            + eAPI * tip + eACI * tijk + eAMI * tim)
    return viscous + ebc - eRHS                     # (i, j, k)


def fc(shape=(n, n, n)):
    return cp.empty(shape, dtype=cp.float64, order="F")


def main():
    ref_path = sys.argv[1] if len(sys.argv) > 1 else None

    # ---- Initial condition (mpi_subdomain_initialization): linear profile + sine perturbation ----
    X = x_sub.reshape(NG, 1, 1); Y = y_sub.reshape(1, NG, 1); Z = z_sub.reshape(1, 1, NG)
    theta = ((theta_cold - theta_hot) / ly * Y + theta_hot
             + cp.sin(4.0 * np.pi / lx * X) * cp.sin(4.0 * np.pi / lz * Z)
             * cp.sin(4.0 * np.pi / ly * Y)) * cp.ones((NG, NG, NG))
    theta = cp.asfortranarray(theta)
    theta[:, 0, :] = theta_hot                       # bottom wall (j = 0)
    theta[:, n + 1, :] = theta_cold                  # top wall (j = n+1)
    wrap_ghost(theta)

    # ---- One plan per direction (solve_theta_plan_many_cuda); all (n, n, n), n % 8 == 0 ----
    plan_z = PlanManyCuda(n, n, n, tx=8, ty=8)       # systems (i, j), solve along k (cyclic)
    plan_y = PlanManyCuda(n, n, n, tx=8, ty=8)       # systems (i, k), solve along j (non-cyclic)
    plan_x = PlanManyCuda(n, n, n, tx=8, ty=8)       # systems (j, k), solve along i (cyclic)

    # LHS coefficients viewed along the solving (third) axis
    idy_j3, idy_jp3 = inv_dmy[1:n+1].reshape(1, 1, n), inv_dmy[2:n+2].reshape(1, 1, n)
    jep3, jem3 = jpbc[1:n+1].reshape(1, 1, n), jmbc[1:n+1].reshape(1, 1, n)
    idx_i3, idx_ip3 = inv_dmx[1:n+1].reshape(1, 1, n), inv_dmx[2:n+2].reshape(1, 1, n)
    idz_k3, idz_kp3 = inv_dmz[1:n+1].reshape(1, 1, n), inv_dmz[2:n+2].reshape(1, 1, n)

    for step in range(Tmax):
        rhs = build_rhs(theta)                       # (i, j, k)

        # ----- z-direction (cyclic): build_LHSz + solve_cycle + copy_ijk2ijk -----
        am = fc(); am[...] = (-coefz * dt) * idz_k3
        ap = fc(); ap[...] = (-coefz * dt) * idz_kp3
        ac = fc(); ac[...] = coefz * dt * (idz_kp3 + idz_k3) + 1.0
        ad = cp.asfortranarray(rhs * dt)
        plan_z.solve(am, ac, ap, ad, cyclic=True)    # solve(a=lower, b=main, c=upper, d=rhs)
        rhs = ad                                     # copy_ijk2ijk: rhs = ad (i, j, k)

        # ----- y-direction (non-cyclic): build_LHSy + solve + transpose_ikj2ijk -----
        ad = cp.asfortranarray(rhs.transpose(0, 2, 1))   # ad(i, k, j) = rhs(i, j, k)
        am = fc(); am[...] = (-coefy * dt) * idy_j3 * jem3
        ap = fc(); ap[...] = (-coefy * dt) * idy_jp3 * jep3
        ac = fc(); ac[...] = coefy * dt * (idy_jp3 + idy_j3) + 1.0
        plan_y.solve(am, ac, ap, ad, cyclic=False)
        rhs = ad.transpose(0, 2, 1)                  # transpose_ikj2ijk: rhs(i, j, k) = ad(i, k, j)

        # ----- x-direction (cyclic): build_LHSx + solve_cycle + update_theta -----
        ad = cp.asfortranarray(rhs.transpose(1, 2, 0))   # ad(j, k, i) = rhs(i, j, k)
        am = fc(); am[...] = (-coefx * dt) * idx_i3
        ap = fc(); ap[...] = (-coefx * dt) * idx_ip3
        ac = fc(); ac[...] = coefx * dt * (idx_ip3 + idx_i3) + 1.0
        plan_x.solve(am, ac, ap, ad, cyclic=True)
        theta[1:n+1, 1:n+1, 1:n+1] += ad.transpose(2, 0, 1)   # update_theta: theta += ad(j, k, i)

        wrap_ghost(theta)

    for p in (plan_z, plan_y, plan_x):
        p.destroy()

    sol = cp.asnumpy(theta[1:n+1, 1:n+1, 1:n+1])     # interior theta (i, j, k), i = 1..n

    # ---- Comparison with the Fortran reference ----
    if ref_path is None or not os.path.exists(ref_path):
        print(f"[WARN] reference file not found: {ref_path} — comparison skipped.")
        print(f"solved theta: min={sol.min():.6f} max={sol.max():.6f} mean={sol.mean():.6f}")
        return 0

    vals = []
    with open(ref_path) as f:
        for line in f:
            p = line.split()
            if len(p) >= 4:
                vals.append(float(p[3].replace("D", "E")))
    ref = np.array(vals).reshape((n, n, n), order="F")   # (i, j, k), i fastest (do k; do j; do i)

    err = np.max(np.abs(sol - ref))
    rel = err / max(np.max(np.abs(ref)), 1e-300)
    print(f"grid={n}^3  Tmax={Tmax}  Ct={Ct:.6f}  dt={dt}")
    print(f"theta range: solved [{sol.min():.6f}, {sol.max():.6f}]  ref [{ref.min():.6f}, {ref.max():.6f}]")
    print(f"max|python(GPU) - fortran(ref)| = {err:.3e}   (relative {rel:.3e})")
    ok = err < 1e-9
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
