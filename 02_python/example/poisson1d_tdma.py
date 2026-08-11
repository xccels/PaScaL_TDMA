#!/usr/bin/env python3
"""Single-GPU validation of the wrapper library against a SciPy reference.

Solves a batch of identical 1D Poisson tridiagonal systems (a=-1, b=2, c=-1)
on the GPU by calling libpascal_tdma_capi.so directly through ctypes, and
compares the result with scipy.linalg.solve_banded. All coefficient arrays
are device-resident (CuPy).

Usage: python poisson1d_tdma.py
"""
import os
import ctypes
import numpy as np
import cupy as cp
from scipy.linalg import solve_banded

HERE = os.path.dirname(os.path.abspath(__file__))
SO = os.path.join(HERE, "..", "lib", "libpascal_tdma_capi.so")

lib = ctypes.CDLL(SO)
lib.ptdma_cuda_create.restype = None
lib.ptdma_cuda_create.argtypes = [ctypes.POINTER(ctypes.c_int)] + [ctypes.c_int] * 8
lib.ptdma_cuda_solve.restype = None
lib.ptdma_cuda_solve.argtypes = [ctypes.c_int] + [ctypes.c_void_p] * 4
lib.ptdma_cuda_destroy.restype = None
lib.ptdma_cuda_destroy.argtypes = [ctypes.c_int]

# nx_sys, ny_sys = grid of independent systems; nz = system length (solving axis)
nx, ny, nz = 8, 8, 64
tx, ty = 8, 8                       # constraints: nx % tx == 0, ny % ty == 0

# 1D Poisson tridiagonal systems, identical for every (i, j) line; Fortran order
a = cp.full((nx, ny, nz), -1.0, dtype=cp.float64, order="F")
b = cp.full((nx, ny, nz),  2.0, dtype=cp.float64, order="F")
c = cp.full((nx, ny, nz), -1.0, dtype=cp.float64, order="F")
rhs = cp.arange(1, nz + 1, dtype=cp.float64)                 # RHS varying along the solving axis
d = cp.broadcast_to(rhs, (nx, ny, nz)).copy(order="F")       # solve overwrites d with the solution
rhs_host = cp.asnumpy(rhs)

h = ctypes.c_int(0)
lib.ptdma_cuda_create(ctypes.byref(h), nx, ny, nz, 0, 1, 0, tx, ty)   # myrank=0, nprocs=1, comm=0
assert h.value >= 1, "plan create failed (handle < 1)"

lib.ptdma_cuda_solve(
    h,
    ctypes.c_void_p(a.data.ptr), ctypes.c_void_p(b.data.ptr),
    ctypes.c_void_p(c.data.ptr), ctypes.c_void_p(d.data.ptr),
)
cp.cuda.runtime.deviceSynchronize()
lib.ptdma_cuda_destroy(h)

x = cp.asnumpy(d)

# Reference: solving one line is enough since all lines are identical
ab = np.zeros((3, nz))
ab[0, 1:] = -1.0                    # super-diagonal
ab[1, :] = 2.0                      # diagonal
ab[2, :-1] = -1.0                   # sub-diagonal
ref = solve_banded((1, 1), ab, rhs_host)

err = float(np.max(np.abs(x - ref[None, None, :])))
print(f"grid (nx,ny,nz)=({nx},{ny},{nz})  systems={nx*ny}  length={nz}")
print(f"max|gpu - scipy| = {err:.3e}")
print("PASS" if err < 1e-9 else "FAIL")
raise SystemExit(0 if err < 1e-9 else 1)
