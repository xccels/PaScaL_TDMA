"""Python interface to the device-resident CUDA solver of PaScaL_TDMA.

Loads libpascal_tdma_capi.so (ctypes) and exposes the batched tridiagonal
solver through the PlanManyCuda class. Coefficient arrays are CuPy device
arrays; only their device addresses cross the language boundary, so all data
stays on the GPU.

Array convention: dtype float64, Fortran order (order='F'),
shape (nx_sys, ny_sys, nz_row). The tridiagonal (solving) direction is the
last axis (nz_row). Argument order follows the library:
solve(a=lower, b=main, c=upper, d=rhs); d is overwritten with the solution.
"""
import ctypes
import os

import cupy as cp

_SO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib", "libpascal_tdma_capi.so")
_lib = ctypes.CDLL(_SO)

_lib.ptdma_cuda_create.restype = None
_lib.ptdma_cuda_create.argtypes = [ctypes.POINTER(ctypes.c_int)] + [ctypes.c_int] * 8
for _n in ("ptdma_cuda_solve", "ptdma_cuda_solve_cycle"):
    _f = getattr(_lib, _n)
    _f.restype = None
    _f.argtypes = [ctypes.c_int] + [ctypes.c_void_p] * 4
_lib.ptdma_cuda_destroy.restype = None
_lib.ptdma_cuda_destroy.argtypes = [ctypes.c_int]


def _devptr(a, shape):
    """Validate a CuPy array against the plan and return its device address."""
    if a.dtype != cp.float64:
        raise TypeError(f"array must be float64, got {a.dtype}")
    if a.shape != shape:
        raise ValueError(f"array shape {a.shape} != plan shape {shape}")
    if not a.flags.f_contiguous:
        raise ValueError("array must be Fortran-contiguous (order='F')")
    return ctypes.c_void_p(a.data.ptr)


class PlanManyCuda:
    """Plan for many tridiagonal systems on the GPU (PaScaL_TDMA_plan_many_cuda).

    shape = (nx_sys, ny_sys, nz_row); constraints: nx_sys % tx == 0, ny_sys % ty == 0.
    Single GPU: nprocs=1, comm=0. Multi GPU: pass the rank, size, and Fortran
    handle (mpi4py COMM.py2f()) of the 1D communicator along the solving direction.
    """

    def __init__(self, nx_sys, ny_sys, nz_row, tx=8, ty=8, myrank=0, nprocs=1, comm=0):
        if nx_sys % tx or ny_sys % ty:
            raise ValueError(f"nx_sys({nx_sys})%tx({tx})!=0 or ny_sys({ny_sys})%ty({ty})!=0")
        h = ctypes.c_int(0)
        _lib.ptdma_cuda_create(ctypes.byref(h), nx_sys, ny_sys, nz_row,
                               myrank, nprocs, comm, tx, ty)
        if h.value < 1:
            raise RuntimeError("ptdma_cuda_create failed (no free plan slot?)")
        self._h = h.value
        self.shape = (nx_sys, ny_sys, nz_row)

    def solve(self, a, b, c, d, cyclic=False):
        """Solve in place; d is overwritten with the solution (c is also modified, a/b preserved)."""
        if self._h is None:
            raise RuntimeError("plan already destroyed")
        fn = _lib.ptdma_cuda_solve_cycle if cyclic else _lib.ptdma_cuda_solve
        fn(self._h, _devptr(a, self.shape), _devptr(b, self.shape),
           _devptr(c, self.shape), _devptr(d, self.shape))

    def destroy(self):
        if getattr(self, "_h", None) is not None:
            _lib.ptdma_cuda_destroy(self._h)
            self._h = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.destroy()
