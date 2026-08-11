!=======================================================================================================================
!> @file        pascal_tdma_capi.f90
!> @brief       This file contains a C-callable wrapper module around the CUDA solver of PaScaL_TDMA.
!> @details     The wrapper exposes PaScaL_TDMA_cuda through ISO_C_BINDING so that the device-resident
!>              GPU solver can be driven from other languages, e.g. Python (ctypes) with CuPy device
!>              arrays and CUDA-aware mpi4py.
!>              - The Fortran derived-type plan (ptdma_plan_many_cuda) cannot cross the C boundary,
!>                so plans are kept in a module-level table and referenced by an integer handle.
!>              - Coefficient arrays are passed as raw device addresses typed as integer(c_intptr_t).
!>                A Fortran device pointer of shape (nx_sys, ny_sys, nz_row) is rebuilt from the address
!>                via c_devptr and c_f_pointer, and the original solver is called. Data stays on the GPU.
!>              - Expected array layout: shape (nx_sys, ny_sys, nz_row), column-major (Fortran order).
!>                The tridiagonal (solving) direction is the third axis (nz_row).
!>
!> @author
!>              - Jungwoo Kim (yasandy@yonsei.ac.kr), School of Mathematics and Computing (Computational Science and Engineering), Yonsei University
!>
!> @date        July 2026
!> @version     2.0
!> @par         License
!>              This project is released under the terms of the MIT License (see LICENSE file).
!=======================================================================================================================

!>
!> @brief       Module wrapping PaScaL_TDMA_cuda behind a C-compatible interface.
!> @details     Plans are stored in a fixed-size slot table. A handle in [1, MAXPLAN] identifies a slot;
!>              0 or negative handles indicate failure. The array sizes given at plan creation are stored
!>              in the slot and reused to rebuild device pointers in the solve routines.
!>
module pascal_tdma_capi

    use iso_c_binding
    use cudafor
    use PaScaL_TDMA_cuda
    implicit none

    integer, parameter :: MAXPLAN = 16

    type :: plan_slot
        logical :: used = .false.
        integer :: nx = 0, ny = 0, nz = 0          ! nx_sys, ny_sys, nz_row (for c_f_pointer shape)
        type(ptdma_plan_many_cuda) :: p
    end type plan_slot

    type(plan_slot), save :: slots(MAXPLAN)

contains

    !>
    !> @brief       Create a plan for many tridiagonal systems and return its handle.
    !> @details     For single-GPU use pass myrank=0, nprocs=1, comm=0 (comm is unused when nprocs==1).
    !>              For multi-GPU use pass the rank, size, and Fortran handle (e.g. mpi4py COMM.py2f())
    !>              of the 1D communicator along the solving direction.
    !> @param       handle      Plan handle in [1, MAXPLAN] on success, -1 on failure (no free slot)
    !> @param       nx_sys      Number of systems in the first axis (constraint: nx_sys % tx == 0)
    !> @param       ny_sys      Number of systems in the second axis (constraint: ny_sys % ty == 0)
    !> @param       nz_row      Length of each tridiagonal system (third axis)
    !> @param       myrank      Rank ID in the communicator
    !> @param       nprocs      Number of processes in the communicator
    !> @param       comm        Fortran handle of the communicator
    !> @param       tx          Thread-block size in the first axis
    !> @param       ty          Thread-block size in the second axis
    !>
    subroutine ptdma_cuda_create(handle, nx_sys, ny_sys, nz_row, myrank, nprocs, comm, tx, ty) &
                                 bind(C, name="ptdma_cuda_create")
        integer(c_int), intent(out) :: handle
        integer(c_int), value :: nx_sys, ny_sys, nz_row
        integer(c_int), value :: myrank, nprocs, comm, tx, ty
        type(dim3) :: thr
        integer :: h

        handle = -1
        do h = 1, MAXPLAN
            if (.not. slots(h)%used) then
                thr = dim3(tx, ty, 1)
                call PaScaL_TDMA_plan_many_create_cuda(slots(h)%p, nx_sys, ny_sys, nz_row, &
                                                       myrank, nprocs, comm, thr)
                slots(h)%nx = nx_sys
                slots(h)%ny = ny_sys
                slots(h)%nz = nz_row
                slots(h)%used = .true.
                handle = h
                return
            end if
        end do
    end subroutine ptdma_cuda_create

    !>
    !> @brief       Solve many tridiagonal systems of equations.
    !> @details     On return d holds the solution (c is also overwritten; a and b are unchanged).
    !> @param       handle      Plan handle returned by ptdma_cuda_create
    !> @param       a_addr      Device address of the lower diagonal, (nx, ny, nz) Fortran-order array
    !> @param       b_addr      Device address of the main diagonal
    !> @param       c_addr      Device address of the upper diagonal
    !> @param       d_addr      Device address of the right-hand side / solution
    !>
    subroutine ptdma_cuda_solve(handle, a_addr, b_addr, c_addr, d_addr) &
                                bind(C, name="ptdma_cuda_solve")
        integer(c_int), value :: handle
        integer(c_intptr_t), value :: a_addr, b_addr, c_addr, d_addr
        type(c_devptr) :: ca, cb, cc, cd
        double precision, device, pointer :: ap(:,:,:), bp(:,:,:), cp(:,:,:), dp(:,:,:)
        integer :: nx, ny, nz

        if (handle < 1 .or. handle > MAXPLAN) return
        if (.not. slots(handle)%used) return
        nx = slots(handle)%nx; ny = slots(handle)%ny; nz = slots(handle)%nz

        ca = transfer(a_addr, ca); call c_f_pointer(ca, ap, [nx, ny, nz])
        cb = transfer(b_addr, cb); call c_f_pointer(cb, bp, [nx, ny, nz])
        cc = transfer(c_addr, cc); call c_f_pointer(cc, cp, [nx, ny, nz])
        cd = transfer(d_addr, cd); call c_f_pointer(cd, dp, [nx, ny, nz])

        call PaScaL_TDMA_many_solve_cuda(slots(handle)%p, ap, bp, cp, dp)
    end subroutine ptdma_cuda_solve

    !>
    !> @brief       Solve many cyclic tridiagonal systems of equations.
    !> @details     Same argument convention as ptdma_cuda_solve.
    !> @param       handle      Plan handle returned by ptdma_cuda_create
    !> @param       a_addr      Device address of the lower diagonal, (nx, ny, nz) Fortran-order array
    !> @param       b_addr      Device address of the main diagonal
    !> @param       c_addr      Device address of the upper diagonal
    !> @param       d_addr      Device address of the right-hand side / solution
    !>
    subroutine ptdma_cuda_solve_cycle(handle, a_addr, b_addr, c_addr, d_addr) &
                                      bind(C, name="ptdma_cuda_solve_cycle")
        integer(c_int), value :: handle
        integer(c_intptr_t), value :: a_addr, b_addr, c_addr, d_addr
        type(c_devptr) :: ca, cb, cc, cd
        double precision, device, pointer :: ap(:,:,:), bp(:,:,:), cp(:,:,:), dp(:,:,:)
        integer :: nx, ny, nz

        if (handle < 1 .or. handle > MAXPLAN) return
        if (.not. slots(handle)%used) return
        nx = slots(handle)%nx; ny = slots(handle)%ny; nz = slots(handle)%nz

        ca = transfer(a_addr, ca); call c_f_pointer(ca, ap, [nx, ny, nz])
        cb = transfer(b_addr, cb); call c_f_pointer(cb, bp, [nx, ny, nz])
        cc = transfer(c_addr, cc); call c_f_pointer(cc, cp, [nx, ny, nz])
        cd = transfer(d_addr, cd); call c_f_pointer(cd, dp, [nx, ny, nz])

        call PaScaL_TDMA_many_solve_cycle_cuda(slots(handle)%p, ap, bp, cp, dp)
    end subroutine ptdma_cuda_solve_cycle

    !>
    !> @brief       Destroy a plan and free its slot.
    !> @param       handle      Plan handle returned by ptdma_cuda_create
    !>
    subroutine ptdma_cuda_destroy(handle) bind(C, name="ptdma_cuda_destroy")
        integer(c_int), value :: handle
        if (handle < 1 .or. handle > MAXPLAN) return
        if (.not. slots(handle)%used) return
        call PaScaL_TDMA_plan_many_destroy_cuda(slots(handle)%p)
        slots(handle)%used = .false.
    end subroutine ptdma_cuda_destroy

end module pascal_tdma_capi
