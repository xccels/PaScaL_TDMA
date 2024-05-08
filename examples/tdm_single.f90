program main

    use mpi
    use PaScaL_TDMA

    implicit none

    integer, parameter  :: N = 100
    integer :: nprocs, myrank, ierr
    integer :: i, n_sub
    integer :: para_range_n
    logical :: is_root = .false.

    double precision, allocatable, dimension(:) :: d, x
    double precision, allocatable, dimension(:) :: a_sub, b_sub, c_sub, d_sub, x_sub
    integer, allocatable, dimension(:) :: cnt, disp
    type(ptdma_plan_single) :: px_single   ! Plan for a single tridiagonal system of equations

    call MPI_Init(ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)

    if (myrank.eq.0) then
        is_root = .true.
    endif

    n_sub = para_range_n(1, N, nprocs, myrank)

    call build_cnt_disp_array()

    if (is_root) then
        allocate (d(N)); d(:) = 0
        allocate (x(N)); x(:) = 0
        call build_global_coeff_array()
    else
        allocate ( d(0) )
        allocate ( x(0) )
    endif

    call build_local_coeff_array()

    ! Scatter rhs vector
    call MPI_Scatterv(d, cnt, disp, MPI_DOUBLE_PRECISION, d_sub, n_sub, &
                      MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Solve equation
    call PaScaL_TDMA_plan_single_create(px_single, myrank, nprocs, MPI_COMM_WORLD, 0)
    call PaScaL_TDMA_single_solve(px_single, a_sub, b_sub, c_sub, d_sub, n_sub)
    call PaScaL_TDMA_plan_single_destroy(px_single)

    ! Gather solution and evaluate norm2
    call MPI_Gatherv(d_sub, n_sub, MPI_DOUBLE_PRECISION, d, cnt, disp, &
                      MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    if (is_root) then
        print *, "Avg. norm2 ( norm2 / N )= ", norm2(d - x) / N
    endif

    call dealloc_all()

    call MPI_Finalize(ierr)

contains

    subroutine build_cnt_disp_array

        allocate ( cnt(nprocs) ); cnt(:) = 0
        allocate ( disp(nprocs) ); disp(:) = 0
    
        ! Build cnt and disp array
        call MPI_Gather(n_sub, 1, MPI_INTEGER, cnt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

        if (myrank.eq.0) then
            disp(1) = 0
            do i = 2, size(cnt)
                disp(i) = disp(i - 1) + cnt(i - 1)
            end do
        endif

    end subroutine build_cnt_disp_array

    subroutine build_local_coeff_array

        allocate (a_sub(n_sub)); a_sub(:) = 1
        allocate (b_sub(n_sub)); b_sub(:) = 2
        allocate (c_sub(n_sub)); c_sub(:) = 1
        allocate (d_sub(n_sub)); d_sub(:) = 0
        allocate (x_sub(n_sub)); x_sub(:) = 0

    end subroutine build_local_coeff_array

    subroutine build_global_coeff_array

        double precision, allocatable, dimension(:) :: a, b, c

        allocate (a(N)); a(:) = 1
        allocate (b(N)); b(:) = 2
        allocate (c(N)); c(:) = 1
    
        ! Generate random x vector and rhs vector in rank 0
        call random_number(x)
        d(1) = b(1) * x(1) + c(1) * x(2)
        do i = 2, N-1
            d(i) = a(i) * x(i - 1) + b(i) * x(i) + c(i) * x(i + 1)
        enddo
        d(N) = a(N) * x(N - 1) + b(N) * x(N)
        deallocate (a, b, c)

    end subroutine build_global_coeff_array

    subroutine dealloc_all

        deallocate (d, x)
        deallocate (cnt, disp)
        deallocate (a_sub, b_sub, c_sub, d_sub, x_sub)

    end subroutine dealloc_all

end program main

integer function para_range_n(n1, n2, nprocs, myrank) result(n)

    implicit none

    integer, intent(in)     :: n1, n2, nprocs, myrank
    integer :: remainder

    n = int((n2 - n1 + 1) / nprocs)
    remainder = mod(n2 - n1 + 1, nprocs)
    if (remainder > myrank) then
        n = n + 1
    endif 

end function para_range_n