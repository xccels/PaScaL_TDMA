program main

    use mpi
    use PaScaL_TDMA

    implicit none

    integer, parameter  :: N = 16
    integer :: nprocs, myrank, ierr
    integer :: i, n_sub
    integer :: para_range_n
    logical :: is_root = .false.
    
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
    call build_local_coeff_array()

    ! Solve equation
    call PaScaL_TDMA_plan_single_create(px_single, myrank, nprocs, MPI_COMM_WORLD, 0)
    call PaScaL_TDMA_single_solve(px_single, a_sub, b_sub, c_sub, x_sub, n_sub)
    call PaScaL_TDMA_plan_single_destroy(px_single)

    ! Evaluate norm2 for solution check
    call check_norm_error()

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

        integer :: m
        integer, allocatable, dimension(:)  :: seed, old
    
        allocate (a_sub(n_sub)); a_sub(:) = -1.0d0
        allocate (b_sub(n_sub)); b_sub(:) =  2.0d0
        allocate (c_sub(n_sub)); c_sub(:) = -1.0d0
        allocate (d_sub(n_sub)); d_sub(:) =  0.0d0
        allocate (x_sub(n_sub)); x_sub(:) =  0.0d0

    ! Generate random rhs
        call random_seed (SIZE = m)
        allocate(seed(m), old(m))   
        seed = myrank * 20240101
        call random_seed (PUT = seed)   ! Sets user seed
        call random_seed (GET = old)
        call random_number(d_sub)

        deallocate(seed, old)   

        x_sub = d_sub
    
        end subroutine build_local_coeff_array

    subroutine initialize_local_coeff_array_tdm
    
        a_sub(:) = -1.0d0
        b_sub(:) =  2.0d0
        c_sub(:) = -1.0d0

    end subroutine initialize_local_coeff_array_tdm

    subroutine dealloc_all

        deallocate (cnt, disp)
        deallocate (a_sub, b_sub, c_sub, d_sub, x_sub)

    end subroutine dealloc_all

    subroutine check_norm_error()

        double precision, allocatable   :: d_compute(:), x_halo(:)
        double precision    :: local_norm2, global_norm2
        integer             :: next_rank, prev_rank
        integer             :: requests(4)
        integer             :: statuses(MPI_STATUS_SIZE,4)

        allocate( d_compute(n_sub) )
        allocate( x_halo(0:n_sub+1) )

        d_compute       = 0.0d0
        x_halo          = 0.0d0
        x_halo(1:n_sub) = x_sub(1:n_sub)

        if( myrank == 0 ) then
            prev_rank = MPI_PROC_NULL
        else
            prev_rank = myrank - 1
        endif

        if( myrank == nprocs - 1 ) then
            next_rank = MPI_PROC_NULL
        else
            next_rank = myrank + 1
        endif

        call MPI_Irecv( x_halo(n_sub+1), 1, MPI_DOUBLE_PRECISION, next_rank, 100, MPI_COMM_WORLD, requests(1), ierr)
        call MPI_Isend( x_halo(n_sub),   1, MPI_DOUBLE_PRECISION, next_rank, 200, MPI_COMM_WORLD, requests(2), ierr)
        call MPI_Isend( x_halo(1),       1, MPI_DOUBLE_PRECISION, prev_rank, 100, MPI_COMM_WORLD, requests(3), ierr)
        call MPI_Irecv( x_halo(0),       1, MPI_DOUBLE_PRECISION, prev_rank, 200, MPI_COMM_WORLD, requests(4), ierr)

        call MPI_Waitall( 4, requests, statuses, ierr)

        call initialize_local_coeff_array_tdm()

        do i = 1, n_sub
            d_compute(i) = a_sub(i) * x_halo(i-1) + b_sub(i) * x_halo(i) + c_sub(i) * x_halo(i+1) 
        enddo

        local_norm2 = norm2(d_compute - d_sub)
        call MPI_reduce(local_norm2, global_norm2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        
        if (is_root) then
            print *, "Avg. norm2 ( norm2 / N )= ", global_norm2 / N
        endif

        deallocate( x_halo )
        deallocate( d_compute )

    end subroutine check_norm_error

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