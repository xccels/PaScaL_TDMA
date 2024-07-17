program main

    use mpi
    use PaScaL_TDMA
    use mpi_topology_2D

    implicit none

    integer :: Nx = 20, Ny = 30
    integer :: nx_sub, ny_sub
    integer :: nprocs, myrank, ierr
    integer :: para_range_n
    logical :: is_root = .false.

    double precision, allocatable, dimension(:,:) :: a_sub, b_sub, c_sub, d_sub, x_sub_tr, x_sub

    type(ptdma_plan_many) :: px_many, py_many   ! Plan for many tridiagonal systems of equations

    call MPI_Init(ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)

    if (myrank.eq.0) then
        is_root = .true.
    endif

    call MPI_Dims_create(nprocs, 2, np_dim, ierr)

    period(0) = .false.
    period(1) = .false.

    call mpi_topology_make

    nx_sub = para_range_n(1, nx, comm_1d_x%nprocs, comm_1d_x%myrank)
    ny_sub = para_range_n(1, ny, comm_1d_y%nprocs, comm_1d_y%myrank)

    call build_local_coeff_array()

    ! Solve equation in y-direction

    call PaScaL_TDMA_plan_many_create(py_many, nx_sub, comm_1d_y%myrank, comm_1d_y%nprocs, comm_1d_y%mpi_comm)
    call PaScaL_TDMA_many_solve(py_many, a_sub, b_sub, c_sub, x_sub, nx_sub, ny_sub)
    call PaScaL_TDMA_plan_many_destroy(py_many, comm_1d_y%nprocs)

    ! Solve equation in x-direction
    allocate ( x_sub_tr(ny_sub, nx_sub) ); x_sub_tr(:,:) = 0
    x_sub_tr = transpose( x_sub )

    call initialize_local_coeff_array_tdm()
    call PaScaL_TDMA_plan_many_create(px_many, ny_sub, comm_1d_x%myrank, comm_1d_x%nprocs, comm_1d_x%mpi_comm)
    call PaScaL_TDMA_many_solve(px_many, a_sub, b_sub, c_sub, x_sub_tr, ny_sub, nx_sub)
    call PaScaL_TDMA_plan_many_destroy(px_many, comm_1d_x%nprocs)

    x_sub = transpose( x_sub_tr )
    deallocate ( x_sub_tr )

    call check_norm_error()
 
    call dealloc_all()
    call mpi_topology_clean
    call MPI_Finalize(ierr)

contains

    subroutine build_local_coeff_array

        integer :: m
        integer, allocatable, dimension(:)  :: seed, old
    
        allocate ( a_sub(nx_sub, ny_sub) ); a_sub(:,:) = -1
        allocate ( b_sub(nx_sub, ny_sub) ); b_sub(:,:) =  2
        allocate ( c_sub(nx_sub, ny_sub) ); c_sub(:,:) = -1

        allocate ( d_sub(nx_sub, ny_sub) ); d_sub(:,:) = 0
        allocate ( x_sub(nx_sub, ny_sub) ); x_sub(:,:) = 0

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
    
        a_sub(:,:) = -1
        b_sub(:,:) =  2
        c_sub(:,:) = -1

    end subroutine initialize_local_coeff_array_tdm

    subroutine check_norm_error()

        double precision, allocatable   :: d_compute(:,:), x_halo(:,:)
        double precision    :: local_norm2, global_norm2
        double precision    :: send_buf_x0(ny_sub), recv_buf_x0(ny_sub)
        double precision    :: send_buf_x1(ny_sub), recv_buf_x1(ny_sub)
        double precision    :: send_buf_y0(nx_sub), recv_buf_y0(nx_sub)
        double precision    :: send_buf_y1(nx_sub), recv_buf_y1(nx_sub)

        integer             :: i, j
        integer             :: next_rank, prev_rank
        integer             :: requests(4)
        integer             :: statuses(MPI_STATUS_SIZE,4)

        call initialize_local_coeff_array_tdm()

        allocate( d_compute(nx_sub, ny_sub) )
        allocate( x_halo(0:nx_sub+1, 0:ny_sub+1) )

        d_compute       = 0.0d0
        x_halo          = 0.0d0
        x_halo(1:nx_sub, 1:ny_sub) = x_sub(1:nx_sub, 1:ny_sub)

        ! Fill halo cells in x-direction
        if( comm_1d_x%myrank == 0 ) then
            prev_rank = MPI_PROC_NULL
        else
            prev_rank = comm_1d_x%west_rank
        endif

        if( comm_1d_x%myrank == comm_1d_x%nprocs - 1 ) then
            next_rank = MPI_PROC_NULL
        else
            next_rank = comm_1d_x%east_rank
        endif

        do j = 1, ny_sub
            send_buf_x0(j) = x_halo(1,j)
            send_buf_x1(j) = x_halo(nx_sub,j)
        enddo

        recv_buf_x0 = 0.0d0
        recv_buf_x1 = 0.0d0

        call MPI_Irecv( recv_buf_x0, ny_sub, MPI_DOUBLE_PRECISION, next_rank, 100, comm_1d_x%mpi_comm, requests(1), ierr)
        call MPI_Isend( send_buf_x1, ny_sub, MPI_DOUBLE_PRECISION, next_rank, 200, comm_1d_x%mpi_comm, requests(2), ierr)
        call MPI_Isend( send_buf_x0, ny_sub, MPI_DOUBLE_PRECISION, prev_rank, 100, comm_1d_x%mpi_comm, requests(3), ierr)
        call MPI_Irecv( recv_buf_x1, ny_sub, MPI_DOUBLE_PRECISION, prev_rank, 200, comm_1d_x%mpi_comm, requests(4), ierr)
        call MPI_Waitall( 4, requests, statuses, ierr)

        do j = 1, ny_sub
            x_halo(0,j) = recv_buf_x1(j)
            x_halo(nx_sub+1,j) = recv_buf_x0(j)
        enddo

        do j = 1, ny_sub
            do i = 1, nx_sub
                d_compute(i,j)  = a_sub(i,j) * x_halo(i-1,j) &
                                + b_sub(i,j) * x_halo(i,j)   &
                                + c_sub(i,j) * x_halo(i+1,j)
            enddo
        enddo

        x_halo(1:nx_sub, 1:ny_sub) = d_compute(1:nx_sub, 1:ny_sub)

        ! Fill halo cells in y-direction
        if( comm_1d_y%myrank == 0 ) then
            prev_rank = MPI_PROC_NULL
        else
            prev_rank = comm_1d_y%myrank - 1
        endif

        if( comm_1d_y%myrank == comm_1d_y%nprocs - 1 ) then
            next_rank = MPI_PROC_NULL
        else
            next_rank = comm_1d_y%myrank + 1
        endif

        do i = 1, nx_sub
            send_buf_y0(i) = x_halo(i,1)
            send_buf_y1(i) = x_halo(i,ny_sub)
        enddo

        recv_buf_y0 = 0.0d0
        recv_buf_y1 = 0.0d0

        call MPI_Irecv( recv_buf_y0, nx_sub, MPI_DOUBLE_PRECISION, next_rank, 100, comm_1d_y%mpi_comm, requests(1), ierr)
        call MPI_Isend( send_buf_y1, nx_sub, MPI_DOUBLE_PRECISION, next_rank, 200, comm_1d_y%mpi_comm, requests(2), ierr)
        call MPI_Isend( send_buf_y0, nx_sub, MPI_DOUBLE_PRECISION, prev_rank, 100, comm_1d_y%mpi_comm, requests(3), ierr)
        call MPI_Irecv( recv_buf_y1, nx_sub, MPI_DOUBLE_PRECISION, prev_rank, 200, comm_1d_y%mpi_comm, requests(4), ierr)
        call MPI_Waitall( 4, requests, statuses, ierr)

        do i = 1, nx_sub
            x_halo(i,0) = recv_buf_y1(i)
            x_halo(i,ny_sub+1) = recv_buf_y0(i)
        enddo

        do j = 1, ny_sub
            do i = 1, nx_sub
                d_compute(i,j)  = a_sub(i,j) * x_halo(i,j-1) &
                                + b_sub(i,j) * x_halo(i,j)   &
                                + c_sub(i,j) * x_halo(i,j+1)
            enddo
        enddo

        local_norm2 = norm2(d_compute - d_sub)
        call MPI_reduce(local_norm2, global_norm2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        
        if (is_root) then
            print *, "Avg. norm2 ( norm2 / N )= ", global_norm2 / (nx * ny)
        endif

        deallocate( x_halo )
        deallocate( d_compute )

    end subroutine check_norm_error

    subroutine dealloc_all

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