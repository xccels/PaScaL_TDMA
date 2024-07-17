program main

    use mpi
    use PaScaL_TDMA
    use mpi_topology_2D

    implicit none

    integer :: nx = 20, ny = 40, nz = 30
    integer :: nx_sub, ny_sub
    integer :: nprocs, myrank, ierr
    integer :: para_range_n
    integer :: npx, npy
    logical :: is_root = .false.

    double precision, allocatable, dimension(:) :: ax_sub, bx_sub, cx_sub
    double precision, allocatable, dimension(:) :: ay_sub, by_sub, cy_sub
    double precision, allocatable, dimension(:) :: az, bz, cz
    double precision, allocatable, dimension(:,:,:) :: d_sub, x_sub_tr, x_sub

    type(ptdma_plan_many_rhs) :: px_many_rhs, py_many_rhs   ! Plan for many tridiagonal systems of equations

    call MPI_Init(ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)

    if (myrank.eq.0) then
        is_root = .true.
    endif

    call MPI_Dims_create(nprocs, 2, np_dim, ierr)

    period(0) = .false.
    period(1) = .false.

    call mpi_topology_make()

    npx = np_dim(0)
    npy = np_dim(1)

    nx_sub = para_range_n( 1, nx, comm_1d_x%nprocs, comm_1d_x%myrank )
    ny_sub = para_range_n( 1, ny, comm_1d_y%nprocs, comm_1d_y%myrank )

    call build_local_coeff_array()

    ! Solve equation in z-direction
    call tdma_many_rhs(az, bz, cz, x_sub, nx_sub * ny_sub, nz)

    ! Solve equation in y-direction
    allocate ( x_sub_tr(nx_sub, nz, ny_sub) ); x_sub_tr(:,:,:) = 0
    x_sub_tr = reshape (x_sub, shape(x_sub_tr), order = [1, 3, 2] )

    call PaScaL_TDMA_plan_many_rhs_create(py_many_rhs, nz * nx_sub, comm_1d_y%myrank, comm_1d_y%nprocs, comm_1d_y%mpi_comm)
    call PaScaL_TDMA_many_rhs_solve(py_many_rhs, ay_sub, by_sub, cy_sub, x_sub_tr, nz * nx_sub, ny_sub)
    call PaScaL_TDMA_plan_many_rhs_destroy(py_many_rhs, comm_1d_y%nprocs)

    x_sub = reshape (x_sub_tr, shape(x_sub), order = [1, 3, 2] )
    deallocate( x_sub_tr )

    ! Solve equation in x-direction
    allocate ( x_sub_tr(ny_sub, nz, nx_sub) ); x_sub_tr(:,:,:) = 0
    x_sub_tr = reshape (x_sub, shape(x_sub_tr), order = [3, 1, 2] )

    call PaScaL_TDMA_plan_many_rhs_create(px_many_rhs, ny_sub * nz, comm_1d_x%myrank, comm_1d_x%nprocs, comm_1d_x%mpi_comm)
    call PaScaL_TDMA_many_rhs_solve(px_many_rhs, ax_sub, bx_sub, cx_sub, x_sub_tr, ny_sub * nz, nx_sub)
    call PaScaL_TDMA_plan_many_rhs_destroy(px_many_rhs, comm_1d_x%nprocs)

    x_sub = reshape (x_sub_tr, shape(x_sub), order = [2, 3, 1] )
    deallocate( x_sub_tr )

    call check_norm_error()

    call mpi_topology_clean()
    call dealloc_all()

    call MPI_Finalize(ierr)

contains

    subroutine build_local_coeff_array

        integer :: m
        integer, allocatable, dimension(:)  :: seed, old
    
        allocate ( ax_sub(nx_sub) ); ax_sub(:) = -1
        allocate ( bx_sub(nx_sub) ); bx_sub(:) =  2
        allocate ( cx_sub(nx_sub) ); cx_sub(:) = -1

        allocate ( ay_sub(ny_sub) ); ay_sub(:) = -1
        allocate ( by_sub(ny_sub) ); by_sub(:) =  2
        allocate ( cy_sub(ny_sub) ); cy_sub(:) = -1

        allocate ( az(nz) ); az(:) = -1
        allocate ( bz(nz) ); bz(:) =  2
        allocate ( cz(nz) ); cz(:) = -1

        allocate ( d_sub(nx_sub, ny_sub, nz) ); d_sub(:,:,:) = 0
        allocate ( x_sub(nx_sub, ny_sub, nz) ); x_sub(:,:,:) = 0

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
    
        ax_sub(:) = -1
        bx_sub(:) =  2
        cx_sub(:) = -1

        ay_sub(:) = -1
        by_sub(:) =  2
        cy_sub(:) = -1

        az(:) = -1
        bz(:) =  2
        cz(:) = -1
        
    end subroutine initialize_local_coeff_array_tdm

    subroutine check_norm_error()

        double precision, allocatable   :: d_compute(:,:,:), x_halo(:,:,:)
        double precision    :: local_norm2, global_norm2
        double precision    :: send_buf_x0(ny_sub*nz), recv_buf_x0(ny_sub*nz)
        double precision    :: send_buf_x1(ny_sub*nz), recv_buf_x1(ny_sub*nz)
        double precision    :: send_buf_y0(nx_sub*nz), recv_buf_y0(nx_sub*nz)
        double precision    :: send_buf_y1(nx_sub*nz), recv_buf_y1(nx_sub*nz)

        integer             :: i, j, k, stride
        integer             :: next_rank, prev_rank
        integer             :: requests(4)
        integer             :: statuses(MPI_STATUS_SIZE,4)

        call initialize_local_coeff_array_tdm()

        allocate( d_compute(nx_sub, ny_sub, nz) )
        allocate( x_halo(0:nx_sub+1, 0:ny_sub+1, 0:nz+1) )

        d_compute       = 0.0d0
        x_halo          = 0.0d0
        x_halo(1:nx_sub, 1:ny_sub, 1:nz) = x_sub(1:nx_sub, 1:ny_sub, 1:nz)

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

        do k = 1, nz
            stride = (k-1)*ny_sub
            do j = 1, ny_sub
                send_buf_x0(j+stride) = x_halo(1,j,k)
                send_buf_x1(j+stride) = x_halo(nx_sub,j,k)
            enddo
        enddo

        recv_buf_x0 = 0.0d0
        recv_buf_x1 = 0.0d0

        call MPI_Irecv( recv_buf_x0, ny_sub * nz, MPI_DOUBLE_PRECISION, next_rank, 100, comm_1d_x%mpi_comm, requests(1), ierr)
        call MPI_Isend( send_buf_x1, ny_sub * nz, MPI_DOUBLE_PRECISION, next_rank, 200, comm_1d_x%mpi_comm, requests(2), ierr)
        call MPI_Isend( send_buf_x0, ny_sub * nz, MPI_DOUBLE_PRECISION, prev_rank, 100, comm_1d_x%mpi_comm, requests(3), ierr)
        call MPI_Irecv( recv_buf_x1, ny_sub * nz, MPI_DOUBLE_PRECISION, prev_rank, 200, comm_1d_x%mpi_comm, requests(4), ierr)
        call MPI_Waitall( 4, requests, statuses, ierr)

        do k = 1, nz
            stride = (k-1)*ny_sub
            do j = 1, ny_sub
                x_halo(0,j,k) = recv_buf_x1(j+stride)
                x_halo(nx_sub+1,j,k) = recv_buf_x0(j+stride)
            enddo
        enddo

        do k = 1, nz
            do j = 1, ny_sub
                do i = 1, nx_sub
                    d_compute(i,j,k)    = ax_sub(i) * x_halo(i-1,j,k) &
                                        + bx_sub(i) * x_halo(i,j,k)   &
                                        + cx_sub(i) * x_halo(i+1,j,k)
                enddo
            enddo
        enddo

        x_halo(1:nx_sub, 1:ny_sub, 1:nz) = d_compute(1:nx_sub, 1:ny_sub, 1:nz)

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

        do k = 1, nz
            stride = (k-1)*nx_sub
            do i = 1, nx_sub
                send_buf_y0(i+stride) = x_halo(i,1,k)
                send_buf_y1(i+stride) = x_halo(i,ny_sub,k)
            enddo
        enddo

        recv_buf_y0 = 0.0d0
        recv_buf_y1 = 0.0d0

        call MPI_Irecv( recv_buf_y0, nx_sub * nz, MPI_DOUBLE_PRECISION, next_rank, 100, comm_1d_y%mpi_comm, requests(1), ierr)
        call MPI_Isend( send_buf_y1, nx_sub * nz, MPI_DOUBLE_PRECISION, next_rank, 200, comm_1d_y%mpi_comm, requests(2), ierr)
        call MPI_Isend( send_buf_y0, nx_sub * nz, MPI_DOUBLE_PRECISION, prev_rank, 100, comm_1d_y%mpi_comm, requests(3), ierr)
        call MPI_Irecv( recv_buf_y1, nx_sub * nz, MPI_DOUBLE_PRECISION, prev_rank, 200, comm_1d_y%mpi_comm, requests(4), ierr)
        call MPI_Waitall( 4, requests, statuses, ierr)

        do k = 1, nz
            stride = (k-1)*nx_sub
            do i = 1, nx_sub
                x_halo(i,0,k) = recv_buf_y1(i+stride)
                x_halo(i,ny_sub+1,k) = recv_buf_y0(i+stride)
            enddo
        enddo

        do k = 1, nz
            do j = 1, ny_sub
                do i = 1, nx_sub
                    d_compute(i,j,k)    = ay_sub(j) * x_halo(i,j-1,k) &
                                        + by_sub(j) * x_halo(i,j,k)   &
                                        + cy_sub(j) * x_halo(i,j+1,k)
                enddo
            enddo
        enddo

        x_halo(1:nx_sub, 1:ny_sub, 1:nz) = d_compute(1:nx_sub, 1:ny_sub, 1:nz)

        ! d = A_z * z
        do k = 1, nz
            do j = 1, ny_sub
                do i = 1, nx_sub
                    d_compute(i, j, k)  = az(k) * x_halo(i, j, k - 1) &
                                        + bz(k) * x_halo(i, j, k) &
                                        + cz(k) * x_halo(i, j, k + 1)
                enddo
            enddo
        enddo

        local_norm2 = norm2(d_compute - d_sub)
        call MPI_reduce(local_norm2, global_norm2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        
        if (is_root) then
            print *, "Avg. norm2 ( norm2 / N )= ", global_norm2 / (nx*ny*nz)
        endif

        deallocate( x_halo )
        deallocate( d_compute )

    end subroutine check_norm_error

    subroutine dealloc_all

        ! if (is_root) then
        !     deallocate (d, x)
        ! endif
        deallocate (ax_sub, bx_sub, cx_sub)
        deallocate (ay_sub, by_sub, cy_sub)
        deallocate (az, bz, cz)
        deallocate (d_sub)

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
