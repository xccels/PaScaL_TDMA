program main

    use omp_lib
    use mpi
    use PaScaL_TDMA

    implicit none

    integer :: nx = 20, ny = 30, nz = 30
    integer :: ny_sub
    integer :: nprocs, myrank, ierr, errorcode, nthds
    integer :: para_range_n
    integer :: kk
    logical :: is_root = .false.

    double precision, allocatable, dimension(:,:,:) :: d_sub, x_sub
    double precision, allocatable, dimension(:,:) :: ap, bp, cp, dp

    type(ptdma_plan_many_thread_team) :: py_many   ! Plan for many tridiagonal systems of equations

    call MPI_Init(ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)

    if (myrank.eq.0) then
        is_root = .true.
    endif

!$omp parallel
    nthds = omp_get_num_threads()
!$omp end parallel

    if ( mod(nz, nthds).ne.0 ) then
        if (myrank.eq.0) then
            print '(a, i4, a, i6)', "n_threads = ", nthds, ", nz = ", nz
            print '(a)', "n_threads should be divisor of nz."
            call MPI_Abort(MPI_COMM_WORLD, errorcode, ierr)
        endif
    endif

    ny_sub = para_range_n(1, ny, nprocs, myrank)

    call build_local_rhs_array()
    call PaScaL_TDMA_plan_many_create_thread_team(py_many, nx, myrank, nprocs, MPI_COMM_WORLD)

!$omp parallel private(ap, bp, cp, dp) default(shared)
    ! Solve equation in y-direction
    allocate ( ap(nx, ny_sub) );
    allocate ( bp(nx, ny_sub) );
    allocate ( cp(nx, ny_sub) );
    allocate ( dp(nx, ny_sub) );

!$omp do private(kk)
    do kk = 1, nz
        ap(:,:) = -1
        bp(:,:) =  2
        cp(:,:) = -1
        dp(:,:) = x_sub(:, :, kk)
        call PaScaL_TDMA_many_solve_thread_team(py_many, ap, bp, cp, dp, nx, ny_sub)
        x_sub(:, :, kk) = dp(:, :)
    enddo
!$omp end do
    deallocate( ap, bp, cp, dp )
!$omp end parallel

    call PaScaL_TDMA_plan_many_destroy_thread_team(py_many, nprocs)

    call check_norm_error()
    call dealloc_all()
    call MPI_Finalize(ierr)

contains

    subroutine build_local_rhs_array

        integer :: m
        integer, allocatable, dimension(:)  :: seed, old

        allocate ( d_sub(nx, ny_sub, nz) ); d_sub(:,:,:) = 0
        allocate ( x_sub(nx, ny_sub, nz) ); x_sub(:,:,:) = 0

        ! Generate random rhs
        call random_seed (SIZE = m)
        allocate(seed(m), old(m))   
        seed = myrank * 20240101
        call random_seed (PUT = seed)   ! Sets user seed
        call random_seed (GET = old)
        call random_number(d_sub)

        deallocate(seed, old)   

        x_sub = d_sub

    end subroutine build_local_rhs_array

    subroutine check_norm_error()

        double precision, allocatable   :: d_compute(:,:,:), x_halo(:,:,:)
        double precision    :: local_norm2, global_norm2
        double precision    :: send_buf_y0(nx*nz), recv_buf_y0(nx*nz)
        double precision    :: send_buf_y1(nx*nz), recv_buf_y1(nx*nz)

        integer             :: i, j, k, stride
        integer             :: next_rank, prev_rank
        integer             :: requests(4)
        integer             :: statuses(MPI_STATUS_SIZE,4)

        allocate ( ap(nx, ny_sub) ); ap = -1
        allocate ( bp(nx, ny_sub) ); bp =  2
        allocate ( cp(nx, ny_sub) ); cp = -1
        ! call initialize_local_coeff_array_tdm()

        allocate( d_compute(nx, ny_sub, nz) )
        allocate( x_halo(nx, 0:ny_sub+1, nz) )

        d_compute       = 0.0d0
        x_halo          = 0.0d0
        x_halo(1:nx, 1:ny_sub, 1:nz) = x_sub(1:nx, 1:ny_sub, 1:nz)

        ! Fill halo cells in y-direction
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

        do k = 1, nz
            stride = (k-1)*nx
            do i = 1, nx
                send_buf_y0(i+stride) = x_halo(i,1,k)
                send_buf_y1(i+stride) = x_halo(i,ny_sub,k)
            enddo
        enddo

        recv_buf_y0 = 0.0d0
        recv_buf_y1 = 0.0d0

        call MPI_Irecv( recv_buf_y0, nx * nz, MPI_DOUBLE_PRECISION, next_rank, 100, MPI_COMM_WORLD, requests(1), ierr)
        call MPI_Isend( send_buf_y1, nx * nz, MPI_DOUBLE_PRECISION, next_rank, 200, MPI_COMM_WORLD, requests(2), ierr)
        call MPI_Isend( send_buf_y0, nx * nz, MPI_DOUBLE_PRECISION, prev_rank, 100, MPI_COMM_WORLD, requests(3), ierr)
        call MPI_Irecv( recv_buf_y1, nx * nz, MPI_DOUBLE_PRECISION, prev_rank, 200, MPI_COMM_WORLD, requests(4), ierr)
        call MPI_Waitall( 4, requests, statuses, ierr)

        do k = 1, nz
            stride = (k-1)*nx
            do i = 1, nx
                x_halo(i,0,k) = recv_buf_y1(i+stride)
                x_halo(i,ny_sub+1,k) = recv_buf_y0(i+stride)
            enddo
        enddo

        do k = 1, nz
            do j = 1, ny_sub
                do i = 1, nx
                    d_compute(i,j,k)    = ap(i,j) * x_halo(i,j-1,k) &
                                        + bp(i,j) * x_halo(i,j,k)   &
                                        + cp(i,j) * x_halo(i,j+1,k)
                enddo
            enddo
        enddo

        local_norm2 = norm2(d_compute - d_sub)
        call MPI_reduce(local_norm2, global_norm2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        
        if (is_root) then
            print *, "Avg. norm2 ( norm2 / N )= ", global_norm2 / (nx*ny*nz)
        endif

        deallocate( ap, bp, cp )
        deallocate( x_halo )
        deallocate( d_compute )

    end subroutine check_norm_error

    subroutine dealloc_all

        deallocate (d_sub, x_sub)
    
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