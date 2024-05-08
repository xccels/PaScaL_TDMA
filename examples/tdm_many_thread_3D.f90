program main

    use omp_lib
    use mpi
    use PaScaL_TDMA

    implicit none

    integer :: nx = 20, ny = 20, nz = 20
    integer :: ny_sub, n_sub
    integer :: nprocs, myrank, ierr, errorcode, nthds
    integer :: para_range_n
    integer :: kk
    integer, allocatable, dimension(:) :: cnt_y, disp_y, cnt_all, disp_all
    logical :: is_root = .false.

    double precision, allocatable, dimension(:,:,:) :: d, x
    double precision, allocatable, dimension(:,:,:) :: d_sub
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
    n_sub = nx * ny_sub * nz

    call build_comm_info_array()

    ! Generate random x vector and rhs vector in rank 0
    if (is_root) then
        allocate ( d(nx, ny, nz) ); d(:,:,:) = 0
        allocate ( x(nx, ny, nz) ); x(:,:,:) = 0
        call build_global_coeff_array()
    endif

    call distribute_rhs_array()

    call PaScaL_TDMA_plan_many_create_thread_team(py_many, nx, myrank, nprocs, MPI_COMM_WORLD)

!$omp parallel private(ap, bp, cp, dp) default(shared)
    ! Solve equation in y-direction
    allocate ( ap(nx, ny_sub) );
    allocate ( bp(nx, ny_sub) );
    allocate ( cp(nx, ny_sub) );
    allocate ( dp(nx, ny_sub) );

!$omp do private(kk)
    do kk = 1, nz
        ap(:,:) = 1
        bp(:,:) = 2
        cp(:,:) = 1
        dp(:,:) = d_sub(:, :, kk)
        call PaScaL_TDMA_many_solve_thread_team(py_many, ap, bp, cp, dp, nx, ny_sub)
        d_sub(:, :, kk) = dp(:, :)
    enddo
!$omp end do
    deallocate( ap, bp, cp, dp )
!$omp end parallel

    call PaScaL_TDMA_plan_many_destroy_thread_team(py_many, nprocs)

    call collect_solution_array()

    if (is_root) then
        print *, "Avg. norm2 (norm2 / (nx * ny)) = ", norm2(d - x) / nx / ny
    endif

    call dealloc_all()

    call MPI_Finalize(ierr)

contains

    subroutine build_comm_info_array

        integer :: i

        allocate ( cnt_y(nprocs) );   cnt_y(:) = 0
        allocate ( cnt_all(nprocs) );   cnt_all(:) = 0
        allocate ( disp_y(nprocs) );  disp_y(:) = 0
        allocate ( disp_all(nprocs) );  disp_all(:) = 0

        ! Build cnt and disp array
        call MPI_Allgather(ny_sub, 1, MPI_INTEGER, cnt_y, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
        call MPI_Allgather(n_sub,  1, MPI_INTEGER, cnt_all, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

        disp_y(1) = 0
        do i = 2, size(cnt_y)
            disp_y(i) = disp_y(i - 1) + cnt_y(i - 1)
        enddo

        disp_all(1) = 0
        do i = 2, nprocs
            disp_all(i) = disp_all(i - 1) + cnt_all(i - 1)
        enddo

    end subroutine build_comm_info_array

    subroutine build_global_coeff_array

        integer :: i, j, k
        double precision, allocatable, dimension(:,:,:) :: a, b, c

        allocate ( a(nx, ny, nz) ); a(:,:,:) = 1
        allocate ( b(nx, ny, nz) ); b(:,:,:) = 2
        allocate ( c(nx, ny, nz) ); c(:,:,:) = 1
    
        call random_number(x(:,:,:))

        ! d = A_x * x
        do k = 1, nz
            do i = 1, nx
                d(i, 1, k) = b(i, 1, k) * x(i, 1, k) + c(i, 1, k) * x(i, 2, k)
            enddo
            do j = 2, ny - 1
                do i = 1, nx
                    d(i, j, k) = a(i, j, k) * x(i, j - 1, k) + b(i, j, k) * x(i, j, k) + c(i, j, k) * x(i, j + 1, k)
                enddo
            enddo
            do i = 1, nx
                d(i, ny, k) = a(i, ny, k) * x(i, ny - 1, k) + b(i, ny, k) * x(i, ny, k)
            enddo
        enddo 
        deallocate (a, b, c)

    end subroutine build_global_coeff_array

    subroutine distribute_rhs_array

        integer :: i, j, k, rank
        double precision, allocatable, dimension(:)   :: d_blk

        ! Scatter rhs vector
        if (is_root) then
            allocate ( d_blk(nx * ny * nz) ); d_blk(:) = 0
            do rank = 1, nprocs
                do k = 1, nz
                    do j = 1, cnt_y(rank)
                        do i = 1, nx
                            d_blk(i + (j - 1) * nx + (k - 1) * cnt_y(rank) * nx &
                                + disp_all(rank)) &
                                = d(i, j + disp_y(rank), k)
                        enddo
                    enddo
                enddo
            enddo
        else
            allocate ( d_blk(0) )
        endif

        allocate ( d_sub(nx, ny_sub, nz) ); d_sub(:,:,:) = 0
        call MPI_Scatterv(d_blk, cnt_all, disp_all, MPI_DOUBLE_PRECISION, d_sub, n_sub, &
                            MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        deallocate( d_blk )

    end subroutine distribute_rhs_array

    subroutine collect_solution_array

        integer :: i, j, k, rank
        double precision, allocatable, dimension(:)   :: d_blk

        if (is_root) then
            allocate ( d_blk(nx * ny * nz) ); d_blk(:) = 0
        else
            allocate ( d_blk(0) )
        endif

        ! Gather solution and evaluate norm2
        call MPI_Gatherv(d_sub, n_sub, MPI_DOUBLE_PRECISION, d_blk, cnt_all, disp_all, &
                        MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

        if (is_root) then
            do rank = 1, nprocs
                do k = 1, nz
                    do j = 1, cnt_y(rank)
                        do i = 1, nx
                            d(i, j+ disp_y(rank), k) &
                                = d_blk(i + (j - 1) * nx + (k - 1) * cnt_y(rank) * nx &
                                + nx * nz * disp_y(rank))
                        enddo
                    enddo
                enddo
            enddo
        endif

        deallocate( d_blk )

    end subroutine collect_solution_array

    subroutine dealloc_all

        if (is_root) then
            deallocate (d, x)
        endif
        deallocate (d_sub)
        deallocate (cnt_y, disp_y)
    
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