program main

    use mpi
    use PaScaL_TDMA
    use mpi_topology_2D

    implicit none

    integer :: Nx = 100, Ny = 100
    integer :: nx_sub, ny_sub, n_sub
    integer :: nprocs, myrank, ierr
    integer :: npx
    integer :: para_range_n
    integer, allocatable, dimension(:) :: cnt_x, disp_x, cnt_y, disp_y, cnt_all, disp_all
    logical :: is_root = .false.

    double precision, allocatable, dimension(:,:) :: d, x
    double precision, allocatable, dimension(:,:) :: a_sub, b_sub, c_sub, d_sub, d_sub_tr

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

    npx = np_dim(0)

    nx_sub = para_range_n(1, nx, comm_1d_x%nprocs, comm_1d_x%myrank)
    ny_sub = para_range_n(1, ny, comm_1d_y%nprocs, comm_1d_y%myrank)

    n_sub = nx_sub * ny_sub

    call build_comm_info_array()

    ! Generate random x vector and rhs vector in rank 0
    if (is_root) then
        allocate ( d(Nx, Ny) ); d(:,:) = 0
        allocate ( x(Nx, Ny) ); x(:,:) = 0
        call build_global_coeff_array()
    endif

    call distribute_rhs_array()

    ! Solve equation in y-direction
    allocate ( a_sub(nx_sub, ny_sub) ); a_sub(:,:) = 1
    allocate ( b_sub(nx_sub, ny_sub) ); b_sub(:,:) = 2
    allocate ( c_sub(nx_sub, ny_sub) ); c_sub(:,:) = 1

    call PaScaL_TDMA_plan_many_create(py_many, nx_sub, comm_1d_y%myrank, comm_1d_y%nprocs, comm_1d_y%mpi_comm)
    call PaScaL_TDMA_many_solve(py_many, a_sub, b_sub, c_sub, d_sub, nx_sub, ny_sub)
    call PaScaL_TDMA_plan_many_destroy(py_many, comm_1d_y%nprocs)

    ! Solve equation in x-direction
    a_sub(:,:) = 1
    b_sub(:,:) = 2
    c_sub(:,:) = 1

    allocate ( d_sub_tr(ny_sub, nx_sub) ); d_sub_tr(:,:) = 0

    d_sub_tr = transpose( d_sub )

    call PaScaL_TDMA_plan_many_create(px_many, ny_sub, comm_1d_x%myrank, comm_1d_x%nprocs, comm_1d_x%mpi_comm)
    call PaScaL_TDMA_many_solve(px_many, a_sub, b_sub, c_sub, d_sub_tr, ny_sub, nx_sub)
    call PaScaL_TDMA_plan_many_destroy(px_many, comm_1d_x%nprocs)

    d_sub = transpose( d_sub_tr )

    deallocate ( d_sub_tr )

    call collect_solution_array()

    if (is_root) then
        print *, "Avg. norm2 (norm2 / (nx * ny)) = ", norm2(d - x) / nx / ny
    endif

    call dealloc_all()

    call mpi_topology_clean

    call MPI_Finalize(ierr)

contains

    subroutine build_comm_info_array

        integer :: i

        allocate ( cnt_x(np_dim(0)) );  cnt_x(:) = 0
        allocate ( cnt_y(np_dim(1)) );  cnt_y(:) = 0
        allocate ( cnt_all(nprocs) );   cnt_all(:) = 0
        allocate ( disp_x(np_dim(0)) ); disp_x(:) = 0
        allocate ( disp_y(np_dim(1)) ); disp_y(:) = 0
        allocate ( disp_all(nprocs) );  disp_all(:) = 0
    
        ! Build cnt and disp array
        call MPI_Allgather(nx_sub, 1, MPI_INTEGER, cnt_x,   1, MPI_INTEGER, comm_1d_x%mpi_comm, ierr)
        call MPI_Allgather(ny_sub, 1, MPI_INTEGER, cnt_y,   1, MPI_INTEGER, comm_1d_y%mpi_comm, ierr)
        call MPI_Allgather(n_sub,  1, MPI_INTEGER, cnt_all, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    
        disp_x(1) = 0
        do i = 2, size(cnt_x)
            disp_x(i) = disp_x(i - 1) + cnt_x(i - 1)
        enddo
    
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

        integer :: i, j
        double precision, allocatable, dimension(:,:) :: a, b, c, y

        allocate ( a(Nx, Ny) ); a(:,:) = 1
        allocate ( b(Nx, Ny) ); b(:,:) = 2
        allocate ( c(Nx, Ny) ); c(:,:) = 1
        allocate ( y(Nx, Ny) ); y(:,:) = 0
    
        call random_number(x(:,:))

        ! y = A_x * x
        do j = 1, ny
            y(1, j) = b(1, j) * x(1, j) + c(1, j) * x(2, j)
            do i = 2, nx - 1
                y(i, j) = a(i, j) * x(i - 1, j) + b(i, j) * x(i, j) + c(i, j) * x(i + 1, j)
            enddo
            y(nx, j) = a(nx, j) * x(nx - 1, j) + b(nx, j) * x(nx, j)
        enddo

        ! d = A_y * y
        do i = 1, nx
            d(i, 1) = b(i, 1) * y(i, 1) + c(i, 1) * y(i, 2)
        enddo
        do j = 2, ny - 1
            do i = 1, nx
                d(i, j) = a(i, j) * y(i, j - 1) + b(i, j) * y(i, j) + c(i, j) * y(i, j + 1)
            enddo
        enddo
        do i = 1, nx
            d(i, ny) = a(i, ny) * y(i, ny - 1) + b(i, ny) * y(i, ny)
        enddo
        deallocate (a, b, c, y)

    end subroutine build_global_coeff_array

    subroutine distribute_rhs_array

        integer :: i, j, iblk
        double precision, allocatable, dimension(:)   :: d_blk

        ! Scatter rhs vector
        if (is_root) then
            allocate ( d_blk(Nx * Ny) ); d_blk(:) = 0
            do iblk = 1, npx
                do j = 1, ny
                    do i = 1, cnt_x(iblk)
                        d_blk(i + (j - 1) * cnt_x(iblk) + disp_x(iblk) * ny) &
                            = d(i + disp_x(iblk), j)
                    enddo
                enddo
            enddo
        else
            allocate ( d_blk(0) )
        endif

        allocate ( d_sub(nx_sub, ny_sub) ); d_sub(:,:) = 0
        call MPI_Scatterv(d_blk, cnt_all, disp_all, MPI_DOUBLE_PRECISION, d_sub, n_sub, &
                            MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        deallocate( d_blk )

    end subroutine distribute_rhs_array

    subroutine collect_solution_array

        integer :: i, j, iblk
        double precision, allocatable, dimension(:)   :: d_blk

        if (is_root) then
            allocate ( d_blk(Nx * Ny) ); d_blk(:) = 0
        else
            allocate ( d_blk(0) )
        endif

        ! Gather solution and evaluate norm2
        call MPI_Gatherv(d_sub, n_sub, MPI_DOUBLE_PRECISION, d_blk, cnt_all, disp_all, &
                        MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

        if (is_root) then
            do iblk = 1, npx
                do j = 1, ny
                    do i = 1, cnt_x(iblk)
                        d(i + disp_x(iblk), j) = &
                            d_blk(i + (j - 1) * cnt_x(iblk) + disp_x(iblk) * ny)
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
        deallocate (a_sub, b_sub, c_sub, d_sub)
        deallocate (cnt_x, disp_x)
        deallocate (cnt_y, disp_y)
        deallocate (cnt_all, disp_all)
    
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