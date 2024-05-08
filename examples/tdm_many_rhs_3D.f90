program main

    use mpi
    use PaScaL_TDMA
    use mpi_topology_2D

    implicit none

    integer :: nx = 20, ny = 20, nz = 5
    integer :: nx_sub, ny_sub, n_sub
    integer :: nprocs, myrank, ierr
    integer :: para_range_n
    integer :: npx, npy
    logical :: is_root = .false.

    integer, allocatable, dimension(:) :: cnt_x, disp_x
    integer, allocatable, dimension(:) :: cnt_y, disp_y
    integer, allocatable, dimension(:) :: cnt_all, disp_all

    double precision, allocatable, dimension(:,:,:) :: d, x
    double precision, allocatable, dimension(:) :: ax_sub, bx_sub, cx_sub
    double precision, allocatable, dimension(:) :: ay_sub, by_sub, cy_sub
    double precision, allocatable, dimension(:) :: az, bz, cz
    double precision, allocatable, dimension(:,:,:) :: d_sub, d_sub_tr

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

    n_sub = nx_sub * ny_sub * nz

    call build_comm_info_array()

    ! Generate random x vector and rhs vector in rank 0
    if (is_root) then
        allocate ( d(nx, ny, nz) ); d(:,:,:) = 0
        allocate ( x(nx, ny, nz) ); x(:,:,:) = 0
        call build_global_coeff_array()
    endif

    call distribute_rhs_array()

    allocate ( az(nz) ); az(:) = 1
    allocate ( bz(nz) ); bz(:) = 2
    allocate ( cz(nz) ); cz(:) = 1

    call tdma_many_rhs(az, bz, cz, d_sub, nx_sub * ny_sub, nz)

    allocate ( ay_sub(ny_sub) ); ay_sub(:) = 1
    allocate ( by_sub(ny_sub) ); by_sub(:) = 2
    allocate ( cy_sub(ny_sub) ); cy_sub(:) = 1
    allocate ( d_sub_tr(nx_sub, nz, ny_sub) ); d_sub_tr(:,:,:) = 0

    d_sub_tr = reshape (d_sub, shape(d_sub_tr), order = [1, 3, 2] )

    ! Solve equation in y-direction
    call PaScaL_TDMA_plan_many_rhs_create(py_many_rhs, nz * nx_sub, comm_1d_y%myrank, comm_1d_y%nprocs, comm_1d_y%mpi_comm)
    call PaScaL_TDMA_many_rhs_solve(py_many_rhs, ay_sub, by_sub, cy_sub, d_sub_tr, nz * nx_sub, ny_sub)
    call PaScaL_TDMA_plan_many_rhs_destroy(py_many_rhs, comm_1d_y%nprocs)

    d_sub = reshape (d_sub_tr, shape(d_sub), order = [1, 3, 2] )

    deallocate( d_sub_tr )

    ! Solve equation in x-direction
    allocate ( ax_sub(nx_sub) ); ax_sub(:) = 1
    allocate ( bx_sub(nx_sub) ); bx_sub(:) = 2
    allocate ( cx_sub(nx_sub) ); cx_sub(:) = 1

    allocate ( d_sub_tr(ny_sub, nz, nx_sub) ); d_sub_tr(:,:,:) = 0

    d_sub_tr = reshape (d_sub, shape(d_sub_tr), order = [3, 1, 2] )

    call PaScaL_TDMA_plan_many_rhs_create(px_many_rhs, ny_sub * nz, comm_1d_x%myrank, comm_1d_x%nprocs, comm_1d_x%mpi_comm)
    call PaScaL_TDMA_many_rhs_solve(px_many_rhs, ax_sub, bx_sub, cx_sub, d_sub_tr, ny_sub * nz, nx_sub)
    call PaScaL_TDMA_plan_many_rhs_destroy(px_many_rhs, comm_1d_x%nprocs)

    d_sub = reshape (d_sub_tr, shape(d_sub), order = [2, 3, 1] )

    call collect_solution_array()

    if (is_root) then
        print *, "Avg. norm2 (norm2 / (nx * ny)) = ", norm2(d - x) / nx / ny / nz
    endif

    call mpi_topology_clean()
    call dealloc_all()

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

        integer :: i, j, k
        double precision, allocatable, dimension(:,:,:) :: y, z
        double precision, allocatable, dimension(:) :: a_x, b_x, c_x
        double precision, allocatable, dimension(:) :: a_y, b_y, c_y
        double precision, allocatable, dimension(:) :: a_z, b_z, c_z

        allocate ( a_x(nx) ); a_x(:) = 1
        allocate ( b_x(nx) ); b_x(:) = 2
        allocate ( c_x(nx) ); c_x(:) = 1
    
        allocate ( a_y(ny) ); a_y(:) = 1
        allocate ( b_y(ny) ); b_y(:) = 2
        allocate ( c_y(ny) ); c_y(:) = 1
    
        allocate ( a_z(nz) ); a_z(:) = 1
        allocate ( b_z(nz) ); b_z(:) = 2
        allocate ( c_z(nz) ); c_z(:) = 1

        allocate ( y(nx, ny, nz) ); y(:,:,:) = 0
        allocate ( z(nx, ny, nz) ); z(:,:,:) = 0
    
        ! Generate random x vector and rhs vector in rank 0
    
        call random_number(x(:,:,:))
    
        ! y = A_x * x
        do k = 1, nz
            do j = 1, ny
                y(1, j, k) = b_x(1) * x(1, j, k) + c_x(1) * x(2, j, k)
                do i = 2, nx - 1
                    y(i, j, k) = a_x(i) * x(i - 1, j, k) + b_x(i) * x(i, j, k) + c_x(i) * x(i + 1, j, k)
                enddo
                y(nx, j, k) = a_x(nx) * x(nx - 1, j, k) + b_x(nx) * x(nx, j, k)
            enddo
        enddo
    
        ! z = A_y * y
        do k = 1, nz
            do i = 1, nx
                z(i, 1, k) = b_y(1) * y(i, 1, k) + c_y(1) * y(i, 2, k)
            enddo
            do j = 2, ny - 1
                do i = 1, nx
                    z(i, j, k) = a_y(j) * y(i, j - 1, k) + b_y(j) * y(i, j, k) + c_y(j) * y(i, j + 1, k)
                enddo
            enddo
            do i = 1, nx
                z(i, ny, k) = a_y(ny) * y(i, ny - 1, k) + b_y(ny) * y(i, ny, k)
            enddo
        enddo

        ! d = A_z * z
        do j = 1, ny
            do i = 1, nx
                d(i, j, 1) = b_z(1) * z(i, j, 1) + c_z(1) * z(i, j, 2)
            enddo
        enddo
        do k = 2, nz - 1
            do j = 1, ny
                do i = 1, nx
                    d(i, j, k) = a_z(k) * z(i, j, k - 1) + b_z(k) * z(i, j, k) + c_z(k) * z(i, j, k + 1)
                enddo
            enddo
        enddo
        do j = 1, ny
            do i = 1, nx
                d(i, j, nz) = a_z(nz) * z(i, j, nz - 1) + b_z(nz) * z(i, j, nz)
            enddo
        enddo

        deallocate (a_x, b_x, c_x)
        deallocate (a_y, b_y, c_y)
        deallocate (a_z, b_z, c_z)
        deallocate (y, z)
    
    end subroutine build_global_coeff_array

    subroutine dealloc_all

        if (is_root) then
            deallocate (d, x)
        endif
        deallocate (ax_sub, bx_sub, cx_sub)
        deallocate (ay_sub, by_sub, cy_sub)
        deallocate (az, bz, cz)
        deallocate (d_sub)
        deallocate (cnt_x, disp_x)
        deallocate (cnt_y, disp_y)
        deallocate (cnt_all, disp_all)

    end subroutine dealloc_all    

    subroutine distribute_rhs_array

        integer :: i, j, k, rank
        integer :: pos_iproc, pos_jproc, pos_k, pos_j
        integer :: disp_iproc, disp_jproc, cnt_iproc, cnt_jproc
        integer :: iproc, jproc

        double precision, allocatable, dimension(:)   :: d_blk

        if (is_root) then

            allocate ( d_blk(nx * ny * nz) ); d_blk(:) = 0

            do rank = 0, nprocs - 1
                iproc = rank / npy + 1
                jproc = mod(rank, npy) + 1
                disp_iproc = disp_x(iproc)
                disp_jproc = disp_y(jproc)
                cnt_iproc = cnt_x(iproc)
                cnt_jproc = cnt_y(jproc)
                pos_iproc = disp_iproc * ny * nz
                pos_jproc = disp_jproc * cnt_iproc * nz
                    
                do k = 1, nz
                    pos_k = (k - 1) * cnt_iproc * cnt_jproc
                    do j = 1, cnt_jproc
                        pos_j = (j - 1) * cnt_iproc
                        do i = 1, cnt_iproc
                            d_blk(i + pos_j + pos_k + pos_jproc + pos_iproc ) &
                                = d(i + disp_iproc, j + disp_jproc, k)
                        enddo
                    enddo
                enddo
            enddo
        else
            allocate ( d_blk(0) )
        endif

        allocate ( d_sub(nx_sub, ny_sub, nz) ); d_sub(:,:,:) = 0
        ! Scatter rhs vector
        call MPI_Scatterv(d_blk, cnt_all, disp_all, MPI_DOUBLE_PRECISION, d_sub, n_sub, &
                          MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    
        deallocate( d_blk )

    end subroutine distribute_rhs_array

    subroutine collect_solution_array

        integer :: i, j, k, rank
        integer :: pos_iproc, pos_jproc, pos_k, pos_j
        integer :: disp_iproc, disp_jproc, cnt_iproc, cnt_jproc
        integer :: iproc, jproc

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
            do rank = 0, nprocs - 1
                iproc = rank / npy + 1
                jproc = mod(rank, npy) + 1
                disp_iproc = disp_x(iproc)
                disp_jproc = disp_y(jproc)
                cnt_iproc = cnt_x(iproc)
                cnt_jproc = cnt_y(jproc)
                pos_iproc = disp_iproc * ny * nz
                pos_jproc = disp_jproc * cnt_iproc * nz

                do k = 1, nz
                    pos_k = (k - 1) * cnt_iproc * cnt_jproc

                    do j = 1, cnt_jproc
                        pos_j = (j - 1) * cnt_iproc

                        do i = 1, cnt_iproc
                            d(i + disp_iproc, j + disp_jproc, k) &
                                = d_blk( i + pos_j + pos_k + pos_jproc + pos_iproc )
                        enddo
                    enddo
                enddo
            enddo
        endif

        deallocate( d_blk )

    end subroutine collect_solution_array

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
