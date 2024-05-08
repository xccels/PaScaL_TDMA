!======================================================================================================================
!> @file        solve_theta_plan_many_thread_team.f90
!> @brief       This file contains a solver subroutine for the example problem of PaScaL_TDMA.
!> @details     The target example problem is the three-dimensional time-dependent heat conduction problem 
!>              in a unit cube domain applied with the boundary conditions of vertically constant temperature 
!>              and horizontally periodic boundaries.
!> @author      
!>              - Kiha Kim (k-kiha@yonsei.ac.kr), Department of Computational Science & Engineering, Yonsei University
!>              - Ji-Hoon Kang (jhkang@kisti.re.kr), Korea Institute of Science and Technology Information
!>              - Jung-Il Choi (jic@yonsei.ac.kr), Department of Computational Science & Engineering, Yonsei University
!>
!> @date        June 2019
!> @version     1.0
!> @par         Copyright
!>              Copyright (c) 2019 Kiha Kim and Jung-Il choi, Yonsei University and 
!>              Ji-Hoon Kang, Korea Institute of Science and Technology Information, All rights reserved.
!> @par         License     
!>              This project is release under the terms of the MIT License (see LICENSE in )
!======================================================================================================================

!>
!> @brief       An example solver for many tridiagonal systems of equations using PaScaL_TDMA.
!> @details     This subroutine is for many tridiagonal systems of equations.
!>              It solves the the three-dimensional time-dependent heat conduction problem using PaScaL_TDMA.
!>              PaScaL_TDMA plans are created for many tridiagonal systems of equations and
!>              the many tridiagonal systems are solved plane-by-plane.
!> @param       theta       Main 3-D variable to be solved
!>
subroutine solve_theta_plan_many_thread_team(theta)

    use omp_lib
    use mpi
    use global
    use mpi_subdomain
    use mpi_topology
    use PaScaL_TDMA
    use timer

    implicit none
    double precision, dimension(0:nx_sub, 0:ny_sub, 0:nz_sub), intent(inout) :: theta
    
    integer :: myrank, nprocs, ierr
    integer :: time_step        ! Current time step
    double precision :: t_curr  ! Current simulation time

    ! Loop and index variables
    integer :: i,j,k, ii,jj
    integer :: ip, jp, kp
    integer :: im, jm, km
    integer :: jem, jep

    integer :: n_thds, tid

    ! Temporary variables for coefficient computations
    double precision :: dedx1, dedx2, dedy3, dedy4, dedz5, dedz6    ! Derivative terms
    double precision :: viscous_e1, viscous_e2, viscous_e3, viscous ! Viscous terms
    double precision :: ebc_down, ebc_up, ebc                       ! Boundary terms
    double precision :: eAPI, eAMI, eACI                            ! Diffusion treatment terms in x-direction
    double precision :: eAPJ, eAMJ, eACJ                            ! Diffusion treatment terms in y-direction
    double precision :: eAPK, eAMK, eACK                            ! Diffusion treatment terms in z-direction
    double precision :: eRHS                                        ! From eAPI to eACK

    double precision, allocatable, dimension(:, :, :)   :: rhs            ! r.h.s. array
    double precision, allocatable, dimension(:, :)      :: ap, am, ac, ad    ! Coefficient of ridiagonal matrix
    double precision, allocatable, dimension(:, :, :)   :: apt, amt, act, adt    ! Coefficient of ridiagonal matrix
    double precision :: t0, t1, tall

    character(len=64)   :: timer_item(64)

    type(ptdma_plan_many_thread_team) :: px_many, pz_many , py_many  ! Plan for many tridiagonal systems of equations

    call MPI_Comm_rank( MPI_COMM_WORLD, myrank, ierr)
    call MPI_Comm_size( MPI_COMM_WORLD, nprocs, ierr)

    n_thds = omp_get_max_threads()

    timer_item(1)    = 'Total '
    timer_item(2)    = 'RHS0'
    timer_item(3)    = 'z-solver'
    timer_item(4)    = 'y-solver'
    timer_item(5)    = 'x-solver'
    timer_item(6)    = 'Ghostcell update'
    timer_item(7)    = 'Others'

    call timer_init(7,timer_item)

    ! Check thread number
    if( mod((nx_sub-1), n_thds)) then
        print '(a,i5,a,i5)', '[Error] nx_sub-1 should be a multiple of threads. threads = ',n_thds,', nx_sub-1 = ',nx_sub-1
        call MPI_Finalize(ierr)
        stop
    endif

    if( mod((ny_sub-1), n_thds)) then
        print '(a,i5,a,i5)', '[Error] ny_sub-1 should be a multiple of threads. threads = ',n_thds,', ny_sub-1 = ',ny_sub-1
        call MPI_Finalize(ierr)
        stop
    endif

    if( mod((nz_sub-1), n_thds)) then
        print '(a,i5,a,i5)', '[Error] nz_sub-1 should be a multiple of threads. threads = ',n_thds,', nz_sub-1 = ',nz_sub-1
        call MPI_Finalize(ierr)
        stop
    endif

    ! Calculating r.h.s.
    allocate( rhs(0:nx_sub, 0:ny_sub, 0:nz_sub))

    ! Simulation begins
    t_curr = tStart
    dt = dtstart

    call timer_start(1)
    do time_step = 1, Tmax

        t_curr = t_curr + dt
        if(myrank==0) write(*,*) '[Main] Current time step = ', time_step
    
        call timer_stamp0()
!$omp parallel do default(shared) &
!$omp& private(i, ip, im, j, jp, jm, jep, jem, k, kp, km) &
!$omp& private(dedx1, dedx2, dedy3, dedy4, dedz5, dedz6) &
!$omp& private(viscous_e1, viscous_e2, viscous_e3, viscous) &
!$omp& private(ebc_down, ebc_up, ebc) &
!$omp& private(eAPI, eAMI, eACI, eAPJ, eAMJ, eACJ, eAPK, eAMK, eACK, eRHS)
        do k = 1, nz_sub-1
            kp = k+1
            km = k-1

            do j = 1, ny_sub-1
                jp = j + 1
                jm = j - 1
                jep = jpbc_index(j)
                jem = jmbc_index(j)

                do i = 1, nx_sub-1
                    ip = i+1
                    im = i-1

                    ! Diffusion term
                    dedx1 = (  theta(i ,j ,k ) - theta(im,j ,k )  )/dmx_sub(i )
                    dedx2 = (  theta(ip,j ,k ) - theta(i ,j ,k )  )/dmx_sub(ip)  
                    dedy3 = (  theta(i ,j ,k ) - theta(i ,jm,k )  ) /dmy_sub(j )
                    dedy4 = (  theta(i ,jp,k ) - theta(i ,j ,k )  ) /dmy_sub(jp)
                    dedz5 = (  theta(i ,j ,k ) - theta(i ,j ,km)  )/dmz_sub(k )
                    dedz6 = (  theta(i ,j ,kp) - theta(i ,j ,k )  )/dmz_sub(kp)

                    viscous_e1 = 1.d0/dx*(dedx2 - dedx1)
                    viscous_e2 = 1.d0/dy*(dedy4 - dedy3)
                    viscous_e3 = 1.d0/dz*(dedz6 - dedz5)
                    viscous = 0.5d0*Ct*(viscous_e1 + viscous_e2 + viscous_e3) 
                    
                    ! Boundary treatment for the y-direction only
                    ebc_down = 0.5d0*Ct/dy/dmy_sub(j)*thetaBC3_sub(i,k)
                    ebc_up = 0.5d0*Ct/dy/dmy_sub(jp)*thetaBC4_sub(i,k)
                    ebc = dble(1. - jem)*ebc_down + dble(1. - jep)*ebc_up

                    ! Diffusion term from incremental notation in next time step: x-direction
                    eAPI = -0.5d0*Ct/dx/dmx_sub(ip)
                    eAMI = -0.5d0*Ct/dx/dmx_sub(i )
                    eACI =  0.5d0*Ct/dx*( 1.d0/dmx_sub(ip) + 1.d0/dmx_sub(i) )

                    ! Diffusion term from incremental notation in next time step: z-direction
                    eAPK = -0.5d0*Ct/dz/dmz_sub(kp)
                    eAMK = -0.5d0*Ct/dz/dmz_sub(k )
                    eACK =  0.5d0*Ct/dz*( 1.d0/dmz_sub(kp) + 1.d0/dmz_sub(k) )

                    ! Diffusion term from incremental notation in next time step: y-direction
                    eAPJ = -0.5d0*Ct/dy*( 1.d0/dmy_sub(jp) )*dble(jep)
                    eAMJ = -0.5d0*Ct/dy*( 1.d0/dmy_sub(j ) )*dble(jem)
                    eACJ =  0.5d0*Ct/dy*( 1.d0/dmy_sub(jp) + 1.d0/dmy_sub(j) )
                    
                    eRHS = eAPK*theta(i,j,kp) + eACK*theta(i,j,k) + eAMK*theta(i,j,km)      &
                        & + eAPJ*theta(i,jp,k) + eACJ*theta(i,j,k) + eAMJ*theta(i,jm,k)      &
                        & + eAPI*theta(ip,j,k) + eACI*theta(i,j,k) + eAMI*theta(im,j,k)

                    ! r.h.s. term 
                    rhs(i,j,k) = theta(i,j,k)/dt + viscous + ebc      &
                            & - (theta(i,j,k)/dt + eRHS)
                enddo
            enddo
            ! print '(a7,x,4(i5,x))', '[RHS-Team]', k, omp_get_thread_num(), omp_get_num_threads(), omp_get_max_threads()
        enddo
!$omp end parallel do
        call timer_stamp(2)
        t1 = MPI_Wtime()

        call timer_stamp0()
        ! solve in the z-direction.
        ! Create a PaScaL_TDMA plan for the tridiagonal systems.
        call PaScaL_TDMA_plan_many_create_thread_team(pz_many, nx_sub-1, comm_1d_z%myrank, comm_1d_z%nprocs, comm_1d_z%mpi_comm)
        call timer_stamp(7)

        call timer_stamp0()
        ! Build a coefficient matrix for the tridiagonal systems into a 2D array.
!$omp parallel private(ap, am, ac, ad) default(shared)
        allocate( ap(1:nx_sub-1, 1:nz_sub-1), am(1:nx_sub-1, 1:nz_sub-1), ac(1:nx_sub-1, 1:nz_sub-1), ad(1:nx_sub-1, 1:nz_sub-1) )

!$omp do private(i, j, k, kp)
        do j = 1, ny_sub-1
            do k = 1, nz_sub-1
                kp = k+1
                do i = 1, nx_sub-1
                    ap(i,k) = -0.5d0*Ct/dz/dmz_sub(kp)*dt
                    am(i,k) = -0.5d0*Ct/dz/dmz_sub(k )*dt
                    ac(i,k) =  0.5d0*Ct/dz*( 1.d0/dmz_sub(kp) + 1.d0/dmz_sub(k) )*dt + 1.d0
                    
                    ad(i,k) = rhs(i,j,k)*dt
                enddo
            enddo

            ! Solve the tridiagonal systems under the defined plan with periodic boundary conditions.
            call PaScaL_TDMA_many_solve_cycle_thread_team(pz_many, am, ac, ap, ad,nx_sub-1,nz_sub-1)

            ! Return the solution to the r.h.s. plane-by-plane
            do k = 1, nz_sub-1
                rhs(1:nx_sub-1,j,k)=ad(1:nx_sub-1,k)
            enddo
        enddo
!$omp end do

        deallocate( ap, am, ac, ad )
!$omp end parallel
        call timer_stamp(3)

        call timer_stamp0()
        ! Destroy the PaScaL_TDMA plan for the tridiagonal systems.
        call PaScaL_TDMA_plan_many_destroy_thread_team(pz_many,comm_1d_z%nprocs)

        ! solve in the x-direction.
        ! Create a PaScaL_TDMA plan for the tridiagonal systems.
        call PaScaL_TDMA_plan_many_create_thread_team(py_many, nx_sub-1, comm_1d_y%myrank, comm_1d_y%nprocs, comm_1d_y%mpi_comm)
        call timer_stamp(7)

        call timer_stamp0()
        ! Build a coefficient matrix for the tridiagonal systems into a 2D array.
!$omp parallel private(ap, am, ac, ad) default(shared)
        allocate( ap(1:nx_sub-1, 1:ny_sub-1), am(1:nx_sub-1, 1:ny_sub-1), ac(1:nx_sub-1, 1:ny_sub-1), ad(1:nx_sub-1, 1:ny_sub-1) )

!$omp do private(i, j, k, jp, jm, jep, jem)
        do k = 1, nz_sub-1
            do j = 1, ny_sub-1
                jp = j + 1
                jm = j - 1
                jep = jpbc_index(j)
                jem = jmbc_index(j)
                
                do i = 1, nx_sub-1
                    ap(i,j) = -0.5d0*Ct/dy/dmy_sub(jp)*dble(jep)*dt
                    am(i,j) = -0.5d0*Ct/dy/dmy_sub(j )*dble(jem)*dt
                    ac(i,j) =  0.5d0*Ct/dy*( 1.d0/dmy_sub(jp) + 1.d0/dmy_sub(j) )*dt + 1.d0
                    ad(i,j) = rhs(i,j,k)
                end do
            end do

            ! Solve the tridiagonal systems under the defined plan.
            call PaScaL_TDMA_many_solve_thread_team(py_many, am, ac, ap, ad, nx_sub-1, ny_sub-1)

            ! Return the solution to the r.h.s. plane-by-plane.
            do j = 1, ny_sub-1
                rhs(1:nx_sub-1,j,k)=ad(1:nx_sub-1,j)
            enddo
        end do
!$omp end do

        deallocate( ap, am, ac, ad )
!$omp end parallel
        call timer_stamp(4)

        call timer_stamp0()
        ! Destroy the PaScaL_TDMA plan for the tridiagonal systems.
        call PaScaL_TDMA_plan_many_destroy_thread_team(py_many,comm_1d_y%nprocs)

        ! solve in the x-direction.
        ! Create a PaScaL_TDMA plan for the tridiagonal systems.
        call PaScaL_TDMA_plan_many_create_thread_team(px_many, ny_sub-1, comm_1d_x%myrank, comm_1d_x%nprocs, comm_1d_x%mpi_comm)
        call timer_stamp(7)

        call timer_stamp0()
        ! Build a coefficient matrix for the tridiagonal systems into a 2D array.
!$omp parallel private(ap, am, ac, ad) default(shared)
        allocate( ap(1:ny_sub-1, 1:nx_sub-1), am(1:ny_sub-1, 1:nx_sub-1), ac(1:ny_sub-1, 1:nx_sub-1), ad(1:ny_sub-1, 1:nx_sub-1) )

!$omp do private(i, j, k, ip, im)
        do k = 1, nz_sub-1
            do i = 1, nx_sub-1
                do j = 1, ny_sub-1
                    ip = i+1
                    im = i-1

                    ap(j,i) = -0.5d0*Ct/dx/dmx_sub(ip)*dt
                    am(j,i) = -0.5d0*Ct/dx/dmx_sub(i )*dt
                    ac(j,i) =  0.5d0*Ct/dx*( 1.d0/dmx_sub(ip) + 1.d0/dmx_sub(i) )*dt + 1.d0
                    ad(j,i) = rhs(i,j,k)
                enddo
            enddo

            ! Solve the tridiagonal systems under the defined plan with periodic boundary conditions.
            call PaScaL_TDMA_many_solve_cycle_thread_team(px_many, am, ac, ap, ad, ny_sub-1, nx_sub-1)

            ! Return the solution to theta plane-by-plane.
            do i = 1, nx_sub-1
                do j = 1, ny_sub-1
                    theta(i,j,k) = theta(i,j,k) + ad(j,i)
                enddo
            enddo
        enddo
!$omp end do

        deallocate( ap, am, ac, ad )
!$omp end parallel
        call timer_stamp(5)

        call timer_stamp0()
        call PaScaL_TDMA_plan_many_destroy_thread_team(px_many,comm_1d_x%nprocs)
        call timer_stamp(7)

        call timer_stamp0()
        ! Update ghostcells from the solution.
        call mpi_subdomain_ghostcell_update(theta, comm_1d_x, comm_1d_y, comm_1d_z)
        call timer_stamp(6)
    enddo

    call timer_end(1)
    call timer_reduction()
    call timer_output(myrank, nprocs)

    deallocate(rhs)

end subroutine solve_theta_plan_many_thread_team

