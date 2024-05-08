module timer

    use mpi

    implicit none

    double precision, private       :: t_zero, t_curr
    double precision, public        :: t_array(64), t_array_reduce(64)
    integer(kind=4)     :: ntimer
    character(len=64)   :: t_str(64)

    contains

    subroutine  timer_init(n, str)

        integer(kind=4)     :: n
        character(len=64)   :: str(64)

        integer(kind=4)     :: i, ierr

        if(n.gt.64) then
            print *,'[Error] Maximun number of timer is 64'
            call MPI_Finalize(ierr)
            stop
        endif

        ntimer = n

        t_array(:) = 0.0d0
        t_array_reduce(:) = 0.0d0

        t_str(:)    = 'null'

        do i = 1, ntimer
            t_str(i) = str(i)
        enddo

    end subroutine timer_init

    subroutine timer_stamp0

        t_zero = MPI_Wtime()

    end subroutine timer_stamp0

    subroutine timer_stamp(timer_id)

        integer(kind=4), intent(in) :: timer_id

        t_curr = MPI_Wtime()
        t_array(timer_id) = t_array(timer_id) + t_curr - t_zero
        t_zero = t_curr

    end subroutine timer_stamp

    subroutine timer_start(timer_id)

        integer(kind=4), intent(in) :: timer_id

        t_array(timer_id) = MPI_Wtime()

    end subroutine timer_start

    subroutine timer_end(timer_id)

        integer(kind=4), intent(in) :: timer_id

        t_array(timer_id) = MPI_Wtime() - t_array(timer_id)

    end subroutine timer_end

    function timer_elapsed(timer_id) result(t_elapsed)

        integer(kind=4), intent(in) :: timer_id
        double precision    :: t_elapsed

        t_elapsed = MPI_Wtime() - t_array(timer_id)

        return

    end function timer_elapsed

    subroutine timer_reduction

        integer(kind=4) :: ierr

        t_array_reduce(:) = 0.0d0
        call MPI_Reduce(t_array, t_array_reduce, ntimer, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    end subroutine timer_reduction

    subroutine timer_output(myrank, nprocs)

        integer(kind=4)     :: i, myrank, nprocs
        character(len=64)   :: filename

        if(myrank.eq.0) then
            do i = 1, ntimer
                if(trim(t_str(i)).ne.'null') then
                    print '(a,a35,a,i3,a, f16.9)','[Timer] ', adjustl(t_str(i)),' : (',i,') : ',t_array_reduce(i) / nprocs
                endif
            enddo
        endif

        write(filename,'(a,i0.5)') 'timer_info.',myrank

        open(11, file=filename,action='write',form='formatted')
        do i = 1, ntimer
            if(trim(t_str(i)).ne.'null') then
                write(11, '(a,a35,a,i3,a, f16.9)') '[Timer] ', adjustl(t_str(i)),' : (',i,') : ',t_array(i)
            endif
        enddo
        close(11)

    end subroutine timer_output
    
end module timer
