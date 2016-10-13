module mpas_performance

    use iso_c_binding, only : c_int, c_long_long

    interface

        subroutine start_timer(id) bind(c)
            use iso_c_binding, only : c_int
            integer(c_int), value :: id
        end subroutine start_timer
  
        subroutine stop_timer(id, secs, u_secs) bind(c)
            use iso_c_binding, only : c_int
            integer(c_int), value :: id
            integer(c_int) :: secs
            integer(c_int) :: u_secs
        end subroutine stop_timer
  
    end interface
 
#ifdef _PAPI
#include "f90papi.h"

    integer(c_int) :: nevents
    integer(c_int), dimension(:), pointer :: events
    integer(c_long_long), dimension(:), pointer :: papi_vals
    integer(c_long_long) :: usecs_start, usecs_end
#endif

    contains


    subroutine mpas_performance_init()

        implicit none

        integer(c_int) :: stat


#ifdef _PAPI
        call PAPIF_library_init(stat)

        nevents = 3
        allocate(events(nevents))
        allocate(papi_vals(nevents))

#ifdef SINGLE_PRECISION
        events(1) = PAPI_SP_OPS
        events(2) = PAPI_VEC_SP
        events(3) = PAPI_TOT_CYC
#else
        events(1) = PAPI_DP_OPS
        events(2) = PAPI_VEC_DP
        events(3) = PAPI_TOT_CYC
#endif

        call PAPIF_start_counters(events, nevents, stat)
#endif

    end subroutine mpas_performance_init


    subroutine mpas_performance_start()

        implicit none

        integer(c_int) :: stat


#ifdef _PAPI
        call PAPIF_read_counters(papi_vals, nevents, stat)
        call PAPIF_get_real_usec(usecs_start, stat)
#endif

    end subroutine mpas_performance_start


    subroutine mpas_performance_stop(timer_id)

        implicit none

        character(len=*), intent(in) :: timer_id

        integer(c_int) :: stat


#ifdef _PAPI
        call PAPIF_get_real_usec(usecs_end, stat)
        call PAPIF_read_counters(papi_vals, nevents, stat)

        write(0,*) '======== '//timer_id//' ========'
        write(0,*) 'FP ops     : ', papi_vals(1)
        write(0,*) 'Vec FP ops : ', papi_vals(2)
        write(0,*) 'cycles     : ', papi_vals(3)
        write(0,*) 'time, us   : ', usecs_end - usecs_start
        write(0,*) ' '
#endif

    end subroutine mpas_performance_stop


    subroutine mpas_performance_finalize()

        implicit none

        integer(c_int) :: stat


#ifdef _PAPI
        call PAPIF_stop_counters(events, nevents, stat)

        deallocate(events)
        deallocate(papi_vals)
#endif

    end subroutine mpas_performance_finalize

end module mpas_performance
