#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module common_movie_mod
  use control_mod, only : test_case, max_string_len
  use common_io_mod, only : output_start_time, output_end_time, &
       max_output_streams, output_frequency, nf_double, nf_int, &
       max_output_variables
  implicit none
  private
  public ::  varrequired, vartype, varnames, varcnt, vardims, &
	dimnames, maxdims, setvarnames, nextoutputstep



#ifdef _PRIM
  integer, parameter :: varcnt =  30

  integer, parameter :: maxdims =  6

  character*(*), parameter :: varnames(varcnt)=(/'ps         ', &
                                                 'geos       ', &
                                                 'area       ', &
                                                 'cv_lat     ', &
                                                 'cv_lon     ', &
                                                 'corners    ', &
                                                 'faceno     ', &
                                                 'zeta       ', &
                                                 'div        ', &
                                                 'T          ', &
                                                 'Th         ', &
                                                 'u          ', &
                                                 'v          ', &
                                                 'ke         ', &
                                                 'hypervis   ', &
                                                 'Q          ', &
                                                 'Q2         ', &
                                                 'Q3         ', &
                                                 'Q4         ', &
                                                 'geo        ', &
                                                 'omega      ', &
                                                 'lat        ', &
                                                 'lon        ', &
                                                 'lev        ', &
                                                 'ilev       ', &
                                                 'hyam       ', &
                                                 'hybm       ', &
                                                 'hyai       ', &
                                                 'hybi       ', &
                                                 'time       '/)


  integer, parameter :: vardims(maxdims,varcnt) =  reshape( (/ 1,5,0,0,0,0, & ! ps
                                                               1,5,0,0,0,0, & ! geos
                                                               1,0,0,0,0,0, & ! area
                                                               1,2,0,0,0,0, & ! cv_lat
                                                               1,2,0,0,0,0, & ! cv_lon
                                                               6,2,0,0,0,0, & ! subelement_corners
                                                               1,2,5,0,0,0, & ! faceno
                                                               1,2,5,0,0,0, & ! zeta
                                                               1,2,5,0,0,0, & ! div
                                                               1,2,5,0,0,0, & ! T
                                                               1,2,5,0,0,0, & ! Th
                                                               1,2,5,0,0,0, & ! u
                                                               1,2,5,0,0,0, & ! v
                                                               1,2,5,0,0,0, & ! ke
                                                               1,5,0,0,0,0, & ! hypervis
                                                               1,2,5,0,0,0, & ! Q
                                                               1,2,5,0,0,0, & ! Q2
                                                               1,2,5,0,0,0, & ! Q3
                                                               1,2,5,0,0,0, & ! Q4
                                                               1,2,5,0,0,0, & ! geo
                                                               1,2,5,0,0,0, & !omega
                                                               1,0,0,0,0,0, &  ! lat
                                                               1,0,0,0,0,0, &  ! lon
                                                               2,0,0,0,0,0, &  ! lev
                                                               3,0,0,0,0,0, &  ! ilev
                                                               2,0,0,0,0,0, &  !hyam
                                                               2,0,0,0,0,0, &  !hybm
                                                               3,0,0,0,0,0, &  !hyai
                                                               3,0,0,0,0,0, &  !hybi
                                                               5,0,0,0,0,0 /),&  ! time
                                                               shape=(/maxdims,varcnt/))

  integer, parameter :: vartype(varcnt)=(/nf_double, nf_double,nf_double,nf_double,nf_double,&
                                          nf_int,    nf_double,nf_double,nf_double,nf_double,&
                                          nf_double, nf_double,nf_double,nf_double,nf_double,&
                                          nf_double, nf_double,nf_double,nf_double,nf_double,&
                                          nf_double, nf_double,nf_double,nf_double,nf_double,&
                                          nf_double, nf_double,nf_double,nf_double,nf_double/)
  logical, parameter :: varrequired(varcnt)=(/.false.,.false.,.false.,.false.,.false.,&
                                              .false.,.false.,.false.,.false.,.false.,&
                                              .false.,.false.,.false.,.false.,.false.,&
                                              .false.,.false.,.false.,.false.,.false.,&
                                              .false.,.true. ,.true. ,&
                                              .true. ,.true. ,&   ! lev,ilev
                                              .true. ,.true. ,&   ! hy arrays
                                              .true. ,.true. ,&   ! hy arrays
                                              .true./)            ! time
  character*(*),parameter :: dimnames(maxdims)=(/'ncol        ',&
                                                 'lev         ',&
                                                 'ilev        ',&
                                                 'nelem       ',&
                                                 'time        ',&
                                                 'nsubelements'/)  
#else
  integer, parameter :: varcnt = 10
  integer, parameter :: maxdims=4
  character*(*),parameter::dimnames(maxdims)=(/'ncol ','nlev ','nelem','time '/)  
  integer, parameter :: vardims(maxdims,varcnt) =  reshape( (/ 1,4,0,0,  &  !ps
                                                               1,2,4,0,  &  !geop
                                                               1,2,4,0,  &  !u
                                                               1,2,4,0,  &  !v
                                                               1,0,0,0,  &  !lon
                                                               1,0,0,0,  &  !lat
                                                               4,0,0,0,  &  ! time
                                                               1,2,4,0,  &  ! zeta
                                                               1,2,4,0,  &  ! div
                                                               1,0,0,0/),&  ! area
                                                               shape=(/maxdims,varcnt/))
  character*(*),parameter::varnames(varcnt)=(/'ps       ','geop     ','u        ','v        ',&
                                              'lon      ','lat      ','time     ',&
                                              'zeta     ','div      ','area     '/)
  integer, parameter :: vartype(varcnt)=(/nf_double,nf_double,nf_double,nf_double,nf_double,&
       nf_double,nf_double,nf_double,nf_double,nf_double/)
  logical, parameter :: varrequired(varcnt)=(/.false.,.false.,.false.,.false.,&
                                              .true.,.true.,.true.,&
                                              .false.,.false.,.true./)
#endif

  ! end of analysis_nl namelist variables

contains

!
! This gets the default var list for namelist_mod
!
  subroutine setvarnames(nlvarnames)
    character*(*), intent(out) :: nlvarnames(:)
    integer :: lvarcnt
    if (varcnt > max_output_variables) then
       print *,__FILE__,__LINE__,"varcnt > max_output_varnames"
       stop
    endif
    lvarcnt=varcnt
    nlvarnames(1:varcnt) = varnames
    !print *,__FILE__,__LINE__,varcnt, size(nlvarnames),varnames
  end subroutine setvarnames
!
! This function returns the next step number in which an output (either restart or movie) 
! needs to be written.
!
  integer function nextoutputstep(tl)
    use time_mod, only : Timelevel_t, nendstep  
    use control_mod, only : restartfreq
    type(timelevel_t), intent(in) :: tl
    integer :: ios, nstep(max_output_streams)

    nstep(:) = nEndStep
    do ios=1,max_output_streams
       if((output_frequency(ios) .gt. 0)) then
!          dont check start_time, since this will skip desired output when
!                  nstep <   start_time <=  nextoutputstep
!          if ((output_start_time(ios) .le. tl%nstep) .and. &
           if ((output_end_time(ios) .ge. tl%nstep)) then
             nstep(ios)=tl%nstep + output_frequency(ios) - &
                  MODULO(tl%nstep,output_frequency(ios))
          end if
       end if
    end do
    nextoutputstep=minval(nstep)
    if(restartfreq>0) then
       nextoutputstep=min(nextoutputstep,tl%nstep+restartfreq-MODULO(tl%nstep,restartfreq))    
    end if
 end function nextoutputstep
end module common_movie_mod
