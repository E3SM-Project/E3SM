#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module common_movie_mod
  use control_mod, only : columnpackage, test_case, max_string_len
  use common_io_mod, only : output_start_time, output_end_time, &
       max_output_streams, output_frequency, nf_double, nf_int, &
       max_output_variables
  implicit none
  private
  public ::  varrequired, vartype, varnames, varcnt, vardims, &
	dimnames, maxdims, setvarnames, nextoutputstep



#ifdef _PRIM
  integer, parameter :: varcnt =  39

  integer, parameter :: maxdims =  7

  character*(*), parameter :: varnames(varcnt)=(/'ps         ', &
                                                 'geos       ', &
                                                 'area       ', &
                                                 'cv_lat     ', &
                                                 'cv_lon     ', &
                                                 'corners    ', &
                                                 'phys_area  ', &
                                                 'phys_lat   ', &
                                                 'phys_lon   ', &
                                                 'phys_cv_lat', &
                                                 'phys_cv_lon', &
                                                 'faceno     ', &
                                                 'zeta       ', &
                                                 'div        ', &
                                                 'C1         ', &
                                                 'C2         ', &
                                                 'C3         ', &
                                                 'C4         ', &
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


  integer, parameter :: vardims(maxdims,varcnt) =  reshape( (/ 1,5,0,0,0,0,0, & ! ps
                                                               1,5,0,0,0,0,0, & ! geos
                                                               1,0,0,0,0,0,0, & ! area
                                                               1,2,0,0,0,0,0, & ! cv_lat
                                                               1,2,0,0,0,0,0, & ! cv_lon
                                                               7,2,0,0,0,0,0, & ! subelement_corners
                                                               6,0,0,0,0,0,0, & ! phys_area
                                                               6,0,0,0,0,0,0, & ! phys_lat
                                                               6,0,0,0,0,0,0, & ! phys_lon
                                                               6,2,0,0,0,0,0, & ! phys_cv_lat
                                                               6,2,0,0,0,0,0, & ! phys_cv_lon
                                                               1,2,5,0,0,0,0, & ! faceno
                                                               1,2,5,0,0,0,0, & ! zeta
                                                               1,2,5,0,0,0,0, & ! div
                                                               6,2,5,0,0,0,0, & ! C1
                                                               6,2,5,0,0,0,0, & ! C2
                                                               6,2,5,0,0,0,0, & ! C3
                                                               6,2,5,0,0,0,0, & ! C4
                                                               1,2,5,0,0,0,0, & ! T
                                                               1,2,5,0,0,0,0, & ! Th
                                                               1,2,5,0,0,0,0, & ! u
                                                               1,2,5,0,0,0,0, & ! v
                                                               1,2,5,0,0,0,0, & ! ke
                                                               1,5,0,0,0,0,0, & ! hypervis
                                                               1,2,5,0,0,0,0, & ! Q
                                                               1,2,5,0,0,0,0, & ! Q2
                                                               1,2,5,0,0,0,0, & ! Q3
                                                               1,2,5,0,0,0,0, & ! Q4
                                                               1,2,5,0,0,0,0, & ! geo
                                                               1,2,5,0,0,0,0, & !omega
                                                               1,0,0,0,0,0,0, &  ! lat
                                                               1,0,0,0,0,0,0, &  ! lon
                                                               2,0,0,0,0,0,0, &  ! lev
                                                               3,0,0,0,0,0,0, &  ! ilev
                                                               2,0,0,0,0,0,0, &  !hy arrays
                                                               2,0,0,0,0,0,0, &
                                                               3,0,0,0,0,0,0, &  
                                                               3,0,0,0,0,0,0, &  
                                                               5,0,0,0,0,0,0 /),&
                                                               shape=(/maxdims,varcnt/))

  integer, parameter :: vartype(varcnt)=(/nf_double, nf_double,nf_double,nf_double,nf_double,&
                                          nf_int, nf_double,nf_double,nf_double,&
                                          nf_double, nf_double, nf_double,nf_double,nf_double,&
                                          nf_double, nf_double,nf_double,nf_double,&
                                          nf_double, nf_double,nf_double,&
                                          nf_double, nf_double,nf_double,&
                                          nf_double, nf_double,nf_double,&
                                          nf_double, nf_double,nf_double,&
                                          nf_double, nf_double,nf_double,nf_double,&
                                          nf_double, nf_double, nf_double, nf_double,&
                                          nf_double/)
  logical, parameter :: varrequired(varcnt)=(/.false.,.false.,.false.,.false.,.false.,&
                                              .false.,&
                                              .false.,&
                                              .false.,&
                                              .false.,&
                                              .false.,&
                                              .false.,&
                                              .false.,&
                                              .false.,&
                                              .false.,&
                                              .false.,&
                                              .false.,&
                                              .false.,&
                                              .false.,&
                                              .false.,&
                                              .false.,&
                                              .false.,&
                                              .false.,&
                                              .false.,&
                                              .false.,&
                                              .false.,.false.,.false.,.false.,&
                                              .false.,.false.,.true. ,.true. ,&
                                              .true. ,.true. ,&   ! lev,ilev
                                              .true. ,.true. ,&   ! hy arrays
                                              .true. ,.true. ,&   ! hy arrays
                                              .true./)
  character*(*),parameter :: dimnames(maxdims)=(/'ncol        ',&
                                                 'lev         ',&
                                                 'ilev        ',&
                                                 'nelem       ',&
                                                 'time        ',&
                                                 'nphys       ',&
                                                 'nsubelements'/)  
#else
#ifdef _PRIMDG  
  integer, parameter :: varcnt =11
  integer, parameter :: maxdims=4
  character*(*),parameter::dimnames(maxdims)=(/'ncolp','nlev ','nelem','time '/)  
  integer, parameter :: vardims(maxdims,varcnt) =  reshape( (/ 1,4,0,0,  &
                                                               1,2,4,0,  &
                                                               1,2,4,0,  &
                                                               1,2,4,0,  &
                                                               1,0,0,0,  &
                                                               1,0,0,0,  &
                                                               4,0,0,0,	 &
							       1,2,4,0,	 &
							       1,2,4,0,	 &
							       1,2,4,0,	 &
							       1,2,4,0/),&
                                                               shape=(/maxdims,varcnt/))
  character*(*),parameter::varnames(varcnt)=(/'ps   ','geop ','u    ','v    ',&
                                              'lonp ','latp ','time ',	      &
					      'T    ','p    ','zeta ','dgs  '/)
  integer, parameter :: vartype(varcnt)=(/nf_double,nf_double,nf_double,nf_double,&
                                          nf_double,nf_double,nf_double,&
					  nf_double,nf_double,nf_double,nf_double/)
  logical, parameter :: varrequired(varcnt)=(/.false.,.false.,.false.,.false.,&
                                              .true.,.true.,.true.,&
					      .false.,.false.,.false.,.false./)
#else
  integer, parameter :: varcnt = 17
  integer, parameter :: maxdims=5
  character*(*),parameter::dimnames(maxdims)=(/'ncol ','nlev ','nelem','time ','nphys'/)  
  integer, parameter :: vardims(maxdims,varcnt) =  reshape( (/ 1,4,0,0,0,  &  !ps
                                                               1,2,4,0,0,  &  !geop
                                                               1,2,4,0,0,  &  !u
                                                               1,2,4,0,0,  &  !v
                                                               1,0,0,0,0,  &  !lon
                                                               1,0,0,0,0,  &  !lat
                                                               5,0,0,0,0,  &  !phys_lon
                                                               5,0,0,0,0,  &  !phys_lat
                                                               5,0,0,0,0,  &  !phys_area
                                                               4,0,0,0,0,  &  ! time
                                                               5,2,4,0,0,  &  ! c1
                                                               5,2,4,0,0,  &  ! c2
                                                               5,2,4,0,0,  &  ! c3
                                                               5,2,4,0,0,  &  ! c4
                                                               1,2,4,0,0,  &  ! zeta
                                                               1,2,4,0,0,  &  ! div
                                                               1,0,0,0,0/),&  ! area
                                                               shape=(/maxdims,varcnt/))
  character*(*),parameter::varnames(varcnt)=(/'ps       ','geop     ','u        ','v        ',&
                                              'lon      ','lat      ',&
                                              'phys_lon ','phys_lat ','phys_area',&
                                              'time     ','c1       ','c2       ','c3       ','c4       ',&
                                              'zeta     ','div      ','area     '/)
  integer, parameter :: vartype(varcnt)=(/nf_double,nf_double,nf_double,nf_double,nf_double,&
       nf_double,nf_double,nf_double,nf_double,nf_double,nf_double,nf_double,&
       nf_double, nf_double, nf_double, nf_double, nf_double/)
  logical, parameter :: varrequired(varcnt)=(/.false.,.false.,.false.,.false.,&
                                              .true.,.true.,.true.,&
                                              .true.,.true.,.true.,&
                                              .false.,.false.,.false.,&
                                              .false.,.false.,.false.,.true./)

#endif
#endif




  ! end of analysis_nl namelist variables


contains

!
! This gets the default var list for namelist_mod
!
  subroutine setvarnames(nlvarnames)
#ifdef _PRIM
    use aquaplanet_io_mod, only : aq_set_varnames
    ! ---------------------
    use physics_io_mod, only : physics_set_varnames
#endif
    character*(*), intent(out) :: nlvarnames(:)
    integer :: lvarcnt
    if (varcnt > max_output_variables) then
       print *,__FILE__,__LINE__,"varcnt > max_output_varnames"
       stop
    endif
    lvarcnt=varcnt
    nlvarnames(1:varcnt) = varnames
    !print *,__FILE__,__LINE__,varcnt, size(nlvarnames),varnames
#ifdef _PRIM 
    if(test_case.eq.'aquaplanet') then
       call aq_set_varnames(lvarcnt,nlvarnames(lvarcnt+1:))
    end if
    if(columnpackage.ne.'none') then
       call physics_set_varnames(lvarcnt,nlvarnames(lvarcnt+1:))
    end if
#endif
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
