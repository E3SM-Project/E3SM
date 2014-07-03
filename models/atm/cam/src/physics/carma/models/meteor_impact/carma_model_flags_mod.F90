!! This module handles reading the namelist and provides access to some other flags
!! that control a specific CARMA model's behavior.
!!
!! By default the specific CARMA model does not have any unique namelist values. If
!! a CARMA model wishes to have its own namelist, then this file needs to be copied
!! from physics/cam to physics/model/<model_name> and the code needed to read in the
!! namelist values added there. This file will take the place of the one in
!! physics/cam. 
!!
!! It needs to be in its own file to resolve some circular dependencies.
!!
!! @author  Chuck Bardeen
!! @version Mar-2011
module carma_model_flags_mod

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use spmd_utils,     only: masterproc

  ! Flags for integration with CAM Microphysics
  public carma_model_readnl                   ! read the carma model namelist
  

  ! Namelist flags
  !
  ! Create a public definition of any new namelist variables that you wish to have,
  ! and default them to an inital value.
  real(r8), public               :: carma_emis_dust  = 0._r8    !! Total dust emission for the event (kg)
  real(r8), public               :: carma_emis_soot  = 0._r8    !! Total soot emission for the event (kg)
  integer, public                :: carma_emis_startdate = 1    !! start year and day of year (yyyyddd)
  integer, public                :: carma_emis_stopdate = 1     !! stop year and day of year (yyyyddd)
  integer, public                :: carma_emis_starttime = 0    !! start time of day (s)
  integer, public                :: carma_emis_stoptime = 0     !! stop time of day (s)
  real(r8), public               :: carma_emis_minlat = -90.    !! minimum latitude
  real(r8), public               :: carma_emis_maxlat = 90.     !! maximum latitude
  real(r8), public               :: carma_emis_minlon = 0.      !! minimum longitude
  real(r8), public               :: carma_emis_maxlon = 360.    !! maximum longitude
  logical, public                :: carma_fractal_soot = .false. !! fractal Soot

contains


  !! Read the CARMA model runtime options from the namelist
  !!
  !! @author  Chuck Bardeen
  !! @version Mar-2011
  subroutine carma_model_readnl(nlfile)
  
    ! Read carma namelist group.
  
    use abortutils,      only: endrun
    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand
  
    ! args
  
    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input
  
    ! local vars
  
    integer :: unitn, ierr
  
    ! read namelist for CARMA
    namelist /carma_model_nl/ &
      carma_emis_dust, carma_emis_soot, carma_emis_startdate, carma_emis_stopdate, &
      carma_emis_starttime, carma_emis_stoptime, carma_emis_minlat, carma_emis_maxlat, &
      carma_emis_minlon, carma_emis_maxlon, carma_fractal_soot
      
    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'carma_model_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, carma_model_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun('carma_model_readnl: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if
  
#ifdef SPMD
    call mpibcast(carma_emis_dust,       1,  mpir8,   0, mpicom)
    call mpibcast(carma_emis_soot,       1,  mpir8,   0, mpicom)
    call mpibcast(carma_emis_startdate,  1,  mpiint,  0, mpicom)
    call mpibcast(carma_emis_stopdate,   1,  mpiint,  0, mpicom)
    call mpibcast(carma_emis_starttime,  1,  mpiint,  0, mpicom)
    call mpibcast(carma_emis_stoptime,   1,  mpiint,  0, mpicom)
    call mpibcast(carma_emis_minlat,     1,  mpir8,   0, mpicom)
    call mpibcast(carma_emis_maxlat,     1,  mpir8,   0, mpicom)
    call mpibcast(carma_emis_minlon,     1,  mpir8,   0, mpicom)
    call mpibcast(carma_emis_maxlon,     1,  mpir8,   0, mpicom)
    call mpibcast(carma_fractal_soot,    1,  mpilog,  0, mpicom)
#endif
  
  end subroutine carma_model_readnl

end module carma_model_flags_mod
