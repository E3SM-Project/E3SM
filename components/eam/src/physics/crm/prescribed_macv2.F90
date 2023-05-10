!-------------------------------------------------------------------
! manages reading and interpolation of prescribed MACV2 anthropogenic
! aerosols' radiative forcing for HighResMIP
! based on sp_driver_v1.f90 and modified to be ported into CAM physics,
! following cloud_rad_props.F90 and aer_rad_props.F90 as CAM code templates
! technical details are summarized in MACv2SP.pdf
! https://drive.google.com/open?id=1uYxPp5ytWELvLQ4Rm0ryc3iYjsOJFciO
! Koichi Sakaguchi
!-------------------------------------------------------------------
module prescribed_macv2

!USE mo_simple_plumes, ONLY: sp_aop_profile
use shr_kind_mod,     only : r8 => shr_kind_r8, r4 => shr_kind_r4
use cam_abortutils,   only : endrun
use cam_logfile,      only: iulog
use spmd_utils,       only: masterproc, mpi_logical, mpicom, mpir8, iam, masterprocid
use cam_history,    only: addfld, horiz_only, outfld
use shr_const_mod,  only: SHR_CONST_PI
use ppgrid,         only: pcols,pver,begchunk,endchunk,psubcols 
use radconstants,   only: nswbands
use time_manager,  only: get_nstep  !for debug
 
IMPLICIT NONE
private
!save !probably needed to retain the variables throughout the integration?

!public interface
public :: macv2_rad_props_init
public :: prescribed_macv2_readnl
public :: sp_aop_profile
public :: sp_aop_dNovrN

logical, public :: do_macv2sp = .false.  ! whether to do MACv2-SP, default is false
                                        ! if do_macv2sp is specified in the atm_in namelist,
                                        ! this value is overwritten and macv2sp is used
                                        
!saved variables
integer, PARAMETER ::                      &
   nplumes   = 9                          ,& !< Number of plumes
   nfeatures = 2                          ,& !< Number of features per plume
   ntimes    = 52                         ,& !< Number of times resolved per year (52 => weekly resolution)
   nyears    = 251                           !< Number of years of available forcing


real(r8)  ::  &
       plume_lat      (nplumes)               ,& !< latitude of plume center (AOD maximum)
       plume_lon      (nplumes)               ,& !< longitude of plume center (AOD maximum)
       beta_a         (nplumes)               ,& !< parameter a for beta function vertical profile
       beta_b         (nplumes)               ,& !< parameter b for beta function vertical profile
       aod_spmx       (nplumes)               ,& !< anthropogenic AOD maximum at 550 for plumes 
       aod_fmbg       (nplumes)               ,& !< anthropogenic AOD at 550 for fine-mode natural background (idealized to mimic Twomey effect)
       asy550         (nplumes)               ,& !< asymmetry parameter at 550nm for plume
       ssa550         (nplumes)               ,& !< single scattering albedo at 550nm for plume
       angstrom       (nplumes)               ,& !< Angstrom parameter for plume 
       sig_lon_E      (nfeatures,nplumes)     ,& !< Eastward extent of plume feature
       sig_lon_W      (nfeatures,nplumes)     ,& !< Westward extent of plume feature
       sig_lat_E      (nfeatures,nplumes)     ,& !< Southward extent of plume feature
       sig_lat_W      (nfeatures,nplumes)     ,& !< Northward extent of plume feature
       theta          (nfeatures,nplumes)     ,& !< Rotation angle of plume feature
       ftr_weight     (nfeatures,nplumes)     ,& !< Feature weights 
       time_weight    (nfeatures,nplumes)     ,& !< Time weights 
       time_weight_bg (nfeatures,nplumes)     ,& !< as time_weight but for natural background in Twomey effect 
       year_weight    (nyears,nplumes)        ,& !< Yearly weight for plume
       ann_cycle      (nfeatures,ntimes,nplumes) !< annual cycle for plume feature

!set by namelist
character(len=256) :: filename = ' '
character(len=256) :: datapath = ' '

integer :: plume_number, plume_feature, year_fr, years 

logical, PARAMETER :: localdebug = .TRUE. ! was false previously !BEH
!MACv2-specific log file (fort.128); workaround for cesm2beta05 to let all processes to write into a log 
integer, PARAMETER :: MACv2_errunit = 128  

!With the unit number of "iulog", the master pcoc writes into atm.log, and other tasks write into cesm.log. 
!If we use * for the unit number, all processes including master proc will write into cesm.log.
!this is different from our cesm2beta code.
!With iulog, the master proc writes to atm.log, but all other proces do not write into any log files
!With *, the master proc writes to log.0000.out but all other proces do not write into any log files

!use for output variable name
character(len=5), public :: swbandnum(nswbands) =(/'_sw01','_sw02','_sw03','_sw04','_sw05','_sw06','_sw07','_sw08','_sw09','_sw10','_sw11','_sw12','_sw13','_sw14'/)

!==============================================================================
contains
!==============================================================================

subroutine macv2_rad_props_init()
    use netcdf
    use ioFileMod,      only: getfil
    use error_messages, only: handle_ncerr

#if ( defined SPMD )
    use mpishorthand
#endif

    character(len=256) :: macv2file 
    character(len=256) :: locfn
    integer :: ncid, dimid, ierr
    integer :: VarID
    integer :: isw


    if(do_macv2sp) then
        macv2file = trim(datapath) // "/" // trim(filename) !from the namelist

        ! masterproc opens the netcdf file and obtrains variable dimensions
        if(masterproc) then
            ! Determine whether file is on local disk.
            call getfil( trim(macv2file), locfn, 0)
            call handle_ncerr( nf90_open(locfn, NF90_NOWRITE, ncid), 'MACv2 file missing')
            write(iulog,*)' macv2_rad_props_init: reading MACv2 antrhopogenic aerosol forcings from file ',locfn

            !check variable dimensions
            ! plume number
            call handle_ncerr(nf90_inq_dimid( ncid, 'plume_number', dimid), 'getting plume_number dim')
            call handle_ncerr(nf90_inquire_dimension( ncid, dimid, len=plume_number), 'getting n plume_number')
            if (plume_number /= nplumes) call endrun('number of plumes does not match')

            !plume feature
            call handle_ncerr(nf90_inq_dimid( ncid, 'plume_feature', dimid), 'getting plume_feature dim')
            call handle_ncerr(nf90_inquire_dimension( ncid, dimid, len=plume_feature), 'getting n plume_feature')
            if (plume_feature /= nfeatures) call endrun('plume_feature does not match')

            !number of years available 
            call handle_ncerr(nf90_inq_dimid( ncid, 'years', dimid), 'getting years dim')
            call handle_ncerr(nf90_inquire_dimension( ncid, dimid, len=years), 'getting n plume_feature')
            if (years /= nyears) call endrun('nyears does not match')

            !number of times resolved per year (52 => weekly resolution)
            call handle_ncerr(nf90_inq_dimid( ncid, 'year_fr', dimid), 'getting years dim')
            call handle_ncerr(nf90_inquire_dimension( ncid, dimid, len=year_fr), 'getting n plume_feature')
            if (year_fr /= ntimes) call endrun('ntimes does not match')

            if (localdebug) then
                write(iulog,*) 'macv2_rad_props_init (KSA): pcols,pver,begchunk,endchunk,psubcols :', pcols,pver,begchunk,endchunk,psubcols 
            end if
    
    
        endif ! if (masterproc)

        !master proc get variables from the input file
        if(masterproc) then
           call handle_ncerr( nf90_inq_varid(ncid, 'plume_lat', VarID),&
              'macv2 plume_lat get')
           !call handle_ncerr( nf90_get_var(ncid, VarID, plume_lat),&
           !   'macv2 read plume_lat')
           call handle_ncerr( nf90_get_var(ncid, VarID, plume_lat(:), &
              start=(/1/),count=(/nplumes/)),'macv2 read plume_lat')

           call handle_ncerr( nf90_inq_varid(ncid, 'plume_lon', VarID),&
              'macv2 plume_lon get')
           call handle_ncerr( nf90_get_var(ncid, VarID, plume_lon(:), start=(/1/),&
              count=(/nplumes/)),'macv2 read plume_lon')

           call handle_ncerr( nf90_inq_varid(ncid, 'beta_a', VarID),&
              'macv2 beta_a get')
           call handle_ncerr( nf90_get_var(ncid, VarID, beta_a(:), start=(/1/),&
              count=(/nplumes/)),'macv2 read beta_a')

           call handle_ncerr( nf90_inq_varid(ncid, 'beta_b', VarID),&
              'macv2 beta_b get')
           call handle_ncerr( nf90_get_var(ncid, VarID, beta_b(:), start=(/1/),&
              count=(/nplumes/)),'macv2 read beta_b')

           call handle_ncerr( nf90_inq_varid(ncid, 'aod_spmx', VarID),&
              'macv2 aod_spmx get')
           call handle_ncerr( nf90_get_var(ncid, VarID, aod_spmx(:), start=(/1/),&
              count=(/nplumes/)),'macv2 read aod_spmx')

           call handle_ncerr( nf90_inq_varid(ncid, 'aod_fmbg', VarID),&
              'macv2 aod_fmbg get')
           call handle_ncerr( nf90_get_var(ncid, VarID, aod_fmbg(:), start=(/1/),&
              count=(/nplumes/)),'macv2 read aod_fmbg')

           call handle_ncerr( nf90_inq_varid(ncid, 'ssa550', VarID),&
              'macv2 ssa550 get')
           call handle_ncerr( nf90_get_var(ncid, VarID, ssa550(:), start=(/1/),&
              count=(/nplumes/)),'macv2 read ssa550')

           call handle_ncerr( nf90_inq_varid(ncid, 'asy550', VarID),&
              'macv2 asy550 get')
           call handle_ncerr( nf90_get_var(ncid, VarID, asy550(:), start=(/1/),&
              count=(/nplumes/)),'macv2 read asy550')

           call handle_ncerr( nf90_inq_varid(ncid, 'angstrom', VarID),&
              'macv2 angstrom get')
           call handle_ncerr( nf90_get_var(ncid, VarID, angstrom(:), start=(/1/),&
              count=(/nplumes/)),'macv2 read angstrom')

           call handle_ncerr( nf90_inq_varid(ncid, 'sig_lat_W', VarID),&
              'macv2 sig_lat_W get')
           call handle_ncerr( nf90_get_var(ncid, VarID, sig_lat_W(:,:), start=(/1,1/),&
              count=(/nfeatures,nplumes/)),'macv2 read sig_lat_W')

           call handle_ncerr( nf90_inq_varid(ncid, 'sig_lat_E', VarID),&
              'macv2 sig_lat_E get')
           call handle_ncerr( nf90_get_var(ncid, VarID, sig_lat_E(:,:), start=(/1,1/),&
              count=(/nfeatures,nplumes/)),'macv2 read sig_lat_E')

           call handle_ncerr( nf90_inq_varid(ncid, 'sig_lon_E', VarID),&
              'macv2 sig_lon_E get')
           call handle_ncerr( nf90_get_var(ncid, VarID, sig_lon_E(:,:), start=(/1,1/),&
              count=(/nfeatures,nplumes/)),'macv2 read sig_lon_E')

           call handle_ncerr( nf90_inq_varid(ncid, 'sig_lon_W', VarID),&
              'macv2 sig_lon_W get')
           call handle_ncerr( nf90_get_var(ncid, VarID, sig_lon_W(:,:), start=(/1,1/),&
              count=(/nfeatures,nplumes/)),'macv2 read sig_lon_W')

           call handle_ncerr( nf90_inq_varid(ncid, 'theta', VarID),&
              'macv2 theta get')
           call handle_ncerr( nf90_get_var(ncid, VarID, theta(:,:), start=(/1,1/),&
              count=(/nfeatures,nplumes/)),'macv2 read theta')

           call handle_ncerr( nf90_inq_varid(ncid, 'ftr_weight', VarID),&
              'macv2 ftr_weight get')
           call handle_ncerr( nf90_get_var(ncid, VarID, ftr_weight(:,:), start=(/1,1/),&
              count=(/nfeatures,nplumes/)),'macv2 read ftr_weight')

           call handle_ncerr( nf90_inq_varid(ncid, 'year_weight', VarID),&
              'macv2 year_weight get')
           call handle_ncerr( nf90_get_var(ncid, VarID, year_weight(:,:), start=(/1,1/),&
              count=(/nyears,nplumes/)),'macv2 read year_weight')

           call handle_ncerr( nf90_inq_varid(ncid, 'ann_cycle', VarID),&
              'macv2 ann_cycle get')
           call handle_ncerr( nf90_get_var(ncid, VarID, ann_cycle(:,:,:), start=(/1,1,1/), &
              count=(/nfeatures,ntimes,nplumes/)),'macv2 read ann_cycle')

        endif ! if masterproc


        !master proc broadcast the MACv2 variables to all the tasks
#if ( defined SPMD )
            !1D arrays
            call mpibcast(plume_lat, nplumes, mpir8, 0, mpicom, ierr)
            call mpibcast(plume_lon, nplumes, mpir8, 0, mpicom, ierr)
            call mpibcast(beta_a, nplumes, mpir8, 0, mpicom, ierr)
            call mpibcast(beta_b, nplumes, mpir8, 0, mpicom, ierr)
            call mpibcast(aod_spmx, nplumes, mpir8, 0, mpicom, ierr)
            call mpibcast(aod_fmbg, nplumes, mpir8, 0, mpicom, ierr)
            call mpibcast(ssa550, nplumes, mpir8, 0, mpicom, ierr)
            call mpibcast(asy550, nplumes, mpir8, 0, mpicom, ierr)
            call mpibcast(angstrom, nplumes, mpir8, 0, mpicom, ierr)
  
            call mpibcast(sig_lat_W, nfeatures*nplumes, mpir8, 0, mpicom, ierr)
            call mpibcast(sig_lat_E, nfeatures*nplumes, mpir8, 0, mpicom, ierr)
            call mpibcast(sig_lon_E, nfeatures*nplumes, mpir8, 0, mpicom, ierr)
            call mpibcast(sig_lon_W, nfeatures*nplumes, mpir8, 0, mpicom, ierr)
            call mpibcast(theta, nfeatures*nplumes, mpir8, 0, mpicom, ierr)
            call mpibcast(ftr_weight, nfeatures*nplumes, mpir8, 0, mpicom, ierr)
            call mpibcast(year_weight, nyears*nplumes, mpir8, 0, mpicom, ierr)
            call mpibcast(ann_cycle, nfeatures*ntimes*nplumes, mpir8, 0, mpicom, ierr)

#endif

        !define output variables for each wavelength
        !for now, do not indicate sampling sequence of field (see cam_history.F90) as done
        !for the other variables from the radiative transfer scheme

        do isw = 1, nswbands

             call addfld('MACv2_aod_loc_2d'//swbandnum(isw),  horiz_only,   'A', '-', & 
                  'MACv2 anthropogenic aerosol optical depth with its vertical coordinate')
             call addfld('MACv2_aod_2d'//swbandnum(isw),  horiz_only,   'A', '-', &
                  'MACv2 anthropogenic aerosol optical depth after interpolation')
             call addfld('MACv2_ssa', horiz_only, 'A', '-', 'MACv2 SSA')

             call addfld('MACv2_aod'//swbandnum(isw), (/ 'lev' /),  'A', '-',  & 
                  'anthropogenic aerosol optical depth, MACv2')
             call addfld('MACv2_ssa'//swbandnum(isw), (/ 'lev' /),  'A', '-',  &
                  'anthropogenic aerosol single scattering albedo, MACv2')
             call addfld('MACv2_asy'//swbandnum(isw), (/ 'lev' /),  'A', '-',  &
                  'anthropogenic aerosol asymmetry parameter, MACv2')

        end do

        
    end if !(if do_macv2)
    
    
    return


end subroutine macv2_rad_props_init

!-------------------------------------------------------------------
!-------------------------------------------------------------------
subroutine prescribed_macv2_readnl(nlfile)
! - need to edit components/cam/bld/namelist_files/namelist_definition.xml
! to define new namelist variables
    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
#if ( defined SPMD )
    use mpishorthand
#endif

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'prescribed_macv2_readnl' 

    character(len=256) :: prescribed_macv2_file
    character(len=256) :: prescribed_macv2_datapath


    namelist /prescribed_macv2_nl/ &     
      do_macv2sp,      &
      prescribed_macv2_file,      &
      prescribed_macv2_datapath


    ! Read namelist
    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'prescribed_macv2_nl', status=ierr) 
       if (ierr == 0) then
          read(unitn, prescribed_macv2_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

!master proc broadcast the MACv2 namelist variables to all the tasks
#if ( defined SPMD )
    call mpibcast(do_macv2sp,   1,  mpi_logical, 0, mpicom, ierr)
#endif

    !Probably don't need to broadcast the filename variables - master process
    ! reads in the variables from the file and broadcast the variables

    ! Update module variables with user settings.
    filename   = prescribed_macv2_file
    datapath   = prescribed_macv2_datapath

end subroutine prescribed_macv2_readnl


! ------------------------------------------------------------------------------------------------------------------------
! SET_TIME_WEIGHT:  The simple plume model assumes that meteorology constrains plume shape and that only source strength
! influences the amplitude of a plume associated with a given source region.   This routine retrieves the temporal weights
! for the plumes.  Each plume feature has its own temporal weights which varies yearly.  The annual cycle is indexed by
! week in the year and superimposed on the yearly mean value of the weight. 
!
subroutine set_time_weight(year_fr)
!
    ! ---------- 
    !
    real(r8), intent(in) ::  &
         year_fr           !< Fractional Year (1850.0 - 2100.99)

    integer          ::  &
         iyear          ,& !< Integer year values between 1 and 156 (1850-2100) 
         iweek          ,& !< Integer index (between 1 and ntimes); for ntimes=52 this corresponds to weeks (roughly)
         iplume            ! plume number
    !
    ! ---------- 
    !
    iyear = FLOOR(year_fr) - 1849
    iweek = FLOOR((year_fr - FLOOR(year_fr)) * ntimes) + 1

    if ((iweek > ntimes) .OR. (iweek < 1) .OR. (iyear > nyears) .OR. (iyear < 1)) STOP 'Time out of bounds in set_time_weight'
    do iplume=1,nplumes
        time_weight(1,iplume) = year_weight(iyear,iplume) * ann_cycle(1,iweek,iplume)
        time_weight(2,iplume) = year_weight(iyear,iplume) * ann_cycle(2,iweek,iplume)
        time_weight_bg(1,iplume) = ann_cycle(1,iweek,iplume)
        time_weight_bg(2,iplume) = ann_cycle(2,iweek,iplume) 
    end do
    
    return

end subroutine set_time_weight


! ------------------------------------------------------------------------------------------------------------------------
! SP_AOP_PROFILE:  This subroutine calculates the simple plume aerosol and cloud active optical properties based on the
! the simple plume fit to the MPI Aerosol Climatology (Version 2).  It sums over nplumes to provide a profile of aerosol
! optical properties on a host models vertical grid. 

subroutine sp_aop_profile (ncol           ,lambda    ,    &
   oro            ,lon            ,lat            , &
   year_fr        ,z_in           ,aod_prof       ,ssa_prof       , &
   asy_prof       ,lchnk          ,isw)

    !the enhancement factor for cloud droplet number concentration, dNovrN, is now calculated
    !in the separate subroutine sp_aop_dNovrN
    !
    ! ---------- 
    ! used dimsizes of pcols for in/out and variables to be written into history files
    ! other variables, used ncol
    integer, intent(in)        :: ncol                       !< number of columns

    integer, intent(in)        :: lchnk  !to write out variables for debug
    integer, intent(in)        :: isw    !index for sw bands

    ! advice from Balwinder to use (:) for in/out variables, which are defined in the parent subroutine
    ! radiation_tend

    real(r8), intent(in)           :: &
         lambda   ,               & !< wavelength (nm)
         year_fr,                 & !< Fractional Year (1903.0 is the 0Z on the first of January 1903, Gregorian)
         oro(:),               & !< orographic height (m)
         lon(:),               & !< longitude (in radians)
         lat(:),               & !< latitude (in radians)
         z_in (:,:)           !< height above the surface (m)


    real(r8), intent(out)          :: &
         !dNovrN(:)           , & !< anthropogenic increase in cloud drop number concentration (factor)
         !dNovrN is now calculated in a separate subroutine sp_aop_dNovrN below
         aod_prof(:,:) , & !< profile of aerosol optical depth
         ssa_prof(:,:) , & !< profile of single scattering albedo
         asy_prof(:,:)     !< profile of asymmetry parameter

    integer                        :: iplume, icol, k, kin, ik !loop indices
    integer, parameter             :: nlevels = 80  !for local vertical coordinate
    real(r8), parameter            :: zmax = 15000.0_r8 !aod = 0 above this height
    real(r8), parameter            :: radtodeg = 180.0_r8/SHR_CONST_PI

    real(r8)                       ::  &
         eta(ncol,nlevels),        & !< normalized height (by 15 km)
         z_beta(ncol,nlevels),     & !< profile for scaling column optical depth
         prof(ncol,nlevels),       & !< scaled profile (by beta function)
         beta_sum(ncol),           & !< vertical sum of beta function
         ssa(ncol),                & !< single scattering albedo 
         asy(ncol),                & !< asymmetry parameter
         cw_an(ncol),              & !< column weight for simple plume (anthropogenic) AOD at 550 nm
         cw_bg(ncol),              & !< column weight for fine-mode natural background AOD at 550 nm
         caod_sp(ncol),            & !< column simple plume anthropogenic AOD at 550 nm
         caod_bg(ncol),            & !< column fine-mode natural background AOD at 550 nm
         a_plume1,                 & !< gaussian longitude factor for feature 1
         a_plume2,                 & !< gaussian longitude factor for feature 2
         b_plume1,                 & !< gaussian latitude factor for feature 1
         b_plume2,                 & !< gaussian latitude factor for feature 2
         delta_lat,                & !< latitude offset
         delta_lon,                & !< longitude offset
         delta_lon_t,              & !< threshold for maximum longitudinal plume extent used in transition from 360 to 0 degrees
         lon1,                     & !< rotated longitude for feature 1
         lat1,                     & !< rotated latitude for feature 2
         lon2,                     & !< rotated longitude for feature 1
         lat2,                     & !< rotated latitude for feature 2
         f1,                       & !< contribution from feature 1
         f2,                       & !< contribution from feature 2
         f3,                       & !< contribution from feature 1 in natural background of Twomey effect
         f4,                       & !< contribution from feature 2 in natural background of Twomey effect
         aod_550,                  & !< aerosol optical depth at 550nm
         aod_lmd,                  & !< aerosol optical depth at input wavelength
         lfactor                     !< factor to compute wavelength dependence of optical properties
         
         !local variables added by KSA to use a different vertical coordinate in MACv2SP
         !copied and modified from mo_simple_plumes_v1.f90
         !real(r8) :: z(ncol,nlevels)            !< heights by colums (ncolumns,level)
         !real(r8) :: dz(ncol,nlevels)           !< layer thicknesses (ncolumns,level)
         real(r8) :: z(nlevels)            !< mid-point heights by colums (level)
         real(r8) :: dz(nlevels)           !< layer thicknesses (level)
         real(r8) :: iz(nlevels+1)           !< interface height (level+1)
         real(r8) :: zdiff(nlevels)           !< for debug; save z-tz as an array (get error in writing/minloc on z-tz)

         real(r8) :: aod_prof_loc(ncol,nlevels)  !< profile of aerosol optical depth in the local vertical coordinate
         real(r8) :: ssa_prof_loc(ncol,nlevels)  !< profile of single scattering albedo in the local vertical coordinate
         real(r8) :: asy_prof_loc(ncol,nlevels)  !< profile of asymmetry parameter in the local vertical coordinate
         real(r8) :: aod_int_loc(pcols)           !< vertically integrated aerosol optical depth from the local variable
         real(r8) :: aod_int(pcols)               !< same, but from the intermedaite, interpolated cam-level aod
         real(r8) :: tz   !height ASL of CAM vertical level inside the kin & icol loops
         
!         real(r8) :: x3(3), y3(3)       !for 3-point interpolation
!         real(r8) :: x2(2), y2(2), iznext(2)       !for 2-point interpolation
!         integer ::  xind(2)  !for 2-point interpolation
         !decided to use optical values from the MACv2's vertical level that is the closest to 
         !the target CAM model level, instead of interpolation
         
         real(r8) :: dblk   !convert integer index k to double precision
         real(r8) :: interp_f  !factor to make the original and interpolated AOD to be the same
         integer :: nstep      !timestep count
         
         !the following three variables (negct, tzneg, and tznegind) are to detect and record 
         !the grid points with negative height ASL. But not fully utilized yet.
         integer :: negct     !negative grid cell count (only once per column)        
         real(r8) :: tzneg(ncol)  !< save tz if any grid column have negative height ASL (save the lowest level)
         integer  :: tznegind(ncol)  !< to save the indices for negative tz columns


    !
    ! ---------- 
    nstep = get_nstep()   !get model time step number
    
    if (localdebug) then
        write(MACv2_errunit,*) 'sp_aop_profile (KSA): sp_aop_profile started, iam: ', iam
        write(MACv2_errunit,*) 'sp_aop_profile (KSA): nstep: ', nstep
    end if


    !initialize local (pcol) variables -------------------------------------
    ! initialize output variables (following aer_rad_props_sw)
    ! initialize to conditions that would cause failure (for columns beyond ncol)
    aod_prof(:,:) = -100._r8
    ssa_prof(:,:) = -100._r8
    asy_prof(:,:) = -100._r8
    
    aod_int_loc(:) = -100._r8
    aod_int(:) = -100._r8
    
    ! also initialize rest of columns with physical values
    aod_prof(1:ncol,:) = 0._r8
    ssa_prof(1:ncol,:) = 0._r8
    asy_prof(1:ncol,:) = 0._r8

    aod_int_loc(1:ncol) = 0._r8
    aod_int(1:ncol) = 0._r8
      
!   The do loop is faster than (:,:,:) = 0.0, according to Balwinder. But for now
!   follow aer_rad_props_sw which initialize the variables as above
!    DO kin=1,pver
!      DO icol=1,ncol
!        aod_prof(icol,kin) = 0.0
!        ssa_prof(icol,kin) = 0.0
!        asy_prof(icol,kin) = 0.0
!      END DO
!    END DO
        
    !initialize local 1D height arrays, again lazy way
    iz(:) = 0._r8
    z(:) = 0._r8
    dz(:) = 0._r8   
    zdiff(:) = 0._r8 
    
    ! initialize local 2D variables
    DO k=1,nlevels
      DO icol=1,ncol
        aod_prof_loc(icol,k) = 0.0_r8
        ssa_prof_loc(icol,k) = 0.0_r8
        asy_prof_loc(icol,k) = 0.0_r8
        eta(icol,k) = 0.0_r8
        z_beta(icol,k) = 0.0_r8
        prof(icol,k) = 0.0_r8
        
      END DO
    END DO

    ! initialize 1D column arrays
    DO icol=1,ncol
      caod_sp(icol)  = 0.0_r8
      caod_bg(icol)  = 0.02_r8   !not sure why 0.02, but without this dNovrN does not work (KSA)
      beta_sum(icol) = 0.0_r8
      ssa(icol) = 0.0_r8
      asy(icol) = 0.0_r8
      cw_an(icol) = 0.0_r8
      cw_bg(icol) = 0.0_r8
      caod_sp(icol) = 0.0_r8
      caod_bg(icol) = 0.0_r8
      caod_bg(icol) = 0.0_r8
      tzneg(icol) = 0.0_r8
      tznegind(icol) = 0
      
    END DO
   
    ! initialize smaller arrays  
    !these are for linear interpolation between vertical coordinates
    !x2(:) = 0.0_r8
    !y2(:) = 0.0_r8
    !iznext(:) = 0.0_r8
    !xind(:) = 0
    
    ! initialize scalars
    a_plume1 = 0.0_r8
    a_plume2 = 0.0_r8
    b_plume1 = 0.0_r8
    b_plume2 = 0.0_r8
    delta_lat = 0.0_r8
    delta_lon = 0.0_r8
    delta_lon_t = 0.0_r8
    lon1 = 0.0_r8
    lat1 = 0.0_r8
    lon2 = 0.0_r8
    lat2 = 0.0_r8
    f1 = 0.0_r8
    f2 = 0.0_r8
    f3 = 0.0_r8
    f4 = 0.0_r8
    aod_550 = 0.0_r8
    aod_lmd = 0.0_r8
    lfactor = 0.0_r8
    tz = 0.0_r8
    dblk = 0.0_r8
    interp_f = 0.0_r8
    negct = 0
    
    ik = 1

    !------------------------------------------------------ initialization
    
    ! get time weights
    !
    call set_time_weight(year_fr)


!B    if ((nstep .eq. 0) .AND. (masterproc) .AND. (localdebug)) then
    if ((masterproc) .AND. (localdebug)) then
        write(MACv2_errunit,*) 'sp_aop_profile (KSA): year_fr:', year_fr
        do iplume=1,nplumes
            write(MACv2_errunit,*) 'sp_aop_profile (KSA): iplume:', iplume        
            write(MACv2_errunit,*) 'sp_aop_profile (KSA): time_weight (1,iplume):', time_weight(1,iplume)
            write(MACv2_errunit,*) 'sp_aop_profile (KSA): time_weight (2,iplume):', time_weight(2,iplume)
        end do
   
        write(MACv2_errunit,*) 'sp_aop_profile (KSA): time_weight (1,:):', time_weight(1,:)
    end if
        
    !define local vertical coordinate
    !original way, leaving dz(1) uninitialized.
    !z(1)  = 0.
    !DO k = 2,nlevels
    !  dz(k) = 50. * exp(0.03*(k-1))
    !  z(k) = z(k-1) + dz(k)
    !END DO
    
    !revised (KSA)
     DO k = 1,nlevels
      dz(k) = 50._r8 * exp(0.03_r8*(k-1))

      dblk = dble(k)
      dz(k) = 50._r8 * dexp(0.03_r8*(dblk-1.0_r8))

      iz(k+1) = iz(k) + dz(k)
      z(k) = iz(k) + 0.5_r8*dz(k)
    END DO
    
!B    if ((nstep .eq. 0) .AND. (masterproc) .AND. (localdebug)) then
    if ((masterproc) .AND. (localdebug)) then
         write(MACv2_errunit,*) 'sp_aop_profile (KSA): z:', z
         write(MACv2_errunit,*) 'sp_aop_profile (KSA): dz:', dz
         write(MACv2_errunit,*) 'sp_aop_profile (KSA): iz:', iz
         write(MACv2_errunit,*) 'sp_aop_profile (KSA): ik:', ik  
    end if
    

    !2D height-related variables
    !z_beta is the Heaviside (step) function H in equation (8)
    DO k=1,nlevels
      DO icol=1,ncol
        
        z_beta(icol,k)   = MERGE(1.0_r8, 0.0_r8, z(k) >= oro(icol))
        eta(icol,k)      = MAX(0.0_r8,MIN(1.0_r8,z(k)/15000._r8))  
        
      END DO
    END DO
    
    
    ! sum contribution from plumes to construct composite profiles of aerosol optical properties
    !
    DO iplume=1,nplumes
      !
      ! calculate vertical distribution function from parameters of beta distribution
      !
      DO icol=1,ncol
        beta_sum(icol) = 0.0_r8
      END DO

      ! prof is the karnel (integrand) of beta function
      ! beta_a = pi, beta_b = qi
      DO k=1,nlevels
        DO icol=1,ncol
          prof(icol,k)   = (eta(icol,k)**(beta_a(iplume)-1._r8) * (1._r8-eta(icol,k))**(beta_b(iplume)-1._r8)) * dz(k)
          
          beta_sum(icol) = beta_sum(icol) + prof(icol,k) !beta function B(p,q),integal
        END DO
      END DO

      ! equation (8) for bi
      ! normalization by beta_sum(icol) before zeroing out topography by z_beta??
      DO k=1,nlevels
        DO icol=1,ncol
          prof(icol,k)   = ( prof(icol,k) / beta_sum(icol) ) * z_beta(icol,k)
        END DO
      END DO
      !
      ! calculate plume weights
      ! Horizontal structure is shapred by the exponential function in eqn (4)
      ! The argument to the exponential function is <xi,R_ij*A_ij^-1*R_ij^-1*xi>,
      ! which is an inner product of xi and the other group of matrices.
      ! Therefore the expression in <> can be written as:
      ! xi^T*R_ij*A_ij^-1*R_ij^-1*xi  (^T denotes transpose)
      ! note that xi^T*R_ij and R_ij^-1*xi result in the same vectors,
      ! being a transpose of each other. This vector is [lon1 lat1] below.
      ! Thus expression inside <> reduces to : [lon1 lat1]*A_ij^-1*[lon1 lat1]^T
      ! that's how the code calculate the exponential funcion (Gaussian shape)
      ! see the slide#2 (KSA)


      DO icol=1,ncol
        !
        ! get plume-center relative spatial parameters for specifying amplitude of plume at given lat and lon
        !
        delta_lat   = lat(icol)*radtodeg - plume_lat(iplume)  !y-yi in the paper
        delta_lon   = lon(icol)*radtodeg - plume_lon(iplume)  !x-xi in the paper

        ! apply threshold for the maximum longitudinal extent used for plume-
        ! centered coordinate
        ! express dela_lon with values smaller than abs(180)
        ! so that we can use delta_lon to express east & west as (>0) & (<0)
        ! as written in eqn (6), (x-xi)>0 or (x-xi)<0
        ! fortan SIGN (a, b): Returns the absolute value of the first argument
        ! times the sign of the second argument.

        delta_lon_t = MERGE (260._r8, 180._r8, iplume == 1)
        delta_lon   = MERGE ( delta_lon - SIGN(360._r8,delta_lon) , delta_lon , ABS(delta_lon) > delta_lon_t)

        ! eqn(6), diagnal elements in the inverse of the covariance matrix, Aij^-1
        ! multipled by 0.5, which is the 1/2 in the argument for exp in eqn (4)
        ! added _r8 (KSA)
        a_plume1  = 0.5_r8 / (MERGE(sig_lon_E(1,iplume), sig_lon_W(1,iplume), delta_lon > 0._r8)**2)
        b_plume1  = 0.5_r8 / (MERGE(sig_lat_E(1,iplume), sig_lat_W(1,iplume), delta_lon > 0._r8)**2)
        a_plume2  = 0.5_r8 / (MERGE(sig_lon_E(2,iplume), sig_lon_W(2,iplume), delta_lon > 0._r8)**2)
        b_plume2  = 0.5_r8 / (MERGE(sig_lat_E(2,iplume), sig_lat_W(2,iplume), delta_lon > 0._r8)**2)
        !
        ! adjust for a plume specific rotation which helps match plume state to climatology.
        ! Eqn (7), inverse of rotation matrix R times offset vector xi
        ! (Rij^-1*xi) and also (xi^T*Rij) (later these variables are
        ! squared to represent dot products of these two vectors)

        ! for feature 1 of each plume
        lon1 =   COS(theta(1,iplume))*(delta_lon) + SIN(theta(1,iplume))*(delta_lat)
        lat1 = - SIN(theta(1,iplume))*(delta_lon) + COS(theta(1,iplume))*(delta_lat)

        ! feature 2 of each plume
        lon2 =   COS(theta(2,iplume))*(delta_lon) + SIN(theta(2,iplume))*(delta_lat)
        lat2 = - SIN(theta(2,iplume))*(delta_lon) + COS(theta(2,iplume))*(delta_lat)
        !
        ! calculate contribution to plume from its different features, to get a column weight for the anthropogenic
        ! (cw_an) and the fine-mode natural background aerosol (cw_bg)
        ! time_weight = wij, after eqn (5) on p4
        ! ftr_weight = fij , after eqn (5) on p4, these two are multiplied
        ! to become a_ij in eqn.4

        ! lon1 and lat1 are squared to represent (xi^T*Rij)*(Rij^-1*xi)
        ! feature 1, anthropogenic
        f1 = time_weight(1,iplume) * ftr_weight(1,iplume) * EXP(-1._r8* (a_plume1 * ((lon1)**2) + (b_plume1 * ((lat1)**2)))) 
        !feature 2, anthropogenic
        f2 = time_weight(2,iplume) * ftr_weight(2,iplume) * EXP(-1._r8* (a_plume2 * ((lon2)**2) + (b_plume2 * ((lat2)**2)))) 
        !feature 1, background
        f3 = time_weight_bg(1,iplume) * ftr_weight(1,iplume) * EXP(-1._r8* (a_plume1 * ((lon1)**2) + (b_plume1 * ((lat1)**2)))) 
        !feature 2, background
        f4 = time_weight_bg(2,iplume) * ftr_weight(2,iplume) * EXP(-1._r8* (a_plume2 * ((lon2)**2) + (b_plume2 * ((lat2)**2))))
        ! f1 + f2 (and f3+f4) equals the summation in eqn4

        ! now multiply AOD at the plume center (anthropogenic and background)
        ! by the horizontal structure obtained above
        cw_an(icol) = f1 * aod_spmx(iplume) + f2 * aod_spmx(iplume)  
        cw_bg(icol) = f3 * aod_fmbg(iplume) + f4 * aod_fmbg(iplume) 
        !
        ! calculate wavelength-dependent scattering properties
        lfactor   = MIN(1.0_r8,700.0_r8/lambda) !inverse of eqn.12

        !single scattering albedo
        !see my slide to verify that this is equivalent to eqn (11)
        ssa(icol) = (ssa550(iplume) * lfactor**4) / ((ssa550(iplume) * lfactor**4) + ((1.0_r8-ssa550(iplume)) * lfactor))
     
        !assymetry parameter
        asy(icol) =  asy550(iplume) * SQRT(lfactor) !eqn 13

      END DO !icol
      !
      ! distribute plume optical properties across its vertical profile weighting by optical depth and scaling for
      ! wavelength using the angstrom parameter. 
      !      
      lfactor = EXP(-angstrom(iplume) * LOG(lambda/550.0_r8)) !eqn 10  !  added _r8 after sigterm error 190831-164930
      !note that lfactor = 1 for lamda = 550 (nm)

      DO k=1,nlevels
        DO icol = 1,ncol
          !extinction coefficient (ei) at each level of each column
          aod_550          = prof(icol,k)     * cw_an(icol)

          !same but for a different wavelength
          aod_lmd          = aod_550          * lfactor

          ! integrate ei across vertical levels:
          ! add column simple plume anthropogenic AOD at 550 nm
          ! from this plume; used for cloud droplet scaling
          caod_sp(icol)    = caod_sp(icol)    + aod_550
          ! same but for the background aerosol optical depth (only used for
          ! cloud droplet scaling, therefore consider only 550nm)
          caod_bg(icol)    = caod_bg(icol)    + prof(icol,k) * cw_bg(icol)

          ! in each grid cell, add contribution from this plume to
          ! the variable saving the total aod, asy, and ssa.
          ! what are these expressions for asy and ssa?

          asy_prof_loc(icol,k) = asy_prof_loc(icol,k) + aod_lmd * ssa(icol) * asy(icol)
          ssa_prof_loc(icol,k) = ssa_prof_loc(icol,k) + aod_lmd * ssa(icol)
          aod_prof_loc(icol,k) = aod_prof_loc(icol,k) + aod_lmd
        END DO
      END DO
    END DO  !iplume
    !

    ! complete optical depth weighting
    ! if single scattering albedo at a given grid cell is too small,
    ! set assymetric parameter = 0
    ! if aerosol optical thickness at a given grid cell is too small,
    ! set single scattering albedo = 1
    ! ssa is multiplied by aod to copy the vertical/horizontal pattern
    ! divide ssa by the integrated aod to normalize the vertical/horizontal pattern
    ! for asy this normalization is done using ssa (before introducing ssa = 1.0 for small aod locations)
    ! for being "small", cime/share/csm_share/shr/shr_flux_mod.F90 has local parameters:
    ! tiny = 1.0e-12_R8 and tiny2 = 1.0e-6_R8
    ! follow tiny2 and use 1.0e-6_R8 instead of TINY(1.) in the original code
        
    DO k=1,nlevels
      DO icol = 1,ncol
        asy_prof_loc(icol,k) = MERGE(asy_prof_loc(icol,k)/ssa_prof_loc(icol,k), 0.0_r8, ssa_prof_loc(icol,k) > 1.0e-6_r8)
        ssa_prof_loc(icol,k) = MERGE(ssa_prof_loc(icol,k)/aod_prof_loc(icol,k), 1.0_r8, aod_prof_loc(icol,k) > 1.0e-6_r8)
      END DO
    END DO
    !

    ! calculate effective radius normalization (divisor) factor
    ! equqtion (15)
    ! Also calculate vertically integrated AOD for scaling vertically interpolated values
    ! following sp_driver_v1.f90, L162 (done as a simple summation)
    DO icol=1,ncol
      !dNovrN(icol) = LOG((1000.0 * (caod_sp(icol) + caod_bg(icol))) + 1.0)/LOG((1000.0 * caod_bg(icol)) + 1.0)
      aod_int_loc(icol) = SUM(aod_prof_loc(icol,:))
    END DO


    !vertical interpolation/translation from zcol = 80 levels to CAM's 32 levels
    !Z3, height above sea level as output, but the original variable 
    !zm is height above the surface
    !MACV2 coordinate is the height above sea level.
    !need to add surface height in doing interpolation in the model
    !code
    
    !The lowest CAM-MPAS level is always greater than ~55m
    
    !aerosol optical thickness (and asy and ssa, which are multiplied by
    !aod) are set to zero above 15,000m. 

    ! do not include write function or endrun subroutine inside the do loops b/c
    ! I/O processes for write is very slow

    DO kin=1,pver
      DO icol = 1,ncol
        
        !tz is the height above the surface for CAM's icol and kin level
        
        tz = z_in(icol,kin) + oro(icol) !need to be the height above sea level 

        if (tz >= zmax) then  !zmax = 15000, as in Stevens et al. 2017

            aod_prof(icol,kin) = 0.0_r8
            ssa_prof(icol,kin) = 1.0_r8
            asy_prof(icol,kin) = 0.0_r8

        else if (tz < 0.0_r8) then !if the model level is below sea level
            if(tznegind(icol) < 1) then
                negct = negct + 1
            end if
            tzneg(icol) = tz
            tznegind(icol) = 1
            
            !use lowest level from the MACv2 profile
            aod_prof(icol,kin) = aod_prof_loc(icol,1)
            ssa_prof(icol,kin) = ssa_prof_loc(icol,1)
            asy_prof(icol,kin) = asy_prof_loc(icol,1)
            
        else
            !find the nearest level in the MACV2 coordinate
                
            zdiff = z - tz    
             
            ik = minloc(abs(zdiff),DIM=1) 
             
            !instead of vertical interpolation, take the closest level's value
            aod_prof(icol,kin) = aod_prof_loc(icol,ik)
            ssa_prof(icol,kin) = ssa_prof_loc(icol,ik)
            asy_prof(icol,kin) = asy_prof_loc(icol,ik)

            !a bit of older code using linear interpolation, for record
            !AOD
            !y2 = aod_prof_loc(icol,xind)
            !aod_prof(icol,kin) = intrpf_2pt(tz,x2,y2)
            !SSA
            !y2 = ssa_prof_loc(icol,xind)
            !ssa_prof(icol,kin) = intrpf_2pt(tz,x2,y2)
            !ASY
            !y2 = asy_prof_loc(icol,xind)
            !asy_prof(icol,kin) = intrpf_2pt(tz,x2,y2)            
               
                
        end if !tz if conditioning

        !vertically integrate following sp_driver_v1.f90, L162

        aod_int(icol) = aod_int(icol) + aod_prof(icol,kin) 


      END DO
    END DO
    
 
    if((negct > 0) .AND. (localdebug)) then
        write(MACv2_errunit,*) 'sp_aop_profile (KSA): iam:', iam  
        write(MACv2_errunit,*) 'sp_aop_profile (KSA): tz is negative at some columns: negct =  ', negct
        write(MACv2_errunit,*) 'sp_aop_profile (KSA): tzneg(1:ncol): ', tzneg(1:ncol)

        !call endrun('sp_aop_profile (KSA): height above the surface is negative: check topography input and/or Z3')           
    end if


    !scale to get the same optical depth - probably not needed for linear interpolation
    !do not apply to ssa and asy
    DO icol=1, ncol
        !2D integrated value as simple sum, following sp_driver_v1.f90
        if (aod_int(icol) > 1.0e-6_r8 ) then
            interp_f = aod_int_loc(icol)/aod_int(icol)
            aod_prof(icol,:) = aod_prof(icol,:)*interp_f
        else
            aod_prof(icol,:) = 0.0_r8
        end if
    END DO

    !fix unphysical values from interpolation
    DO kin=1,pver
      DO icol = 1,ncol
          if(aod_prof(icol,kin) < 0.0_r8) then
             aod_prof(icol,kin) = 0.0_r8
          end if
          
          if(ssa_prof(icol,kin) < 0.0_r8) then
             ssa_prof(icol,kin) = 0.0_r8
          end if
          
          if(ssa_prof(icol,kin) > 1.0_r8) then
             ssa_prof(icol,kin) = 1.0_r8
          end if
          
          if(asy_prof(icol,kin) > 1.0_r8) then
             asy_prof(icol,kin) = 1.0_r8
          end if
          
          if(asy_prof(icol,kin) < -1.0_r8) then
             asy_prof(icol,kin) = -1.0_r8
          end if
          
      END DO
    END DO


    !update the integrated AOD after scaling and negative value correction
    DO icol=1,ncol
      aod_int(icol) = 0.0_r8 
      aod_int(icol) = SUM(aod_prof(icol,:))
    END DO

    !call outfld('MACv2_lat',  lat*radtodeg,  pcols, lchnk)
    !call outfld('MACv2_lon',  lon*radtodeg,  pcols, lchnk)
    call outfld('MACv2_aod_loc_2d'//swbandnum(isw),  aod_int_loc,  pcols, lchnk)
    call outfld('MACv2_aod_2d'//swbandnum(isw),  aod_int,  pcols, lchnk)
    call outfld('MACv2_ssa', ssa, pcols, lchnk)


    return
    
end subroutine sp_aop_profile
  
!================================================================================================


subroutine sp_aop_dNovrN (ncol    ,lambda    ,    &
   oro            ,lon            ,lat            , &
   year_fr        ,z_in           ,dNovrN)


    ! return dNovrN to cloud microphysics for indirect effect of aerosol
    ! ---------- 
  
    integer, intent(in)        :: ncol                       !< number of columns


    real(r8), intent(in)           :: &
         lambda   ,               & !< wavelength (nm)
         year_fr,                 & !< Fractional Year (1903.0 is the 0Z on the first of January 1903, Gregorian)
         oro(:),               & !< orographic height (m)
         lon(:),               & !< longitude (in radians)
         lat(:),               & !< latitude (in radians)
         z_in (:,:)           !< height above the surface (m)
        ! advice from Balwinder to use (:) for in/out variables

        ! oro(ncol),               & !< orographic height (m)
        ! lon(ncol),               & !< longitude (in radians)
        ! lat(ncol),               & !< latitude (in radians)
        ! z_in (ncol,pver)           !< height above the surface (m)

    real(r8), intent(out)          :: &
         !dNovrN(ncol)             !< anthropogenic increase in cloud drop number concentration (factor)
         dNovrN(:)             !< anthropogenic increase in cloud drop number concentration (factor)


    integer                        :: iplume, icol, k, kin, ik
    integer, parameter             :: nlevels = 80  !for local vertical coordinate
    real(r8), parameter            :: zmax = 15000.0 !aod = 0 above this height
    real(r8), parameter            :: radtodeg = 180.0_r8/SHR_CONST_PI

    real(r8)                       ::  &
         eta(ncol,nlevels),        & !< normalized height (by 15 km)
         z_beta(ncol,nlevels),     & !< profile for scaling column optical depth
         prof(ncol,nlevels),       & !< scaled profile (by beta function)
         beta_sum(ncol),           & !< vertical sum of beta function
!         ssa(ncol),                & !< single scattering albedo 
!         asy(ncol),                & !< asymmetry parameter
         cw_an(ncol),              & !< column weight for simple plume (anthropogenic) AOD at 550 nm
         cw_bg(ncol),              & !< column weight for fine-mode natural background AOD at 550 nm
         caod_sp(ncol),            & !< column simple plume anthropogenic AOD at 550 nm
         caod_bg(ncol),            & !< column fine-mode natural background AOD at 550 nm
         a_plume1,                 & !< gaussian longitude factor for feature 1
         a_plume2,                 & !< gaussian longitude factor for feature 2
         b_plume1,                 & !< gaussian latitude factor for feature 1
         b_plume2,                 & !< gaussian latitude factor for feature 2
         delta_lat,                & !< latitude offset
         delta_lon,                & !< longitude offset
         delta_lon_t,              & !< threshold for maximum longitudinal plume extent used in transition from 360 to 0 degrees
         lon1,                     & !< rotated longitude for feature 1
         lat1,                     & !< rotated latitude for feature 2
         lon2,                     & !< rotated longitude for feature 1
         lat2,                     & !< rotated latitude for feature 2
         f1,                       & !< contribution from feature 1
         f2,                       & !< contribution from feature 2
         f3,                       & !< contribution from feature 1 in natural background of Twomey effect
         f4,                       & !< contribution from feature 2 in natural background of Twomey effect
         aod_550,                  & !< aerosol optical depth at 550nm
         aod_lmd,                  & !< aerosol optical depth at input wavelength
         lfactor                     !< factor to compute wavelength dependence of optical properties
         
         !local vertical coordinate: copied from sp_driver_v1.f90
         !real(r8) :: z(ncol,nlevels)            !< heights by colums (ncolumns,level)
         !real(r8) :: dz(ncol,nlevels)           !< layer thicknesses (ncolumns,level)
         real(r8) :: z(nlevels)            !< mid-point heights by colums (level)
         real(r8) :: dz(nlevels)           !< layer thicknesses (level)
         real(r8) :: iz(nlevels+1)           !< interface height (level+1)
         real(r8) :: dblk   !convert integer index k to double precision (debug)
         integer :: nstep      !timestep count

    !
    ! ---------- 
    nstep = get_nstep()   !get model time step number

!B    if ((nstep .eq. 0) .AND. (localdebug)) then
    if (localdebug) then
        write(MACv2_errunit,*) 'sp_aop_dNovrN (KSA): sp_aop_dNovrN started, iam: ', iam
    end if

    !initialize local 1D height arrays, lazy way
    iz(:) = 0._r8
    z(:) = 0._r8
    dz(:) = 0._r8   
    
    ! initialize dNovrN to conditions that would cause failure (for columns beyond ncol but <= pcol)
    dNovrN(:) = -100._r8
    
    ! initialize local 2D variables
    DO k=1,nlevels
      DO icol=1,ncol
        eta(icol,k) = 0.0_r8
        z_beta(icol,k) = 0.0_r8
        prof(icol,k) = 0.0_r8
      END DO
    END DO
    
    ! initialize 1D column arrays
    DO icol=1,ncol
      dNovrN(icol)   = 1.0_r8
      caod_sp(icol)  = 0.0_r8
      caod_bg(icol)  = 0.0_r8      
      beta_sum(icol) = 0.0_r8
      cw_an(icol) = 0.0_r8
      cw_bg(icol) = 0.0_r8
      caod_sp(icol) = 0.0_r8
      caod_bg(icol) = 0.02_r8  !don't change! dNovrN won't work! (KSA)
      
    END DO
    
    
    ! initialize scalars
    a_plume1 = 0.0_r8
    a_plume2 = 0.0_r8
    b_plume1 = 0.0_r8
    b_plume2 = 0.0_r8
    delta_lat = 0.0_r8
    delta_lon = 0.0_r8
    delta_lon_t = 0.0_r8
    lon1 = 0.0_r8
    lat1 = 0.0_r8
    lon2 = 0.0_r8
    lat2 = 0.0_r8
    f1 = 0.0_r8
    f2 = 0.0_r8
    f3 = 0.0_r8
    f4 = 0.0_r8
    aod_550 = 0.0_r8
    aod_lmd = 0.0_r8
    lfactor = 0.0_r8
    dblk = 0.0_r8
        
    ik = 1
    
    ! get time weights
    !
    call set_time_weight(year_fr)


    !define the height array used locally in this subroutine
    !original way, leaving dz(1) uninitialized.
    !z(1)  = 0.
    !DO k = 2,nlevels
    !  dz(k) = 50. * exp(0.03*(k-1))
    !  z(k) = z(k-1) + dz(k)
    !END DO

    DO k = 1,nlevels
      dblk = dble(k)
      dz(k) = 50._r8 * dexp(0.03_r8*(dblk-1.0_r8))

      iz(k+1) = iz(k) + dz(k)
      z(k) = iz(k) + 0.5_r8*dz(k)
    END DO
    
    
    !z_beta is the Heaviside (step) function H in equation (8)
    DO k=1,nlevels
      DO icol=1,ncol
        z_beta(icol,k)   = MERGE(1.0_r8, 0.0_r8, z(k) >= oro(icol))
        eta(icol,k)      = MAX(0.0_r8,MIN(1.0_r8,z(k)/15000._r8))  
      END DO
    END DO

    ! sum contribution from plumes to construct composite profiles of aerosol optical properties
    !
    DO iplume=1,nplumes
      !
      ! calculate vertical distribution function from parameters of beta distribution
      !
      DO icol=1,ncol
        beta_sum(icol) = 0.0_r8
      END DO

      ! prof is the karnel (integrand) of beta function
      ! beta_a = pi, beta_b = qi
      DO k=1,nlevels
        DO icol=1,ncol
          prof(icol,k)   = (eta(icol,k)**(beta_a(iplume)-1._r8) * (1._r8-eta(icol,k))**(beta_b(iplume)-1._r8)) * dz(k)
          beta_sum(icol) = beta_sum(icol) + prof(icol,k) !beta function B(p,q),integal
        END DO
      END DO

      ! equation (8) for bi
      ! normalization by beta_sum(icol) before zeroing out topography by z_beta??
      DO k=1,nlevels
        DO icol=1,ncol
          prof(icol,k)   = ( prof(icol,k) / beta_sum(icol) ) * z_beta(icol,k)
        END DO
      END DO
      !
      ! calculate plume weights
      ! Horizontal structure is shapred by the exponential function in eqn (4)
      ! The argument to the exponential function is <xi,R_ij*A_ij^-1*R_ij^-1*xi>,
      ! which is an inner product of xi and the other group of matrices.
      ! Therefore the expression in <> can be written as:
      ! xi^T*R_ij*A_ij^-1*R_ij^-1*xi  (^T denotes transpose)
      ! note that xi^T*R_ij and R_ij^-1*xi result in the same vectors,
      ! being a transpose of each other. This vector is [lon1 lat1] below.
      ! Thus expression inside <> reduces to : [lon1 lat1]*A_ij^-1*[lon1 lat1]^T
      ! that's how the code calculate the exponential funcion (Gaussian shape)
      ! see the slide#2 (KSA)


      DO icol=1,ncol
        !
        ! get plume-center relative spatial parameters for specifying amplitude of plume at given lat and lon
        !
        delta_lat   = lat(icol)*radtodeg - plume_lat(iplume)  !y-yi in the paper
        delta_lon   = lon(icol)*radtodeg - plume_lon(iplume)  !x-xi in the paper

        !delta_lat   = lat(icol) - plume_lat(iplume)  !y-yi in the paper
        !delta_lon   = lon(icol) - plume_lon(iplume)  !x-xi in the paper
        ! apply threshold for the maximum longitudinal extent used for plume-
        ! centered coordinate
        ! express dela_lon with values smaller than abs(180)
        ! so that we can use delta_lon to express east & west as (>0) & (<0)
        ! as written in eqn (6), (x-xi)>0 or (x-xi)<0
        ! fortan SIGN (a, b): Returns the absolute value of the first argument
        ! times the sign of the second argument.

        delta_lon_t = MERGE (260._r8, 180._r8, iplume == 1)
        delta_lon   = MERGE ( delta_lon - SIGN(360._r8,delta_lon) , delta_lon , ABS(delta_lon) > delta_lon_t)

        ! eqn(6), diagnal elements in the inverse of the covariance matrix, Aij^-1
        ! multipled by 0.5, which is the 1/2 in the argument for exp in eqn (4)
        a_plume1  = 0.5_r8 / (MERGE(sig_lon_E(1,iplume), sig_lon_W(1,iplume), delta_lon > 0._r8)**2)
        b_plume1  = 0.5_r8 / (MERGE(sig_lat_E(1,iplume), sig_lat_W(1,iplume), delta_lon > 0._r8)**2)
        a_plume2  = 0.5_r8 / (MERGE(sig_lon_E(2,iplume), sig_lon_W(2,iplume), delta_lon > 0._r8)**2)
        b_plume2  = 0.5_r8 / (MERGE(sig_lat_E(2,iplume), sig_lat_W(2,iplume), delta_lon > 0._r8)**2)
        !
        ! adjust for a plume specific rotation which helps match plume state to climatology.
        ! Eqn (7), inverse of rotation matrix R times offset vector xi
        ! (Rij^-1*xi) and also (xi^T*Rij) (later these variables are
        ! squared to represent dot products of these two vectors)

        ! for feature 1 of each plume
        lon1 =   COS(theta(1,iplume))*(delta_lon) + SIN(theta(1,iplume))*(delta_lat)
        lat1 = - SIN(theta(1,iplume))*(delta_lon) + COS(theta(1,iplume))*(delta_lat)

        ! feature 2 of each plume
        lon2 =   COS(theta(2,iplume))*(delta_lon) + SIN(theta(2,iplume))*(delta_lat)
        lat2 = - SIN(theta(2,iplume))*(delta_lon) + COS(theta(2,iplume))*(delta_lat)
        !
        ! calculate contribution to plume from its different features, to get a column weight for the anthropogenic
        ! (cw_an) and the fine-mode natural background aerosol (cw_bg)
        ! time_weight = wij, after eqn (5) on p4
        ! ftr_weight = fij , after eqn (5) on p4, these two are multiplied
        ! to become a_ij in eqn.4

        ! lon1 and lat1 are squared to represent (xi^T*Rij)*(Rij^-1*xi)
        ! feature 1, anthropogenic
        f1 = time_weight(1,iplume) * ftr_weight(1,iplume) * EXP(-1._r8* (a_plume1 * ((lon1)**2) + (b_plume1 * ((lat1)**2)))) 
        !feature 2, anthropogenic
        f2 = time_weight(2,iplume) * ftr_weight(2,iplume) * EXP(-1._r8* (a_plume2 * ((lon2)**2) + (b_plume2 * ((lat2)**2)))) 
        !feature 1, background
        f3 = time_weight_bg(1,iplume) * ftr_weight(1,iplume) * EXP(-1._r8* (a_plume1 * ((lon1)**2) + (b_plume1 * ((lat1)**2)))) 
        !feature 2, background
        f4 = time_weight_bg(2,iplume) * ftr_weight(2,iplume) * EXP(-1._r8* (a_plume2 * ((lon2)**2) + (b_plume2 * ((lat2)**2))))
        ! f1 + f2 (and f3+f4) equals the summation in eqn4

        ! now multiply AOD at the plume center (anthropogenic and background)
        ! by the horizontal structure obtained above
        cw_an(icol) = f1 * aod_spmx(iplume) + f2 * aod_spmx(iplume)  
        cw_bg(icol) = f3 * aod_fmbg(iplume) + f4 * aod_fmbg(iplume) 
        !
        ! calculate wavelength-dependent scattering properties
        !
        lfactor   = MIN(1.0_r8,700.0_r8/lambda) !inverse of eqn.12

      END DO !icol
      !
      ! distribute plume optical properties across its vertical profile weighting by optical depth and scaling for
      ! wavelength using the angstrom parameter. 
      !      
      lfactor = EXP(-angstrom(iplume) * LOG(lambda/550.0_r8)) !eqn 10  !  added _r8 after sigterm error 190831-164930
      !note that lfactor = 1 for lamda = 550 (nm)

      DO k=1,nlevels
        DO icol = 1,ncol
          !extinction coefficient (ei) at each level of each column
          aod_550          = prof(icol,k)     * cw_an(icol)

          !same but for a different wavelength
          aod_lmd          = aod_550          * lfactor

          ! integrate ei across vertical levels:
          ! add column simple plume anthropogenic AOD at 550 nm
          ! from this plume used for cloud droplet scaling
          caod_sp(icol)    = caod_sp(icol)    + aod_550
          ! same but for the background aerosol optical depth (only used for
          ! cloud droplet scaling, therefore consider only 550nm)
          caod_bg(icol)    = caod_bg(icol)    + prof(icol,k) * cw_bg(icol)

        END DO
      END DO
    END DO  !iplume
    !

    ! calculate effective radius normalization (divisor) factor
    ! equqtion (15)
    
    if (masterproc .AND. localdebug) then
          write(MACv2_errunit,*) 'sp_aop_dNovrN (KSA): dNovrN(:) before', dNovrN(:)
    end if
    
    DO icol=1,ncol
      dNovrN(icol) = LOG((1000.0_r8 * (caod_sp(icol) + caod_bg(icol))) + 1.0_r8)/LOG((1000.0_r8 * caod_bg(icol)) + 1.0_r8)
    END DO
   
    if (masterproc .AND. localdebug) then
          write(MACv2_errunit,*) 'sp_aop_dNovrN (KSA): dNovrN(:) after', dNovrN(:)
          write(MACv2_errunit,*) 'sp_aop_dNovrN (KSA): finished'
    end if


    return
end subroutine sp_aop_dNovrN

!================================================================================================


real(r8) function intrpf(xi, x, y)
! Function to interpolate between data points (based on intrpf by Numerical Methods for Physics
! by Alejandro Garcia) using Lagrange polynomial (quadratic)
! Inputs
!   x    Vector of x coordinates of data points (3 values)
!   y    Vector of y coordinates of data points (3 values)
!   xi   The x value where interpolation is computed

   real(r8), intent(in) :: xi
   real(r8), intent(in) :: x(3), y(3)

   intrpf = (xi-x(2))*(xi-x(3))/((x(1)-x(2))*(x(1)-x(3)))*y(1) &
             + (xi-x(1))*(xi-x(3))/((x(2)-x(1))*(x(2)-x(3)))*y(2) &
             + (xi-x(1))*(xi-x(2))/((x(3)-x(1))*(x(3)-x(2)))*y(3)

end function intrpf

!================================================================================================

real(r8) function intrpf_2pt(xi, x, y)
! Function to interpolate between two data points 
! Inputs
!   x    Vector of x coordinates of data points (2 values)
!   y    Vector of y coordinates of data points (2 values)
!   xi   The x value where interpolation is computed

   real(r8), intent(in) :: xi
   real(r8), intent(in) :: x(2), y(2)

   intrpf_2pt = y(1) + (xi - x(1))*(y(2) - y(1))/(x(2) - x(1))

end function intrpf_2pt

!==============================================================================

end module prescribed_macv2

