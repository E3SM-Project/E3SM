module restFileMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Reads from or writes to/ the CLM restart file.
  !
  ! !USES:
  use shr_kind_mod         , only : r8 => shr_kind_r8
  use decompMod            , only : bounds_type
  use spmdMod              , only : masterproc, mpicom
  use abortutils           , only : endrun
  use shr_log_mod          , only : errMsg => shr_log_errMsg
  use clm_time_manager     , only : timemgr_restart_io, get_nstep
  use subgridRestMod       , only : SubgridRest
  use accumulMod           , only : accumulRest
  use histFileMod          , only : hist_restart_ncd
  use clm_varpar           , only : crop_prog
  use clm_varctl           , only : use_cn, use_c13, use_c14, use_lch4, use_cndv, use_ed, use_betr
  use clm_varctl           , only : create_glacier_mec_landunit, iulog 
  use clm_varcon           , only : c13ratio, c14ratio
  use clm_varcon           , only : nameg, namel, namec, namep, nameCohort
  use ch4Mod               , only : ch4_type
  use EDRestVectorMod      , only : EDRest
  use EDBioType            , only : EDbio_type
  use CNCarbonFluxType     , only : carbonflux_type
  use CNCarbonStateType    , only : carbonstate_type
  use CNDVType             , only : dgvs_type
  use CNStateType          , only : cnstate_type
  use CNNitrogenFluxType   , only : nitrogenflux_type
  use CNNitrogenStateType  , only : nitrogenstate_type

  use PhosphorusFluxType     , only : phosphorusflux_type
  use PhosphorusStateType    , only : phosphorusstate_type

  use AerosolType          , only : aerosol_type
  use CanopyStateType      , only : canopystate_type
  use EnergyFluxType       , only : energyflux_type
  use FrictionVelocityType , only : frictionvel_type
  use LakeStateType        , only : lakestate_type
  use PhotosynthesisType   , only : photosyns_type
  use SoilHydrologyType    , only : soilhydrology_type  
  use SoilStateType        , only : soilstate_type
  use SolarAbsorbedType    , only : solarabs_type
  use SurfaceAlbedoType    , only : surfalb_type
  use TemperatureType      , only : temperature_type
  use WaterfluxType        , only : waterflux_type
  use WaterstateType       , only : waterstate_type
  use atm2lndType          , only : atm2lnd_type
  use lnd2atmType          , only : lnd2atm_type
  use glc2lndMod           , only : glc2lnd_type
  use lnd2glcMod           , only : lnd2glc_type
  use BeTRTracerType       , only : BeTRTracer_Type    
  use TracerStateType      , only : TracerState_type
  use TracerFluxType       , only : TracerFlux_Type
  use tracercoefftype      , only : tracercoeff_type
  use ncdio_pio            , only : file_desc_t, ncd_pio_createfile, ncd_pio_openfile, ncd_global
  use ncdio_pio            , only : ncd_pio_closefile, ncd_defdim, ncd_putatt, ncd_enddef, check_dim
  use ncdio_pio            , only : check_att, ncd_getatt
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: restFile_read
  public :: restFile_write
  public :: restFile_open
  public :: restFile_close
  public :: restFile_getfile
  public :: restFile_filename        ! Sets restart filename
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: restFile_read_pfile     
  private :: restFile_write_pfile       ! Writes restart pointer file
  private :: restFile_closeRestart      ! Close restart file and write restart pointer file
  private :: restFile_dimset
  private :: restFile_add_ilun_metadata ! Add global metadata defining landunit types
  private :: restFile_add_icol_metadata ! Add global metadata defining column types
  private :: restFile_add_ipft_metadata ! Add global metadata defining pft types
  private :: restFile_dimcheck
  private :: restFile_enddef
  private :: restFile_check_consistency   ! Perform consistency checks on the restart file
  private :: restFile_read_consistency_nl ! Read namelist associated with consistency checks
  private :: restFile_check_fsurdat       ! Check consistency of fsurdat on the restart file
  private :: restFile_check_year          ! Check consistency of year on the restart file
  !
  ! !PRIVATE TYPES: None
  private
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine restFile_write( bounds, file,                                            &
       atm2lnd_vars, aerosol_vars, canopystate_vars, cnstate_vars,                    &
       carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars, carbonflux_vars, &
       ch4_vars, dgvs_vars, energyflux_vars, frictionvel_vars, lakestate_vars,        &
       nitrogenstate_vars, nitrogenflux_vars, photosyns_vars, soilhydrology_vars,     &
       soilstate_vars, solarabs_vars, surfalb_vars, temperature_vars,                 &
       waterflux_vars, waterstate_vars, EDbio_vars,                                   &
       phosphorusstate_vars, phosphorusflux_vars,                                     &
       betrtracer_vars, tracerstate_vars, tracerflux_vars, tracercoeff_vars,          &
       rdate, noptr)
    !
    ! !DESCRIPTION:
    ! Define/write CLM restart file.
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds          
    character(len=*)         , intent(in)    :: file             ! output netcdf restart file
    type(atm2lnd_type)       , intent(in)    :: atm2lnd_vars
    type(aerosol_type)       , intent(in)    :: aerosol_vars
    type(canopystate_type)   , intent(inout) :: canopystate_vars ! due to EDrest call
    type(cnstate_type)       , intent(inout) :: cnstate_vars
    type(carbonstate_type)   , intent(inout) :: carbonstate_vars
    type(carbonstate_type)   , intent(in)    :: c13_carbonstate_vars
    type(carbonstate_type)   , intent(in)    :: c14_carbonstate_vars
    type(carbonflux_type)    , intent(inout) :: carbonflux_vars
    type(ch4_type)           , intent(in)    :: ch4_vars
    type(dgvs_type)          , intent(in)    :: dgvs_vars
    type(energyflux_type)    , intent(in)    :: energyflux_vars
    type(frictionvel_type)   , intent(in)    :: frictionvel_vars
    type(lakestate_type)     , intent(in)    :: lakestate_vars
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    type(nitrogenflux_type)  , intent(in)    :: nitrogenflux_vars
    type(photosyns_type)     , intent(in)    :: photosyns_vars
    type(soilhydrology_type) , intent(in)    :: soilhydrology_vars
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    type(solarabs_type)      , intent(in)    :: solarabs_vars
    type(surfalb_type)       , intent(in)    :: surfalb_vars
    type(temperature_type)   , intent(in)    :: temperature_vars
    type(waterstate_type)    , intent(inout) :: waterstate_vars  ! due to EDrest call
    type(waterflux_type)     , intent(in)    :: waterflux_vars
    type(EDbio_type)         , intent(inout) :: EDbio_vars       ! due to EDrest call
    type(phosphorusstate_type),intent(inout)    :: phosphorusstate_vars
    type(phosphorusflux_type) ,intent(in)     :: phosphorusflux_vars
    type(tracerstate_type)   , intent(inout) :: tracerstate_vars ! due to Betrrest call
    type(BeTRTracer_Type)    , intent(in)    :: betrtracer_vars
    type(tracerflux_type)    , intent(inout) :: tracerflux_vars
    type(tracercoeff_type)   , intent(inout) :: tracercoeff_vars
    character(len=*)         , intent(in), optional :: rdate     ! restart file time stamp for name
    logical                  , intent(in), optional :: noptr     ! if should NOT write to the restart pointer file
    !
    ! !LOCAL VARIABLES:
    type(file_desc_t) :: ncid ! netcdf id
    integer :: i       ! index
    logical :: ptrfile ! write out the restart pointer file
    !-----------------------------------------------------------------------

    if ( present(noptr) )then
       ptrfile = .not. noptr
    else
       ptrfile = .true.
    end if

    ! --------------------------------------------
    ! Open restart file
    ! --------------------------------------------

    call restFile_open( flag='write', file=file, ncid=ncid )

    ! --------------------------------------------
    ! Define dimensions and variables
    ! --------------------------------------------

    call restFile_dimset ( ncid )

    ! Define restart file variables

    call timemgr_restart_io(ncid, flag='define')

    call SubgridRest(bounds, ncid, flag='define' )

    call accumulRest( ncid, flag='define' )

    call atm2lnd_vars%restart (bounds, ncid, flag='define')

    call canopystate_vars%restart (bounds, ncid, flag='define')

    call energyflux_vars%restart (bounds, ncid, flag='define')

    call frictionvel_vars% restart (bounds, ncid, flag='define')

    call lakestate_vars%restart (bounds, ncid, flag='define')

    call photosyns_vars%restart (bounds, ncid, flag='define')

    call soilhydrology_vars%restart (bounds, ncid, flag='define')

    call soilstate_vars%restart (bounds, ncid, flag='define')

    call solarabs_vars%restart (bounds, ncid, flag='define')

    call temperature_vars%restart (bounds, ncid, flag='define')

    call waterflux_vars%restart (bounds, ncid, flag='define')

    call waterstate_vars%restart (bounds, ncid, flag='define', &
         watsat_col=soilstate_vars%watsat_col(bounds%begc:bounds%endc,:)) 

    call aerosol_vars%restart (bounds, ncid,  flag='define', &
         h2osoi_ice_col=waterstate_vars%h2osoi_ice_col(bounds%begc:bounds%endc,:), &
         h2osoi_liq_col=waterstate_vars%h2osoi_liq_col(bounds%begc:bounds%endc,:))

    call surfalb_vars%restart (bounds, ncid, flag='define', &
         tlai_patch=canopystate_vars%tlai_patch(bounds%begp:bounds%endp), &
         tsai_patch=canopystate_vars%tsai_patch(bounds%begp:bounds%endp))

    if (use_lch4) then
       call ch4_vars%restart(bounds, ncid, flag='define')
    end if

    if (use_cn) then

       call cnstate_vars%Restart(bounds, ncid, flag='define')

       call carbonstate_vars%restart(bounds, ncid, flag='define', carbon_type='c12', &
               cnstate_vars=cnstate_vars)
       if (use_c13) then
          call c13_carbonstate_vars%restart(bounds, ncid, flag='define', carbon_type='c13', &
               c12_carbonstate_vars=carbonstate_vars, cnstate_vars=cnstate_vars)
       end if
       if (use_c14) then
          call c14_carbonstate_vars%restart(bounds, ncid, flag='define', carbon_type='c14', &
               c12_carbonstate_vars=carbonstate_vars, cnstate_vars=cnstate_vars)
       end if

       call carbonflux_vars%restart(bounds, ncid, flag='define')
       call nitrogenflux_vars%Restart(bounds, ncid, flag='define')
       call nitrogenstate_vars%Restart(bounds, ncid, flag='define', cnstate_vars=cnstate_vars)

       call phosphorusflux_vars%Restart(bounds, ncid, flag='define')
       call phosphorusstate_vars%Restart(bounds, ncid, flag='define', cnstate_vars=cnstate_vars)

    end if

    if (use_cndv) then
       call dgvs_vars%Restart(bounds, ncid, flag='define')
    end if

    if (use_ed) then
       call EDRest(bounds,  ncid, flag='define', &
            waterstate_vars=waterstate_vars, &
            canopystate_vars=canopystate_vars, &
            EDbio_vars=EDbio_vars, &
            carbonstate_vars=carbonstate_vars, &
            nitrogenstate_vars=nitrogenstate_vars, &
            carbonflux_vars=carbonflux_vars) 
    end if

    if (use_betr) then
       call tracerstate_vars%Restart(bounds, ncid, flag='define', betrtracer_vars=betrtracer_vars)
       call tracerflux_vars%Restart( bounds, ncid, flag='define', betrtracer_vars=betrtracer_vars)
       call tracercoeff_vars%Restart(bounds, ncid, flag='define', betrtracer_vars=betrtracer_vars)
    endif

    if (present(rdate)) then 
       call hist_restart_ncd (bounds, ncid, flag='define', rdate=rdate )
    end if

    call restFile_enddef( ncid )

    ! --------------------------------------------
    ! Write restart file variables
    ! --------------------------------------------
    
    call timemgr_restart_io( ncid, flag='write' )

    call SubgridRest(bounds, ncid, flag='write' )

    call accumulRest( ncid, flag='write' )

    call atm2lnd_vars%restart (bounds, ncid, flag='write')

    call canopystate_vars%restart (bounds, ncid, flag='write')

    call energyflux_vars%restart (bounds, ncid, flag='write')

    call frictionvel_vars% restart (bounds, ncid, flag='write')

    call lakestate_vars%restart (bounds, ncid, flag='write')

    call photosyns_vars%restart (bounds, ncid, flag='write')

    call soilhydrology_vars%restart (bounds, ncid, flag='write')

    call soilstate_vars%restart (bounds, ncid, flag='write')

    call solarabs_vars%restart (bounds, ncid, flag='write')

    call temperature_vars%restart (bounds, ncid, flag='write')

    call waterflux_vars%restart (bounds, ncid, flag='write')

    call waterstate_vars%restart (bounds, ncid, flag='write',  &
         watsat_col=soilstate_vars%watsat_col(bounds%begc:bounds%endc,:) )

    call aerosol_vars%restart (bounds, ncid,  flag='write', &
         h2osoi_ice_col=waterstate_vars%h2osoi_ice_col(bounds%begc:bounds%endc,:), &
         h2osoi_liq_col=waterstate_vars%h2osoi_liq_col(bounds%begc:bounds%endc,:) )

    call surfalb_vars%restart (bounds, ncid, flag='write',  &
         tlai_patch=canopystate_vars%tlai_patch(bounds%begp:bounds%endp), &
         tsai_patch=canopystate_vars%tsai_patch(bounds%begp:bounds%endp))

    if (use_lch4) then
       call ch4_vars%restart(  bounds, ncid, flag='write' )
    end if

    if (use_cn) then
       call cnstate_vars%Restart(bounds, ncid, flag='write')
       call carbonstate_vars%restart(bounds, ncid, flag='write', &
            carbon_type='c12', cnstate_vars=cnstate_vars)
       if (use_c13) then
          call c13_carbonstate_vars%restart(bounds, ncid, flag='write', &
               c12_carbonstate_vars=carbonstate_vars, carbon_type='c13', &
	       cnstate_vars=cnstate_vars)
       end if
       if (use_c14) then
          call c14_carbonstate_vars%restart(bounds, ncid, flag='write', &
               c12_carbonstate_vars=carbonstate_vars, carbon_type='c14', &
	       cnstate_vars=cnstate_vars )
       end if

       call carbonflux_vars%restart(bounds, ncid, flag='write')

       call nitrogenflux_vars%Restart(bounds, ncid, flag='write')
       call nitrogenstate_vars%Restart(bounds, ncid, flag='write', cnstate_vars=cnstate_vars)

       call phosphorusflux_vars%Restart(bounds, ncid, flag='write')
       call phosphorusstate_vars%Restart(bounds, ncid, flag='write', cnstate_vars=cnstate_vars)

    end if

    if (use_cndv) then
       call dgvs_vars%Restart(bounds, ncid, flag='write')
    end if

    if (use_ed) then
       call EDRest( bounds, ncid, flag='write', & 
            waterstate_vars=waterstate_vars, &
            canopystate_vars=canopystate_vars, &
            EDbio_vars=EDbio_vars, & 
            carbonstate_vars=carbonstate_vars, &
            nitrogenstate_vars=nitrogenstate_vars, &
            carbonflux_vars=carbonflux_vars) 
    end if

    if (use_betr) then
       call tracerstate_vars%Restart(bounds, ncid, flag='write', betrtracer_vars=betrtracer_vars)
       call tracerflux_vars%Restart( bounds, ncid, flag='write', betrtracer_vars=betrtracer_vars)
       call tracercoeff_vars%Restart(bounds, ncid, flag='write', betrtracer_vars=betrtracer_vars)
    endif

    call hist_restart_ncd (bounds, ncid, flag='write' )

    ! --------------------------------------------
    ! Close restart file and write restart pointer file
    ! --------------------------------------------
    
    call restFile_close( ncid )
    call restFile_closeRestart( file )
    
    ! Write restart pointer file
    
    if ( ptrfile ) call restFile_write_pfile( file )
    
    ! Write out diagnostic info

    if (masterproc) then
       write(iulog,*) 'Successfully wrote out restart data at nstep = ',get_nstep()
       write(iulog,'(72a1)') ("-",i=1,60)
    end if
    
  end subroutine restFile_write

  !-----------------------------------------------------------------------
  subroutine restFile_read( bounds, file,                                             &
       atm2lnd_vars, aerosol_vars, canopystate_vars, cnstate_vars,                    &
       carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars, carbonflux_vars, &
       ch4_vars, dgvs_vars, energyflux_vars, frictionvel_vars, lakestate_vars,        &
       nitrogenstate_vars, nitrogenflux_vars, photosyns_vars, soilhydrology_vars,     &
       soilstate_vars, solarabs_vars, surfalb_vars, temperature_vars,                 &
       waterflux_vars, waterstate_vars, EDbio_vars,                                   &
       phosphorusstate_vars,phosphorusflux_vars,                                      &
       betrtracer_vars, tracerstate_vars, tracerflux_vars, tracercoeff_vars)
    !
    ! !DESCRIPTION:
    ! Read a CLM restart file.
    !
    ! !USES:
    use subgridRestMod   , only : SubgridRest, subgridRest_read_cleanup
    use accumulMod       , only : accumulRest
    use histFileMod      , only : hist_restart_ncd
    !
    ! !ARGUMENTS:
    character(len=*)         , intent(in)    :: file  ! output netcdf restart file
    type(bounds_type)        , intent(in)    :: bounds  
    type(atm2lnd_type)       , intent(inout) :: atm2lnd_vars
    type(aerosol_type)       , intent(inout) :: aerosol_vars
    type(canopystate_type)   , intent(inout) :: canopystate_vars
    type(cnstate_type)       , intent(inout) :: cnstate_vars
    type(carbonstate_type)   , intent(inout) :: carbonstate_vars
    type(carbonstate_type)   , intent(inout) :: c13_carbonstate_vars
    type(carbonstate_type)   , intent(inout) :: c14_carbonstate_vars
    type(carbonflux_type)    , intent(inout) :: carbonflux_vars
    type(ch4_type)           , intent(inout) :: ch4_vars
    type(dgvs_type)          , intent(inout) :: dgvs_vars
    type(energyflux_type)    , intent(inout) :: energyflux_vars
    type(frictionvel_type)   , intent(inout) :: frictionvel_vars
    type(lakestate_type)     , intent(inout) :: lakestate_vars
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
    type(photosyns_type)     , intent(inout) :: photosyns_vars
    type(soilhydrology_type) , intent(inout) :: soilhydrology_vars
    type(soilstate_type)     , intent(inout) :: soilstate_vars
    type(solarabs_type)      , intent(inout) :: solarabs_vars
    type(temperature_type)   , intent(inout) :: temperature_vars
    type(surfalb_type)       , intent(inout) :: surfalb_vars
    type(waterstate_type)    , intent(inout) :: waterstate_vars
    type(waterflux_type)     , intent(inout) :: waterflux_vars
    type(EDbio_type)         , intent(inout) :: EDbio_vars
    type(phosphorusstate_type) , intent(inout) :: phosphorusstate_vars
    type(phosphorusflux_type)  , intent(inout) :: phosphorusflux_vars
    type(tracerstate_type)   , intent(inout) :: tracerstate_vars ! due to Betrrest call
    type(BeTRTracer_Type)    , intent(in)    :: betrtracer_vars
    type(tracerflux_type)    , intent(inout) :: tracerflux_vars
    type(tracercoeff_type)   , intent(inout) :: tracercoeff_vars
    !
    ! !LOCAL VARIABLES:
    type(file_desc_t) :: ncid ! netcdf id
    integer :: i              ! index
    !-----------------------------------------------------------------------

    ! Open file

    call restFile_open( flag='read', file=file, ncid=ncid )

    ! Read file

    call restFile_dimcheck( ncid )

    call SubgridRest(bounds, ncid, flag='read')

    call accumulRest( ncid, flag='read' )

    call atm2lnd_vars%restart (bounds, ncid, flag='read')

    call canopystate_vars%restart (bounds, ncid, flag='read')

    call energyflux_vars%restart (bounds, ncid, flag='read')

    call frictionvel_vars% restart (bounds, ncid, flag='read')

    call lakestate_vars%restart (bounds, ncid, flag='read')

    call photosyns_vars%restart (bounds, ncid, flag='read')

    call soilhydrology_vars%restart (bounds, ncid, flag='read')

    call soilstate_vars%restart (bounds, ncid, flag='read')

    call solarabs_vars%restart (bounds, ncid, flag='read')

    call temperature_vars%restart (bounds, ncid, flag='read')

    call waterflux_vars%restart (bounds, ncid, flag='read')

    call waterstate_vars%restart (bounds, ncid,  flag='read', &
         watsat_col=soilstate_vars%watsat_col(bounds%begc:bounds%endc,:) )

    ! The following write statement is needed for some tests to pass with pgi 13.9 on
    ! yellowstone, presumably due to a compiler bug. Without this, the following two
    ! tests were failing:
    ! ERS_D_Mmpi-serial.1x1_brazil.ICLM45CNED.yellowstone_pgi.clm-edTest
    ! ERS_Lm3.1x1_smallvilleIA.ICLM45BGCCROP.yellowstone_pgi
    write(iulog,*) 'about to call aerosol_vars%restart: ', ubound(waterstate_vars%h2osoi_ice_col)

    call aerosol_vars%restart (bounds, ncid, flag='read', &
         h2osoi_ice_col=waterstate_vars%h2osoi_ice_col(bounds%begc:bounds%endc,:), &
         h2osoi_liq_col=waterstate_vars%h2osoi_liq_col(bounds%begc:bounds%endc,:) ) 

    call surfalb_vars%restart (bounds, ncid,  flag='read', &
         tlai_patch=canopystate_vars%tlai_patch(bounds%begp:bounds%endp), &
         tsai_patch=canopystate_vars%tsai_patch(bounds%begp:bounds%endp))

    if (use_lch4) then
       call ch4_vars%restart(  bounds, ncid, flag='read' )
    end if

    if (use_cn) then
       call cnstate_vars%Restart(bounds, ncid, flag='read')
       call carbonstate_vars%restart(bounds, ncid, flag='read', &
            carbon_type='c12', cnstate_vars=cnstate_vars)
       if (use_c13) then
          call c13_carbonstate_vars%restart(bounds, ncid, flag='read', &
               c12_carbonstate_vars=carbonstate_vars, carbon_type='c13', &
	       cnstate_vars=cnstate_vars)
       end if
       if (use_c14) then
          call c14_carbonstate_vars%restart(bounds, ncid, flag='read', &
               c12_carbonstate_vars=carbonstate_vars, carbon_type='c14', &
	       cnstate_vars=cnstate_vars)
       end if

       call carbonflux_vars%restart(bounds, ncid, flag='read')

       call nitrogenflux_vars%Restart(bounds, ncid, flag='read')
       call nitrogenstate_vars%Restart(bounds, ncid, flag='read', cnstate_vars=cnstate_vars)

       call phosphorusflux_vars%Restart(bounds, ncid, flag='read')
       call phosphorusstate_vars%Restart(bounds, ncid, flag='read', cnstate_vars=cnstate_vars)

    end if

    if (use_cndv) then
       call dgvs_vars%Restart(bounds, ncid, flag='read')
    end if

    if (use_ed) then
       call EDRest( bounds, ncid, flag='read', & 
            waterstate_vars=waterstate_vars, &
            canopystate_vars=canopystate_vars, &
            EDbio_vars=EDbio_vars, & 
            carbonstate_vars=carbonstate_vars, &
            nitrogenstate_vars=nitrogenstate_vars, &
            carbonflux_vars=carbonflux_vars) 
    end if

    if (use_betr) then
      call tracerstate_vars%Restart(bounds, ncid, flag='read',betrtracer_vars=betrtracer_vars)
      call tracerflux_vars%Restart( bounds, ncid, flag='read',betrtracer_vars=betrtracer_vars)
      call tracercoeff_vars%Restart(bounds, ncid, flag='read', betrtracer_vars=betrtracer_vars)
    endif
        
    call hist_restart_ncd (bounds, ncid, flag='read')

    ! Do error checking on file
    
    call restFile_check_consistency(bounds, ncid)

    ! Close file 

    call subgridRest_read_cleanup
    call restFile_close( ncid )

    ! Write out diagnostic info

    if (masterproc) then
       write(iulog,'(72a1)') ("-",i=1,60)
       write(iulog,*) 'Successfully read restart data for restart run'
       write(iulog,*)
    end if

  end subroutine restFile_read

  !-----------------------------------------------------------------------
  subroutine restFile_getfile( file, path )
    !
    ! !DESCRIPTION:
    ! Determine and obtain netcdf restart file
    !
    ! !USES:
    use clm_varctl, only : caseid, nrevsn, nsrest, brnch_retain_casename
    use clm_varctl, only : nsrContinue, nsrBranch
    use fileutils , only : getfil
    !
    ! !ARGUMENTS:
    character(len=*), intent(out) :: file  ! name of netcdf restart file
    character(len=*), intent(out) :: path  ! full pathname of netcdf restart file
    !
    ! !LOCAL VARIABLES:
    integer :: status                      ! return status
    integer :: length                      ! temporary          
    character(len=256) :: ftest,ctest      ! temporaries
    !-----------------------------------------------------------------------

    ! Continue run:
    ! Restart file pathname is read restart pointer file 

    if (nsrest==nsrContinue) then
       call restFile_read_pfile( path )
       call getfil( path, file, 0 )
    end if

    ! Branch run: 
    ! Restart file pathname is obtained from namelist "nrevsn"
    ! Check case name consistency (case name must be different for branch run, 
    ! unless namelist specification states otherwise)

    if (nsrest==nsrBranch) then
       length = len_trim(nrevsn)
       if (nrevsn(length-2:length) == '.nc') then
          path = trim(nrevsn) 
       else
          path = trim(nrevsn) // '.nc'
       end if
       call getfil( path, file, 0 )

       ! tcraig, adding xx. and .clm2 makes this more robust
       ctest = 'xx.'//trim(caseid)//'.clm2'
       ftest = 'xx.'//trim(file)
       status = index(trim(ftest),trim(ctest))
       if (status /= 0 .and. .not.(brnch_retain_casename)) then
          if (masterproc) then
             write(iulog,*) 'Must change case name on branch run if ',&
                  'brnch_retain_casename namelist is not set'
             write(iulog,*) 'previous case filename= ',trim(file),&
                  ' current case = ',trim(caseid), &
                  ' ctest = ',trim(ctest), &
                  ' ftest = ',trim(ftest)
          end if
          call endrun(msg=errMsg(__FILE__, __LINE__)) 
       end if
    end if

  end subroutine restFile_getfile

  !-----------------------------------------------------------------------
  subroutine restFile_read_pfile( pnamer )
    !
    ! !DESCRIPTION:
    ! Setup restart file and perform necessary consistency checks
    !
    ! !USES:
    use fileutils , only : opnfil, getavu, relavu
    use clm_varctl, only : rpntfil, rpntdir, inst_suffix
    !
    ! !ARGUMENTS:
    character(len=*), intent(out) :: pnamer ! full path of restart file
    !
    ! !LOCAL VARIABLES:
    !EOP
    integer :: i                  ! indices
    integer :: nio                ! restart unit
    integer :: status             ! substring check status
    character(len=256) :: locfn   ! Restart pointer file name
    !-----------------------------------------------------------------------

    ! Obtain the restart file from the restart pointer file. 
    ! For restart runs, the restart pointer file contains the full pathname 
    ! of the restart file. For branch runs, the namelist variable 
    ! [nrevsn] contains the full pathname of the restart file. 
    ! New history files are always created for branch runs.

    if (masterproc) then
       write(iulog,*) 'Reading restart pointer file....'
    endif

    nio = getavu()
    locfn = trim(rpntdir) //'/'// trim(rpntfil)//trim(inst_suffix)
    call opnfil (locfn, nio, 'f')
    read (nio,'(a256)') pnamer
    call relavu (nio)

    if (masterproc) then
       write(iulog,*) 'Reading restart data.....'
       write(iulog,'(72a1)') ("-",i=1,60)
    end if

  end subroutine restFile_read_pfile

  !-----------------------------------------------------------------------
  subroutine restFile_closeRestart( file )
    !
    ! !DESCRIPTION:
    ! Close restart file and write restart pointer file if
    ! in write mode, otherwise just close restart file if in read mode
    !
    ! !USES:
    use clm_time_manager, only : is_last_step
    !
    ! !ARGUMENTS:
    character(len=*) , intent(in) :: file  ! local output filename
    !
    ! !CALLED FROM:
    ! subroutine restart in this module
    !
    ! !REVISION HISTORY:
    ! Author: Mariana Vertenstein
    !
    !
    ! !LOCAL VARIABLES:
    !EOP
    integer :: i                   !index
    !-----------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'Successfully wrote local restart file ',trim(file)
       write(iulog,'(72a1)') ("-",i=1,60)
       write(iulog,*)
    end if

  end subroutine restFile_closeRestart

  !-----------------------------------------------------------------------
  subroutine restFile_write_pfile( fnamer )
    !
    ! !DESCRIPTION:
    ! Open restart pointer file. Write names of current netcdf restart file.
    !
    ! !USES:
    use clm_varctl, only : rpntdir, rpntfil, inst_suffix
    use fileutils , only : relavu
    use fileutils , only : getavu, opnfil
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: fnamer
    !
    ! !LOCAL VARIABLES:
    integer :: m                    ! index
    integer :: nio                  ! restart pointer file
    character(len=256) :: filename  ! local file name
    !-----------------------------------------------------------------------

    if (masterproc) then
       nio = getavu()
       filename= trim(rpntdir) //'/'// trim(rpntfil)//trim(inst_suffix)
       call opnfil( filename, nio, 'f' )

       write(nio,'(a)') fnamer
       call relavu( nio )
       write(iulog,*)'Successfully wrote local restart pointer file'
    end if

  end subroutine restFile_write_pfile

  !-----------------------------------------------------------------------
  subroutine restFile_open( flag, file, ncid )

    use clm_time_manager, only : get_nstep

    character(len=*),  intent(in) :: flag ! flag to specify read or write
    character(len=*),  intent(in) :: file ! filename
    type(file_desc_t), intent(out):: ncid ! netcdf id

    integer :: omode                              ! netCDF dummy variable
    character(len= 32) :: subname='restFile_open' ! subroutine name

    if (flag == 'write') then

       ! Create new netCDF file (in define mode) and set fill mode
       ! to "no fill" to optimize performance

       if (masterproc) then	
          write(iulog,*)
          write(iulog,*)'restFile_open: writing restart dataset at ',&
               trim(file), ' at nstep = ',get_nstep()
          write(iulog,*)
       end if
       call ncd_pio_createfile(ncid, trim(file))

    else if (flag == 'read') then

       ! Open netcdf restart file

       if (masterproc) then
          write(iulog,*) 'Reading restart dataset'
       end if
       call ncd_pio_openfile (ncid, trim(file), 0)

    end if

  end subroutine restFile_open

  !-----------------------------------------------------------------------
  character(len=256) function restFile_filename( rdate )
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use clm_varctl, only : caseid, inst_suffix
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: rdate   ! input date for restart file name 
    !-----------------------------------------------------------------------

    restFile_filename = "./"//trim(caseid)//".clm2"//trim(inst_suffix)//&
         ".r."//trim(rdate)//".nc"
    if (masterproc) then
       write(iulog,*)'writing restart file ',trim(restFile_filename),' for model date = ',rdate
    end if

  end function restFile_filename

  !------------------------------------------------------------------------
  subroutine restFile_dimset( ncid )
    !
    ! !DESCRIPTION:
    ! Read/Write initial data from/to netCDF instantaneous initial data file
    !
    ! !USES:
    use clm_time_manager , only : get_nstep
    use clm_varctl       , only : caseid, ctitle, version, username, hostname, fsurdat
    use clm_varctl       , only : flanduse_timeseries, conventions, source
    use clm_varpar       , only : numrad, nlevlak, nlevsno, nlevgrnd, nlevurb, nlevcan, nlevtrc_full
    use clm_varpar       , only : cft_lb, cft_ub, maxpatch_glcmec
    use decompMod        , only : get_proc_global
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: ncid
    !
    ! !LOCAL VARIABLES:
    integer :: dimid               ! netCDF dimension id
    integer :: numg                ! total number of gridcells across all processors
    integer :: numl                ! total number of landunits across all processors
    integer :: numc                ! total number of columns across all processors
    integer :: nump                ! total number of pfts across all processors
    integer :: numCohort           ! total number of cohorts across all processors
    integer :: ier                 ! error status
    integer :: strlen_dimid        ! string dimension id
    character(len=  8) :: curdate  ! current date
    character(len=  8) :: curtime  ! current time
    character(len=256) :: str
    character(len= 32) :: subname='restFile_dimset' ! subroutine name
    !------------------------------------------------------------------------

    call get_proc_global(ng=numg, nl=numl, nc=numc, np=nump, nCohorts=numCohort)

    ! Define dimensions

    call ncd_defdim(ncid , nameg      , numg           ,  dimid)
    call ncd_defdim(ncid , namel      , numl           ,  dimid)
    call ncd_defdim(ncid , namec      , numc           ,  dimid)
    call ncd_defdim(ncid , namep      , nump           ,  dimid)
    call ncd_defdim(ncid , nameCohort , numCohort      ,  dimid)

    call ncd_defdim(ncid , 'levgrnd' , nlevgrnd       ,  dimid)
    call ncd_defdim(ncid , 'levurb'  , nlevurb        ,  dimid)
    call ncd_defdim(ncid , 'levlak'  , nlevlak        ,  dimid)
    call ncd_defdim(ncid , 'levsno'  , nlevsno        ,  dimid)
    call ncd_defdim(ncid , 'levsno1' , nlevsno+1      ,  dimid)
    call ncd_defdim(ncid , 'levtot'  , nlevsno+nlevgrnd, dimid)
    call ncd_defdim(ncid , 'numrad'  , numrad         ,  dimid)
    call ncd_defdim(ncid , 'levcan'  , nlevcan        ,  dimid)
    call ncd_defdim(ncid , 'string_length', 64        ,  dimid)
    call ncd_defdim(ncid , 'levtrc'  , nlevtrc_full   ,  dimid)    
    if (create_glacier_mec_landunit) then
       call ncd_defdim(ncid , 'glc_nec', maxpatch_glcmec, dimid)
    end if

    ! Define global attributes

    call ncd_putatt(ncid, NCD_GLOBAL, 'Conventions', trim(conventions))
    call getdatetime(curdate, curtime)
    str = 'created on ' // curdate // ' ' // curtime
    call ncd_putatt(ncid, NCD_GLOBAL, 'history' , trim(str))
    call ncd_putatt(ncid, NCD_GLOBAL, 'username', trim(username))
    call ncd_putatt(ncid, NCD_GLOBAL, 'host'    , trim(hostname))
    call ncd_putatt(ncid, NCD_GLOBAL, 'version' , trim(version))
    call ncd_putatt(ncid, NCD_GLOBAL, 'source'  , trim(source))
    str = '$Id: restFileMod.F90 41292 2012-10-26 13:51:45Z erik $'
    call ncd_putatt(ncid, NCD_GLOBAL, 'revision_id'    , trim(str))
    call ncd_putatt(ncid, NCD_GLOBAL, 'case_title'     , trim(ctitle))
    call ncd_putatt(ncid, NCD_GLOBAL, 'case_id'        , trim(caseid))
    call ncd_putatt(ncid, NCD_GLOBAL, 'surface_dataset', trim(fsurdat))
    call ncd_putatt(ncid, NCD_GLOBAL, 'flanduse_timeseries', trim(flanduse_timeseries))
    call ncd_putatt(ncid, NCD_GLOBAL, 'title', 'CLM Restart information')
    if (create_glacier_mec_landunit) then
       call ncd_putatt(ncid, ncd_global, 'created_glacier_mec_landunits', 'true')
    else
       call ncd_putatt(ncid, ncd_global, 'created_glacier_mec_landunits', 'false')
    end if

    call restFile_add_ipft_metadata(ncid)
    call restFile_add_icol_metadata(ncid)
    call restFile_add_ilun_metadata(ncid)

  end subroutine restFile_dimset

  !-----------------------------------------------------------------------
  subroutine restFile_add_ilun_metadata(ncid)
    !
    ! !DESCRIPTION:
    ! Add global metadata defining landunit types
    !
    ! !USES:
    use landunit_varcon, only : max_lunit, landunit_names, landunit_name_length
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: ncid ! local file id
    !
    ! !LOCAL VARIABLES:
    integer :: ltype  ! landunit type
    character(len=*), parameter :: att_prefix = 'ilun_'  ! prefix for attributes
    character(len=len(att_prefix)+landunit_name_length) :: attname ! attribute name

    character(len=*), parameter :: subname = 'restFile_add_ilun_metadata'
    !-----------------------------------------------------------------------
    
    do ltype = 1, max_lunit
       attname = att_prefix // landunit_names(ltype)
       call ncd_putatt(ncid, ncd_global, attname, ltype)
    end do

  end subroutine restFile_add_ilun_metadata

  !-----------------------------------------------------------------------
  subroutine restFile_add_icol_metadata(ncid)
    !
    ! !DESCRIPTION:
    ! Add global metadata defining column types
    !
    ! !USES:
    use column_varcon, only : icol_roof, icol_sunwall, icol_shadewall, icol_road_imperv, icol_road_perv
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: ncid ! local file id
    !
    ! !LOCAL VARIABLES:
    character(len=*), parameter :: att_prefix = 'icol_'  ! prefix for attributes

    character(len=*), parameter :: subname = 'restFile_add_icol_metadata'
    !-----------------------------------------------------------------------
    
    ! Unlike ilun and ipft, the column names currently do not exist in column_varcon.
    ! This is partly because of the trickiness of encoding column values for crop &
    ! icemec.

    call ncd_putatt(ncid, ncd_global, att_prefix // 'vegetated_or_bare_soil', 1) 
    call ncd_putatt(ncid, ncd_global, att_prefix // 'crop'                  , 2) 
    call ncd_putatt(ncid, ncd_global, att_prefix // 'crop_noncompete'       , '2*100+m, m=cft_lb,cft_ub')
    call ncd_putatt(ncid, ncd_global, att_prefix // 'landice'               , 3) 
    call ncd_putatt(ncid, ncd_global, att_prefix // 'landice_multiple_elevation_classes', '4*100+m, m=1,glcnec')  
    call ncd_putatt(ncid, ncd_global, att_prefix // 'deep_lake'             , 5) 
    call ncd_putatt(ncid, ncd_global, att_prefix // 'wetland'               , 6) 
    call ncd_putatt(ncid, ncd_global, att_prefix // 'urban_roof'            , icol_roof)
    call ncd_putatt(ncid, ncd_global, att_prefix // 'urban_sunwall'         , icol_sunwall)
    call ncd_putatt(ncid, ncd_global, att_prefix // 'urban_shadewall'       , icol_shadewall)
    call ncd_putatt(ncid, ncd_global, att_prefix // 'urban_impervious_road' , icol_road_imperv)
    call ncd_putatt(ncid, ncd_global, att_prefix // 'urban_pervious_road'   , icol_road_perv)

  end subroutine restFile_add_icol_metadata

  !-----------------------------------------------------------------------
  subroutine restFile_add_ipft_metadata(ncid)
    !
    ! !DESCRIPTION:
    ! Add global metadata defining pft types
    !
    ! !USES:
    use clm_varpar, only : natpft_lb, mxpft, cft_lb, cft_ub
    use pftvarcon , only : pftname_len, pftname
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: ncid ! local file id
    !
    ! !LOCAL VARIABLES:
    integer :: ptype  ! pft type
    character(len=*), parameter :: att_prefix = 'ipft_'   ! prefix for attributes
    character(len=len(att_prefix)+pftname_len) :: attname ! attribute name

    character(len=*), parameter :: subname = 'restFile_add_ipft_metadata'
    !-----------------------------------------------------------------------
    
    do ptype = natpft_lb, mxpft
       attname = att_prefix // pftname(ptype)
       call ncd_putatt(ncid, ncd_global, attname, ptype)
    end do

    call ncd_putatt(ncid, ncd_global, 'cft_lb', cft_lb)
    call ncd_putatt(ncid, ncd_global, 'cft_ub', cft_ub)

  end subroutine restFile_add_ipft_metadata

  !-----------------------------------------------------------------------
  subroutine restFile_dimcheck( ncid )
    !
    ! !DESCRIPTION:
    ! Check dimensions of restart file
    !
    ! !USES:
    use decompMod,  only : get_proc_global
    use clm_varpar, only : nlevsno, nlevlak, nlevgrnd, nlevurb
    use clm_varctl, only : single_column, nsrest, nsrStartup
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: ncid
    !
    ! !LOCAL VARIABLES:
    integer :: numg      ! total number of gridcells across all processors
    integer :: numl      ! total number of landunits across all processors
    integer :: numc      ! total number of columns across all processors
    integer :: nump      ! total number of pfts across all processors
    integer :: numCohort ! total number of cohorts across all processors
    character(len=32) :: subname='restFile_dimcheck' ! subroutine name
    !-----------------------------------------------------------------------

    ! Get relevant sizes

    if ( .not. single_column .or. nsrest /= nsrStartup )then
       call get_proc_global(ng=numg, nl=numl, nc=numc, np=nump, nCohorts=numCohort)
       call check_dim(ncid, nameg, numg)
       call check_dim(ncid, namel, numl)
       call check_dim(ncid, namec, numc)
       call check_dim(ncid, namep, nump)
       if ( use_ed ) call check_dim(ncid, nameCohort  , numCohort)
    end if
    call check_dim(ncid, 'levsno'  , nlevsno)
    call check_dim(ncid, 'levgrnd' , nlevgrnd)
    call check_dim(ncid, 'levurb'  , nlevurb)
    call check_dim(ncid, 'levlak'  , nlevlak) 

  end subroutine restFile_dimcheck

  !-----------------------------------------------------------------------
  subroutine restFile_enddef( ncid )
    !
    ! !DESCRIPTION:
    ! Read a CLM restart file.
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: ncid
    !-----------------------------------------------------------------------

    call ncd_enddef(ncid)

  end subroutine restFile_enddef

  !-----------------------------------------------------------------------
  subroutine restFile_close( ncid )
    !
    ! !DESCRIPTION:
    ! Read a CLM restart file.
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: ncid
    !
    ! !LOCAL VARIABLES:
    character(len=32) :: subname='restFile_close' ! subroutine name
    !-----------------------------------------------------------------------

    call ncd_pio_closefile(ncid)

  end subroutine restFile_close

  !-----------------------------------------------------------------------
  subroutine restFile_check_consistency(bounds, ncid)
    !
    ! !DESCRIPTION:
    ! Perform some consistency checks on the restart file
    !
    ! !USES:
    use subgridRestMod, only : subgridRest_check_consistency
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in)    :: bounds  ! bounds
    type(file_desc_t), intent(inout) :: ncid    ! netcdf id
    !
    ! !LOCAL VARIABLES:
    logical :: check_finidat_fsurdat_consistency ! whether to check consistency between fsurdat on finidat file and current fsurdat
    logical :: check_finidat_year_consistency    ! whether to check consistency between year on finidat file and current year
    logical :: check_finidat_pct_consistency     ! whether to check consistency between pct_pft on finidat file and surface dataset
    
    character(len=*), parameter :: subname = 'restFile_check_consistency'
    !-----------------------------------------------------------------------
    
    call restFile_read_consistency_nl( &
         check_finidat_fsurdat_consistency, &
         check_finidat_year_consistency, &
         check_finidat_pct_consistency)

    if (check_finidat_fsurdat_consistency) then
       call restFile_check_fsurdat(ncid)
    end if

    if (check_finidat_year_consistency) then
       call restFile_check_year(ncid)
    end if

    if (check_finidat_pct_consistency) then
       call subgridRest_check_consistency(bounds)
    end if

  end subroutine restFile_check_consistency

  !-----------------------------------------------------------------------
  subroutine restFile_read_consistency_nl( &
       check_finidat_fsurdat_consistency, &
       check_finidat_year_consistency, &
       check_finidat_pct_consistency)

    !
    ! !DESCRIPTION:
    ! Read namelist settings related to finidat consistency checks
    !
    ! !USES:
    use fileutils      , only : getavu, relavu
    use clm_nlUtilsMod , only : find_nlgroup_name
    use controlMod     , only : NLFilename
    use shr_mpi_mod    , only : shr_mpi_bcast
    !
    ! !ARGUMENTS:
    logical, intent(out) :: check_finidat_fsurdat_consistency
    logical, intent(out) :: check_finidat_year_consistency
    logical, intent(out) :: check_finidat_pct_consistency
    !
    ! !LOCAL VARIABLES:
    integer :: nu_nml    ! unit for namelist file
    integer :: nml_error ! namelist i/o error flag

    character(len=*), parameter :: subname = 'restFile_read_consistency_nl'
    !-----------------------------------------------------------------------

    namelist /finidat_consistency_checks/ &
         check_finidat_fsurdat_consistency, &
         check_finidat_year_consistency, &
         check_finidat_pct_consistency

    ! Set default namelist values
    check_finidat_fsurdat_consistency = .true.
    check_finidat_year_consistency = .true.
    check_finidat_pct_consistency = .true.

    ! Read namelist
    if (masterproc) then
       nu_nml = getavu()
       open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
       call find_nlgroup_name(nu_nml, 'finidat_consistency_checks', status=nml_error)
       if (nml_error == 0) then
          read(nu_nml, nml=finidat_consistency_checks,iostat=nml_error)
          if (nml_error /= 0) then
             call endrun(msg='ERROR reading finidat_consistency_checks namelist'//errMsg(__FILE__, __LINE__))
          end if
       end if
       close(nu_nml)
       call relavu( nu_nml )
    endif

    call shr_mpi_bcast (check_finidat_fsurdat_consistency, mpicom)
    call shr_mpi_bcast (check_finidat_year_consistency, mpicom)
    call shr_mpi_bcast (check_finidat_pct_consistency, mpicom)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) 'finidat_consistency_checks settings:'
       write(iulog,nml=finidat_consistency_checks)
       write(iulog,*) ' '
    end if

  end subroutine restFile_read_consistency_nl

  !-----------------------------------------------------------------------
  subroutine restFile_check_fsurdat(ncid)
    !
    ! !DESCRIPTION:
    ! Check consistency of the fsurdat value on the restart file and the current fsurdat
    !
    ! !USES:
    use fileutils  , only : get_filename
    use clm_varctl , only : fname_len, fsurdat, flanduse_timeseries
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: ncid    ! netcdf id
    !
    ! !LOCAL VARIABLES:
    character(len=fname_len) :: fsurdat_rest  ! fsurdat from the restart file (includes full path)
    character(len=fname_len) :: filename_cur  ! current fsurdat file name
    character(len=fname_len) :: filename_rest ! fsurdat file name from restart file (does NOT include full path)
    
    character(len=*), parameter :: subname = 'restFile_check_fsurdat'
    !-----------------------------------------------------------------------
    
    ! Only do this check for a transient run. The problem with doing this check for a non-
    ! transient run is the transition from transient to non-transient: It is legitimate to
    ! run with an 1850 surface dataset and a pftdyn file, then use the restart file from
    ! that run to start a present-day (non-transient) run, which would use a 2000 surface
    ! dataset.
    if (flanduse_timeseries /= ' ') then
       call ncd_getatt(ncid, NCD_GLOBAL, 'surface_dataset', fsurdat_rest)

       ! Compare file names, ignoring path
       filename_cur = get_filename(fsurdat)
       filename_rest = get_filename(fsurdat_rest)

       if (filename_rest /= filename_cur) then
          if (masterproc) then
             write(iulog,*) 'ERROR: Initial conditions file (finidat) was generated from a different surface dataset'
             write(iulog,*) 'than the one being used for the current simulation (fsurdat).'
             write(iulog,*) 'Current fsurdat: ', trim(filename_cur)
             write(iulog,*) 'Surface dataset used to generate initial conditions file: ', trim(filename_rest)
             write(iulog,*)
             write(iulog,*) 'Possible solutions to this problem:'
             write(iulog,*) '(1) Make sure you are using the correct surface dataset and initial conditions file'
             write(iulog,*) '(2) If you generated the surface dataset and/or initial conditions file yourself,'
             write(iulog,*) '    then you may need to manually change the surface_dataset global attribute on the'
             write(iulog,*) '    initial conditions file (e.g., using ncatted)'
             write(iulog,*) '(3) If you are confident that you are using the correct surface dataset and initial conditions file,'
             write(iulog,*) '    yet are still experiencing this error, then you can bypass this check by setting:'
             write(iulog,*) '      check_finidat_fsurdat_consistency = .false.'
             write(iulog,*) '    in user_nl_clm'
             write(iulog,*) ' '
          end if
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
    end if

  end subroutine restFile_check_fsurdat

  !-----------------------------------------------------------------------
  subroutine restFile_check_year(ncid)
    !
    ! !DESCRIPTION:
    ! Make sure year on the restart file is consistent with the current model year
    !
    ! !USES:
    use clm_time_manager, only : get_curr_date, get_rest_date
    use clm_varctl      , only : fname_len, flanduse_timeseries
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: ncid    ! netcdf id
    !
    ! !LOCAL VARIABLES:
    logical                  :: att_found                ! whether the attribute was found on the netcdf file
    character(len=fname_len) :: flanduse_timeseries_rest ! flanduse_timeseries from the restart file
    integer                  :: year                     ! current model year
    integer                  :: mon                      ! current model month
    integer                  :: day                      ! current model day of month
    integer                  :: tod                      ! current model time of day
    integer                  :: rest_year                ! year from restart file

    character(len=*), parameter :: subname = 'restFile_check_year'
    !-----------------------------------------------------------------------
    
    ! Only do this check for a transient run
    if (flanduse_timeseries /= ' ') then
       ! Determine if the restart file was generated from a transient run; if so, we will
       ! do this consistency check. For backwards compatibility, we allow for the
       ! possibility that the flanduse_timeseries attribute was not on the restart file;
       ! in that case, we act as if the restart file was generated from a non-transient
       ! run, thus skipping this check.
       call check_att(ncid, NCD_GLOBAL, 'flanduse_timeseries', att_found)
       if (att_found) then
          call ncd_getatt(ncid, NCD_GLOBAL, 'flanduse_timeseries', flanduse_timeseries_rest)
       else
          write(iulog,*) ' '
          write(iulog,*) subname//' WARNING: flanduse_timeseries attribute not found on restart file'
          write(iulog,*) 'Assuming that the restart file was generated from a non-transient run,'
          write(iulog,*) 'and thus skipping the year check'
          write(iulog,*) ' '

          flanduse_timeseries_rest = ' '
       end if
       
       ! If the restart file was generated from a transient run, then confirm that the
       ! year of the restart file matches the current model year.
       if (flanduse_timeseries_rest /= ' ') then
          call get_curr_date(year, mon, day, tod)
          call get_rest_date(ncid, rest_year)
          if (year /= rest_year) then
             if (masterproc) then
                write(iulog,*) 'ERROR: Current model year does not match year on initial conditions file (finidat)'
                write(iulog,*) 'Current year: ', year
                write(iulog,*) 'Year on initial conditions file: ', rest_year
                write(iulog,*) ' '
                write(iulog,*) 'This match is a requirement when both:'
                write(iulog,*) '(a) The current run is a transient run, and'
                write(iulog,*) '(b) The initial conditions file was generated from a transient run'
                write(iulog,*) ' '
                write(iulog,*) 'Possible solutions to this problem:'
                write(iulog,*) '(1) Make sure RUN_STARTDATE is set correctly'
                write(iulog,*) '(2) Make sure you are using the correct initial conditions file (finidat)'
                write(iulog,*) '(3) If you are confident that you are using the correct start date and initial conditions file,'
                write(iulog,*) '    yet are still experiencing this error, then you can bypass this check by setting:'
                write(iulog,*) '      check_finidat_year_consistency = .false.'
                write(iulog,*) '    in user_nl_clm'
                write(iulog,*) ' '
             end if
             call endrun(msg=errMsg(__FILE__, __LINE__))
          end if  ! year /= rest_year
       end if  ! flanduse_timeseries_rest /= ' '
    end if  ! fpftdyn /= ' '

  end subroutine restFile_check_year



end module restFileMod



