module clm_initializeMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clm_initializeMod
!
! !DESCRIPTION:
! Performs land model initialization
!
! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use spmdMod         , only : masterproc
  use shr_sys_mod     , only : shr_sys_flush
  use abortutils      , only : endrun
  use clm_varctl      , only : nsrest, nsrStartup, nsrContinue, nsrBranch, &
                               create_glacier_mec_landunit, iulog, &
                               use_cn, use_cndv
  use clm_varsur      , only : wtxy, vegxy, topoxy
  use perf_mod        , only : t_startf, t_stopf
  use ncdio_pio
  use mct_mod

! !PUBLIC TYPES:
  implicit none
  save

  private    ! By default everything is private

! !PUBLIC MEMBER FUNCTIONS:
  public :: initialize1  ! Phase one initialization
  public :: initialize2  ! Phase two initialization
!
! !REVISION HISTORY:
! Created by Gordon Bonan, Sam Levis and Mariana Vertenstein
!
!
! !PRIVATE MEMBER FUNCTIONS:
  private header         ! echo version numbers
  private do_restread    ! read a restart file
!-----------------------------------------------------------------------
! !PRIVATE DATA MEMBERS: None

!EOP
!-----------------------------------------------------------------------
contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: initialize1
!
! !INTERFACE:
  subroutine initialize1( )
!
! !DESCRIPTION:
! Land model initialization.
! o Initializes run control variables via the [clm_inparm] namelist.
! o Reads surface data on model grid.
! o Defines the multiple plant types and fraction areas for each surface type.
! o Builds the appropriate subgrid <-> grid mapping indices and weights.
! o Set up parallel processing.
! o Initializes time constant variables.
! o Reads restart data for a restart or branch run.
! o Reads initial data and initializes the time variant variables for an initial run.
! o Initializes history file output.
! o Initializes river routing model.
! o Initializes accumulation variables.
!
! !USES:
    use clmtypeInitMod  , only : initClmtype
    use clm_varpar      , only : maxpatch, clm_varpar_init
    use clm_varctl      , only : fsurdat, fatmlndfrc, flndtopo, fglcmask, noland 
    use pftvarcon       , only : pftconrd
    use decompInitMod   , only : decompInit_lnd, decompInit_glcp
    use decompMod       , only : get_proc_bounds
    use domainMod       , only : domain_check, ldomain, domain_init
    use surfrdMod       , only : surfrd_get_globmask, surfrd_get_grid, surfrd_get_topo, &
                                 surfrd_get_data 
    use controlMod      , only : control_init, control_print, nlfilename
    use UrbanInputMod   , only : UrbanInput
    use ncdio_pio       , only : ncd_pio_init
    use clm_atmlnd      , only : init_atm2lnd_type, init_lnd2atm_type, clm_a2l, clm_l2a
    use clm_glclnd      , only : init_glc2lnd_type, init_lnd2glc_type, clm_x2s, clm_s2x
    use initGridCellsMod, only : initGridCells
!
! !ARGUMENTS:
!
! !REVISION HISTORY:
! Created by Gordon Bonan, Sam Levis and Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: ier                   ! error status
    integer  :: i,j,n,k               ! loop indices
    integer  :: nl                    ! gdc and glo lnd indices
    integer  :: ns, ni, nj            ! global grid sizes
    logical  :: isgrid2d              ! true => global grid is regular lat/lon
    integer  :: begp, endp            ! clump beg and ending pft indices
    integer  :: begc, endc            ! clump beg and ending column indices
    integer  :: begl, endl            ! clump beg and ending landunit indices
    integer  :: begg, endg            ! clump beg and ending gridcell indices
    integer ,pointer  :: amask(:)     ! global land mask
    character(len=32) :: subname = 'initialize1' ! subroutine name
!-----------------------------------------------------------------------

    call t_startf('clm_init1')

    ! ------------------------------------------------------------------------
    ! Initialize run control variables, timestep
    ! ------------------------------------------------------------------------

    call header()

    if (masterproc) then
       write(iulog,*) 'Attempting to initialize the land model .....'
       write(iulog,*)
       call shr_sys_flush(iulog)
    endif

    call control_init()
    call clm_varpar_init()
    call ncd_pio_init()

    if (masterproc) call control_print()

    ! ------------------------------------------------------------------------
    ! Read in global land grid and land mask (amask)- needed to set decomposition
    ! ------------------------------------------------------------------------

    ! global memory for amask is allocate in surfrd_get_glomask - must be
    ! deallocated below
    if (masterproc) then
       write(iulog,*) 'Attempting to read global land mask from ',trim(fatmlndfrc)
       call shr_sys_flush(iulog)
    endif
    call surfrd_get_globmask(filename=fatmlndfrc, mask=amask, ni=ni, nj=nj)

    ! Exit early if no valid land points
    if ( all(amask == 0) )then
       if (masterproc) write(iulog,*) trim(subname)//': no valid land points do NOT run clm'
       noland = .true.
       return
    end if

    ! Determine clm decomposition

    call decompInit_lnd(ni, nj, amask)
    deallocate(amask)

    ! Get grid and land fraction (set ldomain)

    if (masterproc) then
       write(iulog,*) 'Attempting to read ldomain from ',trim(fatmlndfrc)
       call shr_sys_flush(iulog)
    endif
    if (create_glacier_mec_landunit) then
       call surfrd_get_grid(ldomain, fatmlndfrc, fglcmask)
    else
       call surfrd_get_grid(ldomain, fatmlndfrc)
    endif
    if (masterproc) then
       call domain_check(ldomain)
    endif
    ldomain%mask = 1  !!! TODO - is this needed?

    ! Get topo if appropriate (set ldomain%topo)

    if (flndtopo /= " ") then
       if (masterproc) then
          write(iulog,*) 'Attempting to read atm topo from ',trim(flndtopo)
          call shr_sys_flush(iulog)
       endif
       call surfrd_get_topo(ldomain, flndtopo)  
    endif

    ! Initialize urban model input (initialize urbinp data structure)

    call UrbanInput(mode='initialize')

    ! Allocate surface grid dynamic memory (for wtxy and vegxy arrays)
    ! Allocate additional dynamic memory for glacier_mec topo and thickness

    call get_proc_bounds(begg, endg)
    allocate (vegxy(begg:endg,maxpatch), wtxy(begg:endg,maxpatch), stat=ier)   
    if (create_glacier_mec_landunit) then
       allocate (topoxy(begg:endg,maxpatch), stat=ier)
    else
       allocate (topoxy(1,1), stat=ier)
    endif
    if (ier /= 0) then
       write(iulog,*)'initialize allocation error'; call endrun()
    endif

    ! Read list of PFTs and their corresponding parameter values
    ! Independent of model resolution, Needs to stay before surfrd_get_data

    call pftconrd()

    ! Read surface dataset and set up vegetation type [vegxy] and 
    ! weight [wtxy] arrays for [maxpatch] subgrid patches.
    
    call surfrd_get_data(ldomain, fsurdat)

    ! Determine decomposition of subgrid scale landunits, columns, pfts

    if (create_glacier_mec_landunit) then
       call decompInit_glcp (ns, ni, nj, ldomain%glcmask)
    else
       call decompInit_glcp (ns, ni, nj)
    endif

    ! Allocate memory and initialize values of clmtype data structures

    call initClmtype()

    ! Initialize atm->lnd, lnd->atm, glc->lnd and lnd->glc data structures

    call init_atm2lnd_type(begg, endg, clm_a2l)
    call init_lnd2atm_type(begg, endg, clm_l2a)
    if (create_glacier_mec_landunit) then
       call init_glc2lnd_type(begg, endg, clm_x2s)
       call init_lnd2glc_type(begg, endg, clm_s2x)
    endif

    ! Build hierarchy and topological info for derived types

    call initGridCells()

    ! Deallocate surface grid dynamic memory (for wtxy and vegxy arrays)

    deallocate (vegxy, wtxy, topoxy)

    call t_stopf('clm_init1')

  end subroutine initialize1

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: initialize2
!
! !INTERFACE:
  subroutine initialize2( )
!
! !DESCRIPTION:
! Land model initialization.
! o Initializes run control variables via the [clm_inparm] namelist.
! o Reads surface data on model grid.
! o Defines the multiple plant types and fraction areas for each surface type.
! o Builds the appropriate subgrid <-> grid mapping indices and weights.
! o Set up parallel processing.
! o Initializes time constant variables.
! o Reads restart data for a restart or branch run.
! o Reads initial data and initializes the time variant variables for an initial run.
! o Initializes history file output.
! o Initializes river routing model.
! o Initializes accumulation variables.
!
! !USES:
    use clm_atmlnd      , only : clm_map2gcell_minimal
    use clm_glclnd      , only : update_clm_s2x
    use clm_varctl      , only : finidat, fpftdyn
    use decompMod       , only : get_proc_clumps, get_proc_bounds
    use filterMod       , only : allocFilters, setFilters
    use histFldsMod     , only : hist_initFlds
    use histFileMod     , only : hist_htapes_build, htapes_fieldlist
    use restFileMod     , only : restFile_getfile, &
                                 restFile_open, restFile_close, restFile_read 
    use accFldsMod      , only : initAccFlds, initAccClmtype
    use mkarbinitMod    , only : mkarbinit
    use pftdynMod       , only : pftdyn_init, pftdyn_interp
    use ndepStreamMod    , only : ndep_init, ndep_interp
    use CNEcosystemDynMod, only : CNEcosystemDynInit
    use pftdynMod             , only : pftwt_init
    use CNDVEcosystemDyniniMod, only : CNDVEcosystemDynini
    use STATICEcosysDynMod , only : EcosystemDynini, readAnnualVegetation
    use STATICEcosysDynMod , only : interpMonthlyVeg
    use DustMod         , only : Dustini
    use clm_time_manager, only : get_curr_date, get_nstep, advance_timestep, &
                                 timemgr_init, timemgr_restart_io, timemgr_restart
    use clm_time_manager, only : get_step_size, get_curr_calday
    use fileutils       , only : getfil
    use UrbanMod        , only : UrbanClumpInit
    use UrbanInitMod    , only : UrbanInitTimeConst, UrbanInitTimeVar, UrbanInitAero 
    use UrbanInputMod   , only : UrbanInput
    use seq_drydep_mod  , only : n_drydep, drydep_method, DD_XLND
    use shr_orb_mod        , only : shr_orb_decl
    use initSurfAlbMod     , only : initSurfAlb, do_initsurfalb 
    use clm_varorb         , only : eccen, mvelpp, lambm0, obliqr
    use VOCEmissionMod  , only : VOCEmission_init


! !Arguments    
    implicit none
!
! !REVISION HISTORY:
! Created by Gordon Bonan, Sam Levis and Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: nl,na,nag             ! indices
    integer  :: i,j,k                 ! indices
    integer  :: yr                    ! current year (0, ...)
    integer  :: mon                   ! current month (1 -> 12)
    integer  :: day                   ! current day (1 -> 31)
    integer  :: ncsec                 ! current time of day [seconds]
    integer  :: nc                    ! clump index
    integer  :: nclumps               ! number of clumps on this processor
    integer  :: begp, endp            ! clump beg and ending pft indices
    integer  :: begc, endc            ! clump beg and ending column indices
    integer  :: begl, endl            ! clump beg and ending landunit indices
    integer  :: begg, endg            ! clump beg and ending gridcell indices
    character(len=256) :: fnamer      ! name of netcdf restart file 
    character(len=256) :: pnamer      ! full pathname of netcdf restart file
    type(file_desc_t)  :: ncid        ! netcdf id
    real(r8) :: dtime                 ! time step increment (sec)
    integer  :: nstep                 ! model time step
    real(r8) :: calday                ! calendar day for nstep
    real(r8) :: caldaym1              ! calendar day for nstep-1
    real(r8) :: declin                ! solar declination angle in radians for nstep
    real(r8) :: declinm1              ! solar declination angle in radians for nstep-1
    real(r8) :: eccf                  ! earth orbit eccentricity factor
    character(len=32) :: subname = 'initialize2' ! subroutine name
!----------------------------------------------------------------------

    call t_startf('clm_init2')

    ! ------------------------------------------------------------------------
    ! Initialize time constant variables 
    ! ------------------------------------------------------------------------

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp )

    ! Initialize Ecosystem Dynamics 

    call t_startf('init_ecosys')
    if (use_cndv) then
       call CNDVEcosystemDynini()
    else if (.not. use_cn) then
       call EcosystemDynini()
    end if

    if (use_cn .or. use_cndv) then
       ! --------------------------------------------------------------
       ! Initialize CLMSP ecosystem dynamics when drydeposition is used
       ! so that estimates of monthly differences in LAI can be computed
       ! --------------------------------------------------------------
       if ( n_drydep > 0 .and. drydep_method == DD_XLND )then
          call EcosystemDynini()
       end if
    end if
    call t_stopf('init_ecosys')

    ! Initialize dust emissions model 

    call t_startf('init_dust')
    call Dustini()
    call t_stopf('init_dust')
    
    ! Initialize MEGAN emissions model 

    call VOCEmission_init( )

    ! ------------------------------------------------------------------------
    ! Initialize time constant urban variables
    ! ------------------------------------------------------------------------

    call t_startf('init_io1')
    call UrbanInitTimeConst()
    call iniTimeConst()

    ! ------------------------------------------------------------------------
    ! Obtain restart file if appropriate
    ! ------------------------------------------------------------------------

    if (do_restread()) then
       call restFile_getfile( file=fnamer, path=pnamer )
    end if

    ! ------------------------------------------------------------------------
    ! Initialize master history list. 
    ! ------------------------------------------------------------------------
    call t_startf('hist_initFlds')

    call hist_initFlds()
    ! On restart process the history namelist. Later the namelist from the restart file
    ! will be used. But, this allows some basic checking to make sure you didn't
    ! try to change the history namelist on restart.
    if (nsrest == nsrContinue ) call htapes_fieldlist()

    call t_stopf('hist_initFlds')
    ! ------------------------------------------------------------------------
    ! Initialize time manager
    ! ------------------------------------------------------------------------

    if (nsrest == nsrStartup) then  
       call timemgr_init()
    else
       call restFile_open( flag='read', file=fnamer, ncid=ncid )
       call timemgr_restart_io( ncid=ncid, flag='read' )
       call restFile_close( ncid=ncid )
       call timemgr_restart()
    end if
    call t_stopf('init_io1')

    ! ------------------------------------------------------------------------
    ! Initialize CN Ecosystem Dynamics (must be after time-manager initialization)
    ! ------------------------------------------------------------------------
    if (use_cn .or. use_cndv) then
       call CNEcosystemDynInit( begc, endc, begp, endp )
    end if

    ! ------------------------------------------------------------------------
    ! Initialize accumulated fields
    ! ------------------------------------------------------------------------

    ! Initialize accumulator fields to be time accumulated for various purposes.
    ! The time manager needs to be initialized before this called is made, since
    ! the step size is needed.

    call t_startf('init_accflds')
    call initAccFlds()
    call t_stopf('init_accflds')

    ! ------------------------------------------------------------------------
    ! Set arbitrary initial conditions for time varying fields 
    ! used in coupled carbon-nitrogen code
    ! ------------------------------------------------------------------------
    
    if (use_cn) then
       call t_startf('init_cninitim')
       if (nsrest == nsrStartup) then
          call CNiniTimeVar()
       end if
       call t_stopf('init_cninitim')
    end if

    ! ------------------------------------------------------------------------
    ! Initialization of dynamic pft weights
    ! ------------------------------------------------------------------------

    ! Determine correct pft weights (interpolate pftdyn dataset if initial run)
    ! Otherwise these are read in for a restart run

    if (use_cndv) then
       call pftwt_init()
    else
       if (fpftdyn /= ' ') then
          call t_startf('init_pftdyn')
          call pftdyn_init()
          call pftdyn_interp( )
          call t_stopf('init_pftdyn')
       end if
    end if

    ! ------------------------------------------------------------------------
    ! Read restart/initial info
    ! ------------------------------------------------------------------------

    ! No weight related information can be contained in the routines,  
    ! "mkarbinit, inicfile and restFile". 

    call t_startf('init_io2')
    if (do_restread()) then
       if (masterproc) write(iulog,*)'reading restart file ',fnamer
       call restFile_read( fnamer )
    else if (nsrest == nsrStartup .and. finidat == ' ') then
       call mkarbinit()
       call UrbanInitTimeVar( )
    end if
    call t_stopf('init_io2')

    ! ------------------------------------------------------------------------
    ! Initialize nitrogen deposition
    ! ------------------------------------------------------------------------

    if (use_cn) then
       call t_startf('init_ndep')
       call ndep_init()
       call ndep_interp()
       call t_stopf('init_ndep')
    end if
    
    ! ------------------------------------------------------------------------
    ! Initialization of model parameterizations that are needed after
    ! restart file is read in
    ! ------------------------------------------------------------------------

    ! ------------------------------------------------------------------------
    ! Initialize history and accumator buffers
    ! ------------------------------------------------------------------------

    call t_startf('init_hist1')

    ! Initialize active history fields. This is only done if not a restart run. 
    ! If a restart run, then this information has already been obtained from the 
    ! restart data read above. Note that routine hist_htapes_build needs time manager 
    ! information, so this call must be made after the restart information has been read.

    if (nsrest == nsrStartup .or. nsrest == nsrBranch) call hist_htapes_build()

    ! Initialize clmtype variables that are obtained from accumulated fields.
    ! This routine is called in an initial run at nstep=0
    ! This routine is also always called for a restart run and must 
    ! therefore be called after the restart file is read in

    call initAccClmtype()

    call t_stopf('init_hist1')

    ! --------------------------------------------------------------
    ! Note - everything below this point needs updated weights
    ! --------------------------------------------------------------

    ! Initialize filters
    
    call t_startf('init_filters')

    call allocFilters()
    nclumps = get_proc_clumps()
!$OMP PARALLEL DO PRIVATE (nc)
    do nc = 1, nclumps
       call setFilters(nc)
    end do
!$OMP END PARALLEL DO

    call t_stopf('init_filters')

    ! Calculate urban "town" roughness length and displacement 
    ! height for urban landunits

    call UrbanInitAero()

    ! Initialize urban radiation model - this uses urbinp data structure

    call UrbanClumpInit()

    ! Finalize urban model initialization
    
    call UrbanInput(mode='finalize')

    !
    ! Even if CN is on, and dry-deposition is active, read CLMSP annual vegetation to get estimates of monthly LAI
    !
    if ( n_drydep > 0 .and. drydep_method == DD_XLND )then
       call readAnnualVegetation()
    end if

    ! End initialization

    call t_startf('init_wlog')
    if (masterproc) then
       write(iulog,*) 'Successfully initialized the land model'
       if (nsrest == nsrStartup) then
          write(iulog,*) 'begin initial run at: '
       else
          write(iulog,*) 'begin continuation run at:'
       end if
       call get_curr_date(yr, mon, day, ncsec)
       write(iulog,*) '   nstep= ',get_nstep(), ' year= ',yr,' month= ',mon,&
            ' day= ',day,' seconds= ',ncsec
       write(iulog,*)
       write(iulog,'(72a1)') ("*",i=1,60)
       write(iulog,*)
    endif
    call t_stopf('init_wlog')

    if (get_nstep() == 0 .or. nsrest == nsrStartup) then
       ! Initialize albedos (correct pft filters are needed)

       if (finidat == ' ' .or. do_initsurfalb) then
          call t_startf('init_orb')
          calday = get_curr_calday()
          call t_startf('init_orbd1')
          call shr_orb_decl( calday, eccen, mvelpp, lambm0, obliqr, declin, eccf )
          call t_stopf('init_orbd1')
          
          dtime = get_step_size()
          caldaym1 = get_curr_calday(offset=-int(dtime))
          call t_startf('init_orbd2')
          call shr_orb_decl( caldaym1, eccen, mvelpp, lambm0, obliqr, declinm1, eccf )
          call t_stopf('init_orbd2')
          
          call t_startf('init_orbSA')
          call initSurfAlb( calday, declin, declinm1 )
          call t_stopf('init_orbSA')
          call t_stopf('init_orb')
       else if ( n_drydep > 0 .and. drydep_method == DD_XLND )then
          ! Call interpMonthlyVeg for dry-deposition so that mlaidiff will be calculated
          ! This needs to be done even if CN or CNDV is on!
          call interpMonthlyVeg()
       end if

       ! Determine gridcell averaged properties to send to atm

       call t_startf('init_map2gc')
       call clm_map2gcell_minimal()
       call t_stopf('init_map2gc')

    end if

    ! Initialize sno export state
    if (create_glacier_mec_landunit) then
       call t_startf('init_create_s2x')
       call update_clm_s2x(init=.true.)
       call t_stopf('init_create_s2x')
    end if

    call t_stopf('clm_init2')

  end subroutine initialize2

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: header
!
! !INTERFACE:
  subroutine header()
!
! !DESCRIPTION:
! Echo and save model version number
!
! !USES:
    use clm_varctl  , only : version
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine initialize in this module
!
! !REVISION HISTORY:
! Created by Gordon Bonan
!
!EOP
!-----------------------------------------------------------------------

    if ( masterproc )then
      write(iulog,*) trim(version)
      write(iulog,*)
    end if

  end subroutine header

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_restread
!
! !INTERFACE:
  logical function do_restread( )
!
! !DESCRIPTION:
! Determine if restart file will be read
!
! !USES:
    use clm_varctl, only : finidat
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine initialize in this module
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

    do_restread = .false.
    if (nsrest == nsrStartup .and. finidat /= ' ') then
       do_restread = .true.
    end if
    if (nsrest == nsrContinue .or. nsrest == nsrBranch) then
       do_restread = .true.
    end if
  end function do_restread
  
end module clm_initializeMod
