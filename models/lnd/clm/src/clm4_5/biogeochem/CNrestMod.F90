module CNrestMod

#if (defined CN)
!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: CNrestMod
! 
! !DESCRIPTION: 
! Read/Write to/from CN info to CLM restart file. 
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use spmdMod     , only : masterproc
  use clm_varctl  , only : iulog, override_bgc_restart_mismatch_dump
  use abortutils  , only : endrun
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: CNrest

! !PUBLIC DATA MEMBERS:
!
! !REVISION HISTORY:
! 11/05/03: Module created by Peter Thornton
!F. Li and S. Levis (11/06/12)
!EOP
!----------------------------------------------------------------------- 

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNRest
!
! !INTERFACE:
  subroutine CNRest ( ncid, flag )
!
! !DESCRIPTION: 
! Read/write CN restart data
!
! !USES:
    use clmtype
    use clm_atmlnd, only : clm_a2l
    use clm_varpar, only : numrad, ndecomp_pools, nlevdecomp
    use decompMod , only : get_proc_bounds
    use clm_time_manager, only : is_restart, get_nstep
    use clm_varcon, only : nlevgrnd
    use ncdio_pio
    use clm_varctl, only : use_c13, use_c14, spinup_state
    use clm_varcon, only : c13ratio, c14ratio, spval
    use nanMod    , only : nan
    use shr_const_mod,only : SHR_CONST_PDB
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t)  :: ncid   ! netcdf id
    character(len=*), intent(in) :: flag   !'read' or 'write'
!
! !CALLED FROM:
! subroutine restart in module restFileMod
!
! !REVISION HISTORY:
! Author: Peter Thornton
!
!
!Other local variables
    real(r8) :: c3_del13c     ! typical del13C for C3 photosynthesis (permil, relative to PDB)
    real(r8) :: c4_del13c     ! typical del13C for C4 photosynthesis (permil, relative to PDB)
    real(r8) :: c3_r1         ! isotope ratio (13c/12c) for C3 photosynthesis
    real(r8) :: c4_r1         ! isotope ratio (13c/12c) for C4 photosynthesis
    real(r8) :: c3_r2         ! isotope ratio (13c/[12c+13c]) for C3 photosynthesis
    real(r8) :: c4_r2         ! isotope ratio (13c/[12c+13c]) for C4 photosynthesis
!   real(r8), pointer :: rc13_annsum_npp(:)
!   real(r8), pointer :: rc13_cannsum_npp(:)
   type(pft_cstate_type), pointer :: pcisos
   type(pft_cstate_type), pointer :: pcbulks
! !LOCAL VARIABLES:
!EOP
    integer :: c,p,j,k,i,l           ! indices 
    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices 
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
    real(r8):: m            ! multiplier for the exit_spinup code
    logical :: readvar      ! determine if variable is on initial file
    character(len=128) :: varname         ! temporary
    type(gridcell_type), pointer :: gptr  ! pointer to gridcell derived subtype
    type(landunit_type), pointer :: lptr  ! pointer to landunit derived subtype
    type(column_type)  , pointer :: cptr  ! pointer to column derived subtype
    type(pft_type)     , pointer :: pptr  ! pointer to pft derived subtype
    integer , pointer :: iptemp(:) ! pointer to memory to be allocated
    integer :: ier                 ! error status
    real(r8), pointer :: ptr1d(:), ptr2d(:,:) !temporary arrays for slicing larger arrays
    integer  :: nstep                    ! time step number
    integer  :: restart_file_spinup_state  ! spinup state as read from restart file, for determining whether to enter or exit spinup mode.
    logical :: exit_spinup = .false.
    logical :: enter_spinup = .false.
    integer :: decomp_cascade_state, restart_file_decomp_cascade_state      ! flags for comparing the model and restart decomposition cascades
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    if ( use_c13 ) then
       pcisos => clm3%g%l%c%p%pc13s
       pcbulks => clm3%g%l%c%p%pcs
    endif
    ! Set pointers into derived type

    gptr => clm3%g
    lptr => clm3%g%l
    cptr => clm3%g%l%c
    pptr => clm3%g%l%c%p


    ! Determine necessary subgrid bounds

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    if ( use_c13 ) then
       c3_del13c = -28._r8
       c4_del13c = -13._r8
       c3_r1 = SHR_CONST_PDB + ((c3_del13c*SHR_CONST_PDB)/1000._r8)
       c3_r2 = c3_r1/(1._r8 + c3_r1)
       c4_r1 = SHR_CONST_PDB + ((c4_del13c*SHR_CONST_PDB)/1000._r8)
       c4_r2 = c4_r1/(1._r8 + c4_r1)
    endif

    !--------------------------------
    ! pft ecophysiological variables 
    !--------------------------------
    
    ! dormant_flag
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='dormant_flag', xtype=ncd_double,  &
            dim1name='pft',long_name='dormancy flag',units='unitless' )
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='dormant_flag', data=pptr%pepv%dormant_flag, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! days_active
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='days_active', xtype=ncd_double,  &
            dim1name='pft',long_name='number of days since last dormancy',units='days' )
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='days_active', data=pptr%pepv%days_active, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! onset_flag
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='onset_flag', xtype=ncd_double,  &
            dim1name='pft',long_name='flag if critical growing degree-day sum is exceeded',units='unitless' )
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='onset_flag', data=pptr%pepv%onset_flag, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! onset_counter
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='onset_counter', xtype=ncd_double,  &
            dim1name='pft',long_name='onset days counter',units='sec' )
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='onset_counter', data=pptr%pepv%onset_counter, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! onset_gddflag
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='onset_gddflag', xtype=ncd_double,  &
            dim1name='pft',long_name='onset flag for growing degree day sum',units='' )
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='onset_gddflag', data=pptr%pepv%onset_gddflag, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! onset_fdd
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='onset_fdd', xtype=ncd_double,  &
            dim1name='pft',long_name='onset freezing degree days counter',units='days' )
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='onset_fdd', data=pptr%pepv%onset_fdd, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! onset_gdd
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='onset_gdd', xtype=ncd_double,  &
            dim1name='pft',long_name='onset growing degree days',units='days' )
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='onset_gdd', data=pptr%pepv%onset_gdd, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! onset_swi
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='onset_swi', xtype=ncd_double,  &
            dim1name='pft',long_name='onset soil water index',units='days' )
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='onset_swi', data=pptr%pepv%onset_swi, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! offset_flag
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='offset_flag', xtype=ncd_double,  &
            dim1name='pft',long_name='offset flag',units='unitless' )
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='offset_flag', data=pptr%pepv%offset_flag, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! offset_counter
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='offset_counter', xtype=ncd_double,  &
            dim1name='pft',long_name='offset days counter',units='sec' )
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='offset_counter', data=pptr%pepv%offset_counter, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! offset_fdd
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='offset_fdd', xtype=ncd_double,  &
            dim1name='pft',long_name='offset freezing degree days counter',units='days' )
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='offset_fdd', data=pptr%pepv%offset_fdd, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! offset_swi
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='offset_swi', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='offset_swi', data=pptr%pepv%offset_swi, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

#if (defined CROP)
    ! fert_counter
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='fert_counter', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='fert_counter', data=pptr%pepv%fert_counter, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

   ! fert
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='fert', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='fert', data=pptr%pnf%fert, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if
#endif

    ! lgsf
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='lgsf', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='lgsf', data=pptr%pepv%lgsf, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! bglfr
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='bglfr', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='bglfr', data=pptr%pepv%bglfr, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! bgtr
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='bgtr', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='bgtr', data=pptr%pepv%bgtr, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! dayl
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='dayl', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='dayl', data=pptr%pepv%dayl, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! prev_dayl
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='prev_dayl', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='prev_dayl', data=pptr%pepv%prev_dayl, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! annavg_t2m
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='annavg_t2m', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='annavg_t2m', data=pptr%pepv%annavg_t2m, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! tempavg_t2m
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='tempavg_t2m', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='tempavg_t2m', data=pptr%pepv%tempavg_t2m, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! gpp
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='gpp_pepv', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='gpp_pepv', data=pptr%pepv%gpp, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! availc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='availc', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='availc', data=pptr%pepv%availc, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! xsmrpool_recover
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='xsmrpool_recover', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='xsmrpool_recover', data=pptr%pepv%xsmrpool_recover, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    if ( use_c13 ) then
       ! xsmrpool_c13ratio
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='xsmrpool_c13ratio', xtype=ncd_double,  &
               dim1name='pft',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='xsmrpool_c13ratio', data=pptr%pepv%xsmrpool_c13ratio, &
               dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if
    endif

    ! alloc_pnow
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='alloc_pnow', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='alloc_pnow', data=pptr%pepv%alloc_pnow, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! c_allometry
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='c_allometry', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='c_allometry', data=pptr%pepv%c_allometry, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! n_allometry
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='n_allometry', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='n_allometry', data=pptr%pepv%n_allometry, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! plant_ndemand
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='plant_ndemand', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='plant_ndemand', data=pptr%pepv%plant_ndemand, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! tempsum_potential_gpp
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='tempsum_potential_gpp', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='tempsum_potential_gpp', data=pptr%pepv%tempsum_potential_gpp, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    !annsum_potential_gpp 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='annsum_potential_gpp', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='annsum_potential_gpp', data=pptr%pepv%annsum_potential_gpp, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! tempmax_retransn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='tempmax_retransn', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='tempmax_retransn', data=pptr%pepv%tempmax_retransn, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! annmax_retransn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='annmax_retransn', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='annmax_retransn', data=pptr%pepv%annmax_retransn, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! avail_retransn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='avail_retransn', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='avail_retransn', data=pptr%pepv%avail_retransn, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
        if (is_restart()) call endrun
       end if 
    end if

    ! plant_nalloc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='plant_nalloc', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='plant_nalloc', data=pptr%pepv%plant_nalloc, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! plant_calloc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='plant_calloc', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='plant_calloc', data=pptr%pepv%plant_calloc, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! excess_cflux
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='excess_cflux', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='excess_cflux', data=pptr%pepv%excess_cflux, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! downreg
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='downreg', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='downreg', data=pptr%pepv%downreg, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! prev_leafc_to_litter
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='prev_leafc_to_litter', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='prev_leafc_to_litter', data=pptr%pepv%prev_leafc_to_litter, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! prev_frootc_to_litter
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='prev_frootc_to_litter', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='prev_frootc_to_litter', data=pptr%pepv%prev_frootc_to_litter, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! tempsum_npp
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='tempsum_npp', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='tempsum_npp', data=pptr%pepv%tempsum_npp, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! annsum_npp
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='annsum_npp', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='annsum_npp', data=pptr%pepv%annsum_npp, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    if ( use_c13 ) then
       ! rc13_canair
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='rc13_canair', xtype=ncd_double,  &
               dim1name='pft',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='rc13_canair', data=pptr%pepv%rc13_canair, &
               dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if
       
       ! rc13_psnsun
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='rc13_psnsun', xtype=ncd_double,  &
               dim1name='pft',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='rc13_psnsun', data=pptr%pepv%rc13_psnsun, &
               dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if
       
       ! rc13_psnsha
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='rc13_psnsha', xtype=ncd_double,  &
               dim1name='pft',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='rc13_psnsha', data=pptr%pepv%rc13_psnsha, &
               dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if
    endif

#if (defined CROP)
    ! grain_flag
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='grain_flag', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='grain_flag', data=pptr%pepv%grain_flag, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if
#endif

    !--------------------------------
    ! pft carbon state variables 
    !--------------------------------

    ! leafc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafc', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='leafc', data=pptr%pcs%leafc, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! leafc_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafc_storage', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='leafc_storage', data=pptr%pcs%leafc_storage, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! leafc_xfer
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafc_xfer', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='leafc_xfer', data=pptr%pcs%leafc_xfer, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! frootc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='frootc', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='frootc', data=pptr%pcs%frootc, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! frootc_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='frootc_storage', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='frootc_storage', data=pptr%pcs%frootc_storage, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    !frootc_xfer 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='frootc_xfer', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='frootc_xfer', data=pptr%pcs%frootc_xfer, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livestemc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livestemc', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='livestemc', data=pptr%pcs%livestemc, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livestemc_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livestemc_storage', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='livestemc_storage', data=pptr%pcs%livestemc_storage, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livestemc_xfer
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livestemc_xfer', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='livestemc_xfer', data=pptr%pcs%livestemc_xfer, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! deadstemc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadstemc', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='deadstemc', data=pptr%pcs%deadstemc, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! deadstemc_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadstemc_storage', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='deadstemc_storage', data=pptr%pcs%deadstemc_storage, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! deadstemc_xfer
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadstemc_xfer', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='deadstemc_xfer', data=pptr%pcs%deadstemc_xfer, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livecrootc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livecrootc', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='livecrootc', data=pptr%pcs%livecrootc, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livecrootc_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livecrootc_storage', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='livecrootc_storage', data=pptr%pcs%livecrootc_storage, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livecrootc_xfer
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livecrootc_xfer', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='livecrootc_xfer', data=pptr%pcs%livecrootc_xfer, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! deadcrootc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadcrootc', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='deadcrootc', data=pptr%pcs%deadcrootc, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! deadcrootc_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadcrootc_storage', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='deadcrootc_storage', data=pptr%pcs%deadcrootc_storage, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! deadcrootc_xfer
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadcrootc_xfer', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='deadcrootc_xfer', data=pptr%pcs%deadcrootc_xfer, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! gresp_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='gresp_storage', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='gresp_storage', data=pptr%pcs%gresp_storage, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! gresp_xfer
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='gresp_xfer', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='gresp_xfer', data=pptr%pcs%gresp_xfer, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='cpool', data=pptr%pcs%cpool, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! xsmrpool
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='xsmrpool', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='xsmrpool', data=pptr%pcs%xsmrpool, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! pft_ctrunc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='pft_ctrunc', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='pft_ctrunc', data=pptr%pcs%pft_ctrunc, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! totvegc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='totvegc', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='totvegc', data=pptr%pcs%totvegc, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if


    ! ! leaf_prof
    ! call cnrest_addfld_decomp(ncid=ncid, varname='leaf_prof', longname='', units='', flag=flag, data_rl=pptr%pps%leaf_prof, readvar=readvar)

    ! ! froot_prof
    ! call cnrest_addfld_decomp(ncid=ncid, varname='froot_prof', longname='', units='', flag=flag, data_rl=pptr%pps%froot_prof, readvar=readvar)

    ! ! croot_prof
    ! call cnrest_addfld_decomp(ncid=ncid, varname='croot_prof', longname='', units='', flag=flag, data_rl=pptr%pps%croot_prof, readvar=readvar)

    ! ! stem_prof
    ! call cnrest_addfld_decomp(ncid=ncid, varname='stem_prof', longname='', units='', flag=flag, data_rl=pptr%pps%stem_prof, readvar=readvar)


     if ( use_c13 ) then
        !--------------------------------
        ! C13 pft carbon state variables 
        !--------------------------------
        
        ! leafc
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='leafc_13', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='leafc_13', data=pptr%pc13s%leafc, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc13s%leafc with atmospheric c13 value'
              do i = begp,endp
                 if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1._r8) then
                    pcisos%leafc(i) = pcbulks%leafc(i) * c3_r2
                 else
                    pcisos%leafc(i) = pcbulks%leafc(i) * c4_r2
                 endif
!                 if (pptr%pcs%leafc(i) .ne. spval .and. pptr%pcs%leafc(i) .ne. nan ) then
!                    pptr%pc13s%leafc(i) = pptr%pcs%leafc(i) * c13ratio
!                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! leafc_storage
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='leafc_storage_13', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='leafc_storage_13', data=pptr%pc13s%leafc_storage, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc13s%leafc_storage with atmospheric c13 value'
              do i = begp,endp
                 if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1._r8) then
                    pcisos%leafc_storage(i) = pcbulks%leafc_storage(i) * c3_r2
                 else
                    pcisos%leafc_storage(i) = pcbulks%leafc_storage(i) * c4_r2
                 endif
!                 if (pptr%pcs%leafc_storage(i) .ne. spval .and. pptr%pcs%leafc_storage(i) .ne. nan ) then
!                    pptr%pc13s%leafc_storage(i) = pptr%pcs%leafc_storage(i) * c13ratio
!                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! leafc_xfer
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='leafc_xfer_13', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='leafc_xfer_13', data=pptr%pc13s%leafc_xfer, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc13s%leafc_xfer with atmospheric c13 value'
              do i = begp,endp
                 if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1._r8) then
                    pcisos%leafc_xfer(i) = pcbulks%leafc_xfer(i) * c3_r2
                 else
                    pcisos%leafc_xfer(i) = pcbulks%leafc_xfer(i) * c4_r2
                 endif
!                 if (pptr%pcs%leafc_xfer(i) .ne. spval .and. pptr%pcs%leafc_xfer(i) .ne. nan ) then
!                    pptr%pc13s%leafc_xfer(i) = pptr%pcs%leafc_xfer(i) * c13ratio
!                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! frootc
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='frootc_13', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='frootc_13', data=pptr%pc13s%frootc, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc13s%frootc with atmospheric c13 value'
              do i = begp,endp
                 if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1._r8) then
                    pcisos%frootc(i) = pcbulks%frootc(i) * c3_r2
                 else
                    pcisos%frootc(i) = pcbulks%frootc(i) * c4_r2
                 endif
!                 if (pptr%pcs%frootc(i) .ne. spval .and. pptr%pcs%frootc(i) .ne. nan ) then
!                    pptr%pc13s%frootc(i) = pptr%pcs%frootc(i) * c13ratio
!                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! frootc_storage
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='frootc_storage_13', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='frootc_storage_13', data=pptr%pc13s%frootc_storage, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc13s%frootc_storage with atmospheric c13 value'
              do i = begp,endp
                 if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1._r8) then
                    pcisos%frootc_storage(i) = pcbulks%frootc_storage(i) * c3_r2
                 else
                    pcisos%frootc_storage(i) = pcbulks%frootc_storage(i) * c4_r2
                 endif
!                 if (pptr%pcs%frootc_storage(i) .ne. spval .and. pptr%pcs%frootc_storage(i) .ne. nan ) then
!                    pptr%pc13s%frootc_storage(i) = pptr%pcs%frootc_storage(i) * c13ratio
!                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        !frootc_xfer 
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='frootc_xfer_13', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='frootc_xfer_13', data=pptr%pc13s%frootc_xfer, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc13s%frootc_xfer with atmospheric c13 value'
              do i = begp,endp
                 if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1._r8) then
                    pcisos%frootc_xfer(i) = pcbulks%frootc_xfer(i) * c3_r2
                 else
                    pcisos%frootc_xfer(i) = pcbulks%frootc_xfer(i) * c4_r2
                 endif
!                 if (pptr%pcs%frootc_xfer(i) .ne. spval .and. pptr%pcs%frootc_xfer(i) .ne. nan ) then
!                    pptr%pc13s%frootc_xfer(i) = pptr%pcs%frootc_xfer(i) * c13ratio
!                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! livestemc
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='livestemc_13', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='livestemc_13', data=pptr%pc13s%livestemc, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc13s%livestemc with atmospheric c13 value'
              do i = begp,endp
                 if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1._r8) then
                    pcisos%livestemc(i) = pcbulks%livestemc(i) * c3_r2
                 else
                    pcisos%livestemc(i) = pcbulks%livestemc(i) * c4_r2
                 endif
!                 if (pptr%pcs%livestemc(i) .ne. spval .and. pptr%pcs%livestemc(i) .ne. nan ) then
!                    pptr%pc13s%livestemc(i) = pptr%pcs%livestemc(i) * c13ratio
!                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! livestemc_storage
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='livestemc_storage_13', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='livestemc_storage_13', data=pptr%pc13s%livestemc_storage, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc13s%livestemc_storage with atmospheric c13 value'
              do i = begp,endp
                 if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1._r8) then
                    pcisos%livestemc_storage(i) = pcbulks%livestemc_storage(i) * c3_r2
                 else
                    pcisos%livestemc_storage(i) = pcbulks%livestemc_storage(i) * c4_r2
                 endif
!                if (pptr%pcs%livestemc_storage(i) .ne. spval .and. pptr%pcs%livestemc_storage(i) .ne. nan ) then
!                    pptr%pc13s%livestemc_storage(i) = pptr%pcs%livestemc_storage(i) * c13ratio
!                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! livestemc_xfer
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='livestemc_xfer_13', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='livestemc_xfer_13', data=pptr%pc13s%livestemc_xfer, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc13s%livestemc_xfer with atmospheric c13 value'
              do i = begp,endp
                 if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1._r8) then
                    pcisos%livestemc_xfer(i) = pcbulks%livestemc_xfer(i) * c3_r2
                 else
                    pcisos%livestemc_xfer(i) = pcbulks%livestemc_xfer(i) * c4_r2
                 endif
!                 if (pptr%pcs%livestemc_xfer(i) .ne. spval .and. pptr%pcs%livestemc_xfer(i) .ne. nan ) then
!                    pptr%pc13s%livestemc_xfer(i) = pptr%pcs%livestemc_xfer(i) * c13ratio
!                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! deadstemc
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='deadstemc_13', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='deadstemc_13', data=pptr%pc13s%deadstemc, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc13s%deadstemc with atmospheric c13 value'
              do i = begp,endp
                 if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1._r8) then
                    pcisos%deadstemc(i) = pcbulks%deadstemc(i) * c3_r2
                 else
                    pcisos%deadstemc(i) = pcbulks%deadstemc(i) * c4_r2
                 endif 
!                 if (pptr%pcs%deadstemc(i) .ne. spval .and. pptr%pcs%deadstemc(i) .ne. nan ) then
!                    pptr%pc13s%deadstemc(i) = pptr%pcs%deadstemc(i) * c13ratio
!                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! deadstemc_storage
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='deadstemc_storage_13', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='deadstemc_storage_13', data=pptr%pc13s%deadstemc_storage, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc13s%deadstemc_storage with atmospheric c13 value'
              do i = begp,endp
                 if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1._r8) then
                    pcisos%deadstemc_storage(i) = pcbulks%deadstemc_storage(i) * c3_r2
                 else
                    pcisos%deadstemc_storage(i) = pcbulks%deadstemc_storage(i) * c4_r2
                 endif
!                 if (pptr%pcs%deadstemc_storage(i) .ne. spval .and. pptr%pcs%deadstemc_storage(i) .ne. nan ) then
!                    pptr%pc13s%deadstemc_storage(i) = pptr%pcs%deadstemc_storage(i) * c13ratio
!                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! deadstemc_xfer
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='deadstemc_xfer_13', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='deadstemc_xfer_13', data=pptr%pc13s%deadstemc_xfer, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc13s%deadstemc_xfer with atmospheric c13 value'
              do i = begp,endp
                 if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1._r8) then
                    pcisos%deadstemc_xfer(i) = pcbulks%deadstemc_xfer(i) * c3_r2
                 else
                    pcisos%deadstemc_xfer(i) = pcbulks%deadstemc_xfer(i) * c4_r2
                 endif
!                 if (pptr%pcs%deadstemc_xfer(i) .ne. spval .and. pptr%pcs%deadstemc_xfer(i) .ne. nan ) then
!                    pptr%pc13s%deadstemc_xfer(i) = pptr%pcs%deadstemc_xfer(i) * c13ratio
!                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! livecrootc
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='livecrootc_13', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='livecrootc_13', data=pptr%pc13s%livecrootc, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc13s%livecrootc with atmospheric c13 value'
              do i = begp,endp
                 if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1._r8) then
                    pcisos%livecrootc(i) = pcbulks%livecrootc(i) * c3_r2
                 else
                    pcisos%livecrootc(i) = pcbulks%livecrootc(i) * c4_r2
                 endif
!                 if (pptr%pcs%livecrootc(i) .ne. spval .and. pptr%pcs%livecrootc(i) .ne. nan ) then
!                    pptr%pc13s%livecrootc(i) = pptr%pcs%livecrootc(i) * c13ratio
!                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! livecrootc_storage
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='livecrootc_storage_13', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='livecrootc_storage_13', data=pptr%pc13s%livecrootc_storage, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc13s%livecrootc_storage with atmospheric c13 value'
              do i = begp,endp
                 if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1._r8) then
                    pcisos%livecrootc_storage(i) = pcbulks%livecrootc_storage(i) * c3_r2
                 else
                    pcisos%livecrootc_storage(i) = pcbulks%livecrootc_storage(i) * c4_r2
                 endif
!                 if (pptr%pcs%livecrootc_storage(i) .ne. spval .and. pptr%pcs%livecrootc_storage(i) .ne. nan ) then
!                    pptr%pc13s%livecrootc_storage(i) = pptr%pcs%livecrootc_storage(i) * c13ratio
!                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! livecrootc_xfer
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='livecrootc_xfer_13', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='livecrootc_xfer_13', data=pptr%pc13s%livecrootc_xfer, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc13s%livecrootc_xfer with atmospheric c13 value'
              do i = begp,endp
                 if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1._r8) then
                    pcisos%livecrootc_xfer(i) = pcbulks%livecrootc_xfer(i) * c3_r2
                 else
                    pcisos%livecrootc_xfer(i) = pcbulks%livecrootc_xfer(i) * c4_r2
                 endif
!                 if (pptr%pcs%livecrootc_xfer(i) .ne. spval .and. pptr%pcs%livecrootc_xfer(i) .ne. nan ) then
!                    pptr%pc13s%livecrootc_xfer(i) = pptr%pcs%livecrootc_xfer(i) * c13ratio
!                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! deadcrootc
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='deadcrootc_13', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='deadcrootc_13', data=pptr%pc13s%deadcrootc, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc13s%deadcrootc with atmospheric c13 value'
              do i = begp,endp
                 if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1._r8) then
                    pcisos%deadcrootc(i) = pcbulks%deadcrootc(i) * c3_r2
                 else
                    pcisos%deadcrootc(i) = pcbulks%deadcrootc(i) * c4_r2
                 endif
!                 if (pptr%pcs%deadcrootc(i) .ne. spval .and. pptr%pcs%deadcrootc(i) .ne. nan ) then
!                    pptr%pc13s%deadcrootc(i) = pptr%pcs%deadcrootc(i) * c13ratio
!                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! deadcrootc_storage
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='deadcrootc_storage_13', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='deadcrootc_storage_13', data=pptr%pc13s%deadcrootc_storage, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc13s%deadcrootc_storage with atmospheric c13 value'
              do i = begp,endp
                 if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1._r8) then
                    pcisos%deadcrootc_storage(i) = pcbulks%deadcrootc_storage(i) * c3_r2
                 else
                    pcisos%deadcrootc_storage(i) = pcbulks%deadcrootc_storage(i) * c4_r2
                 endif
!                 if (pptr%pcs%deadcrootc_storage(i) .ne. spval .and. pptr%pcs%deadcrootc_storage(i) .ne. nan ) then
!                    pptr%pc13s%deadcrootc_storage(i) = pptr%pcs%deadcrootc_storage(i) * c13ratio
!                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! deadcrootc_xfer
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='deadcrootc_xfer_13', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='deadcrootc_xfer_13', data=pptr%pc13s%deadcrootc_xfer, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc13s%deadcrootc_xfer with atmospheric c13 value'
              do i = begp,endp
                 if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1._r8) then
                    pcisos%deadcrootc_xfer(i) = pcbulks%deadcrootc_xfer(i) * c3_r2
                 else
                    pcisos%deadcrootc_xfer(i) = pcbulks%deadcrootc_xfer(i) * c4_r2
                 endif
!                 if (pptr%pcs%deadcrootc_xfer(i) .ne. spval .and. pptr%pcs%deadcrootc_xfer(i) .ne. nan ) then
!                    pptr%pc13s%deadcrootc_xfer(i) = pptr%pcs%deadcrootc_xfer(i) * c13ratio
!                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! gresp_storage
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='gresp_storage_13', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='gresp_storage_13', data=pptr%pc13s%gresp_storage, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc13s%gresp_storage with atmospheric c13 value'
              do i = begp,endp
                 if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1._r8) then
                    pcisos%gresp_storage(i) = pcbulks%gresp_storage(i) * c3_r2
                 else
                    pcisos%gresp_storage(i) = pcbulks%gresp_storage(i) * c4_r2
                 endif
!                 if (pptr%pcs%gresp_storage(i) .ne. spval .and. pptr%pcs%gresp_storage(i) .ne. nan ) then
!                    pptr%pc13s%gresp_storage(i) = pptr%pcs%gresp_storage(i) * c13ratio
!                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! gresp_xfer
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='gresp_xfer_13', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='gresp_xfer_13', data=pptr%pc13s%gresp_xfer, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc13s%gresp_xfer with atmospheric c13 value'
              do i = begp,endp
                 if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1._r8) then
                    pcisos%gresp_xfer(i) = pcbulks%gresp_xfer(i) * c3_r2
                 else
                    pcisos%gresp_xfer(i) = pcbulks%gresp_xfer(i) * c4_r2
                 endif
!                 if (pptr%pcs%gresp_xfer(i) .ne. spval .and. pptr%pcs%gresp_xfer(i) .ne. nan ) then
!                    pptr%pc13s%gresp_xfer(i) = pptr%pcs%gresp_xfer(i) * c13ratio
!                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! cpool
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='cpool_13', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='cpool_13', data=pptr%pc13s%cpool, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc13s%cpool with atmospheric c13 value'
              do i = begp,endp
                 if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1._r8) then
                    pcisos%cpool(i) = pcbulks%cpool(i) * c3_r2
                 else
                    pcisos%cpool(i) = pcbulks%cpool(i) * c4_r2
                 endif
!                 if (pptr%pcs%cpool(i) .ne. spval .and. pptr%pcs%cpool(i) .ne. nan ) then
!                    pptr%pc13s%cpool(i) = pptr%pcs%cpool(i) * c13ratio
!                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! xsmrpool
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='xsmrpool_13', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='xsmrpool_13', data=pptr%pc13s%xsmrpool, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc13s%xsmrpool with atmospheric c13 value'
              do i = begp,endp
                 if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1._r8) then
                    pcisos%xsmrpool(i) = pcbulks%xsmrpool(i) * c3_r2
                 else
                    pcisos%xsmrpool(i) = pcbulks%xsmrpool(i) * c4_r2
                 endif
!                 if (pptr%pcs%xsmrpool(i) .ne. spval .and. pptr%pcs%xsmrpool(i) .ne. nan ) then
!                    pptr%pc13s%xsmrpool(i) = pptr%pcs%xsmrpool(i) * c13ratio
!                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! pft_ctrunc
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='pft_ctrunc_13', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='pft_ctrunc_13', data=pptr%pc13s%pft_ctrunc, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc13s%pft_ctrunc with atmospheric c13 value'
              do i = begp,endp
                 if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1._r8) then
                    pcisos%pft_ctrunc(i) = pcbulks%pft_ctrunc(i) * c3_r2
                 else
                    pcisos%pft_ctrunc(i) = pcbulks%pft_ctrunc(i) * c4_r2
                 endif
!                 if (pptr%pcs%pft_ctrunc(i) .ne. spval .and. pptr%pcs%pft_ctrunc(i) .ne. nan ) then
!                    pptr%pc13s%pft_ctrunc(i) = pptr%pcs%pft_ctrunc(i) * c13ratio
!                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! totvegc
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='totvegc_13', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='totvegc_13', data=pptr%pc13s%totvegc, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc13s%totvegc with atmospheric c13 value'
              do i = begp,endp
                 if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1._r8) then
                    pcisos%totvegc(i) = pcbulks%totvegc(i) * c3_r2
                 else
                    pcisos%totvegc(i) = pcbulks%totvegc(i) * c4_r2
                 endif
!                 if (pptr%pcs%totvegc(i) .ne. spval .and. pptr%pcs%totvegc(i) .ne. nan ) then
!                    pptr%pc13s%totvegc(i) = pptr%pcs%totvegc(i) * c13ratio
!                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
     endif
     
     if ( use_c14 ) then
        !--------------------------------
        ! C14 pft carbon state variables 
        !--------------------------------
        
        ! leafc
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='leafc_14', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='leafc_14', data=pptr%pc14s%leafc, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc14s%leafc with atmospheric c14 value'
              do i = begp,endp
                 if (pptr%pcs%leafc(i) .ne. spval .and. pptr%pcs%leafc(i) .ne. nan ) then
                    pptr%pc14s%leafc(i) = pptr%pcs%leafc(i) * c14ratio
                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        
        ! leafc_storage
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='leafc_storage_14', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='leafc_storage_14', data=pptr%pc14s%leafc_storage, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc14s%leafc_storage with atmospheric c14 value'
              do i = begp,endp
                 if (pptr%pcs%leafc_storage(i) .ne. spval .and. pptr%pcs%leafc_storage(i) .ne. nan ) then
                    pptr%pc14s%leafc_storage(i) = pptr%pcs%leafc_storage(i) * c14ratio
                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! leafc_xfer
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='leafc_xfer_14', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='leafc_xfer_14', data=pptr%pc14s%leafc_xfer, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc14s%leafc_xfer with atmospheric c14 value'
              do i = begp,endp
                 if (pptr%pcs%leafc_xfer(i) .ne. spval .and. pptr%pcs%leafc_xfer(i) .ne. nan ) then
                    pptr%pc14s%leafc_xfer(i) = pptr%pcs%leafc_xfer(i) * c14ratio
                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! frootc
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='frootc_14', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='frootc_14', data=pptr%pc14s%frootc, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc14s%frootc with atmospheric c14 value'
              do i = begp,endp
                 if (pptr%pcs%frootc(i) .ne. spval .and. pptr%pcs%frootc(i) .ne. nan ) then
                    pptr%pc14s%frootc(i) = pptr%pcs%frootc(i) * c14ratio
                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! frootc_storage
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='frootc_storage_14', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='frootc_storage_14', data=pptr%pc14s%frootc_storage, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc14s%frootc_storage with atmospheric c14 value'
              do i = begp,endp
                 if (pptr%pcs%frootc_storage(i) .ne. spval .and. pptr%pcs%frootc_storage(i) .ne. nan ) then
                    pptr%pc14s%frootc_storage(i) = pptr%pcs%frootc_storage(i) * c14ratio
                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        !frootc_xfer 
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='frootc_xfer_14', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='frootc_xfer_14', data=pptr%pc14s%frootc_xfer, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc14s%frootc_xfer with atmospheric c14 value'
              do i = begp,endp
                 if (pptr%pcs%frootc_xfer(i) .ne. spval .and. pptr%pcs%frootc_xfer(i) .ne. nan ) then
                    pptr%pc14s%frootc_xfer(i) = pptr%pcs%frootc_xfer(i) * c14ratio
                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! livestemc
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='livestemc_14', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='livestemc_14', data=pptr%pc14s%livestemc, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc14s%livestemc with atmospheric c14 value'
              do i = begp,endp
                 if (pptr%pcs%livestemc(i) .ne. spval .and. pptr%pcs%livestemc(i) .ne. nan ) then
                    pptr%pc14s%livestemc(i) = pptr%pcs%livestemc(i) * c14ratio
                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! livestemc_storage
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='livestemc_storage_14', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='livestemc_storage_14', data=pptr%pc14s%livestemc_storage, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc14s%livestemc_storage with atmospheric c14 value'
              do i = begp,endp
                 if (pptr%pcs%livestemc_storage(i) .ne. spval .and. pptr%pcs%livestemc_storage(i) .ne. nan ) then
                    pptr%pc14s%livestemc_storage(i) = pptr%pcs%livestemc_storage(i) * c14ratio
                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! livestemc_xfer
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='livestemc_xfer_14', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='livestemc_xfer_14', data=pptr%pc14s%livestemc_xfer, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc14s%livestemc_xfer with atmospheric c14 value'
              do i = begp,endp
                 if (pptr%pcs%livestemc_xfer(i) .ne. spval .and. pptr%pcs%livestemc_xfer(i) .ne. nan ) then
                    pptr%pc14s%livestemc_xfer(i) = pptr%pcs%livestemc_xfer(i) * c14ratio
                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! deadstemc
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='deadstemc_14', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='deadstemc_14', data=pptr%pc14s%deadstemc, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc14s%deadstemc with atmospheric c14 value'
              do i = begp,endp
                 if (pptr%pcs%deadstemc(i) .ne. spval .and. pptr%pcs%deadstemc(i) .ne. nan ) then
                    pptr%pc14s%deadstemc(i) = pptr%pcs%deadstemc(i) * c14ratio
                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! deadstemc_storage
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='deadstemc_storage_14', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='deadstemc_storage_14', data=pptr%pc14s%deadstemc_storage, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc14s%deadstemc_storage with atmospheric c14 value'
              do i = begp,endp
                 if (pptr%pcs%deadstemc_storage(i) .ne. spval .and. pptr%pcs%deadstemc_storage(i) .ne. nan ) then
                    pptr%pc14s%deadstemc_storage(i) = pptr%pcs%deadstemc_storage(i) * c14ratio
                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! deadstemc_xfer
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='deadstemc_xfer_14', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='deadstemc_xfer_14', data=pptr%pc14s%deadstemc_xfer, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc14s%deadstemc_xfer with atmospheric c14 value'
              do i = begp,endp
                 if (pptr%pcs%deadstemc_xfer(i) .ne. spval .and. pptr%pcs%deadstemc_xfer(i) .ne. nan ) then
                    pptr%pc14s%deadstemc_xfer(i) = pptr%pcs%deadstemc_xfer(i) * c14ratio
                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! livecrootc
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='livecrootc_14', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='livecrootc_14', data=pptr%pc14s%livecrootc, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc14s%livecrootc with atmospheric c14 value'
              do i = begp,endp
                 if (pptr%pcs%livecrootc(i) .ne. spval .and. pptr%pcs%livecrootc(i) .ne. nan ) then
                    pptr%pc14s%livecrootc(i) = pptr%pcs%livecrootc(i) * c14ratio
                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! livecrootc_storage
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='livecrootc_storage_14', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='livecrootc_storage_14', data=pptr%pc14s%livecrootc_storage, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc14s%livecrootc_storage with atmospheric c14 value'
              do i = begp,endp
                 if (pptr%pcs%livecrootc_storage(i) .ne. spval .and. pptr%pcs%livecrootc_storage(i) .ne. nan ) then
                    pptr%pc14s%livecrootc_storage(i) = pptr%pcs%livecrootc_storage(i) * c14ratio
                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! livecrootc_xfer
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='livecrootc_xfer_14', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='livecrootc_xfer_14', data=pptr%pc14s%livecrootc_xfer, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc14s%livecrootc_xfer with atmospheric c14 value'
              do i = begp,endp
                 if (pptr%pcs%livecrootc_xfer(i) .ne. spval .and. pptr%pcs%livecrootc_xfer(i) .ne. nan ) then
                    pptr%pc14s%livecrootc_xfer(i) = pptr%pcs%livecrootc_xfer(i) * c14ratio
                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! deadcrootc
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='deadcrootc_14', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='deadcrootc_14', data=pptr%pc14s%deadcrootc, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc14s%deadcrootc with atmospheric c14 value'
              do i = begp,endp
                 if (pptr%pcs%deadcrootc(i) .ne. spval .and. pptr%pcs%deadcrootc(i) .ne. nan ) then
                    pptr%pc14s%deadcrootc(i) = pptr%pcs%deadcrootc(i) * c14ratio
                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! deadcrootc_storage
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='deadcrootc_storage_14', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='deadcrootc_storage_14', data=pptr%pc14s%deadcrootc_storage, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc14s%deadcrootc_storage with atmospheric c14 value'
              do i = begp,endp
                 if (pptr%pcs%deadcrootc_storage(i) .ne. spval .and. pptr%pcs%deadcrootc_storage(i) .ne. nan ) then
                    pptr%pc14s%deadcrootc_storage(i) = pptr%pcs%deadcrootc_storage(i) * c14ratio
                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! deadcrootc_xfer
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='deadcrootc_xfer_14', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='deadcrootc_xfer_14', data=pptr%pc14s%deadcrootc_xfer, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc14s%deadcrootc_xfer with atmospheric c14 value'
              do i = begp,endp
                 if (pptr%pcs%deadcrootc_xfer(i) .ne. spval .and. pptr%pcs%deadcrootc_xfer(i) .ne. nan ) then
                    pptr%pc14s%deadcrootc_xfer(i) = pptr%pcs%deadcrootc_xfer(i) * c14ratio
                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! gresp_storage
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='gresp_storage_14', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='gresp_storage_14', data=pptr%pc14s%gresp_storage, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc14s%gresp_storage with atmospheric c14 value'
              do i = begp,endp
                 if (pptr%pcs%gresp_storage(i) .ne. spval .and. pptr%pcs%gresp_storage(i) .ne. nan ) then
                    pptr%pc14s%gresp_storage(i) = pptr%pcs%gresp_storage(i) * c14ratio
                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! gresp_xfer
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='gresp_xfer_14', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='gresp_xfer_14', data=pptr%pc14s%gresp_xfer, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc14s%gresp_xfer with atmospheric c14 value'
              do i = begp,endp
                 if (pptr%pcs%gresp_xfer(i) .ne. spval .and. pptr%pcs%gresp_xfer(i) .ne. nan ) then
                    pptr%pc14s%gresp_xfer(i) = pptr%pcs%gresp_xfer(i) * c14ratio
                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! cpool
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='cpool_14', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='cpool_14', data=pptr%pc14s%cpool, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc14s%cpool with atmospheric c14 value'
              do i = begp,endp
                 if (pptr%pcs%cpool(i) .ne. spval .and. pptr%pcs%cpool(i) .ne. nan ) then
                    pptr%pc14s%cpool(i) = pptr%pcs%cpool(i) * c14ratio
                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! xsmrpool
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='xsmrpool_14', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='xsmrpool_14', data=pptr%pc14s%xsmrpool, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc14s%xsmrpool with atmospheric c14 value'
              do i = begp,endp
                 if (pptr%pcs%xsmrpool(i) .ne. spval .and. pptr%pcs%xsmrpool(i) .ne. nan ) then
                    pptr%pc14s%xsmrpool(i) = pptr%pcs%xsmrpool(i) * c14ratio
                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! pft_ctrunc
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='pft_ctrunc_14', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='pft_ctrunc_14', data=pptr%pc14s%pft_ctrunc, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc14s%pft_ctrunc with atmospheric c14 value'
              do i = begp,endp
                 if (pptr%pcs%pft_ctrunc(i) .ne. spval .and. pptr%pcs%pft_ctrunc(i) .ne. nan ) then
                    pptr%pc14s%pft_ctrunc(i) = pptr%pcs%pft_ctrunc(i) * c14ratio
                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! totvegc
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='totvegc_14', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='totvegc_14', data=pptr%pc14s%totvegc, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing pptr%pc14s%totvegc with atmospheric c14 value'
              do i = begp,endp
                 if (pptr%pcs%totvegc(i) .ne. spval .and. pptr%pcs%totvegc(i) .ne. nan ) then
                    pptr%pc14s%totvegc(i) = pptr%pcs%totvegc(i) * c14ratio
                 endif
              end do
              if (is_restart()) call endrun
           end if
        end if
        
        ! rc14_atm
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='rc14_atm', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='rc14_atm', data=pptr%pepv%rc14_atm, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
           if (flag=='read' .and. .not. readvar) then
              write(iulog,*) 'initializing cptr%cc14s%rc14_atm with atmospheric c14 value'
              do i = begp,endp
                 pptr%pepv%rc14_atm(i) = c14ratio
              end do
              if (is_restart()) call endrun
           end if
        end if

     endif

    !--------------------------------
    ! pft nitrogen state variables
    !--------------------------------
    
    ! leafn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafn', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='leafn', data=pptr%pns%leafn, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! leafn_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafn_storage', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='leafn_storage', data=pptr%pns%leafn_storage, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! leafn_xfer
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafn_xfer', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='leafn_xfer', data=pptr%pns%leafn_xfer, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! frootn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='frootn', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='frootn', data=pptr%pns%frootn, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! frootn_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='frootn_storage', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='frootn_storage', data=pptr%pns%frootn_storage, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! frootn_xfer
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='frootn_xfer', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='frootn_xfer', data=pptr%pns%frootn_xfer, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livestemn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livestemn', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='livestemn', data=pptr%pns%livestemn, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livestemn_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livestemn_storage', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='livestemn_storage', data=pptr%pns%livestemn_storage, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livestemn_xfer
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livestemn_xfer', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='livestemn_xfer', data=pptr%pns%livestemn_xfer, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! deadstemn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadstemn', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='deadstemn', data=pptr%pns%deadstemn, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    !deadstemn_storage 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadstemn_storage', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='deadstemn_storage', data=pptr%pns%deadstemn_storage, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    !deadstemn_xfer 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadstemn_xfer', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='deadstemn_xfer', data=pptr%pns%deadstemn_xfer, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livecrootn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livecrootn', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='livecrootn', data=pptr%pns%livecrootn, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livecrootn_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livecrootn_storage', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='livecrootn_storage', data=pptr%pns%livecrootn_storage, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    !livecrootn_xfer 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livecrootn_xfer', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='livecrootn_xfer', data=pptr%pns%livecrootn_xfer, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! deadcrootn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadcrootn', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='deadcrootn', data=pptr%pns%deadcrootn, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! deadcrootn_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadcrootn_storage', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='deadcrootn_storage', data=pptr%pns%deadcrootn_storage, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! deadcrootn_xfer
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadcrootn_xfer', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='deadcrootn_xfer', data=pptr%pns%deadcrootn_xfer, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    !retransn 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='retransn', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='retransn', data=pptr%pns%retransn, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! npool
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='npool', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='npool', data=pptr%pns%npool, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! pft_ntrunc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='pft_ntrunc', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='pft_ntrunc', data=pptr%pns%pft_ntrunc, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    !--------------------------------
    ! column physical state variables
    !--------------------------------
    
    ! decl
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='decl', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='decl', data=cptr%cps%decl, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! fpi
    call cnrest_addfld_decomp(ncid=ncid, varname='fpi',                        &
                              longname='fraction of potential immobilization', &
                              units='unitless', flag=flag,                     &
                              data_rl=cptr%cps%fpi_vr, readvar=readvar)

#ifdef VERTSOILC
    ! som_adv_coef
    call cnrest_addfld_decomp(ncid=ncid, varname='som_adv_coef',   &
                      longname='SOM advective flux', units='m/s',  &
                      flag=flag, fill_value=spval,                 &
                      data_rl=cptr%cps%som_adv_coef, readvar=readvar)

    ! som_diffus_coef
    call cnrest_addfld_decomp(ncid=ncid, varname='som_diffus_coef',          &
                      longname='SOM diffusivity due to bio/cryo-turbation',  &
                      units='m^2/s',  flag=flag, fill_value=spval,           &
                      data_rl=cptr%cps%som_diffus_coef, readvar=readvar)
#endif

! #ifdef NITRIF_DENITRIF
!     ! tmean_monthly_max
!     call cnrest_addfld_decomp(ncid=ncid, varname='tmean_monthly_max', longname='', units='', flag=flag, data_rl=cptr%cps%tmean_monthly_max, readvar=readvar)

!     ! tmean_monthly
!     call cnrest_addfld_decomp(ncid=ncid, varname='tmean_monthly', longname='', units='', flag=flag, data_rl=cptr%cps%tmean_monthly, readvar=readvar)
! #endif
    ! fpg
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='fpg', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='fpg', data=cptr%cps%fpg, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! annsum_counter
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='annsum_counter', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='annsum_counter', data=cptr%cps%annsum_counter, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! cannsum_npp
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cannsum_npp', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='cannsum_npp', data=cptr%cps%cannsum_npp, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! col_lag_npp
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='col_lag_npp', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='col_lag_npp', data=cptr%cps%col_lag_npp, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! cannavg_t2m
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cannavg_t2m', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='cannavg_t2m', data=cptr%cps%cannavg_t2m, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

  ! for fire model changed by F. Li and S. Levis
  !  burndate
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='burndate', xtype=ncd_int,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='burndate', data=pptr%pps%burndate, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
      
   !lfc
     if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='lfc', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='lfc', data=cptr%cps%lfc, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if  
   
    !wf
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='wf', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='wf', data=cptr%cps%wf, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    
    !btran2
      if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='btran2', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='btran2', data=pptr%pps%btran2, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if


     
   
    !--------------------------------
    ! column carbon state variables
    !--------------------------------

    do k = 1, ndecomp_pools
       ptr2d => cptr%ccs%decomp_cpools_vr(:,:,k)
       varname=trim(decomp_cascade_con%decomp_pool_name_restart(k))//'c'
       call cnrest_addfld_decomp(ncid=ncid, varname=varname, longname='', units='', flag=flag, data_rl=ptr2d, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          call endrun( 'ERROR:: '//trim(varname)//' is required on an initialization dataset' )
       end if
    end do

    ! col_ctrunc
    varname = 'col_ctrunc'
    call cnrest_addfld_decomp(ncid=ncid, varname=varname, longname='', units='', flag=flag, data_rl=cptr%ccs%col_ctrunc_vr, readvar=readvar)
    if (flag=='read' .and. .not. readvar) then
       call endrun( 'ERROR:: '//trim(varname)//' is required on an initialization dataset' )
    end if

    ! ! nfixation_prof
    ! call cnrest_addfld_decomp(ncid=ncid, varname='nfixation_prof', longname='', units='', flag=flag, data_rl=cptr%cps%nfixation_prof, readvar=readvar)

    ! ! ndep_prof
    ! call cnrest_addfld_decomp(ncid=ncid, varname='ndep_prof', longname='', units='', flag=flag, data_rl=cptr%cps%ndep_prof, readvar=readvar)

    ! altmax
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='altmax', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='altmax', data=cptr%cps%altmax, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! altmax_lastyear
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='altmax_lastyear', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='altmax_lastyear', data=cptr%cps%altmax_lastyear, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! altmax_indx
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='altmax_indx', xtype=ncd_int,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='altmax_indx', data=cptr%cps%altmax_indx, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! altmax_lastyear_indx
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='altmax_lastyear_indx', xtype=ncd_int,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='altmax_lastyear_indx', data=cptr%cps%altmax_lastyear_indx, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if



    ! seedc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='seedc', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='seedc', data=cptr%ccs%seedc, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
        
    ! totlitc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='totlitc', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='totlitc', data=cptr%ccs%totlitc, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! totcolc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='totcolc', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='totcolc', data=cptr%ccs%totcolc, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! prod10c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='prod10c', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='prod10c', data=cptr%ccs%prod10c, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! prod100c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='prod100c', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='prod100c', data=cptr%ccs%prod100c, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! totsomc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='totsomc', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='totsomc', data=cptr%ccs%totsomc, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if


    if ( use_c13 ) then
       !--------------------------------
       ! C13 column carbon state variables
       !--------------------------------
       
       do k = 1, ndecomp_pools
          ptr2d => cptr%cc13s%decomp_cpools_vr(:,:,k)
          varname=trim(decomp_cascade_con%decomp_pool_name_restart(k))//'c_13'
          call cnrest_addfld_decomp(ncid=ncid, varname=varname, longname='', units='', flag=flag, data_rl=ptr2d, readvar=readvar)
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing cptr%cc13s%decomp_cpools_vr with atmospheric c13 value for: '//varname
             do i = begc,endc
                do j = 1, nlevdecomp
                   if (cptr%ccs%decomp_cpools_vr(i,j,k) .ne. spval .and. cptr%ccs%decomp_cpools_vr(i,j,k) .ne. nan ) then
                      cptr%cc13s%decomp_cpools_vr(i,j,k) = cptr%ccs%decomp_cpools_vr(i,j,k) * c3_r2
                   endif
                end do
             end do
          end if
       end do
       
       ! seedc
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='seedc_13', xtype=ncd_double,  &
               dim1name='column',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='seedc_13', data=cptr%cc13s%seedc, &
               dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             do i = begc,endc
                if (cptr%ccs%seedc(i) .ne. spval .and. cptr%ccs%seedc(i) .ne. nan ) then
                   cptr%cc13s%seedc(i) = cptr%ccs%seedc(i) * c3_r2
                endif
             end do
             if (is_restart()) call endrun
          end if
       end if
       
       ! col_ctrunc_13
       call cnrest_addfld_decomp(ncid=ncid, varname='col_ctrunc_13_vr', longname='', units='', flag=flag, data_rl=cptr%cc13s%col_ctrunc_vr, readvar=readvar)

       ! ! col_ctrunc
       ! if (flag == 'define') then
       !    call ncd_defvar(ncid=ncid, varname='col_ctrunc_13', xtype=ncd_double,  &
       !         dim1name='column',long_name='',units='')
       ! else if (flag == 'read' .or. flag == 'write') then
       !    call ncd_io(varname='col_ctrunc_13', data=cptr%cc13s%col_ctrunc, &
       !         dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       !    if (flag=='read' .and. .not. readvar) then
       !       if (is_restart()) call endrun
       !    end if
       ! end if
       
       ! totlitc
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='totlitc_13', xtype=ncd_double,  &
               dim1name='column',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='totlitc_13', data=cptr%cc13s%totlitc, &
               dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             do i = begc,endc
                if (cptr%ccs%totlitc(i) .ne. spval .and. cptr%ccs%totlitc(i) .ne. nan ) then
                   cptr%cc13s%totlitc(i) = cptr%ccs%totlitc(i) * c3_r2
                endif
             end do
             if (is_restart()) call endrun
          end if
       end if
       
       ! totcolc
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='totcolc_13', xtype=ncd_double,  &
               dim1name='column',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='totcolc_13', data=cptr%cc13s%totcolc, &
               dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             do i = begc,endc
                if (cptr%ccs%totcolc(i) .ne. spval .and. cptr%ccs%totcolc(i) .ne. nan ) then
                   cptr%cc13s%totcolc(i) = cptr%ccs%totcolc(i) * c3_r2
                endif
             end do
             if (is_restart()) call endrun
          end if
       end if
       
       ! prod10c
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='prod10c_13', xtype=ncd_double,  &
               dim1name='column',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='prod10c_13', data=cptr%cc13s%prod10c, &
               dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             do i = begc,endc
                if (cptr%ccs%prod10c(i) .ne. spval .and. cptr%ccs%prod10c(i) .ne. nan ) then
                   cptr%cc13s%prod10c(i) = cptr%ccs%prod10c(i) * c3_r2
                endif
             end do
             if (is_restart()) call endrun
          end if
       end if
       
       ! prod100c
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='prod100c_13', xtype=ncd_double,  &
               dim1name='column',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='prod100c_13', data=cptr%cc13s%prod100c, &
               dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             do i = begc,endc
                if (cptr%ccs%prod100c(i) .ne. spval .and. cptr%ccs%prod100c(i) .ne. nan ) then
                   cptr%cc13s%prod100c(i) = cptr%ccs%prod100c(i) * c3_r2
                endif
             end do
             if (is_restart()) call endrun
          end if
       end if
    endif

    if ( use_c14 ) then
       !--------------------------------
       ! C14 column carbon state variables
       !--------------------------------
       
       do k = 1, ndecomp_pools
          ptr2d => cptr%cc14s%decomp_cpools_vr(:,:,k)
          varname=trim(decomp_cascade_con%decomp_pool_name_restart(k))//'c_14'
          call cnrest_addfld_decomp(ncid=ncid, varname=varname, longname='', units='', flag=flag, data_rl=ptr2d, readvar=readvar)
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing cptr%cc14s%decomp_cpools_vr with atmospheric c14 value for: '//varname
             do i = begc,endc
                do j = 1, nlevdecomp
                   if (cptr%ccs%decomp_cpools_vr(i,j,k) .ne. spval .and. cptr%ccs%decomp_cpools_vr(i,j,k) .ne. nan ) then
                      cptr%cc14s%decomp_cpools_vr(i,j,k) = cptr%ccs%decomp_cpools_vr(i,j,k) * c14ratio
                   endif
                end do
             end do
          end if
       end do
       
       ! seedc
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='seedc_14', xtype=ncd_double,  &
               dim1name='column',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='seedc_14', data=cptr%cc14s%seedc, &
               dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing cptr%cc14s%seedc with atmospheric c14 value'
             do i = begc,endc
                if (cptr%ccs%seedc(i) .ne. spval .and. cptr%ccs%seedc(i) .ne. nan ) then
                   cptr%cc14s%seedc(i) = cptr%ccs%seedc(i) * c14ratio
                endif
             end do
             if (is_restart()) call endrun
          end if
       end if
       
       ! col_ctrunc_c14
       call cnrest_addfld_decomp(ncid=ncid, varname='col_ctrunc_14_vr', longname='', units='', flag=flag, data_rl=cptr%cc14s%col_ctrunc_vr, readvar=readvar)

       ! ! col_ctrunc
       ! if (flag == 'define') then
       !    call ncd_defvar(ncid=ncid, varname='col_ctrunc_14', xtype=ncd_double,  &
       !         dim1name='column',long_name='',units='')
       ! else if (flag == 'read' .or. flag == 'write') then
       !    call ncd_io(varname='col_ctrunc_14', data=cptr%cc14s%col_ctrunc, &
       !         dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       !    if (flag=='read' .and. .not. readvar) then
       !       write(iulog,*) 'initializing cptr%cc14s%col_ctrunc with atmospheric c14 value'
       !       do i = begc,endc
       !          if (cptr%ccs%col_ctrunc(i) .ne. spval .and. cptr%ccs%col_ctrunc(i) .ne. nan ) then
       !             cptr%cc14s%col_ctrunc(i) = cptr%ccs%col_ctrunc(i) * c14ratio
       !          endif
       !       end do
       !       if (is_restart()) call endrun
       !    end if
       ! end if
       
       ! totlitc
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='totlitc_14', xtype=ncd_double,  &
               dim1name='column',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='totlitc_14', data=cptr%cc14s%totlitc, &
               dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing cptr%cc14s%totlitc with atmospheric c14 value'
             do i = begc,endc
                if (cptr%ccs%totlitc(i) .ne. spval .and. cptr%ccs%totlitc(i) .ne. nan ) then
                   cptr%cc14s%totlitc(i) = cptr%ccs%totlitc(i) * c14ratio
                endif
             end do
             if (is_restart()) call endrun
          end if
       end if
       
       ! totcolc
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='totcolc_14', xtype=ncd_double,  &
               dim1name='column',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='totcolc_14', data=cptr%cc14s%totcolc, &
               dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing cptr%cc14s%totcolc with atmospheric c14 value'
             do i = begc,endc
                if (cptr%ccs%totcolc(i) .ne. spval .and. cptr%ccs%totcolc(i) .ne. nan ) then
                   cptr%cc14s%totcolc(i) = cptr%ccs%totcolc(i) * c14ratio
                endif
             end do
             if (is_restart()) call endrun
          end if
       end if
       
       ! prod10c
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='prod10c_14', xtype=ncd_double,  &
               dim1name='column',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='prod10c_14', data=cptr%cc14s%prod10c, &
               dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing cptr%cc14s%prod10c with atmospheric c14 value'
             do i = begc,endc
                if (cptr%ccs%prod10c(i) .ne. spval .and. cptr%ccs%prod10c(i) .ne. nan ) then
                   cptr%cc14s%prod10c(i) = cptr%ccs%prod10c(i) * c14ratio
                endif
             end do
             if (is_restart()) call endrun
          end if
       end if
       
       ! prod100c
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='prod100c_14', xtype=ncd_double,  &
               dim1name='column',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='prod100c_14', data=cptr%cc14s%prod100c, &
               dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing cptr%cc14s%prod100c with atmospheric c14 value'
             do i = begc,endc
                if (cptr%ccs%prod100c(i) .ne. spval .and. cptr%ccs%prod100c(i) .ne. nan ) then
                   cptr%cc14s%prod100c(i) = cptr%ccs%prod100c(i) * c14ratio
                endif
             end do
             if (is_restart()) call endrun
          end if
       end if
    endif
       
    !--------------------------------
    ! column nitrogen state variables
    !--------------------------------
    
    ! sminn
    varname='sminn'
    call cnrest_addfld_decomp(ncid=ncid, varname=varname, longname='', units='', flag=flag, data_rl=cptr%cns%sminn_vr, readvar=readvar)
    if (flag=='read' .and. .not. readvar) then
       call endrun( 'ERROR:: '//trim(varname)//' is required on an initialization dataset' )
    end if
    
    ! decomposing N pools
    do k = 1, ndecomp_pools
       ptr2d => cptr%cns%decomp_npools_vr(:,:,k)
       varname=trim(decomp_cascade_con%decomp_pool_name_restart(k))//'n'
       call cnrest_addfld_decomp(ncid=ncid, varname=varname, longname='', units='', flag=flag, data_rl=ptr2d, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          call endrun( 'ERROR:: '//trim(varname)//' is required on an initialization dataset' )
       end if
    end do
    
    
    ! col_ntrunc
    call cnrest_addfld_decomp(ncid=ncid, varname='col_ntrunc', longname='', units='', flag=flag, data_rl=cptr%cns%col_ntrunc_vr, readvar=readvar)

#ifdef NITRIF_DENITRIF
    ! f_nit_vr
    call cnrest_addfld_decomp(ncid=ncid, varname='f_nit_vr', &
                              longname='soil nitrification flux', &
                              units='gN/m3/s', &
                              flag=flag, data_rl=cptr%cnf%f_nit_vr, readvar=readvar)

    ! pot_f_nit_vr
    call cnrest_addfld_decomp(ncid=ncid, varname='pot_f_nit_vr', &
                              longname='potential soil nitrification flux', &
                              units='gN/m3/s', &
                              flag=flag, data_rl=cptr%cnf%pot_f_nit_vr, readvar=readvar)

    ! smin_no3_vr
    varname = 'smin_no3'
    call cnrest_addfld_decomp(ncid=ncid, varname=varname, longname='', units='', flag=flag, data_rl=cptr%cns%smin_no3_vr, readvar=readvar)
    if (flag=='read' .and. .not. readvar) then
       call endrun( 'ERROR:: '//trim(varname)//' is required on an initialization dataset' )
    end if

    ! smin_nh4
    varname = 'smin_nh4'
    call cnrest_addfld_decomp(ncid=ncid, varname=varname, longname='', units='', flag=flag, data_rl=cptr%cns%smin_nh4_vr, readvar=readvar)
    if (flag=='read' .and. .not. readvar) then
       call endrun( 'ERROR:: '//trim(varname)//' is required on an initialization dataset' )
    end if

#endif

    ! totcoln
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='totcoln', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='totcoln', data=cptr%cns%totcoln, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! seedn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='seedn', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='seedn', data=cptr%cns%seedn, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! prod10n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='prod10n', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='prod10n', data=cptr%cns%prod10n, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! prod100n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='prod100n', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='prod100n', data=cptr%cns%prod100n, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! decomp_cascade_state  
    ! the purpose of this is to check to make sure the bgc used matches what the restart file was generated with.  
    ! add info about the SOM decomposition cascade
#ifdef CENTURY_DECOMP
    decomp_cascade_state = 1
#else
    decomp_cascade_state = 0
#endif
    ! add info about the nitrification / denitrification state
#ifdef NITRIF_DENITRIF
    decomp_cascade_state = decomp_cascade_state + 10
#endif
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='decomp_cascade_state', xtype=ncd_int,  &
            long_name='BGC of the model that wrote this restart file: 1s column: 0 = CLM-CN cascade, 1 = Century cascade;' &
            // ' 10s column: 0 = CLM-CN denitrification, 10 = Century denitrification',units='')
    else if (flag == 'read') then
       call ncd_io(varname='decomp_cascade_state', data=restart_file_decomp_cascade_state, &
            ncid=ncid, flag=flag, readvar=readvar) 
       if ( .not. readvar) then
  	  !!! assume, for sake of backwards compatibility, that if decomp_cascade_state is not in the restart file, then the current model state is the same as the prior model state
          restart_file_decomp_cascade_state = decomp_cascade_state
          if ( masterproc ) write(iulog,*) ' CNRest: WARNING!  Restart file does not contain info on decomp_cascade_state used to generate the restart file.  '
          if ( masterproc ) write(iulog,*) '   Assuming the same as current setting: ', decomp_cascade_state
       end if	
    else if (flag == 'write') then
       call ncd_io(varname='decomp_cascade_state', data=decomp_cascade_state, &
            ncid=ncid, flag=flag, readvar=readvar) 
    end if
    if ( flag == 'read' .and. decomp_cascade_state .ne. restart_file_decomp_cascade_state ) then
       if ( masterproc ) then
           write(iulog,*) 'CNRest: ERROR--the decomposition cascade differs between the current model state and the model that wrote the restart file. '
           write(iulog,*) 'This means that the model will be horribly out of equilibrium until after a lengthy spinup. '
           write(iulog,*) 'Stopping here since this is probably an error in configuring the run. If you really wish to proceed, '
           write(iulog,*) 'then override by setting override_bgc_restart_mismatch_dump to .true. in the namelist'
           if ( .not. override_bgc_restart_mismatch_dump ) then
              call endrun( ' CNRest: Stopping. Decomposition cascade mismatch error.')
           endif
        endif
    endif

    ! spinup_state
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='spinup_state', xtype=ncd_int,  &
            long_name='Spinup state of the model that wrote this restart file: 0 = normal model mode, 1 = AD spinup',units='')
    else if (flag == 'read') then
       call ncd_io(varname='spinup_state', data=restart_file_spinup_state, &
            ncid=ncid, flag=flag, readvar=readvar) 
       if ( .not. readvar) then
  	  !!! assume, for sake of backwards compatibility, that if spinup_state is not in the restart file, then the current model state is the same as the prior model state
          restart_file_spinup_state = spinup_state
          if ( masterproc ) write(iulog,*) ' CNRest: WARNING!  Restart file does not contain info on spinup state used to generate the restart file. '
          if ( masterproc ) write(iulog,*) '   Assuming the same as current setting: ', spinup_state
       end if	
    else if (flag == 'write') then
       call ncd_io(varname='spinup_state', data=spinup_state, &
            ncid=ncid, flag=flag, readvar=readvar) 
    end if

    ! now compare the model and restart file spinup states, and either take the model into spinup mode or out of it if they are not identical
    ! taking model out of spinup mode requires multiplying each decomposing pool by the associated AD factor.
    ! putting model into spinup mode requires dividing each decomposing pool by the associated AD factor.
    nstep = get_nstep()  ! only allow this to occur on first timestep of model run.

    if (flag == 'read' .and. spinup_state .ne. restart_file_spinup_state ) then
       if (spinup_state .eq. 0 .and. restart_file_spinup_state .eq. 1 ) then
          if ( masterproc ) write(iulog,*) ' CNRest: taking SOM pools out of AD spinup mode'
          exit_spinup = .true.
       else if (spinup_state .eq. 1 .and. restart_file_spinup_state .eq. 0 ) then
          if ( masterproc ) write(iulog,*) ' CNRest: taking SOM pools into AD spinup mode'
          enter_spinup = .true.
       else
          call endrun(' CNRest: error in entering/exiting spinup.  spinup_state != restart_file_spinup_state, but do not know what to do')
       end if
       if (nstep .ge. 2) then
          call endrun(' CNRest: error in entering/exiting spinup. this should occur only when nstep = 1 ')
       endif
       do k = 1, ndecomp_pools
          if ( exit_spinup ) then
             m = decomp_cascade_con%spinup_factor(k)
          else if ( enter_spinup ) then
             m = 1. / decomp_cascade_con%spinup_factor(k)
          end if
          do c = begc, endc
             do j = 1, nlevdecomp
                clm3%g%l%c%ccs%decomp_cpools_vr(c,j,k) = clm3%g%l%c%ccs%decomp_cpools_vr(c,j,k) * m
                
                if ( use_c13 ) then
                   clm3%g%l%c%cc13s%decomp_cpools_vr(c,j,k) = clm3%g%l%c%cc13s%decomp_cpools_vr(c,j,k) * m
                endif
                
                if ( use_c14 ) then
                   clm3%g%l%c%cc14s%decomp_cpools_vr(c,j,k) = clm3%g%l%c%cc14s%decomp_cpools_vr(c,j,k) * m
                endif
                
                clm3%g%l%c%cns%decomp_npools_vr(c,j,k) = clm3%g%l%c%cns%decomp_npools_vr(c,j,k) * m
             end do
          end do
       end do
    end if

    if ( .not. is_restart() .and. nstep .eq. 1 ) then
       do i = begp, endp
          if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1._r8) then
             pcisos%grainc(i) = pcbulks%grainc(i) * c3_r2
             pcisos%grainc_storage(i) = pcbulks%grainc_storage(i) * c3_r2
             pcisos%grainc_xfer(i) = pcbulks%grainc_xfer(i) * c3_r2
             pcisos%dispvegc(i) = pcbulks%dispvegc(i) * c3_r2
             pcisos%storvegc(i) = pcbulks%storvegc(i) * c3_r2
             pcisos%totvegc(i) = pcbulks%totvegc(i) * c3_r2
             pcisos%totpftc(i) = pcbulks%totpftc(i) * c3_r2
             pcisos%woodc(i) = pcbulks%woodc(i) * c3_r2
          else
             pcisos%grainc(i) = pcbulks%grainc(i) * c4_r2
             pcisos%grainc_storage(i) = pcbulks%grainc_storage(i) * c4_r2
             pcisos%grainc_xfer(i) = pcbulks%grainc_xfer(i) * c4_r2
             pcisos%dispvegc(i) = pcbulks%dispvegc(i) * c4_r2
             pcisos%storvegc(i) = pcbulks%storvegc(i) * c4_r2
             pcisos%totvegc(i) = pcbulks%totvegc(i) * c4_r2
             pcisos%totpftc(i) = pcbulks%totpftc(i) * c4_r2
             pcisos%woodc(i) = pcbulks%woodc(i) * c4_r2
          end if
       end do
    end if

#if (defined CNDV)
    ! pft type dgvm physical state - crownarea
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CROWNAREA', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CROWNAREA', data=pptr%pdgvs%crownarea, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! tempsum_litfall
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='tempsum_litfall', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='tempsum_litfall', data=pptr%pepv%tempsum_litfall, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! annsum_litfall
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='annsum_litfall', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='annsum_litfall', data=pptr%pepv%annsum_litfall, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! nind
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='nind', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='nind', data=pptr%pdgvs%nind, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! fpcgrid
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='fpcgrid', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='fpcgrid', data=pptr%pdgvs%fpcgrid, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! fpcgridold
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='fpcgridold', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='fpcgridold', data=pptr%pdgvs%fpcgridold, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! gridcell type dgvm physical state - tmomin20
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='TMOMIN20', xtype=ncd_double,  &
            dim1name='gridcell',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='TMOMIN20', data=gptr%gdgvs%tmomin20, &
            dim1name=nameg, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! gridcell type dgvm physical state - agdd20
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='AGDD20', xtype=ncd_double,  &
            dim1name='gridcell',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='AGDD20', data=gptr%gdgvs%agdd20, &
            dim1name=nameg, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! pft type dgvm physical state - t_mo_min
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_MO_MIN', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='T_MO_MIN', data=pptr%pdgvs%t_mo_min, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! present
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='present', xtype=ncd_int,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       allocate (iptemp(begp:endp), stat=ier)
       if (ier /= 0) then
          call endrun('CNrest: allocation error ')
       end if
       if (flag == 'write') then
          do p = begp,endp
             iptemp(p) = 0
             if (pptr%pdgvs%present(p)) iptemp(p) = 1
          end do
       end if
       call ncd_io(varname='present', data=iptemp, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read') then
          if (.not. readvar) then
             if (is_restart()) call endrun
          else
             do p = begp,endp
                pptr%pdgvs%present(p) = .false.
                if (iptemp(p) == 1) pptr%pdgvs%present(p) = .true.
             end do
          end if
       end if
       deallocate (iptemp)
    end if

    ! leafcmax
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafcmax', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='leafcmax', data=pptr%pcs%leafcmax, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! heatstress
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='heatstress', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='heatstress', data=pptr%pdgvs%heatstress, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! greffic
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='greffic', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='greffic', data=pptr%pdgvs%greffic, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if
#endif

  end subroutine CNRest



  subroutine cnrest_addfld_decomp(ncid, varname, longname, units, flag, data_rl, &
                                  fill_value, readvar)

!
! !DESCRIPTION: 
! Read/write CN restart data, for vertical decomp grid that can be set to have length = 1 or nlevgrnd
!
! !USES:
    use clm_varcon, only : nlevgrnd
    use ncdio_pio
    use clmtype, only: namec
    use clm_time_manager, only : is_restart
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t)  :: ncid   ! netcdf id
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: longname
    character(len=*), intent(in) :: units
    character(len=*), intent(in) :: flag   !'read' or 'write'
    real(r8), optional :: fill_value
    real(r8), optional, pointer :: data_rl(:,:)
    logical, optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    real(r8), pointer :: ptr1d(:)

    
#ifdef VERTSOILC

       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname=trim(varname)//'_vr', xtype=ncd_double,  &
               dim1name='column',dim2name='levgrnd', switchdim=.true., &
               long_name=longname,units=units)
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname=trim(varname)//'_vr', data=data_rl, &
               dim1name=namec,switchdim=.true., ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if

#else
       !! nlevdecomp = 1; so treat as 1D variable
       ptr1d => data_rl(:,1)
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname=trim(varname), xtype=ncd_double,  &
               dim1name='column',long_name=longname,units=units)
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname=trim(varname), data=ptr1d, &
               dim1name=namec,ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if
#endif

    
  end subroutine cnrest_addfld_decomp
  
#endif

end module CNrestMod

