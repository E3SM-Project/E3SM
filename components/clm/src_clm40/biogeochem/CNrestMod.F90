module CNrestMod

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
  use abortutils  , only : endrun
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: CNrest
!
! !REVISION HISTORY:
! 11/05/03: Module created by Peter Thornton
!
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
    use clm_atmlnd      , only : clm_a2l
    use clm_varpar      , only : numrad
    use clm_varctl      , only : use_c13, use_cndv, use_exit_spinup
    use decompMod       , only : get_proc_bounds
    use clm_time_manager, only : is_restart
    use ncdio_pio
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
! !LOCAL VARIABLES:
!EOP
    integer :: c,p,j        ! indices 
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
!-----------------------------------------------------------------------

    ! Set pointers into derived type

    gptr => grc
    lptr => lun
    cptr => col
    pptr => pft

    ! Determine necessary subgrid bounds

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    !--------------------------------
    ! pft ecophysiological variables 
    !--------------------------------
    
    ! dormant_flag
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='dormant_flag', xtype=ncd_double,  &
            dim1name='pft',long_name='dormancy flag',units='unitless' )
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='dormant_flag', data=pepv%dormant_flag, &
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
       call ncd_io(varname='days_active', data=pepv%days_active, &
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
       call ncd_io(varname='onset_flag', data=pepv%onset_flag, &
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
       call ncd_io(varname='onset_counter', data=pepv%onset_counter, &
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
       call ncd_io(varname='onset_gddflag', data=pepv%onset_gddflag, &
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
       call ncd_io(varname='onset_fdd', data=pepv%onset_fdd, &
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
       call ncd_io(varname='onset_gdd', data=pepv%onset_gdd, &
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
       call ncd_io(varname='onset_swi', data=pepv%onset_swi, &
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
       call ncd_io(varname='offset_flag', data=pepv%offset_flag, &
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
       call ncd_io(varname='offset_counter', data=pepv%offset_counter, &
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
       call ncd_io(varname='offset_fdd', data=pepv%offset_fdd, &
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
       call ncd_io(varname='offset_swi', data=pepv%offset_swi, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! lgsf
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='lgsf', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='lgsf', data=pepv%lgsf, &
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
       call ncd_io(varname='bglfr', data=pepv%bglfr, &
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
       call ncd_io(varname='bgtr', data=pepv%bgtr, &
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
       call ncd_io(varname='dayl', data=pepv%dayl, &
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
       call ncd_io(varname='prev_dayl', data=pepv%prev_dayl, &
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
       call ncd_io(varname='annavg_t2m', data=pepv%annavg_t2m, &
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
       call ncd_io(varname='tempavg_t2m', data=pepv%tempavg_t2m, &
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
       call ncd_io(varname='gpp_pepv', data=pepv%gpp, &
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
       call ncd_io(varname='availc', data=pepv%availc, &
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
       call ncd_io(varname='xsmrpool_recover', data=pepv%xsmrpool_recover, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if


    if (use_c13) then
       ! xsmrpool_c13ratio
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='xsmrpool_c13ratio', xtype=ncd_double,  &
               dim1name='pft',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='xsmrpool_c13ratio', data=pepv%xsmrpool_c13ratio, &
               dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if
    end if

    ! alloc_pnow
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='alloc_pnow', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='alloc_pnow', data=pepv%alloc_pnow, &
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
       call ncd_io(varname='c_allometry', data=pepv%c_allometry, &
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
       call ncd_io(varname='n_allometry', data=pepv%n_allometry, &
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
       call ncd_io(varname='plant_ndemand', data=pepv%plant_ndemand, &
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
       call ncd_io(varname='tempsum_potential_gpp', data=pepv%tempsum_potential_gpp, &
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
       call ncd_io(varname='annsum_potential_gpp', data=pepv%annsum_potential_gpp, &
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
       call ncd_io(varname='tempmax_retransn', data=pepv%tempmax_retransn, &
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
       call ncd_io(varname='annmax_retransn', data=pepv%annmax_retransn, &
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
       call ncd_io(varname='avail_retransn', data=pepv%avail_retransn, &
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
       call ncd_io(varname='plant_nalloc', data=pepv%plant_nalloc, &
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
       call ncd_io(varname='plant_calloc', data=pepv%plant_calloc, &
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
       call ncd_io(varname='excess_cflux', data=pepv%excess_cflux, &
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
       call ncd_io(varname='downreg', data=pepv%downreg, &
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
       call ncd_io(varname='prev_leafc_to_litter', data=pepv%prev_leafc_to_litter, &
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
       call ncd_io(varname='prev_frootc_to_litter', data=pepv%prev_frootc_to_litter, &
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
       call ncd_io(varname='tempsum_npp', data=pepv%tempsum_npp, &
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
       call ncd_io(varname='annsum_npp', data=pepv%annsum_npp, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    if (use_c13) then
       ! rc13_canair
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='rc13_canair', xtype=ncd_double,  &
               dim1name='pft',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='rc13_canair', data=pepv%rc13_canair, &
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
          call ncd_io(varname='rc13_psnsun', data=pepv%rc13_psnsun, &
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
          call ncd_io(varname='rc13_psnsha', data=pepv%rc13_psnsha, &
               dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if
    end if

    !--------------------------------
    ! pft carbon state variables 
    !--------------------------------

    ! leafc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafc', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='leafc', data=pcs%leafc, &
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
       call ncd_io(varname='leafc_storage', data=pcs%leafc_storage, &
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
       call ncd_io(varname='leafc_xfer', data=pcs%leafc_xfer, &
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
       call ncd_io(varname='frootc', data=pcs%frootc, &
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
       call ncd_io(varname='frootc_storage', data=pcs%frootc_storage, &
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
       call ncd_io(varname='frootc_xfer', data=pcs%frootc_xfer, &
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
       call ncd_io(varname='livestemc', data=pcs%livestemc, &
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
       call ncd_io(varname='livestemc_storage', data=pcs%livestemc_storage, &
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
       call ncd_io(varname='livestemc_xfer', data=pcs%livestemc_xfer, &
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
       call ncd_io(varname='deadstemc', data=pcs%deadstemc, &
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
       call ncd_io(varname='deadstemc_storage', data=pcs%deadstemc_storage, &
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
       call ncd_io(varname='deadstemc_xfer', data=pcs%deadstemc_xfer, &
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
       call ncd_io(varname='livecrootc', data=pcs%livecrootc, &
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
       call ncd_io(varname='livecrootc_storage', data=pcs%livecrootc_storage, &
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
       call ncd_io(varname='livecrootc_xfer', data=pcs%livecrootc_xfer, &
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
       call ncd_io(varname='deadcrootc', data=pcs%deadcrootc, &
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
       call ncd_io(varname='deadcrootc_storage', data=pcs%deadcrootc_storage, &
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
       call ncd_io(varname='deadcrootc_xfer', data=pcs%deadcrootc_xfer, &
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
       call ncd_io(varname='gresp_storage', data=pcs%gresp_storage, &
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
       call ncd_io(varname='gresp_xfer', data=pcs%gresp_xfer, &
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
       call ncd_io(varname='cpool', data=pcs%cpool, &
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
       call ncd_io(varname='xsmrpool', data=pcs%xsmrpool, &
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
       call ncd_io(varname='pft_ctrunc', data=pcs%pft_ctrunc, &
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
       call ncd_io(varname='totvegc', data=pcs%totvegc, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    if (use_c13) then
       !--------------------------------
       ! C13 pft carbon state variables 
       !--------------------------------
       
       ! leafc
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='leafc_13', xtype=ncd_double,  &
               dim1name='pft',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='leafc_13', data=pc13s%leafc, &
               dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if

       ! leafc_storage
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='leafc_storage_13', xtype=ncd_double,  &
               dim1name='pft',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='leafc_storage_13', data=pc13s%leafc_storage, &
               dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if

       ! leafc_xfer
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='leafc_xfer_13', xtype=ncd_double,  &
               dim1name='pft',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='leafc_xfer_13', data=pc13s%leafc_xfer, &
               dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if

       ! frootc
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='frootc_13', xtype=ncd_double,  &
               dim1name='pft',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='frootc_13', data=pc13s%frootc, &
               dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if

       ! frootc_storage
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='frootc_storage_13', xtype=ncd_double,  &
               dim1name='pft',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='frootc_storage_13', data=pc13s%frootc_storage, &
               dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if

       !frootc_xfer 
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='frootc_xfer_13', xtype=ncd_double,  &
               dim1name='pft',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='frootc_xfer_13', data=pc13s%frootc_xfer, &
               dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if

       ! livestemc
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='livestemc_13', xtype=ncd_double,  &
               dim1name='pft',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='livestemc_13', data=pc13s%livestemc, &
               dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if

       ! livestemc_storage
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='livestemc_storage_13', xtype=ncd_double,  &
               dim1name='pft',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='livestemc_storage_13', data=pc13s%livestemc_storage, &
               dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if

       ! livestemc_xfer
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='livestemc_xfer_13', xtype=ncd_double,  &
               dim1name='pft',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='livestemc_xfer_13', data=pc13s%livestemc_xfer, &
               dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if

       ! deadstemc
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='deadstemc_13', xtype=ncd_double,  &
               dim1name='pft',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='deadstemc_13', data=pc13s%deadstemc, &
               dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if

       ! deadstemc_storage
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='deadstemc_storage_13', xtype=ncd_double,  &
               dim1name='pft',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='deadstemc_storage_13', data=pc13s%deadstemc_storage, &
               dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if

       ! deadstemc_xfer
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='deadstemc_xfer_13', xtype=ncd_double,  &
               dim1name='pft',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='deadstemc_xfer_13', data=pc13s%deadstemc_xfer, &
               dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if

       ! livecrootc
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='livecrootc_13', xtype=ncd_double,  &
               dim1name='pft',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='livecrootc_13', data=pc13s%livecrootc, &
               dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if

       ! livecrootc_storage
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='livecrootc_storage_13', xtype=ncd_double,  &
               dim1name='pft',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='livecrootc_storage_13', data=pc13s%livecrootc_storage, &
               dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if

       ! livecrootc_xfer
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='livecrootc_xfer_13', xtype=ncd_double,  &
               dim1name='pft',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='livecrootc_xfer_13', data=pc13s%livecrootc_xfer, &
               dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if

       ! deadcrootc
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='deadcrootc_13', xtype=ncd_double,  &
               dim1name='pft',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='deadcrootc_13', data=pc13s%deadcrootc, &
               dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if

       ! deadcrootc_storage
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='deadcrootc_storage_13', xtype=ncd_double,  &
               dim1name='pft',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='deadcrootc_storage_13', data=pc13s%deadcrootc_storage, &
               dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if

       ! deadcrootc_xfer
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='deadcrootc_xfer_13', xtype=ncd_double,  &
               dim1name='pft',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='deadcrootc_xfer_13', data=pc13s%deadcrootc_xfer, &
               dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if

       ! gresp_storage
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='gresp_storage_13', xtype=ncd_double,  &
               dim1name='pft',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='gresp_storage_13', data=pc13s%gresp_storage, &
               dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if

       ! gresp_xfer
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='gresp_xfer_13', xtype=ncd_double,  &
               dim1name='pft',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='gresp_xfer_13', data=pc13s%gresp_xfer, &
               dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if

       ! cpool
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='cpool_13', xtype=ncd_double,  &
               dim1name='pft',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='cpool_13', data=pc13s%cpool, &
               dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if

       ! xsmrpool
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='xsmrpool_13', xtype=ncd_double,  &
               dim1name='pft',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='xsmrpool_13', data=pc13s%xsmrpool, &
               dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if

       ! pft_ctrunc
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='pft_ctrunc_13', xtype=ncd_double,  &
               dim1name='pft',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='pft_ctrunc_13', data=pc13s%pft_ctrunc, &
               dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if

       ! totvegc
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='totvegc_13', xtype=ncd_double,  &
               dim1name='pft',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='totvegc_13', data=pc13s%totvegc, &
               dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if
    end if

    !--------------------------------
    ! pft nitrogen state variables
    !--------------------------------
    
    ! leafn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafn', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='leafn', data=pns%leafn, &
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
       call ncd_io(varname='leafn_storage', data=pns%leafn_storage, &
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
       call ncd_io(varname='leafn_xfer', data=pns%leafn_xfer, &
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
       call ncd_io(varname='frootn', data=pns%frootn, &
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
       call ncd_io(varname='frootn_storage', data=pns%frootn_storage, &
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
       call ncd_io(varname='frootn_xfer', data=pns%frootn_xfer, &
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
       call ncd_io(varname='livestemn', data=pns%livestemn, &
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
       call ncd_io(varname='livestemn_storage', data=pns%livestemn_storage, &
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
       call ncd_io(varname='livestemn_xfer', data=pns%livestemn_xfer, &
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
       call ncd_io(varname='deadstemn', data=pns%deadstemn, &
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
       call ncd_io(varname='deadstemn_storage', data=pns%deadstemn_storage, &
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
       call ncd_io(varname='deadstemn_xfer', data=pns%deadstemn_xfer, &
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
       call ncd_io(varname='livecrootn', data=pns%livecrootn, &
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
       call ncd_io(varname='livecrootn_storage', data=pns%livecrootn_storage, &
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
       call ncd_io(varname='livecrootn_xfer', data=pns%livecrootn_xfer, &
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
       call ncd_io(varname='deadcrootn', data=pns%deadcrootn, &
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
       call ncd_io(varname='deadcrootn_storage', data=pns%deadcrootn_storage, &
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
       call ncd_io(varname='deadcrootn_xfer', data=pns%deadcrootn_xfer, &
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
       call ncd_io(varname='retransn', data=pns%retransn, &
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
       call ncd_io(varname='npool', data=pns%npool, &
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
       call ncd_io(varname='pft_ntrunc', data=pns%pft_ntrunc, &
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
       call ncd_io(varname='decl', data=cps%decl, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! fpi
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='fpi', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='fpi', data=cps%fpi, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! fpg
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='fpg', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='fpg', data=cps%fpg, &
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
       call ncd_io(varname='annsum_counter', data=cps%annsum_counter, &
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
       call ncd_io(varname='cannsum_npp', data=cps%cannsum_npp, &
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
       call ncd_io(varname='cannavg_t2m', data=cps%cannavg_t2m, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! wf
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='wf', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='wf', data=cps%wf, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! me
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='me', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='me', data=cps%me, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! fire_prob
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='fire_prob', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='fire_prob', data=cps%fire_prob, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! mean_fire_prob
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='mean_fire_prob', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='mean_fire_prob', data=cps%mean_fire_prob, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! fireseasonl
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='fireseasonl', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='fireseasonl', data=cps%fireseasonl, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! farea_burned
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='farea_burned', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='farea_burned', data=cps%farea_burned, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! ann_farea_burned
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='ann_farea_burned', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='ann_farea_burned', data=cps%ann_farea_burned, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    !--------------------------------
    ! column carbon state variables
    !--------------------------------

    ! cwdc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cwdc', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='cwdc', data=ccs%cwdc, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='litr1c', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='litr1c', data=ccs%litr1c, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    !litr2c 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='litr2c', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='litr2c', data=ccs%litr2c, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! litr3c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='litr3c', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='litr3c', data=ccs%litr3c, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    !soil1c 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soil1c', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='soil1c', data=ccs%soil1c, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! soil2c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soil2c', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='soil2c', data=ccs%soil2c, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! soil3c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soil3c', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='soil3c', data=ccs%soil3c, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! soil4c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soil4c', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='soil4c', data=ccs%soil4c, &
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
       call ncd_io(varname='seedc', data=ccs%seedc, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! col_ctrunc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='col_ctrunc', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='col_ctrunc', data=ccs%col_ctrunc, &
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
       call ncd_io(varname='totlitc', data=ccs%totlitc, &
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
       call ncd_io(varname='totcolc', data=ccs%totcolc, &
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
       call ncd_io(varname='prod10c', data=ccs%prod10c, &
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
       call ncd_io(varname='prod100c', data=ccs%prod100c, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    if (use_c13) then
       !--------------------------------
       ! C13 column carbon state variables
       !--------------------------------
       
       ! cwdc
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='cwdc_13', xtype=ncd_double,  &
               dim1name='column',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='cwdc_13', data=cc13s%cwdc, &
               dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if

       ! litr1c
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='litr1c_13', xtype=ncd_double,  &
               dim1name='column',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='litr1c_13', data=cc13s%litr1c, &
               dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if

       !litr2c 
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='litr2c_13', xtype=ncd_double,  &
               dim1name='column',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='litr2c_13', data=cc13s%litr2c, &
               dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if

       ! litr3c
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='litr3c_13', xtype=ncd_double,  &
               dim1name='column',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='litr3c_13', data=cc13s%litr3c, &
               dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if

       !soil1c 
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='soil1c_13', xtype=ncd_double,  &
               dim1name='column',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='soil1c_13', data=cc13s%soil1c, &
               dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if

       ! soil2c
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='soil2c_13', xtype=ncd_double,  &
               dim1name='column',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='soil2c_13', data=cc13s%soil2c, &
               dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if

       ! soil3c
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='soil3c_13', xtype=ncd_double,  &
               dim1name='column',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='soil3c_13', data=cc13s%soil3c, &
               dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if

       ! soil4c
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='soil4c_13', xtype=ncd_double,  &
               dim1name='column',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='soil4c_13', data=cc13s%soil4c, &
               dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if

       ! seedc
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='seedc_13', xtype=ncd_double,  &
               dim1name='column',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='seedc_13', data=cc13s%seedc, &
               dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if

       ! col_ctrunc
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='col_ctrunc_13', xtype=ncd_double,  &
               dim1name='column',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='col_ctrunc_13', data=cc13s%col_ctrunc, &
               dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if

       ! totlitc
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='totlitc_13', xtype=ncd_double,  &
               dim1name='column',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='totlitc_13', data=cc13s%totlitc, &
               dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if

       ! totcolc
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='totcolc_13', xtype=ncd_double,  &
               dim1name='column',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='totcolc_13', data=cc13s%totcolc, &
               dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if

       ! prod10c
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='prod10c_13', xtype=ncd_double,  &
               dim1name='column',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='prod10c_13', data=cc13s%prod10c, &
               dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if

       ! prod100c
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='prod100c_13', xtype=ncd_double,  &
               dim1name='column',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='prod100c_13', data=cc13s%prod100c, &
               dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if
    end if
    
    !--------------------------------
    ! column nitrogen state variables
    !--------------------------------
    
    ! cwdn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cwdn', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='cwdn', data=cns%cwdn, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    !litr1n 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='litr1n', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='litr1n', data=cns%litr1n, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! litr2n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='litr2n', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='litr2n', data=cns%litr2n, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! litr3n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='litr3n', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='litr3n', data=cns%litr3n, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! soil1n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soil1n', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='soil1n', data=cns%soil1n, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! soil2n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soil2n', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='soil2n', data=cns%soil2n, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! soil3n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soil3n', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='soil3n', data=cns%soil3n, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! soil4n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soil4n', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='soil4n', data=cns%soil4n, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! sminn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sminn', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='sminn', data=cns%sminn, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! col_ntrunc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='col_ntrunc', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='col_ntrunc', data=cns%col_ntrunc, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! totcoln
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='totcoln', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='totcoln', data=cns%totcoln, &
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
       call ncd_io(varname='seedn', data=cns%seedn, &
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
       call ncd_io(varname='prod10n', data=cns%prod10n, &
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
       call ncd_io(varname='prod100n', data=cns%prod100n, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    if (use_exit_spinup) then
       if (flag == 'read') then
          m = 20._r8
          do c = begc, endc
             ccs%soil1c(c) = ccs%soil1c(c) * m
             ccs%soil2c(c) = ccs%soil2c(c) * m
             ccs%soil3c(c) = ccs%soil3c(c) * m
             ccs%soil4c(c) = ccs%soil4c(c) * m
             if (use_c13) then
                ! adding code for 13C, 12/25/05, PET 
                cc13s%soil1c(c) = cc13s%soil1c(c) * m
                cc13s%soil2c(c) = cc13s%soil2c(c) * m
                cc13s%soil3c(c) = cc13s%soil3c(c) * m
                cc13s%soil4c(c) = cc13s%soil4c(c) * m
             end if
             cns%soil1n(c) = cns%soil1n(c) * m
             cns%soil2n(c) = cns%soil2n(c) * m
             cns%soil3n(c) = cns%soil3n(c) * m
             cns%soil4n(c) = cns%soil4n(c) * m
          end do
       end if
    end if

    if (use_cndv) then
       ! pft type dgvm physical state - crownarea
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='CROWNAREA', xtype=ncd_double,  &
               dim1name='pft',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='CROWNAREA', data=pdgvs%crownarea, &
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
          call ncd_io(varname='tempsum_litfall', data=pepv%tempsum_litfall, &
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
          call ncd_io(varname='annsum_litfall', data=pepv%annsum_litfall, &
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
          call ncd_io(varname='nind', data=pdgvs%nind, &
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
          call ncd_io(varname='fpcgrid', data=pdgvs%fpcgrid, &
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
          call ncd_io(varname='fpcgridold', data=pdgvs%fpcgridold, &
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
          call ncd_io(varname='TMOMIN20', data=gdgvs%tmomin20, &
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
          call ncd_io(varname='AGDD20', data=gdgvs%agdd20, &
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
          call ncd_io(varname='T_MO_MIN', data=pdgvs%t_mo_min, &
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
                if (pdgvs%present(p)) iptemp(p) = 1
             end do
          end if
          call ncd_io(varname='present', data=iptemp, &
               dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
          if (flag=='read') then
             if (.not. readvar) then
                if (is_restart()) call endrun
             else
                do p = begp,endp
                   pdgvs%present(p) = .false.
                   if (iptemp(p) == 1) pdgvs%present(p) = .true.
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
          call ncd_io(varname='leafcmax', data=pcs%leafcmax, &
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
          call ncd_io(varname='heatstress', data=pdgvs%heatstress, &
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
          call ncd_io(varname='greffic', data=pdgvs%greffic, &
               dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun
          end if
       end if
    end if

  end subroutine CNRest

end module CNrestMod

