!===============================================================================
! Modal Aerosol Model
!===============================================================================
module aero_model
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use constituents,   only: pcnst, cnst_name, cnst_get_ind
  use ppgrid,         only: pcols, pver, pverp
  use abortutils,     only: endrun
  use cam_logfile,    only: iulog
  use perf_mod,       only: t_startf, t_stopf
  use camsrfexch,     only: cam_in_t, cam_out_t
  use aerodep_flx,    only: aerodep_flx_prescribed
  use physics_types,  only: physics_state, physics_ptend, physics_ptend_init
  use physics_buffer, only: physics_buffer_desc
  use physics_buffer, only: pbuf_get_field, pbuf_get_index, pbuf_set_field
  use physconst,      only: gravit, rair, rhoh2o
  use spmd_utils,     only: masterproc

  use cam_history,    only: outfld, fieldname_len
  use chem_mods,      only: gas_pcnst, adv_mass
  use mo_tracname,    only: solsym

  use modal_aero_data,only: cnst_name_cw
  use modal_aero_data,only: ntot_amode, modename_amode

  implicit none
  private

  public :: aero_model_readnl
  public :: aero_model_register
  public :: aero_model_init
  public :: aero_model_gasaerexch ! create, grow, change, and shrink aerosols.
  public :: aero_model_drydep     ! aerosol dry deposition and sediment
  public :: aero_model_wetdep     ! aerosol wet removal
  public :: aero_model_emissions  ! aerosol emissions
  public :: aero_model_surfarea   ! aerosol surface area for chemistry

 ! Misc private data 

  ! number of modes
  integer :: nmodes
  integer :: pblh_idx            = 0
  integer :: dgnum_idx           = 0
  integer :: dgnumwet_idx        = 0
  integer :: rate1_cw2pr_st_idx  = 0  

  integer :: wetdens_ap_idx      = 0
  integer :: qaerwat_idx         = 0

  integer :: fracis_idx          = 0
  integer :: prain_idx           = 0
  integer :: nevapr_idx          = 0
  integer :: rprddp_idx          = 0 
  integer :: rprdsh_idx          = 0 

  ! variables for table lookup of aerosol impaction/interception scavenging rates
  integer, parameter :: nimptblgrow_mind=-7, nimptblgrow_maxd=12
  real(r8) :: dlndg_nimptblgrow
  real(r8) :: scavimptblnum(nimptblgrow_mind:nimptblgrow_maxd, ntot_amode)
  real(r8) :: scavimptblvol(nimptblgrow_mind:nimptblgrow_maxd, ntot_amode)

  ! for aero_model_surfarea called from mo_usrrxt
  integer :: aitken_idx = -1
  integer, dimension(ntot_amode) :: num_idx = -1
  integer :: index_tot_mass(ntot_amode,10) = -1
  integer :: index_chm_mass(ntot_amode,10) = -1

  integer :: ndx_h2so4
  character(len=fieldname_len) :: dgnum_name(ntot_amode)

  ! Namelist variables
  character(len=16) :: wetdep_list(pcnst) = ' '
  character(len=16) :: drydep_list(pcnst) = ' '
  real(r8)          :: sol_facti_cloud_borne = 1._r8

  integer :: ndrydep = 0
  integer,allocatable :: drydep_indices(:)
  integer :: nwetdep = 0
  integer,allocatable :: wetdep_indices(:)
  logical :: drydep_lq(pcnst)
  logical :: wetdep_lq(pcnst)

contains
  
  !=============================================================================
  ! reads aerosol namelist options
  !=============================================================================
  subroutine aero_model_readnl(nlfile)

    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'aero_model_readnl'

    ! Namelist variables
    character(len=16) :: aer_wetdep_list(pcnst) = ' '
    character(len=16) :: aer_drydep_list(pcnst) = ' '

    namelist /aerosol_nl/ aer_wetdep_list, aer_drydep_list, sol_facti_cloud_borne

    !-----------------------------------------------------------------------------

    ! Read namelist
    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'aerosol_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, aerosol_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    ! Broadcast namelist variables
    call mpibcast(aer_wetdep_list,   len(aer_wetdep_list(1))*pcnst, mpichar, 0, mpicom)
    call mpibcast(aer_drydep_list,   len(aer_drydep_list(1))*pcnst, mpichar, 0, mpicom)
    call mpibcast(sol_facti_cloud_borne, 1,                         mpir8,   0, mpicom)
#endif

    wetdep_list = aer_wetdep_list
    drydep_list = aer_drydep_list

  end subroutine aero_model_readnl

  !=============================================================================
  !=============================================================================
  subroutine aero_model_register
    use modal_aero_initialize_data, only : modal_aero_register

    call modal_aero_register()

  end subroutine aero_model_register

  !=============================================================================
  !=============================================================================
  subroutine aero_model_init( pbuf2d )

    use mo_chem_utls,    only: get_inv_ndx
    use cam_history,     only: addfld, add_default, phys_decomp
    use phys_control,    only: phys_getopts
    use mo_chem_utls,    only: get_rxt_ndx, get_spc_ndx
    use modal_aero_data, only: cnst_name_cw
    use modal_aero_initialize_data, only: modal_aero_initialize
    use rad_constituents,           only: rad_cnst_get_info
    use dust_model,      only: dust_init, dust_names, dust_active, dust_nbin, dust_nnum
    use seasalt_model,   only: seasalt_init, seasalt_names, seasalt_active,seasalt_nbin
    use drydep_mod,      only: inidrydep
    use wetdep,          only: wetdep_init

    ! args
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    ! local vars
    character(len=*), parameter :: subrname = 'aero_model_init'
    integer :: m, n, id
    character(len=20) :: dummy

    logical  :: history_aerosol ! Output MAM or SECT aerosol tendencies

    integer :: l
    character(len=6) :: test_name
    character(len=64) :: errmes

    character(len=2)  :: unit_basename  ! Units 'kg' or '1' 

    dgnum_idx      = pbuf_get_index('DGNUM')
    dgnumwet_idx   = pbuf_get_index('DGNUMWET')
    
    call phys_getopts( history_aerosol_out=history_aerosol )

    call rad_cnst_get_info(0, nmodes=nmodes)

    call modal_aero_initialize(pbuf2d)
    call modal_aero_bcscavcoef_init()

    call dust_init()
    call seasalt_init()
    call wetdep_init()

    fracis_idx      = pbuf_get_index('FRACIS') 
    prain_idx       = pbuf_get_index('PRAIN')  
    nevapr_idx      = pbuf_get_index('NEVAPR') 
    rprddp_idx      = pbuf_get_index('RPRDDP')  
    rprdsh_idx      = pbuf_get_index('RPRDSH')  

    nwetdep = 0
    ndrydep = 0

    count_species: do m = 1,pcnst
       if ( len_trim(wetdep_list(m)) /= 0 ) then
          nwetdep = nwetdep+1
       endif
       if ( len_trim(drydep_list(m)) /= 0 ) then
          ndrydep = ndrydep+1
       endif
    enddo count_species
    
    if (nwetdep>0) &
         allocate(wetdep_indices(nwetdep))
    if (ndrydep>0) &
         allocate(drydep_indices(ndrydep))

    do m = 1,ndrydep
       call cnst_get_ind ( drydep_list(m), id, abort=.false. )
       if (id>0) then
          drydep_indices(m) = id
       else
          call endrun(subrname//': invalid drydep species: '//trim(drydep_list(m)) )
       endif

       if (masterproc) then
          write(iulog,*) subrname//': '//drydep_list(m)//' will have drydep applied'
       endif
    enddo
    do m = 1,nwetdep
       call cnst_get_ind ( wetdep_list(m), id, abort=.false. )
       if (id>0) then
          wetdep_indices(m) = id
       else
          call endrun(subrname//': invalid wetdep species: '//trim(wetdep_list(m)) )
       endif
       
       if (masterproc) then
          write(iulog,*) subrname//': '//wetdep_list(m)//' will have wet removal'
       endif
    enddo

    if (ndrydep>0) then

       call inidrydep(rair, gravit)

       dummy = 'RAM1'
       call addfld (dummy,'frac ',1, 'A','RAM1',phys_decomp)
       if ( history_aerosol ) then  
          call add_default (dummy, 1, ' ')
       endif
       dummy = 'airFV'
       call addfld (dummy,'frac ',1, 'A','FV',phys_decomp)
       if ( history_aerosol ) then  
          call add_default (dummy, 1, ' ')
       endif

    endif

    if (dust_active) then
       ! emissions diagnostics ....

       do m = 1, dust_nbin+dust_nnum
          dummy = trim(dust_names(m)) // 'SF'
          call addfld (dummy,'kg/m2/s ',1, 'A',trim(dust_names(m))//' dust surface emission',phys_decomp)
          if (history_aerosol) then
             call add_default (dummy, 1, ' ')
          endif
       enddo

       dummy = 'DSTSFMBL'
       call addfld (dummy,'kg/m2/s',1, 'A','Mobilization flux at surface',phys_decomp)
       if (history_aerosol) then
          call add_default (dummy, 1, ' ')
       endif

       dummy = 'LND_MBL'
       call addfld (dummy,'frac ',1, 'A','Soil erodibility factor',phys_decomp)
       if (history_aerosol) then
          call add_default (dummy, 1, ' ')
       endif

    endif

    if (seasalt_active) then
       
       dummy = 'SSTSFMBL'
       call addfld (dummy,'kg/m2/s',1, 'A','Mobilization flux at surface',phys_decomp)
       if (history_aerosol) then
          call add_default (dummy, 1, ' ')
       endif

       do m = 1, seasalt_nbin
          dummy = trim(seasalt_names(m)) // 'SF'
          call addfld (dummy,'kg/m2/s ',1, 'A',trim(seasalt_names(m))//' seasalt surface emission',phys_decomp)
          if (history_aerosol) then
             call add_default (dummy, 1, ' ')
          endif
       enddo

    endif

    
    ! set flags for drydep tendencies
    drydep_lq(:) = .false.
    do m=1,ndrydep 
       id = drydep_indices(m)
       drydep_lq(id) =  .true.
    enddo

    ! set flags for wetdep tendencies
    wetdep_lq(:) = .false.
    do m=1,nwetdep
       id = wetdep_indices(m)
       wetdep_lq(id) = .true.
    enddo

    wetdens_ap_idx = pbuf_get_index('WETDENS_AP')
    qaerwat_idx    = pbuf_get_index('QAERWAT')
    pblh_idx       = pbuf_get_index('pblh')

    rate1_cw2pr_st_idx  = pbuf_get_index('RATE1_CW2PR_ST') 
    call pbuf_set_field(pbuf2d, rate1_cw2pr_st_idx, 0.0_r8)

    do m = 1,ndrydep
       
       ! units 
       if (drydep_list(m)(1:3) == 'num') then
          unit_basename = ' 1'
       else
          unit_basename = 'kg'  
       endif

       call addfld (trim(drydep_list(m))//'DDF',unit_basename//'/m2/s ',   1, 'A', &
            trim(drydep_list(m))//' dry deposition flux at bottom (grav + turb)',phys_decomp)
       call addfld (trim(drydep_list(m))//'TBF',unit_basename//'/m2/s',   1, 'A', &
            trim(drydep_list(m))//' turbulent dry deposition flux',phys_decomp)
       call addfld (trim(drydep_list(m))//'GVF',unit_basename//'/m2/s ',   1, 'A', &
            trim(drydep_list(m))//' gravitational dry deposition flux',phys_decomp)
       call addfld (trim(drydep_list(m))//'DTQ',unit_basename//'/kg/s ',pver, 'A', &
            trim(drydep_list(m))//' dry deposition',phys_decomp)
       call addfld (trim(drydep_list(m))//'DDV','m/s     ',pver, 'A', &
            trim(drydep_list(m))//' deposition velocity',phys_decomp)

       if ( history_aerosol ) then 
          call add_default (trim(drydep_list(m))//'DDF', 1, ' ')
          call add_default (trim(drydep_list(m))//'TBF', 1, ' ')
          call add_default (trim(drydep_list(m))//'GVF', 1, ' ')
       endif

    enddo

    do m = 1,nwetdep
       
       ! units 
       if (wetdep_list(m)(1:3) == 'num') then
          unit_basename = ' 1'
       else
          unit_basename = 'kg'  
       endif

       call addfld (trim(wetdep_list(m))//'SFWET',unit_basename//'/m2/s ', &
            1,  'A','Wet deposition flux at surface',phys_decomp)
       call addfld (trim(wetdep_list(m))//'SFSIC',unit_basename//'/m2/s ', &
            1,  'A','Wet deposition flux (incloud, convective) at surface',phys_decomp)
       call addfld (trim(wetdep_list(m))//'SFSIS',unit_basename//'/m2/s ', &
            1,  'A','Wet deposition flux (incloud, stratiform) at surface',phys_decomp)
       call addfld (trim(wetdep_list(m))//'SFSBC',unit_basename//'/m2/s ', &
            1,  'A','Wet deposition flux (belowcloud, convective) at surface',phys_decomp)
       call addfld (trim(wetdep_list(m))//'SFSBS',unit_basename//'/m2/s ', &
            1,  'A','Wet deposition flux (belowcloud, stratiform) at surface',phys_decomp)
       call addfld (trim(wetdep_list(m))//'WET',unit_basename//'/kg/s ',pver, 'A','wet deposition tendency',phys_decomp)
       call addfld (trim(wetdep_list(m))//'SIC',unit_basename//'/kg/s ',pver, 'A', &
            trim(wetdep_list(m))//' ic wet deposition',phys_decomp)
       call addfld (trim(wetdep_list(m))//'SIS',unit_basename//'/kg/s ',pver, 'A', &
            trim(wetdep_list(m))//' is wet deposition',phys_decomp)
       call addfld (trim(wetdep_list(m))//'SBC',unit_basename//'/kg/s ',pver, 'A', &
            trim(wetdep_list(m))//' bc wet deposition',phys_decomp)
       call addfld (trim(wetdep_list(m))//'SBS',unit_basename//'/kg/s ',pver, 'A', &
            trim(wetdep_list(m))//' bs wet deposition',phys_decomp)
       
       if ( history_aerosol ) then          
          call add_default (trim(wetdep_list(m))//'SFWET', 1, ' ')
          call add_default (trim(wetdep_list(m))//'SFSIC', 1, ' ')
          call add_default (trim(wetdep_list(m))//'SFSIS', 1, ' ')
          call add_default (trim(wetdep_list(m))//'SFSBC', 1, ' ')
          call add_default (trim(wetdep_list(m))//'SFSBS', 1, ' ')
       endif

    enddo

    do m = 1,gas_pcnst

       if  ( solsym(m)(1:3) == 'num') then
          unit_basename = ' 1'  ! Units 'kg' or '1' 
       else
          unit_basename = 'kg'  ! Units 'kg' or '1' 
       end if

       call addfld( 'GS_'//trim(solsym(m)), unit_basename//'/m2/s ',1,  'A', &
                    trim(solsym(m))//' gas chemistry/wet removal (for gas species)', phys_decomp)
       call addfld( 'AQ_'//trim(solsym(m)), unit_basename//'/m2/s ',1,  'A', &
                    trim(solsym(m))//' aqueous chemistry (for gas species)', phys_decomp)
       if ( history_aerosol ) then 
          call add_default( 'GS_'//trim(solsym(m)), 1, ' ')
          call add_default( 'AQ_'//trim(solsym(m)), 1, ' ')
       endif
    enddo
    do n = 1,pcnst
       if( .not. (cnst_name_cw(n) == ' ') ) then

          if (cnst_name_cw(n)(1:3) == 'num') then
             unit_basename = ' 1'
          else
             unit_basename = 'kg'  
          endif

          call addfld( cnst_name_cw(n),                unit_basename//'/kg ', pver, 'A', &
               trim(cnst_name_cw(n))//' in cloud water',phys_decomp)
          call addfld (trim(cnst_name_cw(n))//'SFWET', unit_basename//'/m2/s ',1,  'A', &
               trim(cnst_name_cw(n))//' wet deposition flux at surface',phys_decomp)
          call addfld (trim(cnst_name_cw(n))//'SFSIC', unit_basename//'/m2/s ',1,  'A', &
               trim(cnst_name_cw(n))//' wet deposition flux (incloud, convective) at surface',phys_decomp)
          call addfld (trim(cnst_name_cw(n))//'SFSIS', unit_basename//'/m2/s ',1,  'A', &
               trim(cnst_name_cw(n))//' wet deposition flux (incloud, stratiform) at surface',phys_decomp)
          call addfld (trim(cnst_name_cw(n))//'SFSBC', unit_basename//'/m2/s ',1,  'A', &
               trim(cnst_name_cw(n))//' wet deposition flux (belowcloud, convective) at surface',phys_decomp)
          call addfld (trim(cnst_name_cw(n))//'SFSBS', unit_basename//'/m2/s ',1,  'A', &
               trim(cnst_name_cw(n))//' wet deposition flux (belowcloud, stratiform) at surface',phys_decomp)
          call addfld (trim(cnst_name_cw(n))//'DDF',   unit_basename//'/m2/s ',   1, 'A', &
               trim(cnst_name_cw(n))//' dry deposition flux at bottom (grav + turb)',phys_decomp)
          call addfld (trim(cnst_name_cw(n))//'TBF',   unit_basename//'/m2/s ',   1, 'A', &
               trim(cnst_name_cw(n))//' turbulent dry deposition flux',phys_decomp)
          call addfld (trim(cnst_name_cw(n))//'GVF',   unit_basename//'/m2/s ',   1, 'A', &
               trim(cnst_name_cw(n))//' gravitational dry deposition flux',phys_decomp)     

          if ( history_aerosol ) then 
             call add_default( cnst_name_cw(n), 1, ' ' )
             call add_default (trim(cnst_name_cw(n))//'GVF', 1, ' ')
             call add_default (trim(cnst_name_cw(n))//'SFWET', 1, ' ') 
             call add_default (trim(cnst_name_cw(n))//'TBF', 1, ' ')
             call add_default (trim(cnst_name_cw(n))//'DDF', 1, ' ')
             call add_default (trim(cnst_name_cw(n))//'SFSBS', 1, ' ')      
             call add_default (trim(cnst_name_cw(n))//'SFSIC', 1, ' ')
             call add_default (trim(cnst_name_cw(n))//'SFSBC', 1, ' ')
             call add_default (trim(cnst_name_cw(n))//'SFSIS', 1, ' ')
          endif
       endif
    enddo
    do n=1,ntot_amode
       dgnum_name(n) = ' '
       write(dgnum_name(n),fmt='(a,i1)') 'dgnumwet',n
       call addfld( dgnum_name(n), 'm', pver, 'I', 'Aerosol mode wet diameter', phys_decomp )
       if ( history_aerosol ) then 
          call add_default( dgnum_name(n), 1, ' ' )
       endif
    end do

    ndx_h2so4 = get_spc_ndx('H2SO4')

    ! for aero_model_surfarea called from mo_usrrxt
    do l=1,ntot_amode
       if ( trim(modename_amode(l)) == 'aitken' ) then
          aitken_idx = l
       end if
       test_name = ' '
       write(test_name,fmt='(a5,i1)') 'num_a',l
       num_idx(l) = get_spc_ndx( trim(test_name) )
       if (num_idx(l) < 0) then
          write(errmes,fmt='(a,i1)') 'usrrxt_inti: cannot find MAM num_idx ',l
          write(iulog,*) errmes
          call endrun(errmes)
       endif
    end do
    dgnumwet_idx = pbuf_get_index('DGNUMWET')
    if ( aitken_idx < 0 ) then
       errmes = 'usrrxt_inti: cannot find aitken_idx'
       call endrun(errmes)
    end if

    !
    ! define indeces associated with the various names (defined in
    ! chemistry/modal_aero/modal_aero_initialize_data.F90)
    !
#if ( defined MODAL_AERO_3MODE )
    !
    ! accumulation mode #1
    !
    index_tot_mass(1,1) = get_spc_ndx('so4_a1')
    index_tot_mass(1,2) = get_spc_ndx('pom_a1')
    index_tot_mass(1,3) = get_spc_ndx('soa_a1')
    index_tot_mass(1,4) = get_spc_ndx('bc_a1' )
    index_tot_mass(1,5) = get_spc_ndx('dst_a1')
    index_tot_mass(1,6) = get_spc_ndx('ncl_a1')
    index_chm_mass(1,1) = get_spc_ndx('so4_a1')
    index_chm_mass(1,2) = get_spc_ndx('soa_a1')
    index_chm_mass(1,3) = get_spc_ndx('bc_a1' )
    !
    ! aitken mode
    !
    index_tot_mass(2,1) = get_spc_ndx('so4_a2')
    index_tot_mass(2,2) = get_spc_ndx('soa_a2')
    index_tot_mass(2,3) = get_spc_ndx('ncl_a2')
    index_chm_mass(2,1) = get_spc_ndx('so4_a2')
    index_chm_mass(2,2) = get_spc_ndx('soa_a2')
    !
    ! coarse mode
    !
    index_tot_mass(3,1) = get_spc_ndx('dst_a3')
    index_tot_mass(3,2) = get_spc_ndx('ncl_a3')
    index_tot_mass(3,3) = get_spc_ndx('so4_a3')
    index_chm_mass(3,1) = get_spc_ndx('so4_a3')
    !
#endif
#if ( defined MODAL_AERO_7MODE )
    !
    ! accumulation mode #1
    !
    index_tot_mass(1,1) = get_spc_ndx('so4_a1')
    index_tot_mass(1,2) = get_spc_ndx('nh4_a1')
    index_tot_mass(1,3) = get_spc_ndx('pom_a1')
    index_tot_mass(1,4) = get_spc_ndx('soa_a1')
    index_tot_mass(1,5) = get_spc_ndx('bc_a1' )
    index_tot_mass(1,6) = get_spc_ndx('ncl_a1')
    index_chm_mass(1,1) = get_spc_ndx('so4_a1')
    index_chm_mass(1,2) = get_spc_ndx('nh4_a1')
    index_chm_mass(1,3) = get_spc_ndx('soa_a1')
    index_chm_mass(1,4) = get_spc_ndx('bc_a1' )
    !
    ! aitken mode
    !
    index_tot_mass(2,1) = get_spc_ndx('so4_a2')
    index_tot_mass(2,2) = get_spc_ndx('nh4_a2')
    index_tot_mass(2,3) = get_spc_ndx('soa_a2')
    index_tot_mass(2,4) = get_spc_ndx('ncl_a2')
    index_chm_mass(2,1) = get_spc_ndx('so4_a2')
    index_chm_mass(2,2) = get_spc_ndx('nh4_a2')
    index_chm_mass(2,3) = get_spc_ndx('soa_a2')
    !
    ! primary carbon mode not added 
    !
    ! fine sea salt 
    !
    index_tot_mass(4,1) = get_spc_ndx('so4_a4')
    index_tot_mass(4,2) = get_spc_ndx('nh4_a4')
    index_tot_mass(4,3) = get_spc_ndx('ncl_a4')
    index_chm_mass(4,1) = get_spc_ndx('so4_a4')
    index_chm_mass(4,2) = get_spc_ndx('nh4_a4')
    !
    ! fine soil dust 
    !
    index_tot_mass(5,1) = get_spc_ndx('so4_a5')
    index_tot_mass(5,2) = get_spc_ndx('nh4_a5')
    index_tot_mass(5,3) = get_spc_ndx('dst_a5')
    index_chm_mass(5,1) = get_spc_ndx('so4_a5')
    index_chm_mass(5,2) = get_spc_ndx('nh4_a5')
    !
    ! coarse sea salt 
    !
    index_tot_mass(6,1) = get_spc_ndx('so4_a6')
    index_tot_mass(6,2) = get_spc_ndx('nh4_a6')
    index_tot_mass(6,3) = get_spc_ndx('ncl_a6')
    index_chm_mass(6,1) = get_spc_ndx('so4_a6')
    index_chm_mass(6,2) = get_spc_ndx('nh4_a6')
    !
    ! coarse soil dust 
    !
    index_tot_mass(7,1) = get_spc_ndx('so4_a7')
    index_tot_mass(7,2) = get_spc_ndx('nh4_a7')
    index_tot_mass(7,3) = get_spc_ndx('dst_a7')
    index_chm_mass(7,1) = get_spc_ndx('so4_a7')
    index_chm_mass(7,2) = get_spc_ndx('nh4_a7')
    !
#endif

  end subroutine aero_model_init

  !=============================================================================
  !=============================================================================
  subroutine aero_model_drydep  ( state, pbuf, obklen, ustar, cam_in, dt, cam_out, ptend )

    use dust_sediment_mod, only: dust_sediment_tend
    use drydep_mod,        only: d3ddflux, calcram
    use modal_aero_data,   only: qqcw_get_field
    use modal_aero_data,   only: cnst_name_cw
    use modal_aero_data,   only: alnsg_amode
    use modal_aero_data,   only: sigmag_amode
    use modal_aero_data,   only: nspec_amode
    use modal_aero_data,   only: numptr_amode
    use modal_aero_data,   only: numptrcw_amode
    use modal_aero_data,   only: lmassptr_amode
    use modal_aero_data,   only: lmassptrcw_amode
    use modal_aero_deposition, only: set_srf_drydep

  ! args 
    type(physics_state),    intent(in)    :: state     ! Physics state variables
    real(r8),               intent(in)    :: obklen(:)          
    real(r8),               intent(in)    :: ustar(:)  ! sfc fric vel
    type(cam_in_t), target, intent(in)    :: cam_in    ! import state
    real(r8),               intent(in)    :: dt             ! time step
    type(cam_out_t),        intent(inout) :: cam_out   ! export state
    type(physics_ptend),    intent(out)   :: ptend     ! indivdual parameterization tendencies
    type(physics_buffer_desc),    pointer :: pbuf(:)

  ! local vars
    real(r8), pointer :: landfrac(:) ! land fraction
    real(r8), pointer :: icefrac(:)  ! ice fraction
    real(r8), pointer :: ocnfrac(:)  ! ocean fraction
    real(r8), pointer :: fvin(:)     !
    real(r8), pointer :: ram1in(:)   ! for dry dep velocities from land model for progseasalts

    real(r8) :: fv(pcols)            ! for dry dep velocities, from land modified over ocean & ice
    real(r8) :: ram1(pcols)          ! for dry dep velocities, from land modified over ocean & ice

    integer :: lchnk                   ! chunk identifier
    integer :: ncol                    ! number of atmospheric columns
    integer :: jvlc                    ! index for last dimension of vlc_xxx arrays
    integer :: lphase                  ! index for interstitial / cloudborne aerosol
    integer :: lspec                   ! index for aerosol number / chem-mass / water-mass
    integer :: m                       ! aerosol mode index
    integer :: mm                      ! tracer index
    integer :: i

    real(r8) :: tvs(pcols,pver)
    real(r8) :: rho(pcols,pver)      ! air density in kg/m3
    real(r8) :: sflx(pcols)          ! deposition flux
    real(r8) :: dep_trb(pcols)       !kg/m2/s
    real(r8) :: dep_grv(pcols)       !kg/m2/s (total of grav and trb)
    real(r8) :: pvmzaer(pcols,pverp) ! sedimentation velocity in Pa
    real(r8) :: dqdt_tmp(pcols,pver) ! temporary array to hold tendency for 1 species

    real(r8) :: rad_drop(pcols,pver)
    real(r8) :: dens_drop(pcols,pver)
    real(r8) :: sg_drop(pcols,pver)
    real(r8) :: rad_aer(pcols,pver)
    real(r8) :: dens_aer(pcols,pver)
    real(r8) :: sg_aer(pcols,pver)

    real(r8) :: vlc_dry(pcols,pver,4)     ! dep velocity
    real(r8) :: vlc_grv(pcols,pver,4)     ! dep velocity
    real(r8)::  vlc_trb(pcols,4)          ! dep velocity
    real(r8) :: aerdepdryis(pcols,pcnst)  ! aerosol dry deposition (interstitial)
    real(r8) :: aerdepdrycw(pcols,pcnst)  ! aerosol dry deposition (cloud water)
    real(r8), pointer :: fldcw(:,:)
    real(r8), pointer :: dgncur_awet(:,:,:)
    real(r8), pointer :: wetdens(:,:,:)
    real(r8), pointer :: qaerwat(:,:,:)

    landfrac => cam_in%landfrac(:)
    icefrac  => cam_in%icefrac(:)
    ocnfrac  => cam_in%ocnfrac(:)
    fvin     => cam_in%fv(:)
    ram1in   => cam_in%ram1(:)

    lchnk = state%lchnk
    ncol  = state%ncol

    ! calc ram and fv over ocean and sea ice ...
    call calcram( ncol,landfrac,icefrac,ocnfrac,obklen,&
                  ustar,ram1in,ram1,state%t(:,pver),state%pmid(:,pver),&
                  state%pdel(:,pver),fvin,fv)

    call outfld( 'airFV', fv(:), pcols, lchnk )
    call outfld( 'RAM1', ram1(:), pcols, lchnk )
 
    ! note that tendencies are not only in sfc layer (because of sedimentation)
    ! and that ptend is updated within each subroutine for different species
    
    call physics_ptend_init(ptend, state%psetcols, 'aero_model_drydep', lq=drydep_lq)

    call pbuf_get_field(pbuf, dgnumwet_idx,   dgncur_awet, start=(/1,1,1/), kount=(/pcols,pver,nmodes/) ) 
    call pbuf_get_field(pbuf, wetdens_ap_idx, wetdens,     start=(/1,1,1/), kount=(/pcols,pver,nmodes/) ) 
    call pbuf_get_field(pbuf, qaerwat_idx,    qaerwat,     start=(/1,1,1/), kount=(/pcols,pver,nmodes/) ) 

    tvs(:ncol,:) = state%t(:ncol,:)!*(1+state%q(:ncol,k)
    rho(:ncol,:)=  state%pmid(:ncol,:)/(rair*state%t(:ncol,:))

!
! calc settling/deposition velocities for cloud droplets (and cloud-borne aerosols)
!
! *** mean drop radius should eventually be computed from ndrop and qcldwtr
    rad_drop(:,:) = 5.0e-6_r8
    dens_drop(:,:) = rhoh2o
    sg_drop(:,:) = 1.46_r8
    jvlc = 3
    call modal_aero_depvel_part( ncol,state%t(:,:), state%pmid(:,:), ram1, fv,  &
                     vlc_dry(:,:,jvlc), vlc_trb(:,jvlc), vlc_grv(:,:,jvlc),  &
                     rad_drop(:,:), dens_drop(:,:), sg_drop(:,:), 0, lchnk)
    jvlc = 4
    call modal_aero_depvel_part( ncol,state%t(:,:), state%pmid(:,:), ram1, fv,  &
                     vlc_dry(:,:,jvlc), vlc_trb(:,jvlc), vlc_grv(:,:,jvlc),  &
                     rad_drop(:,:), dens_drop(:,:), sg_drop(:,:), 3, lchnk)



    do m = 1, ntot_amode   ! main loop over aerosol modes

       do lphase = 1, 2   ! loop over interstitial / cloud-borne forms

          if (lphase == 1) then   ! interstial aerosol - calc settling/dep velocities of mode

! rad_aer = volume mean wet radius (m)
! dgncur_awet = geometric mean wet diameter for number distribution (m)
             rad_aer(1:ncol,:) = 0.5_r8*dgncur_awet(1:ncol,:,m)   &
                                 *exp(1.5_r8*(alnsg_amode(m)**2))
! dens_aer(1:ncol,:) = wet density (kg/m3)
             dens_aer(1:ncol,:) = wetdens(1:ncol,:,m)
             sg_aer(1:ncol,:) = sigmag_amode(m)

             jvlc = 1
             call modal_aero_depvel_part( ncol, state%t(:,:), state%pmid(:,:), ram1, fv,  & 
                        vlc_dry(:,:,jvlc), vlc_trb(:,jvlc), vlc_grv(:,:,jvlc),  &
                        rad_aer(:,:), dens_aer(:,:), sg_aer(:,:), 0, lchnk)
             jvlc = 2
             call modal_aero_depvel_part( ncol, state%t(:,:), state%pmid(:,:), ram1, fv,  & 
                        vlc_dry(:,:,jvlc), vlc_trb(:,jvlc), vlc_grv(:,:,jvlc),  &
                        rad_aer(:,:), dens_aer(:,:), sg_aer(:,:), 3, lchnk)
          end if

          do lspec = 0, nspec_amode(m)+1   ! loop over number + constituents + water

             if (lspec == 0) then   ! number
                if (lphase == 1) then
                   mm = numptr_amode(m)
                   jvlc = 1
                else
                   mm = numptrcw_amode(m)
                   jvlc = 3
                endif
             else if (lspec <= nspec_amode(m)) then   ! non-water mass
                if (lphase == 1) then
                   mm = lmassptr_amode(lspec,m)
                   jvlc = 2
                else
                   mm = lmassptrcw_amode(lspec,m)
                   jvlc = 4
                endif
             else   ! water mass
!   bypass dry deposition of aerosol water
                cycle
                if (lphase == 1) then
                   mm = 0
!                  mm = lwaterptr_amode(m)
                   jvlc = 2
                else
                   mm = 0
                   jvlc = 4
                endif
             endif


          if (mm <= 0) cycle

!         if (lphase == 1) then
          if ((lphase == 1) .and. (lspec <= nspec_amode(m))) then
             ptend%lq(mm) = .TRUE.

             ! use pvprogseasalts instead (means making the top level 0)
             pvmzaer(:ncol,1)=0._r8
             pvmzaer(:ncol,2:pverp) = vlc_dry(:ncol,:,jvlc)

             call outfld( trim(cnst_name(mm))//'DDV', pvmzaer(:,2:pverp), pcols, lchnk )

             if(.true.) then ! use phil's method
             !      convert from meters/sec to pascals/sec
             !      pvprogseasalts(:,1) is assumed zero, use density from layer above in conversion
                pvmzaer(:ncol,2:pverp) = pvmzaer(:ncol,2:pverp) * rho(:ncol,:)*gravit

             !      calculate the tendencies and sfc fluxes from the above velocities
                call dust_sediment_tend( &
                     ncol,             dt,       state%pint(:,:), state%pmid, state%pdel, state%t , &
                     state%q(:,:,mm),  pvmzaer,  ptend%q(:,:,mm), sflx  )
             else   !use charlie's method
                call d3ddflux( ncol, vlc_dry(:,:,jvlc), state%q(:,:,mm), state%pmid, &
                               state%pdel, tvs, sflx, ptend%q(:,:,mm), dt )
             endif

             ! apportion dry deposition into turb and gravitational settling for tapes
             do i=1,ncol
                dep_trb(i)=sflx(i)*vlc_trb(i,jvlc)/vlc_dry(i,pver,jvlc)
                dep_grv(i)=sflx(i)*vlc_grv(i,pver,jvlc)/vlc_dry(i,pver,jvlc)
             enddo

             call outfld( trim(cnst_name(mm))//'DDF', sflx, pcols, lchnk)
             call outfld( trim(cnst_name(mm))//'TBF', dep_trb, pcols, lchnk )
             call outfld( trim(cnst_name(mm))//'GVF', dep_grv, pcols, lchnk )
             call outfld( trim(cnst_name(mm))//'DTQ', ptend%q(:,:,mm), pcols, lchnk)
             aerdepdryis(:ncol,mm) = sflx(:ncol)

          else if ((lphase == 1) .and. (lspec == nspec_amode(m)+1)) then  ! aerosol water
             ! use pvprogseasalts instead (means making the top level 0)
             pvmzaer(:ncol,1)=0._r8
             pvmzaer(:ncol,2:pverp) = vlc_dry(:ncol,:,jvlc)

             if(.true.) then ! use phil's method
             !      convert from meters/sec to pascals/sec
             !      pvprogseasalts(:,1) is assumed zero, use density from layer above in conversion
                pvmzaer(:ncol,2:pverp) = pvmzaer(:ncol,2:pverp) * rho(:ncol,:)*gravit

             !      calculate the tendencies and sfc fluxes from the above velocities
                call dust_sediment_tend( &
                     ncol,             dt,       state%pint(:,:), state%pmid, state%pdel, state%t , &
                     qaerwat(:,:,mm),  pvmzaer,  dqdt_tmp(:,:), sflx  )
             else   !use charlie's method
                call d3ddflux( ncol, vlc_dry(:,:,jvlc), qaerwat(:,:,mm), state%pmid, &
                               state%pdel, tvs, sflx, dqdt_tmp(:,:), dt )
             endif

             ! apportion dry deposition into turb and gravitational settling for tapes
             do i=1,ncol
                dep_trb(i)=sflx(i)*vlc_trb(i,jvlc)/vlc_dry(i,pver,jvlc)
                dep_grv(i)=sflx(i)*vlc_grv(i,pver,jvlc)/vlc_dry(i,pver,jvlc)
             enddo

             qaerwat(1:ncol,:,mm) = qaerwat(1:ncol,:,mm) + dqdt_tmp(1:ncol,:) * dt

          else  ! lphase == 2
             ! use pvprogseasalts instead (means making the top level 0)
             pvmzaer(:ncol,1)=0._r8
             pvmzaer(:ncol,2:pverp) = vlc_dry(:ncol,:,jvlc)
             fldcw => qqcw_get_field(pbuf, mm,lchnk)

             if(.true.) then ! use phil's method
             !      convert from meters/sec to pascals/sec
             !      pvprogseasalts(:,1) is assumed zero, use density from layer above in conversion
                pvmzaer(:ncol,2:pverp) = pvmzaer(:ncol,2:pverp) * rho(:ncol,:)*gravit

             !      calculate the tendencies and sfc fluxes from the above velocities
                call dust_sediment_tend( &
                     ncol,             dt,       state%pint(:,:), state%pmid, state%pdel, state%t , &
                     fldcw(:,:),  pvmzaer,  dqdt_tmp(:,:), sflx  )
             else   !use charlie's method
                call d3ddflux( ncol, vlc_dry(:,:,jvlc), fldcw(:,:), state%pmid, &
                               state%pdel, tvs, sflx, dqdt_tmp(:,:), dt )
             endif

             ! apportion dry deposition into turb and gravitational settling for tapes
             do i=1,ncol
                dep_trb(i)=sflx(i)*vlc_trb(i,jvlc)/vlc_dry(i,pver,jvlc)
                dep_grv(i)=sflx(i)*vlc_grv(i,pver,jvlc)/vlc_dry(i,pver,jvlc)
             enddo

             fldcw(1:ncol,:) = fldcw(1:ncol,:) + dqdt_tmp(1:ncol,:) * dt

             call outfld( trim(cnst_name_cw(mm))//'DDF', sflx, pcols, lchnk)
             call outfld( trim(cnst_name_cw(mm))//'TBF', dep_trb, pcols, lchnk )
             call outfld( trim(cnst_name_cw(mm))//'GVF', dep_grv, pcols, lchnk )
             aerdepdrycw(:ncol,mm) = sflx(:ncol)

          endif

          enddo   ! lspec = 0, nspec_amode(m)+1
       enddo   ! lphase = 1, 2
    enddo   ! m = 1, ntot_amode

    ! if the user has specified prescribed aerosol dep fluxes then 
    ! do not set cam_out dep fluxes according to the prognostic aerosols
    if (.not.aerodep_flx_prescribed()) then
       call set_srf_drydep(aerdepdryis, aerdepdrycw, cam_out)
    endif

  endsubroutine aero_model_drydep

  !=============================================================================
  !=============================================================================
  subroutine aero_model_wetdep( state, dt, dlf, cam_out, ptend, pbuf)

    use modal_aero_deposition, only: set_srf_wetdep
    use wetdep,                only: wetdepa_v2, wetdep_inputs_set, wetdep_inputs_t
    use modal_aero_data
    use modal_aero_calcsize,   only: modal_aero_calcsize_sub
    use modal_aero_wateruptake,only: modal_aero_wateruptake_dr


    ! args

    type(physics_state), intent(in)    :: state       ! Physics state variables
    real(r8),            intent(in)    :: dt          ! time step
    real(r8),            intent(in)    :: dlf(:,:)    ! shallow+deep convective detrainment [kg/kg/s]
    type(cam_out_t),     intent(inout) :: cam_out     ! export state
    type(physics_ptend), intent(out)   :: ptend       ! indivdual parameterization tendencies
    type(physics_buffer_desc), pointer :: pbuf(:)

    ! local vars

    integer :: m ! tracer index

    integer :: lchnk ! chunk identifier
    integer :: ncol ! number of atmospheric columns

    real(r8) :: iscavt(pcols, pver)

    integer :: mm
    integer :: i,k

    real(r8) :: icscavt(pcols, pver)
    real(r8) :: isscavt(pcols, pver)
    real(r8) :: bcscavt(pcols, pver)
    real(r8) :: bsscavt(pcols, pver)
    real(r8) :: sol_factb, sol_facti
    real(r8) :: sol_factic(pcols,pver)
    real(r8) :: sol_factbi, sol_factii, sol_factiic

    real(r8) :: sflx(pcols) ! deposition flux

    integer :: jnv ! index for scavcoefnv 3rd dimension
    integer :: lphase ! index for interstitial / cloudborne aerosol
    integer :: lspec ! index for aerosol number / chem-mass / water-mass
    integer :: lcoardust, lcoarnacl ! indices for coarse mode dust and seasalt masses
    real(r8) :: dqdt_tmp(pcols,pver) ! temporary array to hold tendency for 1 species
    real(r8) :: f_act_conv(pcols,pver) ! prescribed aerosol activation fraction for convective cloud ! rce 2010/05/01
    real(r8) :: f_act_conv_coarse(pcols,pver) ! similar but for coarse mode ! rce 2010/05/02
    real(r8) :: f_act_conv_coarse_dust, f_act_conv_coarse_nacl ! rce 2010/05/02
    real(r8) :: fracis_cw(pcols,pver)
    real(r8) :: hygro_sum_old(pcols,pver) ! before removal [sum of (mass*hydro/dens)]
    real(r8) :: hygro_sum_del(pcols,pver) ! removal change to [sum of (mass*hydro/dens)]
    real(r8) :: hygro_sum_old_ik, hygro_sum_new_ik
    real(r8) :: prec(pcols) ! precipitation rate
    real(r8) :: q_tmp(pcols,pver) ! temporary array to hold "most current" mixing ratio for 1 species
    real(r8) :: qqcw_tmp(pcols,pver) ! temporary array to hold qqcw ! rce 2010/05/01
    real(r8) :: scavcoefnv(pcols,pver,0:2) ! Dana and Hales coefficient (/mm) for
                                           ! cloud-borne num & vol (0),
                                           ! interstitial num (1), interstitial vol (2)
    real(r8) :: tmpa, tmpb
    real(r8) :: tmpdust, tmpnacl
    real(r8) :: water_old, water_new ! temporary old/new aerosol water mix-rat
    logical  :: isprx(pcols,pver) ! true if precipation
    real(r8) :: aerdepwetis(pcols,pcnst) ! aerosol wet deposition (interstitial)
    real(r8) :: aerdepwetcw(pcols,pcnst) ! aerosol wet deposition (cloud water)
    real(r8), pointer :: fldcw(:,:)

    real(r8), pointer :: dgnumwet(:,:,:)
    real(r8), pointer :: qaerwat(:,:,:)  ! aerosol water
    real(r8), pointer :: rate1ord_cw2pr_st(:,:)

    real(r8), pointer :: fracis(:,:,:)   ! fraction of transported species that are insoluble

    type(wetdep_inputs_t) :: dep_inputs

    lchnk = state%lchnk
    ncol  = state%ncol

    call physics_ptend_init(ptend, state%psetcols, 'aero_model_wetdep', lq=wetdep_lq)
    
    ! Do calculations of mode radius and water uptake if:
    ! 1) modal aerosols are affecting the climate, or
    ! 2) prognostic modal aerosols are enabled
    
    call t_startf('calcsize')
    ! for prognostic modal aerosols the transfer of mass between aitken and accumulation
    ! modes is done in conjunction with the dry radius calculation
    call modal_aero_calcsize_sub(state, ptend, dt, pbuf)
    call t_stopf('calcsize')

    call t_startf('wateruptake')
    call modal_aero_wateruptake_dr(state, pbuf)
    call t_stopf('wateruptake')

    if (nwetdep<1) return

    call wetdep_inputs_set( state, pbuf, dep_inputs )

    call pbuf_get_field(pbuf, dgnumwet_idx,       dgnumwet, start=(/1,1,1/), kount=(/pcols,pver,nmodes/) )
    call pbuf_get_field(pbuf, qaerwat_idx,        qaerwat,  start=(/1,1,1/), kount=(/pcols,pver,nmodes/) )
    call pbuf_get_field(pbuf, rate1_cw2pr_st_idx, rate1ord_cw2pr_st)
    call pbuf_get_field(pbuf, fracis_idx,         fracis, start=(/1,1,1/), kount=(/pcols, pver, pcnst/) )

    prec(:ncol)=0._r8
    do k=1,pver
       where (prec(:ncol) >= 1.e-7_r8)
          isprx(:ncol,k) = .true.
       elsewhere
          isprx(:ncol,k) = .false.
       endwhere
       prec(:ncol) = prec(:ncol) + (dep_inputs%prain(:ncol,k) + dep_inputs%cmfdqr(:ncol,k) - dep_inputs%evapr(:ncol,k)) &
            *state%pdel(:ncol,k)/gravit
    end do

    ! calculate the mass-weighted sol_factic for coarse mode species
    ! sol_factic_coarse(:,:) = 0.30_r8 ! tuned 1/4
    f_act_conv_coarse(:,:) = 0.60_r8 ! rce 2010/05/02
    f_act_conv_coarse_dust = 0.40_r8 ! rce 2010/05/02
    f_act_conv_coarse_nacl = 0.80_r8 ! rce 2010/05/02
    if (modeptr_coarse > 0) then
       lcoardust = lptr_dust_a_amode(modeptr_coarse)
       lcoarnacl = lptr_nacl_a_amode(modeptr_coarse)
       if ((lcoardust > 0) .and. (lcoarnacl > 0)) then
          do k = 1, pver
             do i = 1, ncol
                tmpdust = max( 0.0_r8, state%q(i,k,lcoardust) + ptend%q(i,k,lcoardust)*dt )
                tmpnacl = max( 0.0_r8, state%q(i,k,lcoarnacl) + ptend%q(i,k,lcoarnacl)*dt )
                if ((tmpdust+tmpnacl) > 1.0e-30_r8) then
                   ! sol_factic_coarse(i,k) = (0.2_r8*tmpdust + 0.4_r8*tmpnacl)/(tmpdust+tmpnacl) ! tuned 1/6
                   f_act_conv_coarse(i,k) = (f_act_conv_coarse_dust*tmpdust &
                        + f_act_conv_coarse_nacl*tmpnacl)/(tmpdust+tmpnacl) ! rce 2010/05/02
                end if
             end do
          end do
       end if
    end if

    scavcoefnv(:,:,0) = 0.0_r8 ! below-cloud scavcoef = 0.0 for cloud-borne species

    do m = 1, ntot_amode ! main loop over aerosol modes

       do lphase = 1, 2 ! loop over interstitial (1) and cloud-borne (2) forms

          ! sol_factb and sol_facti values
          ! sol_factb - currently this is basically a tuning factor
          ! sol_facti & sol_factic - currently has a physical basis, and reflects activation fraction
          !
          ! 2008-mar-07 rce - sol_factb (interstitial) changed from 0.3 to 0.1
          ! - sol_factic (interstitial, dust modes) changed from 1.0 to 0.5
          ! - sol_factic (cloud-borne, pcarb modes) no need to set it to 0.0
          ! because the cloud-borne pcarbon == 0 (no activation)
          !
          ! rce 2010/05/02
          ! prior to this date, sol_factic was used for convective in-cloud wet removal,
          ! and its value reflected a combination of an activation fraction (which varied between modes)
          ! and a tuning factor
          ! from this date forward, two parameters are used for convective in-cloud wet removal
          ! f_act_conv is the activation fraction
          ! note that "non-activation" of aerosol in air entrained into updrafts should
          ! be included here
          ! eventually we might use the activate routine (with w ~= 1 m/s) to calculate
          ! this, but there is still the entrainment issue
          ! sol_factic is strictly a tuning factor
          !
          if (lphase == 1) then ! interstial aerosol
             hygro_sum_old(:,:) = 0.0_r8
             hygro_sum_del(:,:) = 0.0_r8
             call modal_aero_bcscavcoef_get( m, ncol, isprx, dgnumwet, &
                  scavcoefnv(:,:,1), scavcoefnv(:,:,2) )

             sol_factb = 0.1_r8 ! all below-cloud scav ON (0.1 "tuning factor")
             ! sol_factb = 0.03_r8 ! all below-cloud scav ON (0.1 "tuning factor") ! tuned 1/6

             sol_facti = 0.0_r8 ! strat in-cloud scav totally OFF for institial

             sol_factic = 0.4_r8 ! xl 2010/05/20

             if (m == modeptr_pcarbon) then
                ! sol_factic = 0.0_r8 ! conv in-cloud scav OFF (0.0 activation fraction)
                f_act_conv = 0.0_r8 ! rce 2010/05/02
             else if ((m == modeptr_finedust) .or. (m == modeptr_coardust)) then
                ! sol_factic = 0.2_r8 ! conv in-cloud scav ON (0.5 activation fraction) ! tuned 1/4
                f_act_conv = 0.4_r8 ! rce 2010/05/02
             else
                ! sol_factic = 0.4_r8 ! conv in-cloud scav ON (1.0 activation fraction) ! tuned 1/4
                f_act_conv = 0.8_r8 ! rce 2010/05/02
             end if

          else ! cloud-borne aerosol (borne by stratiform cloud drops)

             sol_factb  = 0.0_r8   ! all below-cloud scav OFF (anything cloud-borne is located "in-cloud")
             sol_facti  = sol_facti_cloud_borne   ! strat  in-cloud scav cloud-borne tuning factor
             sol_factic = 0.0_r8   ! conv   in-cloud scav OFF (having this on would mean
                                   !        that conv precip collects strat droplets)
             f_act_conv = 0.0_r8   ! conv   in-cloud scav OFF (having this on would mean

          end if
          !
          ! rce 2010/05/03
          ! wetdepa has 6 "sol_fact" parameters:
          ! sol_facti, sol_factic, sol_factb for liquid cloud
          ! sol_factii, sol_factiic, sol_factbi for ice cloud
          ! the ice cloud parameters are optional, and if not provided, they default to
          ! one of the other sol_fact parameters (see subr. wetdepa about this)
          ! for now, we set the ice cloud parameters equal
          ! to their liquid cloud counterparts
          ! currently the ice parameters are not used in wetdepa as
          ! wetdepa sets "weight" (the ice cloud fraction) to 0.0
          ! if this changes, we will have to give more thought to
          ! the ice cloud parameter values
          !
          sol_factbi = sol_factb
          sol_factii = sol_facti
          sol_factiic = sol_factic(1,1)


          do lspec = 0, nspec_amode(m)+1 ! loop over number + chem constituents + water

             if (lspec == 0) then ! number
                if (lphase == 1) then
                   mm = numptr_amode(m)
                   jnv = 1
                else
                   mm = numptrcw_amode(m)
                   jnv = 0
                endif
             else if (lspec <= nspec_amode(m)) then ! non-water mass
                if (lphase == 1) then
                   mm = lmassptr_amode(lspec,m)
                   jnv = 2
                else
                   mm = lmassptrcw_amode(lspec,m)
                   jnv = 0
                endif
             else ! water mass
                ! bypass wet removal of aerosol water
                cycle
                if (lphase == 1) then
                   mm = 0
                   ! mm = lwaterptr_amode(m)
                   jnv = 2
                else
                   mm = 0
                   jnv = 0
                endif
             endif

             if (mm <= 0) cycle


             ! set f_act_conv for interstitial (lphase=1) coarse mode species
             ! for the convective in-cloud, we conceptually treat the coarse dust and seasalt
             ! as being externally mixed, and apply f_act_conv = f_act_conv_coarse_dust/nacl to dust/seasalt
             ! number and sulfate are conceptually partitioned to the dust and seasalt
             ! on a mass basis, so the f_act_conv for number and sulfate are
             ! mass-weighted averages of the values used for dust/seasalt
             if ((lphase == 1) .and. (m == modeptr_coarse)) then
                ! sol_factic = sol_factic_coarse
                f_act_conv = f_act_conv_coarse ! rce 2010/05/02
                if (lspec > 0) then
                   if (lmassptr_amode(lspec,m) == lptr_dust_a_amode(m)) then
                      ! sol_factic = 0.2_r8 ! tuned 1/4
                      f_act_conv = f_act_conv_coarse_dust ! rce 2010/05/02
                   else if (lmassptr_amode(lspec,m) == lptr_nacl_a_amode(m)) then
                      ! sol_factic = 0.4_r8 ! tuned 1/6
                      f_act_conv = f_act_conv_coarse_nacl ! rce 2010/05/02
                   end if
                end if
             end if


             if ((lphase == 1) .and. (lspec <= nspec_amode(m))) then
                ptend%lq(mm) = .TRUE.
                dqdt_tmp(:,:) = 0.0_r8
                ! q_tmp reflects changes from modal_aero_calcsize and is the "most current" q
                q_tmp(1:ncol,:) = state%q(1:ncol,:,mm) + ptend%q(1:ncol,:,mm)*dt
                fldcw => qqcw_get_field(pbuf, mm,lchnk)

                call wetdepa_v2( state%t, state%pmid, state%q(:,:,1), state%pdel, &
                     dep_inputs%cldt, dep_inputs%cldcu, dep_inputs%cmfdqr, &
                     dep_inputs%evapc, dep_inputs%conicw, dep_inputs%prain, dep_inputs%qme, &
                     dep_inputs%evapr, dep_inputs%totcond, q_tmp, dt, &
                     dqdt_tmp, iscavt, dep_inputs%cldv, dep_inputs%cldvcu, dep_inputs%cldvst, &
                     dlf, fracis(:,:,mm), sol_factb, ncol, &
                     scavcoefnv(:,:,jnv), &
                     is_strat_cloudborne=.false.,  &
                     rate1ord_cw2pr_st=rate1ord_cw2pr_st,  &
                     qqcw=fldcw,  &
                     f_act_conv=f_act_conv, &
                     icscavt=icscavt, isscavt=isscavt, bcscavt=bcscavt, bsscavt=bsscavt, &
                     sol_facti_in=sol_facti, sol_factbi_in=sol_factbi, sol_factii_in=sol_factii, &   ! rce 2010/05/03
                     sol_factic_in=sol_factic, sol_factiic_in=sol_factiic )                          ! rce 2010/05/03

                ptend%q(1:ncol,:,mm) = ptend%q(1:ncol,:,mm) + dqdt_tmp(1:ncol,:)

                call outfld( trim(cnst_name(mm))//'WET', dqdt_tmp(:,:), pcols, lchnk)
                call outfld( trim(cnst_name(mm))//'SIC', icscavt, pcols, lchnk)
                call outfld( trim(cnst_name(mm))//'SIS', isscavt, pcols, lchnk)
                call outfld( trim(cnst_name(mm))//'SBC', bcscavt, pcols, lchnk)
                call outfld( trim(cnst_name(mm))//'SBS', bsscavt, pcols, lchnk)

                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+dqdt_tmp(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
                call outfld( trim(cnst_name(mm))//'SFWET', sflx, pcols, lchnk)
                aerdepwetis(:ncol,mm) = sflx(:ncol)

                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+icscavt(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
                call outfld( trim(cnst_name(mm))//'SFSIC', sflx, pcols, lchnk)
                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+isscavt(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
                call outfld( trim(cnst_name(mm))//'SFSIS', sflx, pcols, lchnk)
                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+bcscavt(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
                call outfld( trim(cnst_name(mm))//'SFSBC', sflx, pcols, lchnk)
                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+bsscavt(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
                call outfld( trim(cnst_name(mm))//'SFSBS', sflx, pcols, lchnk)

                if (lspec > 0) then
                   tmpa = spechygro(lspectype_amode(lspec,m))/ &
                        specdens_amode(lspectype_amode(lspec,m))
                   tmpb = tmpa*dt
                   hygro_sum_old(1:ncol,:) = hygro_sum_old(1:ncol,:) &
                        + tmpa*q_tmp(1:ncol,:)
                   hygro_sum_del(1:ncol,:) = hygro_sum_del(1:ncol,:) &
                        + tmpb*dqdt_tmp(1:ncol,:)
                end if

             else if ((lphase == 1) .and. (lspec == nspec_amode(m)+1)) then
                ! aerosol water -- because of how wetdepa treats evaporation of stratiform
                ! precip, it is not appropriate to apply wetdepa to aerosol water
                ! instead, "hygro_sum" = [sum of (mass*hygro/dens)] is calculated before and
                ! after wet removal, and new water is calculated using
                ! new_water = old_water*min(10,(hygro_sum_new/hygro_sum_old))
                ! the "min(10,...)" is to avoid potential problems when hygro_sum_old ~= 0
                ! also, individual wet removal terms (ic,is,bc,bs) are not output to history
                ! ptend%lq(mm) = .TRUE.
                ! dqdt_tmp(:,:) = 0.0_r8
                do k = 1, pver
                   do i = 1, ncol
                      ! water_old = max( 0.0_r8, state%q(i,k,mm)+ptend%q(i,k,mm)*dt )
                      water_old = max( 0.0_r8, qaerwat(i,k,mm) )
                      hygro_sum_old_ik = max( 0.0_r8, hygro_sum_old(i,k) )
                      hygro_sum_new_ik = max( 0.0_r8, hygro_sum_old_ik+hygro_sum_del(i,k) )
                      if (hygro_sum_new_ik >= 10.0_r8*hygro_sum_old_ik) then
                         water_new = 10.0_r8*water_old
                      else
                         water_new = water_old*(hygro_sum_new_ik/hygro_sum_old_ik)
                      end if
                      ! dqdt_tmp(i,k) = (water_new - water_old)/dt
                      qaerwat(i,k,mm) = water_new
                   end do
                end do

                ! ptend%q(1:ncol,:,mm) = ptend%q(1:ncol,:,mm) + dqdt_tmp(1:ncol,:)

                ! call outfld( trim(cnst_name(mm))

                ! sflx(:)=0._r8
                ! do k=1,pver
                ! do i=1,ncol
                ! sflx(i)=sflx(i)+dqdt_tmp(i,k)*state%pdel(i,k)/gravit
                ! enddo
                ! enddo
                ! call outfld( trim(cnst_name(mm))

             else ! lphase == 2
                dqdt_tmp(:,:) = 0.0_r8
                qqcw_tmp(:,:) = 0.0_r8 ! rce 2010/05/01
                fldcw => qqcw_get_field(pbuf, mm,lchnk)

                call wetdepa_v2(state%t, state%pmid, state%q(:,:,1), state%pdel, &
                     dep_inputs%cldt, dep_inputs%cldcu, dep_inputs%cmfdqr, &
                     dep_inputs%evapc, dep_inputs%conicw, dep_inputs%prain, dep_inputs%qme, &
                     dep_inputs%evapr, dep_inputs%totcond, fldcw, dt, &
                     dqdt_tmp, iscavt, dep_inputs%cldv, dep_inputs%cldvcu, dep_inputs%cldvst, &
                     dlf, fracis_cw, sol_factb, ncol, &
                     scavcoefnv(:,:,jnv), &
                     is_strat_cloudborne=.true.,  &
                     rate1ord_cw2pr_st=rate1ord_cw2pr_st,  &
                     qqcw=qqcw_tmp,  &
                     f_act_conv=f_act_conv, &
                     icscavt=icscavt, isscavt=isscavt, bcscavt=bcscavt, bsscavt=bsscavt, &
                     sol_facti_in=sol_facti, sol_factbi_in=sol_factbi, sol_factii_in=sol_factii, &   ! rce 2010/05/03
                     sol_factic_in=sol_factic, sol_factiic_in=sol_factiic )                          ! rce 2010/05/03

                fldcw(1:ncol,:) = fldcw(1:ncol,:) + dqdt_tmp(1:ncol,:) * dt

                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+dqdt_tmp(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
                call outfld( trim(cnst_name_cw(mm))//'SFWET', sflx, pcols, lchnk)
                aerdepwetcw(:ncol,mm) = sflx(:ncol)

                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+icscavt(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
                call outfld( trim(cnst_name_cw(mm))//'SFSIC', sflx, pcols, lchnk)
                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+isscavt(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
                call outfld( trim(cnst_name_cw(mm))//'SFSIS', sflx, pcols, lchnk)
                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+bcscavt(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
                call outfld( trim(cnst_name_cw(mm))//'SFSBC', sflx, pcols, lchnk)
                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+bsscavt(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
                call outfld( trim(cnst_name_cw(mm))//'SFSBS', sflx, pcols, lchnk)

             endif

          enddo ! lspec = 0, nspec_amode(m)+1
       enddo ! lphase = 1, 2
    enddo ! m = 1, ntot_amode

    ! if the user has specified prescribed aerosol dep fluxes then
    ! do not set cam_out dep fluxes according to the prognostic aerosols
    if (.not.aerodep_flx_prescribed()) then
       call set_srf_wetdep(aerdepwetis, aerdepwetcw, cam_out)
    endif

  endsubroutine aero_model_wetdep

  !-------------------------------------------------------------------------
  ! provides aerosol surface area info for modal aerosols
  ! called from mo_usrrxt
  !-------------------------------------------------------------------------
  subroutine aero_model_surfarea( &
                  mmr, radmean, relhum, pmid, temp, strato_sad, &
                  sulfate, rho, ltrop, het1_ndx, pbuf, ncol, sfc, dm_aer, sad_total )

    use ref_pres,        only : top_lev=>trop_cloud_top_lev
    use modal_aero_data, only : sigmag_amode
    use mo_constants,    only : pi
    use modal_aero_data, only : nspec_amode
    use modal_aero_data, only : alnsg_amode

    ! dummy args
    real(r8), intent(in)    :: pmid(:,:)
    real(r8), intent(in)    :: temp(:,:)
    real(r8), intent(in)    :: mmr(:,:,:)
    real(r8), intent(in)    :: radmean      ! mean radii in cm
    real(r8), intent(in)    :: strato_sad(:,:)
    integer,  intent(in)    :: ncol
    integer,  intent(in)    :: ltrop(:)
    integer,  intent(in)    :: het1_ndx
    real(r8), intent(in)    :: relhum(:,:)
    real(r8), intent(in)    :: rho(:,:) ! total atm density (/cm^3)
    real(r8), intent(in)    :: sulfate(:,:)
    type(physics_buffer_desc), pointer :: pbuf(:)

    real(r8), intent(inout) :: sfc(:,:,:)
    real(r8), intent(inout) :: dm_aer(:,:,:)
    real(r8), intent(inout) :: sad_total(:,:)

    ! local vars

    real(r8), target :: sad_mode(pcols,pver,ntot_amode)
    real(r8) :: rho_air
    real(r8), pointer, dimension(:,:,:) :: dgnumwet
    integer :: l
    integer :: i,k,m 

    real(r8) :: chm_mass,tot_mass

    call pbuf_get_field(pbuf, dgnumwet_idx, dgnumwet )

    !
    ! compute surface aero for each mode; however, at this point we only use Aitken mode (mode 2 in MAM3; how
    ! can we move from hard-wiring this?) as the surface area for chemical reactions.
    !
    sad_mode = 0._r8
    sad_total = 0._r8
    dm_aer = 0._r8
    do i = 1,ncol
       do k = ltrop(i),pver
          rho_air = pmid(i,k)/(temp(i,k)*287.04_r8)
          do l=1,ntot_amode
             !
             ! compute a mass weighting of the number
             !
             tot_mass = 0._r8
             chm_mass = 0._r8
             do m=1,nspec_amode(l)
               if ( index_tot_mass(l,m) > 0 ) &
                    tot_mass = tot_mass + mmr(i,k,index_tot_mass(l,m))
               if ( index_chm_mass(l,m) > 0 ) &
                    chm_mass = chm_mass + mmr(i,k,index_chm_mass(l,m))
             end do
             if ( tot_mass > 0._r8 ) then
               sad_mode(i,k,l) = chm_mass/tot_mass * &
                    mmr(i,k,num_idx(l))*rho_air*pi*dgnumwet(i,k,l)**2*&
                    exp(2*alnsg_amode(l)**2)  ! m^2/m^3
               sad_mode(i,k,l) = 1.e-2_r8 * sad_mode(i,k,l) ! cm^2/cm^3
             else
               sad_mode(i,k,l) = 0._r8
             end if
          end do
          sad_total(i,k) = sum(sad_mode(i,k,:))
          dm_aer(i,k,:) = dgnumwet(i,k,:) * 1.e2_r8 ! convert m to cm

       enddo
    enddo

    sfc(:,:,:) = sad_mode(:,:,:) 

  end subroutine aero_model_surfarea

  !=============================================================================
  !=============================================================================
  subroutine aero_model_gasaerexch( loffset, ncol, lchnk, delt, reaction_rates, &
                                    tfld, pmid, pdel, mbar, relhum, &
                                    zm,  qh2o, cwat, cldfr, cldnum, &
                                    airdens, invariants, del_h2so4_gasprod,  &
                                    vmr0, vmr, pbuf )

    use time_manager,          only : get_nstep
    use modal_aero_coag,       only : modal_aero_coag_sub
    use modal_aero_gasaerexch, only : modal_aero_gasaerexch_sub
    use modal_aero_newnuc,     only : modal_aero_newnuc_sub
    use mo_setsox,             only : setsox, has_sox
    use modal_aero_data,       only : cnst_name_cw, qqcw_get_field

    !-----------------------------------------------------------------------
    !      ... dummy arguments
    !-----------------------------------------------------------------------
    integer,  intent(in) :: loffset                ! offset applied to modal aero "pointers"
    integer,  intent(in) :: ncol                   ! number columns in chunk
    integer,  intent(in) :: lchnk                  ! chunk index
    real(r8), intent(in) :: delt                   ! time step size (sec)
    real(r8), intent(in) :: reaction_rates(:,:,:)  ! reaction rates
    real(r8), intent(in) :: tfld(:,:)              ! temperature (K)
    real(r8), intent(in) :: pmid(:,:)              ! pressure at model levels (Pa)
    real(r8), intent(in) :: pdel(:,:)              ! pressure thickness of levels (Pa)
    real(r8), intent(in) :: mbar(:,:)              ! mean wet atmospheric mass ( amu )
    real(r8), intent(in) :: relhum(:,:)            ! relative humidity
    real(r8), intent(in) :: airdens(:,:)           ! total atms density (molec/cm**3)
    real(r8), intent(in) :: invariants(:,:,:)
    real(r8), intent(in) :: del_h2so4_gasprod(:,:) 
    real(r8), intent(in) :: zm(:,:) 
    real(r8), intent(in) :: qh2o(:,:) 
    real(r8), intent(in) :: cwat(:,:)          ! cloud liquid water content (kg/kg)
    real(r8), intent(in) :: cldfr(:,:) 
    real(r8), intent(in) :: cldnum(:,:)       ! droplet number concentration (#/kg)
    real(r8), intent(in) :: vmr0(:,:,:)       ! initial mixing ratios (before gas-phase chem changes)
    real(r8), intent(inout) :: vmr(:,:,:)         ! mixing ratios ( vmr )
    type(physics_buffer_desc), pointer :: pbuf(:)
    
    ! local vars 
    
    integer :: n, m
    integer :: i,k
    integer :: nstep

    real(r8) :: del_h2so4_aeruptk(ncol,pver)

    real(r8), pointer :: dgnum(:,:,:), dgnumwet(:,:,:), wetdens(:,:,:)
    real(r8), pointer :: pblh(:)                    ! pbl height (m)

    real(r8), dimension(ncol) :: wrk
    character(len=32)         :: name
    real(r8) :: dvmrcwdt(ncol,pver,gas_pcnst)
    real(r8) :: dvmrdt(ncol,pver,gas_pcnst)
    real(r8) :: vmrcw(ncol,pver,gas_pcnst)            ! cloud-borne aerosol (vmr)

    real(r8), pointer :: fldcw(:,:)

    call pbuf_get_field(pbuf, dgnum_idx,      dgnum,  start=(/1,1,1/), kount=(/pcols,pver,ntot_amode/) )
    call pbuf_get_field(pbuf, dgnumwet_idx,   dgnumwet )
    call pbuf_get_field(pbuf, wetdens_ap_idx, wetdens )
    call pbuf_get_field(pbuf, pblh_idx,       pblh)

    do n=1,ntot_amode
       call outfld(dgnum_name(n),dgnumwet(1:ncol,1:pver,n), ncol, lchnk )
    end do

! do gas-aerosol exchange (h2so4, msa, nh3 condensation)

    nstep = get_nstep()

    ! calculate tendency due to gas phase chemistry and processes
    dvmrdt(:ncol,:,:) = (vmr(:ncol,:,:) - vmr0(:ncol,:,:)) / delt
    do m = 1, gas_pcnst
      wrk(:) = 0._r8
      do k = 1,pver
        wrk(:ncol) = wrk(:ncol) + dvmrdt(:ncol,k,m)*adv_mass(m)/mbar(:ncol,k)*pdel(:ncol,k)/gravit
      end do
      name = 'GS_'//trim(solsym(m))
      call outfld( name, wrk(:ncol), ncol, lchnk )
    enddo

!
! Aerosol processes ...
!
    call qqcw2vmr( lchnk, vmrcw, mbar, ncol, loffset, pbuf )

    dvmrdt(:ncol,:,:) = vmr(:ncol,:,:)
    dvmrcwdt(:ncol,:,:) = vmrcw(:ncol,:,:)

  ! aqueous chemistry ...

    if( has_sox ) then
       call setsox(   &
            ncol,     &
            lchnk,    &
            loffset,  &
            delt,     &
            pmid,     &
            pdel,     &
            tfld,     &
            mbar,     &
            cwat,     &
            cldfr,    &
            cldnum,   &
            airdens,  &
            invariants, &
            vmrcw,    &
            vmr       &
            )
    endif

!   Tendency due to aqueous chemistry 
    dvmrdt   = (vmr - dvmrdt) / delt
    dvmrcwdt = (vmrcw - dvmrcwdt) / delt
    do m = 1, gas_pcnst
      wrk(:) = 0._r8
      do k = 1,pver
        wrk(:ncol) = wrk(:ncol) + dvmrdt(:ncol,k,m) * adv_mass(m)/mbar(:ncol,k)*pdel(:ncol,k)/gravit
      end do
      name = 'AQ_'//trim(solsym(m))
      call outfld( name, wrk(:ncol), ncol, lchnk )
    enddo

! do gas-aerosol exchange (h2so4, msa, nh3 condensation)

    if (ndx_h2so4 > 0) then
       del_h2so4_aeruptk(1:ncol,:) = vmr(1:ncol,:,ndx_h2so4)
    else
       del_h2so4_aeruptk(:,:) = 0.0_r8
    endif

    call t_startf('modal_gas-aer_exchng')

    call modal_aero_gasaerexch_sub(                         &
         lchnk,    ncol,     nstep,            &
         loffset,            delt,             &
         tfld,     pmid,     pdel,             &
         vmr,                vmrcw,            &
         dvmrdt,             dvmrcwdt,     &
         dgnum,              dgnumwet     )

    if (ndx_h2so4 > 0) then
       del_h2so4_aeruptk(1:ncol,:) = vmr(1:ncol,:,ndx_h2so4) - del_h2so4_aeruptk(1:ncol,:)
    endif

    call t_stopf('modal_gas-aer_exchng')

    call t_startf('modal_nucl')

    ! do aerosol nucleation (new particle formation)
    call modal_aero_newnuc_sub(                             &
         lchnk,    ncol,     nstep,            &
         loffset,            delt,             &
         tfld,     pmid,     pdel,             &
         zm,       pblh,                       &
         qh2o,     cldfr,                      &
         vmr,                                  &
         del_h2so4_gasprod,  del_h2so4_aeruptk )

    call t_stopf('modal_nucl')

    call t_startf('modal_coag')

    ! do aerosol coagulation
    call modal_aero_coag_sub(                               &
         lchnk,    ncol,     nstep,            &
         loffset,            delt,             &
         tfld,     pmid,     pdel,             &
         vmr,                                  &
         dgnum,              dgnumwet,         &
         wetdens                          )

    call t_stopf('modal_coag')

    call vmr2qqcw( lchnk, vmrcw, mbar, ncol, loffset, pbuf )

    ! diagnostics for cloud-borne aerosols... 
    do n = 1,pcnst
       fldcw => qqcw_get_field(pbuf,n,lchnk,errorhandle=.true.)
       if(associated(fldcw)) then
          call outfld( cnst_name_cw(n), fldcw(:,:), pcols, lchnk )
       endif
    end do

  end subroutine aero_model_gasaerexch

  !=============================================================================
  !=============================================================================
  subroutine aero_model_emissions( state, cam_in )
    use seasalt_model, only: seasalt_emis, seasalt_names, seasalt_indices, seasalt_active,seasalt_nbin
    use dust_model,    only: dust_emis, dust_names, dust_indices, dust_active,dust_nbin, dust_nnum
    use physics_types, only: physics_state

    ! Arguments:

    type(physics_state),    intent(in)    :: state   ! Physics state variables
    type(cam_in_t),         intent(inout) :: cam_in  ! import state

    ! local vars

    integer :: lchnk, ncol
    integer :: m, mm
    real(r8) :: soil_erod_tmp(pcols)
    real(r8) :: sflx(pcols)   ! accumulate over all bins for output
    real(r8) :: u10cubed(pcols)
    real (r8), parameter :: z0=0.0001_r8  ! m roughness length over oceans--from ocean model

    lchnk = state%lchnk
    ncol = state%ncol

    if (dust_active) then

       call dust_emis( ncol, lchnk, cam_in%dstflx, cam_in%cflx, soil_erod_tmp )

       ! some dust emis diagnostics ...
       sflx(:)=0._r8
       do m=1,dust_nbin+dust_nnum
          mm = dust_indices(m)
          if (m<=dust_nbin) sflx(:ncol)=sflx(:ncol)+cam_in%cflx(:ncol,mm)
          call outfld(trim(dust_names(m))//'SF',cam_in%cflx(:,mm),pcols, lchnk)
       enddo
       call outfld('DSTSFMBL',sflx(:),pcols,lchnk)
       call outfld('LND_MBL',soil_erod_tmp(:),pcols, lchnk )
    endif

    if (seasalt_active) then
       u10cubed(:ncol)=sqrt(state%u(:ncol,pver)**2+state%v(:ncol,pver)**2)
       ! move the winds to 10m high from the midpoint of the gridbox:
       ! follows Tie and Seinfeld and Pandis, p.859 with math.

       u10cubed(:ncol)=u10cubed(:ncol)*log(10._r8/z0)/log(state%zm(:ncol,pver)/z0)

       ! we need them to the 3.41 power, according to Gong et al., 1997:
       u10cubed(:ncol)=u10cubed(:ncol)**3.41_r8

       sflx(:)=0._r8

       call seasalt_emis( u10cubed, cam_in%sst, cam_in%ocnfrac, ncol, cam_in%cflx )

       do m=1,seasalt_nbin
          mm = seasalt_indices(m)
          sflx(:ncol)=sflx(:ncol)+cam_in%cflx(:ncol,mm)
          call outfld(trim(seasalt_names(m))//'SF',cam_in%cflx(:,mm),pcols,lchnk)
       enddo
       call outfld('SSTSFMBL',sflx(:),pcols,lchnk)
    endif

  end subroutine aero_model_emissions

  !===============================================================================
  ! private methods


  !===============================================================================
  !===============================================================================
  subroutine modal_aero_bcscavcoef_init
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! Computes lookup table for aerosol impaction/interception scavenging rates
    !
    ! Authors: R. Easter
    !
    !-----------------------------------------------------------------------
    
    use shr_kind_mod,only: r8 => shr_kind_r8
    use modal_aero_data
    use abortutils,  only: endrun

    implicit none


    !   local variables
    integer nnfit_maxd
    parameter (nnfit_maxd=27)

    integer i, jgrow, jdens, jpress, jtemp, ll, mode, nnfit
    integer lunerr

    real(r8) dg0, dg0_cgs, press, &
         rhodryaero, rhowetaero, rhowetaero_cgs, rmserr, &
         scavratenum, scavratevol, sigmag,                &
         temp, wetdiaratio, wetvolratio
    real(r8) aafitnum(1), xxfitnum(1,nnfit_maxd), yyfitnum(nnfit_maxd)
    real(r8) aafitvol(1), xxfitvol(1,nnfit_maxd), yyfitvol(nnfit_maxd)

    
    lunerr = 6
    dlndg_nimptblgrow = log( 1.25_r8 )

    modeloop: do mode = 1, ntot_amode

       sigmag = sigmag_amode(mode)

       ll = lspectype_amode(1,mode)
       rhodryaero = specdens_amode(ll)

       growloop: do jgrow = nimptblgrow_mind, nimptblgrow_maxd

          wetdiaratio = exp( jgrow*dlndg_nimptblgrow )
          dg0 = dgnum_amode(mode)*wetdiaratio

          wetvolratio = exp( jgrow*dlndg_nimptblgrow*3._r8 )
          rhowetaero = 1.0_r8 + (rhodryaero-1.0_r8)/wetvolratio
          rhowetaero = min( rhowetaero, rhodryaero )

          !
          !   compute impaction scavenging rates at 1 temp-press pair and save
          !
          nnfit = 0

          temp = 273.16_r8
          press = 0.75e6_r8   ! dynes/cm2
          rhowetaero = rhodryaero

          dg0_cgs = dg0*1.0e2_r8   ! m to cm
          rhowetaero_cgs = rhowetaero*1.0e-3_r8   ! kg/m3 to g/cm3
          call calc_1_impact_rate( &
               dg0_cgs, sigmag, rhowetaero_cgs, temp, press, &
               scavratenum, scavratevol, lunerr )

          nnfit = nnfit + 1
          if (nnfit .gt. nnfit_maxd) then
             write(lunerr,9110)
             call endrun()
          end if
9110      format( '*** subr. modal_aero_bcscavcoef_init -- nnfit too big' )

          xxfitnum(1,nnfit) = 1._r8
          yyfitnum(nnfit) = log( scavratenum )

          xxfitvol(1,nnfit) = 1._r8
          yyfitvol(nnfit) = log( scavratevol )

5900      continue

          !
          ! skip mlinfit stuff because scav table no longer has dependencies on
          !    air temp, air press, and particle wet density
          ! just load the log( scavrate--- ) values
          !
          !!
          !!   do linear regression
          !!	log(scavrate) = a1 + a2*log(wetdens)
          !!
          !	call mlinft( xxfitnum, yyfitnum, aafitnum, nnfit, 1, 1, rmserr )
          !	call mlinft( xxfitvol, yyfitvol, aafitvol, nnfit, 1, 1, rmserr )
          !
          !	scavimptblnum(jgrow,mode) = aafitnum(1)
          !	scavimptblvol(jgrow,mode) = aafitvol(1)

          scavimptblnum(jgrow,mode) = yyfitnum(1)
          scavimptblvol(jgrow,mode) = yyfitvol(1)

       enddo growloop
    enddo modeloop
    return
  end subroutine modal_aero_bcscavcoef_init

  !===============================================================================
  !===============================================================================
  subroutine modal_aero_depvel_part( ncol, t, pmid, ram1, fv, vlc_dry, vlc_trb, vlc_grv,  &
                                     radius_part, density_part, sig_part, moment, lchnk )

!    calculates surface deposition velocity of particles
!    L. Zhang, S. Gong, J. Padro, and L. Barrie
!    A size-seggregated particle dry deposition scheme for an atmospheric aerosol module
!    Atmospheric Environment, 35, 549-560, 2001.
!
!    Authors: X. Liu

    !
    ! !USES
    !
    use physconst,     only: pi,boltz, gravit, rair
    use mo_drydep,     only: n_land_type, fraction_landuse

    ! !ARGUMENTS:
    !
    implicit none
    !
    real(r8), intent(in) :: t(pcols,pver)       !atm temperature (K)
    real(r8), intent(in) :: pmid(pcols,pver)    !atm pressure (Pa)
    real(r8), intent(in) :: fv(pcols)           !friction velocity (m/s)
    real(r8), intent(in) :: ram1(pcols)         !aerodynamical resistance (s/m)
    real(r8), intent(in) :: radius_part(pcols,pver)    ! mean (volume/number) particle radius (m)
    real(r8), intent(in) :: density_part(pcols,pver)   ! density of particle material (kg/m3)
    real(r8), intent(in) :: sig_part(pcols,pver)       ! geometric standard deviation of particles
    integer,  intent(in) :: moment ! moment of size distribution (0 for number, 2 for surface area, 3 for volume)
    integer,  intent(in) :: ncol
    integer,  intent(in) :: lchnk

    real(r8), intent(out) :: vlc_trb(pcols)       !Turbulent deposn velocity (m/s)
    real(r8), intent(out) :: vlc_grv(pcols,pver)       !grav deposn velocity (m/s)
    real(r8), intent(out) :: vlc_dry(pcols,pver)       !dry deposn velocity (m/s)
    !------------------------------------------------------------------------

    !------------------------------------------------------------------------
    ! Local Variables
    integer  :: m,i,k,ix                !indices
    real(r8) :: rho     !atm density (kg/m**3)
    real(r8) :: vsc_dyn_atm(pcols,pver)   ![kg m-1 s-1] Dynamic viscosity of air
    real(r8) :: vsc_knm_atm(pcols,pver)   ![m2 s-1] Kinematic viscosity of atmosphere
    real(r8) :: shm_nbr       ![frc] Schmidt number
    real(r8) :: stk_nbr       ![frc] Stokes number
    real(r8) :: mfp_atm(pcols,pver)       ![m] Mean free path of air
    real(r8) :: dff_aer       ![m2 s-1] Brownian diffusivity of particle
    real(r8) :: slp_crc(pcols,pver) ![frc] Slip correction factor
    real(r8) :: rss_trb       ![s m-1] Resistance to turbulent deposition
    real(r8) :: rss_lmn       ![s m-1] Quasi-laminar layer resistance
    real(r8) :: brownian      ! collection efficiency for Browning diffusion
    real(r8) :: impaction     ! collection efficiency for impaction
    real(r8) :: interception  ! collection efficiency for interception
    real(r8) :: stickfrac     ! fraction of particles sticking to surface
    real(r8) :: radius_moment(pcols,pver) ! median radius (m) for moment
    real(r8) :: lnsig         ! ln(sig_part)
    real(r8) :: dispersion    ! accounts for influence of size dist dispersion on bulk settling velocity
                              ! assuming radius_part is number mode radius * exp(1.5 ln(sigma))

    integer  :: lt
    real(r8) :: lnd_frc
    real(r8) :: wrk1, wrk2, wrk3

    ! constants
    real(r8) gamma(11)      ! exponent of schmidt number
!   data gamma/0.54d+00,  0.56d+00,  0.57d+00,  0.54d+00,  0.54d+00, &
!              0.56d+00,  0.54d+00,  0.54d+00,  0.54d+00,  0.56d+00, &
!              0.50d+00/
    data gamma/0.56e+00_r8,  0.54e+00_r8,  0.54e+00_r8,  0.56e+00_r8,  0.56e+00_r8, &        
               0.56e+00_r8,  0.50e+00_r8,  0.54e+00_r8,  0.54e+00_r8,  0.54e+00_r8, &
               0.54e+00_r8/
    save gamma

    real(r8) alpha(11)      ! parameter for impaction
!   data alpha/50.00d+00,  0.95d+00,  0.80d+00,  1.20d+00,  1.30d+00, &
!               0.80d+00, 50.00d+00, 50.00d+00,  2.00d+00,  1.50d+00, &
!             100.00d+00/
    data alpha/1.50e+00_r8,   1.20e+00_r8,  1.20e+00_r8,  0.80e+00_r8,  1.00e+00_r8, &
               0.80e+00_r8, 100.00e+00_r8, 50.00e+00_r8,  2.00e+00_r8,  1.20e+00_r8, &
              50.00e+00_r8/
    save alpha

    real(r8) radius_collector(11) ! radius (m) of surface collectors
!   data radius_collector/-1.00d+00,  5.10d-03,  3.50d-03,  3.20d-03, 10.00d-03, &
!                          5.00d-03, -1.00d+00, -1.00d+00, 10.00d-03, 10.00d-03, &
!                         -1.00d+00/
    data radius_collector/10.00e-03_r8,  3.50e-03_r8,  3.50e-03_r8,  5.10e-03_r8,  2.00e-03_r8, &
                           5.00e-03_r8, -1.00e+00_r8, -1.00e+00_r8, 10.00e-03_r8,  3.50e-03_r8, &
                          -1.00e+00_r8/
    save radius_collector

    integer            :: iwet(11) ! flag for wet surface = 1, otherwise = -1
!   data iwet/1,   -1,   -1,   -1,   -1,  &
!            -1,   -1,   -1,    1,   -1,  &
!             1/
    data iwet/-1,  -1,   -1,   -1,   -1,  &
              -1,   1,   -1,    1,   -1,  &
              -1/
    save iwet


    !------------------------------------------------------------------------
    do k=1,pver
       do i=1,ncol

          lnsig = log(sig_part(i,k))
! use a maximum radius of 50 microns when calculating deposition velocity
          radius_moment(i,k) = min(50.0e-6_r8,radius_part(i,k))*   &
                          exp((float(moment)-1.5_r8)*lnsig*lnsig)
          dispersion = exp(2._r8*lnsig*lnsig)

          rho=pmid(i,k)/rair/t(i,k)

          ! Quasi-laminar layer resistance: call rss_lmn_get
          ! Size-independent thermokinetic properties
          vsc_dyn_atm(i,k) = 1.72e-5_r8 * ((t(i,k)/273.0_r8)**1.5_r8) * 393.0_r8 / &
               (t(i,k)+120.0_r8)      ![kg m-1 s-1] RoY94 p. 102
          mfp_atm(i,k) = 2.0_r8 * vsc_dyn_atm(i,k) / &   ![m] SeP97 p. 455
               (pmid(i,k)*sqrt(8.0_r8/(pi*rair*t(i,k))))
          vsc_knm_atm(i,k) = vsc_dyn_atm(i,k) / rho ![m2 s-1] Kinematic viscosity of air

          slp_crc(i,k) = 1.0_r8 + mfp_atm(i,k) * &
                  (1.257_r8+0.4_r8*exp(-1.1_r8*radius_moment(i,k)/(mfp_atm(i,k)))) / &
                  radius_moment(i,k)   ![frc] Slip correction factor SeP97 p. 464
          vlc_grv(i,k) = (4.0_r8/18.0_r8) * radius_moment(i,k)*radius_moment(i,k)*density_part(i,k)* &
                  gravit*slp_crc(i,k) / vsc_dyn_atm(i,k) ![m s-1] Stokes' settling velocity SeP97 p. 466
          vlc_grv(i,k) = vlc_grv(i,k) * dispersion

          vlc_dry(i,k)=vlc_grv(i,k)
       enddo
    enddo
    k=pver  ! only look at bottom level for next part
    do i=1,ncol
       dff_aer = boltz * t(i,k) * slp_crc(i,k) / &    ![m2 s-1]
                 (6.0_r8*pi*vsc_dyn_atm(i,k)*radius_moment(i,k)) !SeP97 p.474
       shm_nbr = vsc_knm_atm(i,k) / dff_aer                        ![frc] SeP97 p.972

       wrk2 = 0._r8
       wrk3 = 0._r8
       do lt = 1,n_land_type
          lnd_frc = fraction_landuse(i,lt,lchnk)
          if ( lnd_frc /= 0._r8 ) then
             brownian = shm_nbr**(-gamma(lt))
             if (radius_collector(lt) > 0.0_r8) then
!       vegetated surface
                stk_nbr = vlc_grv(i,k) * fv(i) / (gravit*radius_collector(lt))
                interception = 2.0_r8*(radius_moment(i,k)/radius_collector(lt))**2.0_r8
             else
!       non-vegetated surface
                stk_nbr = vlc_grv(i,k) * fv(i) * fv(i) / (gravit*vsc_knm_atm(i,k))  ![frc] SeP97 p.965
                interception = 0.0_r8
             endif
             impaction = (stk_nbr/(alpha(lt)+stk_nbr))**2.0_r8   

             if (iwet(lt) > 0) then
                stickfrac = 1.0_r8
             else
                stickfrac = exp(-sqrt(stk_nbr))
                if (stickfrac < 1.0e-10_r8) stickfrac = 1.0e-10_r8
             endif
             rss_lmn = 1.0_r8 / (3.0_r8 * fv(i) * stickfrac * (brownian+interception+impaction))
             rss_trb = ram1(i) + rss_lmn + ram1(i)*rss_lmn*vlc_grv(i,k)

             wrk1 = 1.0_r8 / rss_trb
             wrk2 = wrk2 + lnd_frc*( wrk1 )
             wrk3 = wrk3 + lnd_frc*( wrk1 + vlc_grv(i,k) )
          endif
       enddo  ! n_land_type
       vlc_trb(i) = wrk2
       vlc_dry(i,k) = wrk3
    enddo !ncol

    return
  end subroutine modal_aero_depvel_part

  !===============================================================================
  subroutine modal_aero_bcscavcoef_get( m, ncol, isprx, dgn_awet, scavcoefnum, scavcoefvol )

    use modal_aero_data
    !-----------------------------------------------------------------------
    implicit none

    integer,intent(in) :: m, ncol
    logical,intent(in):: isprx(pcols,pver)
    real(r8), intent(in) :: dgn_awet(pcols,pver,ntot_amode)
    real(r8), intent(out) :: scavcoefnum(pcols,pver), scavcoefvol(pcols,pver)

    integer i, k, jgrow
    real(r8) dumdgratio, xgrow, dumfhi, dumflo, scavimpvol, scavimpnum


    do k = 1, pver
       do i = 1, ncol

          ! do only if no precip
          if ( isprx(i,k) ) then
             !
             ! interpolate table values using log of (actual-wet-size)/(base-dry-size)

             dumdgratio = dgn_awet(i,k,m)/dgnum_amode(m)

             if ((dumdgratio .ge. 0.99_r8) .and. (dumdgratio .le. 1.01_r8)) then
                scavimpvol = scavimptblvol(0,m)
                scavimpnum = scavimptblnum(0,m)
             else
                xgrow = log( dumdgratio ) / dlndg_nimptblgrow
                jgrow = int( xgrow )
                if (xgrow .lt. 0._r8) jgrow = jgrow - 1
                if (jgrow .lt. nimptblgrow_mind) then
                   jgrow = nimptblgrow_mind
                   xgrow = jgrow
                else
                   jgrow = min( jgrow, nimptblgrow_maxd-1 )
                end if

                dumfhi = xgrow - jgrow
                dumflo = 1._r8 - dumfhi

                scavimpvol = dumflo*scavimptblvol(jgrow,m) + &
                     dumfhi*scavimptblvol(jgrow+1,m)
                scavimpnum = dumflo*scavimptblnum(jgrow,m) + &
                     dumfhi*scavimptblnum(jgrow+1,m)

             end if

             ! impaction scavenging removal amount for volume
             scavcoefvol(i,k) = exp( scavimpvol )
             ! impaction scavenging removal amount to number
             scavcoefnum(i,k) = exp( scavimpnum )

             ! scavcoef = impaction scav rate (1/h) for precip = 1 mm/h
             ! scavcoef = impaction scav rate (1/s) for precip = pfx_inrain
             ! (scavcoef/3600) = impaction scav rate (1/s) for precip = 1 mm/h
             ! (pfx_inrain*3600) = in-rain-area precip rate (mm/h)
             ! impactrate = (scavcoef/3600) * (pfx_inrain*3600)
          else
             scavcoefvol(i,k) = 0._r8
             scavcoefnum(i,k) = 0._r8
          end if

       end do
    end do

    return
  end subroutine modal_aero_bcscavcoef_get

  !===============================================================================
	subroutine calc_1_impact_rate(             &
     		dg0, sigmag, rhoaero, temp, press, &
     		scavratenum, scavratevol, lunerr )
   !
   !   routine computes a single impaction scavenging rate
   !	for precipitation rate of 1 mm/h
   !
   !   dg0 = geometric mean diameter of aerosol number size distrib. (cm)
   !   sigmag = geometric standard deviation of size distrib.
   !   rhoaero = density of aerosol particles (g/cm^3)
   !   temp = temperature (K)
   !   press = pressure (dyne/cm^2)
   !   scavratenum = number scavenging rate (1/h)
   !   scavratevol = volume or mass scavenging rate (1/h)
   !   lunerr = logical unit for error message
   !
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mo_constants, only: boltz_cgs, pi, rhowater => rhoh2o_cgs, &
                           gravity => gravity_cgs, rgas => rgas_cgs

   implicit none

   !   subr. parameters
   integer lunerr
   real(r8) dg0, sigmag, rhoaero, temp, press, scavratenum, scavratevol

   !   local variables
   integer nrainsvmax
   parameter (nrainsvmax=50)
   real(r8) rrainsv(nrainsvmax), xnumrainsv(nrainsvmax),&
        vfallrainsv(nrainsvmax)

   integer naerosvmax
   parameter (naerosvmax=51)
   real(r8) aaerosv(naerosvmax), &
     	ynumaerosv(naerosvmax), yvolaerosv(naerosvmax)

   integer i, ja, jr, na, nr
   real(r8) a, aerodiffus, aeromass, ag0, airdynvisc, airkinvisc
   real(r8) anumsum, avolsum, cair, chi
   real(r8) d, dr, dum, dumfuchs, dx
   real(r8) ebrown, eimpact, eintercept, etotal, freepath
   real(r8) precip, precipmmhr, precipsum
   real(r8) r, rainsweepout, reynolds, rhi, rhoair, rlo, rnumsum
   real(r8) scavsumnum, scavsumnumbb
   real(r8) scavsumvol, scavsumvolbb
   real(r8) schmidt, sqrtreynolds, sstar, stokes, sx              
   real(r8) taurelax, vfall, vfallstp
   real(r8) x, xg0, xg3, xhi, xlo, xmuwaterair                     

   
   rlo = .005_r8
   rhi = .250_r8
   dr = 0.005_r8
   nr = 1 + nint( (rhi-rlo)/dr )
   if (nr .gt. nrainsvmax) then
      write(lunerr,9110)
      call endrun()
   end if

9110 format( '*** subr. calc_1_impact_rate -- nr > nrainsvmax' )

   precipmmhr = 1.0_r8
   precip = precipmmhr/36000._r8

   ag0 = dg0/2._r8
   sx = log( sigmag )
   xg0 = log( ag0 )
   xg3 = xg0 + 3._r8*sx*sx

   xlo = xg3 - 4._r8*sx
   xhi = xg3 + 4._r8*sx
   dx = 0.2_r8*sx

   dx = max( 0.2_r8*sx, 0.01_r8 )
   xlo = xg3 - max( 4._r8*sx, 2._r8*dx )
   xhi = xg3 + max( 4._r8*sx, 2._r8*dx )

   na = 1 + nint( (xhi-xlo)/dx )
   if (na .gt. naerosvmax) then
      write(lunerr,9120)
      call endrun()
   end if

9120 format( '*** subr. calc_1_impact_rate -- na > naerosvmax' )

   !   air molar density
   cair = press/(rgas*temp)
   !   air mass density
   rhoair = 28.966_r8*cair
   !   molecular freepath
   freepath = 2.8052e-10_r8/cair
   !   air dynamic viscosity
   airdynvisc = 1.8325e-4_r8 * (416.16_r8/(temp+120._r8)) *    &
        ((temp/296.16_r8)**1.5_r8)
   !   air kinemaic viscosity
   airkinvisc = airdynvisc/rhoair
   !   ratio of water viscosity to air viscosity (from Slinn)
   xmuwaterair = 60.0_r8

   !
   !   compute rain drop number concentrations
   !	rrainsv = raindrop radius (cm)
   !	xnumrainsv = raindrop number concentration (#/cm^3)
   !		(number in the bin, not number density)
   !	vfallrainsv = fall velocity (cm/s)
   !
   precipsum = 0._r8
   do i = 1, nr
      r = rlo + (i-1)*dr
      rrainsv(i) = r
      xnumrainsv(i) = exp( -r/2.7e-2_r8 )

      d = 2._r8*r
      if (d .le. 0.007_r8) then
         vfallstp = 2.88e5_r8 * d**2._r8
      else if (d .le. 0.025_r8) then
         vfallstp = 2.8008e4_r8 * d**1.528_r8
      else if (d .le. 0.1_r8) then
         vfallstp = 4104.9_r8 * d**1.008_r8
      else if (d .le. 0.25_r8) then
         vfallstp = 1812.1_r8 * d**0.638_r8
      else
         vfallstp = 1069.8_r8 * d**0.235_r8
      end if

      vfall = vfallstp * sqrt(1.204e-3_r8/rhoair)
      vfallrainsv(i) = vfall
      precipsum = precipsum + vfall*(r**3)*xnumrainsv(i)
   end do
   precipsum = precipsum*pi*1.333333_r8

   rnumsum = 0._r8
   do i = 1, nr
      xnumrainsv(i) = xnumrainsv(i)*(precip/precipsum)
      rnumsum = rnumsum + xnumrainsv(i)
   end do

   !
   !   compute aerosol concentrations
   !	aaerosv = particle radius (cm)
   !	fnumaerosv = fraction of total number in the bin (--)
   !	fvolaerosv = fraction of total volume in the bin (--)
   !
   anumsum = 0._r8
   avolsum = 0._r8
   do i = 1, na
      x = xlo + (i-1)*dx
      a = exp( x )
      aaerosv(i) = a
      dum = (x - xg0)/sx
      ynumaerosv(i) = exp( -0.5_r8*dum*dum )
      yvolaerosv(i) = ynumaerosv(i)*1.3333_r8*pi*a*a*a
      anumsum = anumsum + ynumaerosv(i)
      avolsum = avolsum + yvolaerosv(i)
   end do

   do i = 1, na
      ynumaerosv(i) = ynumaerosv(i)/anumsum
      yvolaerosv(i) = yvolaerosv(i)/avolsum
   end do


   !
   !   compute scavenging
   !
   scavsumnum = 0._r8
   scavsumvol = 0._r8
   !
   !   outer loop for rain drop radius
   !
   jr_loop: do jr = 1, nr

      r = rrainsv(jr)
      vfall = vfallrainsv(jr)

      reynolds = r * vfall / airkinvisc
      sqrtreynolds = sqrt( reynolds )

      !
      !   inner loop for aerosol particle radius
      !
      scavsumnumbb = 0._r8
      scavsumvolbb = 0._r8

      ja_loop: do ja = 1, na

         a = aaerosv(ja)

         chi = a/r

         dum = freepath/a
         dumfuchs = 1._r8 + 1.246_r8*dum + 0.42_r8*dum*exp(-0.87_r8/dum)
         taurelax = 2._r8*rhoaero*a*a*dumfuchs/(9._r8*rhoair*airkinvisc)

         aeromass = 4._r8*pi*a*a*a*rhoaero/3._r8
         aerodiffus = boltz_cgs*temp*taurelax/aeromass

         schmidt = airkinvisc/aerodiffus
         stokes = vfall*taurelax/r

         ebrown = 4._r8*(1._r8 + 0.4_r8*sqrtreynolds*(schmidt**0.3333333_r8)) /  &
              (reynolds*schmidt)

         dum = (1._r8 + 2._r8*xmuwaterair*chi) /         &
              (1._r8 + xmuwaterair/sqrtreynolds)
         eintercept = 4._r8*chi*(chi + dum)

         dum = log( 1._r8 + reynolds )
         sstar = (1.2_r8 + dum/12._r8) / (1._r8 + dum)
         eimpact = 0._r8
         if (stokes .gt. sstar) then
	    dum = stokes - sstar
	    eimpact = (dum/(dum+0.6666667_r8)) ** 1.5_r8
         end if

         etotal = ebrown + eintercept + eimpact
         etotal = min( etotal, 1.0_r8 )

         rainsweepout = xnumrainsv(jr)*4._r8*pi*r*r*vfall

         scavsumnumbb = scavsumnumbb + rainsweepout*etotal*ynumaerosv(ja)
         scavsumvolbb = scavsumvolbb + rainsweepout*etotal*yvolaerosv(ja)

      enddo ja_loop

      scavsumnum = scavsumnum + scavsumnumbb
      scavsumvol = scavsumvol + scavsumvolbb

   enddo jr_loop

   scavratenum = scavsumnum*3600._r8
   scavratevol = scavsumvol*3600._r8

   return
 end subroutine calc_1_impact_rate
  
  !=============================================================================
  !=============================================================================
  subroutine qqcw2vmr(lchnk, vmr, mbar, ncol, im, pbuf)
    use modal_aero_data, only : qqcw_get_field
    use physics_buffer, only : physics_buffer_desc
    !-----------------------------------------------------------------
    !	... Xfrom from mass to volume mixing ratio
    !-----------------------------------------------------------------

    use chem_mods, only : adv_mass, gas_pcnst

    implicit none

    !-----------------------------------------------------------------
    !	... Dummy args
    !-----------------------------------------------------------------
    integer, intent(in)     :: lchnk, ncol, im
    real(r8), intent(in)    :: mbar(ncol,pver)
    real(r8), intent(inout) :: vmr(ncol,pver,gas_pcnst)
    type(physics_buffer_desc), pointer :: pbuf(:)

    !-----------------------------------------------------------------
    !	... Local variables
    !-----------------------------------------------------------------
    integer :: k, m
    real(r8), pointer :: fldcw(:,:)

    do m=1,gas_pcnst
       if( adv_mass(m) /= 0._r8 ) then
          fldcw => qqcw_get_field(pbuf, m+im,lchnk,errorhandle=.true.)
          if(associated(fldcw)) then
             do k=1,pver
                vmr(:ncol,k,m) = mbar(:ncol,k) * fldcw(:ncol,k) / adv_mass(m)
             end do
          else
             vmr(:,:,m) = 0.0_r8
          end if
       end if
    end do
  end subroutine qqcw2vmr


  !=============================================================================
  !=============================================================================
  subroutine vmr2qqcw( lchnk, vmr, mbar, ncol, im, pbuf )
    !-----------------------------------------------------------------
    !	... Xfrom from volume to mass mixing ratio
    !-----------------------------------------------------------------

    use m_spc_id
    use chem_mods,       only : adv_mass, gas_pcnst
    use modal_aero_data, only : qqcw_get_field
    use physics_buffer,  only : physics_buffer_desc

    implicit none

    !-----------------------------------------------------------------
    !	... Dummy args
    !-----------------------------------------------------------------
    integer, intent(in)     :: lchnk, ncol, im
    real(r8), intent(in)    :: mbar(ncol,pver)
    real(r8), intent(in)    :: vmr(ncol,pver,gas_pcnst)
    type(physics_buffer_desc), pointer :: pbuf(:)

    !-----------------------------------------------------------------
    !	... Local variables
    !-----------------------------------------------------------------
    integer :: k, m
    real(r8), pointer :: fldcw(:,:)
    !-----------------------------------------------------------------
    !	... The non-group species
    !-----------------------------------------------------------------
    do m = 1,gas_pcnst
       fldcw => qqcw_get_field(pbuf, m+im,lchnk,errorhandle=.true.)
       if( adv_mass(m) /= 0._r8 .and. associated(fldcw)) then
          do k = 1,pver
             fldcw(:ncol,k) = adv_mass(m) * vmr(:ncol,k,m) / mbar(:ncol,k)
          end do
       end if
    end do

  end subroutine vmr2qqcw

end module aero_model
