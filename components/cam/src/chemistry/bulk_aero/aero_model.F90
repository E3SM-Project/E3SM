!===============================================================================
! Bulk Aerosol Model
!===============================================================================
module aero_model
  use shr_kind_mod,  only: r8 => shr_kind_r8
  use constituents,  only: pcnst, cnst_name, cnst_get_ind
  use ppgrid,        only: pcols, pver, pverp
  use cam_abortutils,    only: endrun
  use cam_logfile,   only: iulog
  use perf_mod,      only: t_startf, t_stopf
  use camsrfexch,    only: cam_in_t, cam_out_t
  use aerodep_flx,   only: aerodep_flx_prescribed
  use physics_types, only: physics_state, physics_ptend, physics_ptend_init
  use physics_buffer,only: physics_buffer_desc
  use physconst,     only: gravit, rair
  use dust_model,    only: dust_active, dust_names, dust_nbin
  use seasalt_model, only: sslt_active=>seasalt_active, seasalt_names, seasalt_nbin
  use spmd_utils,    only: masterproc
  use physics_buffer,only: pbuf_get_field, pbuf_get_index
  use cam_history,   only: outfld

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

  integer :: so4_ndx, cb2_ndx, oc2_ndx, nit_ndx
  integer :: soa_ndx, soai_ndx, soam_ndx, soab_ndx, soat_ndx, soax_ndx

  ! Namelist variables
  character(len=16) :: wetdep_list(pcnst) = ' '
  character(len=16) :: drydep_list(pcnst) = ' '

  integer :: ndrydep = 0
  integer,allocatable :: drydep_indices(:)
  integer :: nwetdep = 0
  integer,allocatable :: wetdep_indices(:)
  logical :: drydep_lq(pcnst)
  logical :: wetdep_lq(pcnst)

  integer :: fracis_idx = 0

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

    namelist /aerosol_nl/ aer_wetdep_list, aer_drydep_list

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
#endif

    wetdep_list = aer_wetdep_list
    drydep_list = aer_drydep_list

  end subroutine aero_model_readnl

  !=============================================================================
  !=============================================================================
  subroutine aero_model_register
    use mo_setsoa, only : soa_register

    call soa_register()
  end subroutine aero_model_register

  !=============================================================================
  !=============================================================================
  subroutine aero_model_init( pbuf2d )

    use mo_chem_utls,  only: get_inv_ndx, get_spc_ndx
    use cam_history,   only: addfld, horiz_only, add_default
    use phys_control,  only: phys_getopts
    use mo_aerosols,   only: aerosols_inti
    use mo_setsoa,     only: soa_inti
    use dust_model,    only: dust_init
    use seasalt_model, only: seasalt_init
    use drydep_mod,    only: inidrydep
    use wetdep,        only: wetdep_init

    ! args
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    ! local vars
    character(len=12), parameter :: subrname = 'aero_model_init'
    integer :: m, id
    character(len=20) :: dummy
    logical  :: history_aerosol ! Output MAM or SECT aerosol tendencies
    
    call phys_getopts( history_aerosol_out=history_aerosol   )
    call aerosols_inti()
    call soa_inti(pbuf2d)
    call dust_init()
    call seasalt_init()
    call wetdep_init()

    fracis_idx = pbuf_get_index('FRACIS') 

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

    do m = 1,ndrydep
       
       dummy = trim(drydep_list(m)) // 'TB'
       call addfld (dummy,horiz_only, 'A','kg/m2/s',trim(drydep_list(m))//' turbulent dry deposition flux')
       if ( history_aerosol ) then  
          call add_default (dummy, 1, ' ')
       endif
       dummy = trim(drydep_list(m))  // 'GV'
       call addfld (dummy,horiz_only, 'A','kg/m2/s',trim(drydep_list(m)) //' gravitational dry deposition flux')
       if ( history_aerosol ) then  
          call add_default (dummy, 1, ' ')
       endif
       dummy = trim(drydep_list(m))  // 'DD'
       call addfld (dummy,horiz_only, 'A','kg/m2/s',trim(drydep_list(m)) //' dry deposition flux at bottom (grav + turb)')
       if ( history_aerosol ) then  
          call add_default (dummy, 1, ' ')
       endif
       dummy = trim(drydep_list(m)) // 'DT'
       call addfld (dummy,(/ 'lev' /), 'A','kg/kg/s',trim(drydep_list(m))//' dry deposition')
       if ( history_aerosol ) then  
          call add_default (dummy, 1, ' ')
       endif
       dummy = trim(drydep_list(m)) // 'DV'
       call addfld (dummy,(/ 'lev' /), 'A','m/s',trim(drydep_list(m))//' deposition velocity')
       if ( history_aerosol ) then  
          call add_default (dummy, 1, ' ')
       endif

    enddo
    
    if (ndrydep>0) then

       call inidrydep(rair, gravit)

       dummy = 'RAM1'
       call addfld (dummy,horiz_only, 'A','frac','RAM1')
       if ( history_aerosol ) then  
          call add_default (dummy, 1, ' ')
       endif
       dummy = 'airFV'
       call addfld (dummy,horiz_only, 'A','frac','FV')
       if ( history_aerosol ) then  
          call add_default (dummy, 1, ' ')
       endif

       if (sslt_active) then
          dummy = 'SSTSFDRY'
          call addfld (dummy,horiz_only, 'A','kg/m2/s','Sea salt deposition flux at surface')
          if ( history_aerosol ) then  
             call add_default (dummy, 1, ' ')
          endif
       endif
       if (dust_active) then
          dummy = 'DSTSFDRY'
          call addfld (dummy,horiz_only, 'A','kg/m2/s','Dust deposition flux at surface')
          if ( history_aerosol ) then  
             call add_default (dummy, 1, ' ')
          endif
       endif

    endif

    do m = 1,nwetdep

       call addfld (trim(wetdep_list(m))//'SFWET', &
            horiz_only,  'A','kg/m2/s','Wet deposition flux at surface')
       call addfld (trim(wetdep_list(m))//'SFSIC', &
            horiz_only,  'A','kg/m2/s','Wet deposition flux (incloud, convective) at surface')
       call addfld (trim(wetdep_list(m))//'SFSIS', &
            horiz_only,  'A','kg/m2/s','Wet deposition flux (incloud, stratiform) at surface')
       call addfld (trim(wetdep_list(m))//'SFSBC', &
            horiz_only,  'A','kg/m2/s','Wet deposition flux (belowcloud, convective) at surface')
       call addfld (trim(wetdep_list(m))//'SFSBS', &
            horiz_only,  'A','kg/m2/s','Wet deposition flux (belowcloud, stratiform) at surface')
       call addfld (trim(wetdep_list(m))//'WET', &
            (/ 'lev' /), 'A','kg/kg/s','wet deposition tendency')
       call addfld (trim(wetdep_list(m))//'SIC', &
            (/ 'lev' /), 'A','kg/kg/s', trim(wetdep_list(m))//' ic wet deposition')
       call addfld (trim(wetdep_list(m))//'SIS', &
            (/ 'lev' /), 'A','kg/kg/s', trim(wetdep_list(m))//' is wet deposition')
       call addfld (trim(wetdep_list(m))//'SBC', &
            (/ 'lev' /), 'A','kg/kg/s', trim(wetdep_list(m))//' bc wet deposition')
       call addfld (trim(wetdep_list(m))//'SBS', &
            (/ 'lev' /), 'A','kg/kg/s', trim(wetdep_list(m))//' bs wet deposition')
    enddo
    
    if (nwetdep>0) then
       if (sslt_active) then
          dummy = 'SSTSFWET'
          call addfld (dummy,horiz_only, 'A','kg/m2/s','Sea salt wet deposition flux at surface')
          if ( history_aerosol ) then  
             call add_default (dummy, 1, ' ')
          endif
       endif
       if (dust_active) then
          dummy = 'DSTSFWET'
          call addfld (dummy,horiz_only, 'A','kg/m2/s','Dust wet deposition flux at surface')
          if ( history_aerosol ) then  
             call add_default (dummy, 1, ' ')
          endif
       endif
    endif
    
    if (dust_active) then
       ! emissions diagnostics ....

       do m = 1, dust_nbin
          dummy = trim(dust_names(m)) // 'SF'
          call addfld (dummy,horiz_only, 'A','kg/m2/s',trim(dust_names(m))//' dust surface emission')
          if (history_aerosol) then
             call add_default (dummy, 1, ' ')
          endif
       enddo

       dummy = 'DSTSFMBL'
       call addfld (dummy,horiz_only, 'A','kg/m2/s','Mobilization flux at surface')
       if (history_aerosol) then
          call add_default (dummy, 1, ' ')
       endif

       dummy = 'LND_MBL'
       call addfld (dummy,horiz_only, 'A','frac','Soil erodibility factor')
       if (history_aerosol) then
          call add_default (dummy, 1, ' ')
       endif

    endif
    
    if (sslt_active) then

       dummy = 'SSTSFMBL'
       call addfld (dummy,horiz_only, 'A','kg/m2/s','Mobilization flux at surface')
       if (history_aerosol) then
          call add_default (dummy, 1, ' ')
       endif

       do m = 1, seasalt_nbin
          dummy = trim(seasalt_names(m)) // 'SF'
          call addfld (dummy,horiz_only, 'A','kg/m2/s',trim(seasalt_names(m))//' seasalt surface emission')
          if (history_aerosol) then
             call add_default (dummy, 1, ' ')
          endif
       enddo

    endif

    so4_ndx    = get_spc_ndx( 'SO4' )
    soa_ndx    = get_spc_ndx( 'SOA' )
    soai_ndx   = get_spc_ndx( 'SOAI' )
    soam_ndx   = get_spc_ndx( 'SOAM' )
    soab_ndx   = get_spc_ndx( 'SOAB' )
    soat_ndx   = get_spc_ndx( 'SOAT' )
    soax_ndx   = get_spc_ndx( 'SOAX' )
    cb2_ndx    = get_spc_ndx( 'CB2' )
    oc2_ndx    = get_spc_ndx( 'OC2' )
    nit_ndx    = get_spc_ndx( 'NH4NO3' )

  end subroutine aero_model_init

  !=============================================================================
  !=============================================================================
  subroutine aero_model_drydep  ( state, pbuf, obklen, ustar, cam_in, dt, cam_out, ptend )

    use dust_sediment_mod, only: dust_sediment_tend
    use drydep_mod,        only: d3ddflux, calcram
    use dust_model,        only: dust_depvel, dust_nbin, dust_names
    use seasalt_model,     only: sslt_depvel=>seasalt_depvel, sslt_nbin=>seasalt_nbin, sslt_names=>seasalt_names

    ! args 
    type(physics_state),    intent(in)    :: state     ! Physics state variables
    real(r8),               intent(in)    :: obklen(:)          
    real(r8),               intent(in)    :: ustar(:)  ! sfc fric vel
    type(cam_in_t), target, intent(in)    :: cam_in    ! import state
    real(r8),               intent(in)    :: dt        ! time step
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

     ! local decarations

    integer, parameter :: naero = sslt_nbin+dust_nbin
    integer, parameter :: begslt = 1
    integer, parameter :: endslt = sslt_nbin
    integer, parameter :: begdst = sslt_nbin+1
    integer, parameter :: enddst = sslt_nbin+dust_nbin

    integer :: ncol, lchnk 

    character(len=6) :: aeronames(naero) ! = (/ sslt_names, dust_names /)

    real(r8) :: vlc_trb(pcols,naero)    !Turbulent deposn velocity (m/s)
    real(r8) :: vlc_grv(pcols,pver,naero)  !grav deposn velocity (m/s)
    real(r8) :: vlc_dry(pcols,pver,naero)  !dry deposn velocity (m/s)

    real(r8) :: dep_trb(pcols)       !kg/m2/s
    real(r8) :: dep_grv(pcols)       !kg/m2/s (total of grav and trb)

    real(r8) :: tsflx_dst(pcols)
    real(r8) :: tsflx_slt(pcols)
    real(r8) :: pvaeros(pcols,pverp)    ! sedimentation velocity in Pa
    real(r8) :: sflx(pcols)

    real(r8) :: tvs(pcols,pver)
    real(r8) :: rho(pcols,pver)      ! air density in kg/m3

    integer :: m,mm, i, im
    
    if (ndrydep<1) return

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

    aeronames(:sslt_nbin)   = sslt_names(:)
    aeronames(sslt_nbin+1:) = dust_names(:)

    lchnk = state%lchnk
    ncol  = state%ncol

    tvs(:ncol,:) = state%t(:ncol,:)
    rho(:ncol,:) = state%pmid(:ncol,:)/(rair*state%t(:ncol,:))

    ! compute dep velocities for sea salt and dust...
    if (sslt_active) then
       call sslt_depvel( state%t(:,:), state%pmid(:,:), state%q(:,:,1), ram1, fv, ncol, lchnk, &
                         vlc_dry(:,:,begslt:endslt), vlc_trb(:,begslt:endslt), vlc_grv(:,:,begslt:endslt))
    endif
    if (dust_active) then
       call dust_depvel( state%t(:,:), state%pmid(:,:),                 ram1, fv, ncol, &
                         vlc_dry(:,:,begdst:enddst), vlc_trb(:,begdst:enddst), vlc_grv(:,:,begdst:enddst) )
    endif

    tsflx_dst(:)=0._r8
    tsflx_slt(:)=0._r8

    ! do drydep for each of the bins of dust and seasalt
    do m=1,ndrydep

       mm = drydep_indices(m)
       findindex: do im = 1,naero
         if (trim(cnst_name(mm))==trim(aeronames(im))) exit findindex
       enddo findindex

       pvaeros(:ncol,1)=0._r8
       pvaeros(:ncol,2:pverp) = vlc_dry(:ncol,:,im)

       call outfld( trim(cnst_name(mm))//'DV', pvaeros(:,2:pverp), pcols, lchnk )

       if(.true.) then ! use phil's method
          !      convert from meters/sec to pascals/sec
          !      pvaeros(:,1) is assumed zero, use density from layer above in conversion
          pvaeros(:ncol,2:pverp) = pvaeros(:ncol,2:pverp) * rho(:ncol,:)*gravit        

          !      calculate the tendencies and sfc fluxes from the above velocities
          call dust_sediment_tend( &
               ncol,             dt,       state%pint(:,:), state%pmid, state%pdel, state%t , &
               state%q(:,:,mm) , pvaeros  , ptend%q(:,:,mm), sflx  )
       else   !use charlie's method
          call d3ddflux(ncol, vlc_dry(:,:,im), state%q(:,:,mm),state%pmid,state%pdel, tvs,sflx,ptend%q(:,:,mm),dt)
       endif
       ! apportion dry deposition into turb and gravitational settling for tapes
       do i=1,ncol
          dep_trb(i)=sflx(i)*vlc_trb(i,im)/vlc_dry(i,pver,im)
          dep_grv(i)=sflx(i)*vlc_grv(i,pver,im)/vlc_dry(i,pver,im)
       enddo

       if ( any( sslt_names(:)==trim(cnst_name(mm)) ) ) &
            tsflx_slt(:ncol)=tsflx_slt(:ncol)+sflx(:ncol)
       if ( any( dust_names(:)==trim(cnst_name(mm)) ) ) &
            tsflx_dst(:ncol)=tsflx_dst(:ncol)+sflx(:ncol)

       ! if the user has specified prescribed aerosol dep fluxes then 
       ! do not set cam_out dep fluxes according to the prognostic aerosols
       if (.not. aerodep_flx_prescribed()) then
          ! set deposition in export state
          if (im==begdst) then
             cam_out%dstdry1(:ncol) = max(sflx(:ncol), 0._r8)
          elseif(im==begdst+1) then
             cam_out%dstdry2(:ncol) = max(sflx(:ncol), 0._r8)
          elseif(im==begdst+2) then
             cam_out%dstdry3(:ncol) = max(sflx(:ncol), 0._r8)
          elseif(im==begdst+3) then
             cam_out%dstdry4(:ncol) = max(sflx(:ncol), 0._r8)
          endif
       endif

       call outfld( trim(cnst_name(mm))//'DD', sflx, pcols, lchnk)
       call outfld( trim(cnst_name(mm))//'TB', dep_trb, pcols, lchnk )
       call outfld( trim(cnst_name(mm))//'GV', dep_grv, pcols, lchnk )
       call outfld( trim(cnst_name(mm))//'DT', ptend%q(:,:,mm), pcols, lchnk)

    end do
    
    ! output the total dry deposition
    if (sslt_active) then
       call outfld( 'SSTSFDRY', tsflx_slt, pcols, lchnk)
    endif
    if (dust_active) then
       call outfld( 'DSTSFDRY', tsflx_dst, pcols, lchnk)
    endif

  endsubroutine aero_model_drydep

  !=============================================================================
  !=============================================================================
  subroutine aero_model_wetdep( dt, dlf, dlf2, cmfmc2, state, sh_e_ed_ratio,    & !Intent-ins
       mu, md, du, eu, ed, dp, dsubcld, jt, maxg, ideep, lengath, species_class,&
       cam_out,                                                                 & !Intent-inout
       pbuf,                                                                    & !Pointer
       ptend                                                                    ) !Intent-out

    use wetdep,        only : wetdepa_v1, wetdep_inputs_set, wetdep_inputs_t
    use dust_model,    only : dust_names
    use seasalt_model, only : sslt_names=>seasalt_names

    ! args

    type(physics_state), intent(in)    :: state       ! Physics state variables
    real(r8),            intent(in)    :: dt          ! time step
    real(r8),            intent(in)    :: dlf(:,:)    ! shallow+deep convective detrainment [kg/kg/s]
    real(r8),            intent(in)    :: dlf2(:,:)   ! Shal conv cldwtr detrainment (kg/kg/s - grid avg)
    real(r8),            intent(in)    :: cmfmc2(pcols,pverp) ! Shal conv mass flux (kg/m2/s)
    real(r8),            intent(in)    :: sh_e_ed_ratio(pcols,pver)  ! shallow conv [ent/(ent+det)] ratio
                                                ! mu, md, ..., ideep, lengath are all deep conv variables
                                                ! *** AND ARE GATHERED ***
    real(r8),            intent(in)    :: mu(pcols,pver)   ! Updraft mass flux (positive)
    real(r8),            intent(in)    :: md(pcols,pver)   ! Downdraft mass flux (negative)
    real(r8),            intent(in)    :: du(pcols,pver)   ! Mass detrain rate from updraft
    real(r8),            intent(in)    :: eu(pcols,pver)   ! Mass entrain rate into updraft
    real(r8),            intent(in)    :: ed(pcols,pver)   ! Mass entrain rate into downdraft
    ! eu, ed, du are "d(massflux)/dp" and are all positive
    real(r8),            intent(in)    :: dp(pcols,pver)   ! Delta pressure between interfaces
    real(r8),            intent(in)    :: dsubcld(pcols)   ! Delta pressure from cloud base to sfc

    integer,             intent(in)    :: jt(pcols)         ! Index of cloud top for each column
    integer,             intent(in)    :: maxg(pcols)       ! Index of cloud top for each column
    integer,             intent(in)    :: ideep(pcols)      ! Gathering array
    integer,             intent(in)    :: lengath           ! Gathered min lon indices over which to operate
    integer,             intent(in)    :: species_class(:)

    type(cam_out_t),     intent(inout) :: cam_out     ! export state
    type(physics_buffer_desc), pointer :: pbuf(:)

    type(physics_ptend), intent(out)   :: ptend       ! indivdual parameterization tendencies

    ! local vars

    integer  :: ncol                     ! number of atmospheric columns
    integer  :: lchnk                    ! chunk identifier
    integer  :: m,mm, i,k

    real(r8) :: sflx_tot_dst(pcols)
    real(r8) :: sflx_tot_slt(pcols)

    real(r8) :: iscavt(pcols, pver)
    real(r8) :: scavt(pcols, pver)
    real(r8) :: scavcoef(pcols,pver)     ! Dana and Hales coefficient (/mm) (0.1)
    real(r8) :: sflx(pcols)              ! deposition flux

    real(r8) :: icscavt(pcols, pver)
    real(r8) :: isscavt(pcols, pver)
    real(r8) :: bcscavt(pcols, pver)
    real(r8) :: bsscavt(pcols, pver)

    real(r8) :: sol_factb, sol_facti

    real(r8) :: rainmr(pcols,pver)       ! mixing ratio of rain within cloud volume
    real(r8) :: cldv(pcols,pver)         ! cloudy volume undergoing scavenging
    real(r8) :: cldvcu(pcols,pver)       ! Convective precipitation area at the top interface of current layer
    real(r8) :: cldvst(pcols,pver)       ! Stratiform precipitation area at the top interface of current layer
 
    real(r8), pointer :: fracis(:,:,:)   ! fraction of transported species that are insoluble

    type(wetdep_inputs_t) :: dep_inputs  ! obj that contains inputs to wetdepa routine

    if (nwetdep<1) return

    call pbuf_get_field(pbuf, fracis_idx, fracis, start=(/1,1,1/), kount=(/pcols, pver, pcnst/) )

    call physics_ptend_init(ptend, state%psetcols, 'aero_model_wetdep', lq=wetdep_lq)

    call wetdep_inputs_set( state, pbuf, dep_inputs )

    lchnk = state%lchnk
    ncol  = state%ncol

    sflx_tot_dst(:) = 0._r8
    sflx_tot_slt(:) = 0._r8

    do m = 1, nwetdep

       mm = wetdep_indices(m)

       if ( any( dust_names(:)==trim(cnst_name(mm)) ) ) then
          sol_factb = 0.15_r8
          sol_facti = 0.15_r8
       else
          sol_factb = 0.3_r8
          sol_facti = 0.3_r8
       endif

       scavcoef(:ncol,:)=0.1_r8

       call wetdepa_v1( state%t, state%pmid, state%q(:,:,1), state%pdel, &
            dep_inputs%cldt, dep_inputs%cldcu, dep_inputs%cmfdqr, &
            dep_inputs%conicw, dep_inputs%prain, dep_inputs%qme, &
            dep_inputs%evapr, dep_inputs%totcond, state%q(:,:,mm), dt, &
            scavt, iscavt, dep_inputs%cldv, &
            fracis(:,:,mm), sol_factb, ncol, &
            scavcoef, &
            sol_facti_in=sol_facti, &
            icscavt=icscavt, isscavt=isscavt, bcscavt=bcscavt, bsscavt=bsscavt )

       ptend%q(:ncol,:,mm)=scavt(:ncol,:)

       call outfld( trim(cnst_name(mm))//'WET', ptend%q(:,:,mm), pcols, lchnk)
       call outfld( trim(cnst_name(mm))//'SIC', icscavt , pcols, lchnk)
       call outfld( trim(cnst_name(mm))//'SIS', isscavt, pcols, lchnk)
       call outfld( trim(cnst_name(mm))//'SBC', bcscavt, pcols, lchnk)
       call outfld( trim(cnst_name(mm))//'SBS', bsscavt, pcols, lchnk)

       sflx(:)=0._r8

       do k=1,pver
          do i=1,ncol
             sflx(i)=sflx(i)+ptend%q(i,k,mm)*state%pdel(i,k)/gravit
          enddo
       enddo
       call outfld( trim(cnst_name(mm))//'SFWET', sflx, pcols, lchnk)
       
       if ( any( sslt_names(:)==trim(cnst_name(mm)) ) ) &
            sflx_tot_slt(:ncol) = sflx_tot_slt(:ncol) + sflx(:ncol)
       if ( any( dust_names(:)==trim(cnst_name(mm)) ) ) &
            sflx_tot_dst(:ncol) = sflx_tot_dst(:ncol) + sflx(:ncol)

       ! if the user has specified prescribed aerosol dep fluxes then 
       ! do not set cam_out dep fluxes according to the prognostic aerosols
       if (.not.aerodep_flx_prescribed()) then
          ! export deposition fluxes to coupler ??? why "-" sign ???
          if (trim(cnst_name(mm))=='CB2') then
             cam_out%bcphiwet(:) = max(-sflx(:), 0._r8)
          elseif (trim(cnst_name(mm))=='OC2') then
             cam_out%ocphiwet(:) = max(-sflx(:), 0._r8)
          elseif (trim(cnst_name(mm))==trim(dust_names(1))) then
             cam_out%dstwet1(:) = max(-sflx(:), 0._r8)
          elseif (trim(cnst_name(mm))==trim(dust_names(2))) then
             cam_out%dstwet2(:) = max(-sflx(:), 0._r8)
          elseif (trim(cnst_name(mm))==trim(dust_names(3))) then
             cam_out%dstwet3(:) = max(-sflx(:), 0._r8)
          elseif (trim(cnst_name(mm))==trim(dust_names(4))) then
             cam_out%dstwet4(:) = max(-sflx(:), 0._r8)
          endif
       endif

    enddo
    
    if (sslt_active) then
       call outfld( 'SSTSFWET', sflx_tot_slt, pcols, lchnk)
    endif
    if (dust_active) then
       call outfld( 'DSTSFWET', sflx_tot_dst, pcols, lchnk)
    endif

  endsubroutine aero_model_wetdep

  !-------------------------------------------------------------------------
  ! provides aerosol surface area info for sectional aerosols
  ! called from mo_usrrxt
  !-------------------------------------------------------------------------
  subroutine aero_model_surfarea( &
                  mmr, radmean, relhum, pmid, temp, strato_sad, &
                  sulfate,  m, ltrop, het1_ndx, pbuf, ncol, sfc, dm_aer, sad_total )

    use mo_constants, only : pi, avo => avogadro

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
    real(r8), intent(in)    :: m(:,:) ! total atm density (/cm^3)
    real(r8), intent(in)    :: sulfate(:,:)
    type(physics_buffer_desc), pointer :: pbuf(:)

    real(r8), intent(inout) :: sfc(:,:,:)
    real(r8), intent(inout) :: dm_aer(:,:,:)
    real(r8), intent(inout) :: sad_total(:,:)

    ! local vars

    integer  :: i,k
    real(r8) :: rho_air
    real(r8) :: v, n, n_exp, r_rd, r_sd
    real(r8) :: dm_sulf, dm_sulf_wet, log_sd_sulf, sfc_sulf, sfc_nit
    real(r8) :: dm_orgc, dm_orgc_wet, log_sd_orgc, sfc_oc, sfc_soa
    real(r8) :: sfc_soai, sfc_soam, sfc_soab, sfc_soat, sfc_soax
    real(r8) :: dm_bc, dm_bc_wet, log_sd_bc, sfc_bc
    real(r8) :: rxt_sulf, rxt_nit, rxt_oc, rxt_soa
    real(r8) :: c_n2o5, c_ho2, c_no2, c_no3
    real(r8) :: s_exp

    !-----------------------------------------------------------------
    ! 	... parameters for log-normal distribution by number
    ! references:
    !   Chin et al., JAS, 59, 461, 2003
    !   Liao et al., JGR, 108(D1), 4001, 2003
    !   Martin et al., JGR, 108(D3), 4097, 2003
    !-----------------------------------------------------------------
    real(r8), parameter :: rm_sulf  = 6.95e-6_r8        ! mean radius of sulfate particles (cm) (Chin)
    real(r8), parameter :: sd_sulf  = 2.03_r8           ! standard deviation of radius for sulfate (Chin)
    real(r8), parameter :: rho_sulf = 1.7e3_r8          ! density of sulfate aerosols (kg/m3) (Chin) 

    real(r8), parameter :: rm_orgc  = 2.12e-6_r8        ! mean radius of organic carbon particles (cm) (Chin)
    real(r8), parameter :: sd_orgc  = 2.20_r8           ! standard deviation of radius for OC (Chin)
    real(r8), parameter :: rho_orgc = 1.8e3_r8          ! density of OC aerosols (kg/m3) (Chin)

    real(r8), parameter :: rm_bc    = 1.18e-6_r8        ! mean radius of soot/BC particles (cm) (Chin)
    real(r8), parameter :: sd_bc    = 2.00_r8           ! standard deviation of radius for BC (Chin)
    real(r8), parameter :: rho_bc   = 1.0e3_r8          ! density of BC aerosols (kg/m3) (Chin)

    real(r8), parameter :: mw_so4 = 98.e-3_r8     ! so4 molecular wt (kg/mole)

    integer  ::  irh, rh_l, rh_u
    real(r8) ::  factor, rfac_sulf, rfac_oc, rfac_bc, rfac_ss

    !-----------------------------------------------------------------
    ! 	... table for hygroscopic growth effect on radius (Chin et al)
    !           (no growth effect for mineral dust)
    !-----------------------------------------------------------------
    real(r8), dimension(7) :: table_rh, table_rfac_sulf, table_rfac_bc, table_rfac_oc, table_rfac_ss

    data table_rh(1:7)        / 0.0_r8, 0.5_r8, 0.7_r8, 0.8_r8, 0.9_r8, 0.95_r8, 0.99_r8/
    data table_rfac_sulf(1:7) / 1.0_r8, 1.4_r8, 1.5_r8, 1.6_r8, 1.8_r8, 1.9_r8,  2.2_r8/
    data table_rfac_oc(1:7)   / 1.0_r8, 1.2_r8, 1.4_r8, 1.5_r8, 1.6_r8, 1.8_r8,  2.2_r8/
    data table_rfac_bc(1:7)   / 1.0_r8, 1.0_r8, 1.0_r8, 1.2_r8, 1.4_r8, 1.5_r8,  1.9_r8/
    data table_rfac_ss(1:7)   / 1.0_r8, 1.6_r8, 1.8_r8, 2.0_r8, 2.4_r8, 2.9_r8,  4.8_r8/

    !-----------------------------------------------------------------
    ! 	... exponent for calculating number density
    !-----------------------------------------------------------------
    n_exp = exp( -4.5_r8*log(sd_sulf)*log(sd_sulf) )

    dm_sulf = 2._r8 * rm_sulf
    dm_orgc = 2._r8 * rm_orgc
    dm_bc   = 2._r8 * rm_bc

    log_sd_sulf = log(sd_sulf)
    log_sd_orgc = log(sd_orgc)
    log_sd_bc   = log(sd_bc)

    ver_loop: do k = 1,pver
       col_loop: do i = 1,ncol
          !-------------------------------------------------------------------------
          ! 	... air density (kg/m3)
          !-------------------------------------------------------------------------
          rho_air = pmid(i,k)/(temp(i,k)*287.04_r8)
          !-------------------------------------------------------------------------
          !       ... aerosol growth interpolated from M.Chin's table
          !-------------------------------------------------------------------------
          if (relhum(i,k) >= table_rh(7)) then
             rfac_sulf = table_rfac_sulf(7)
             rfac_oc = table_rfac_oc(7)
             rfac_bc = table_rfac_bc(7)
          else
             do irh = 2,7
                if (relhum(i,k) <= table_rh(irh)) then
                   exit
                end if
             end do
             rh_l = irh-1
             rh_u = irh

             factor = (relhum(i,k) - table_rh(rh_l))/(table_rh(rh_u) - table_rh(rh_l))

             rfac_sulf = table_rfac_sulf(rh_l) + factor*(table_rfac_sulf(rh_u) - table_rfac_sulf(rh_l))
             rfac_oc = table_rfac_oc(rh_u) + factor*(table_rfac_oc(rh_u) - table_rfac_oc(rh_l))
             rfac_bc = table_rfac_bc(rh_u) + factor*(table_rfac_bc(rh_u) - table_rfac_bc(rh_l))
          end if

          dm_sulf_wet = dm_sulf * rfac_sulf
          dm_orgc_wet = dm_orgc * rfac_oc
          dm_bc_wet = dm_bc * rfac_bc

          dm_bc_wet   = min(dm_bc_wet  ,50.e-6_r8) ! maximum size is 0.5 micron (Chin)
          dm_orgc_wet = min(dm_orgc_wet,50.e-6_r8) ! maximum size is 0.5 micron (Chin)


          !-------------------------------------------------------------------------
          ! 	... sulfate aerosols
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !       ... use ubvals climatology for stratospheric sulfate surface area density
          !-------------------------------------------------------------------------
          if( k < ltrop(i) ) then
             sfc_sulf = strato_sad(i,k)
             if ( het1_ndx > 0 ) then
                sfc_sulf = 0._r8        ! reaction already taken into account in mo_strato_rates.F90
             end if
          else

             if( so4_ndx > 0 ) then
                !-------------------------------------------------------------------------
                ! convert mass mixing ratio of aerosol to cm3/cm3 (cm^3_aerosol/cm^3_air)
                ! v=volume density (m^3/m^3)
                ! rho_aer=density of aerosol (kg/m^3)
                ! v=m*rho_air/rho_aer   [kg/kg * (kg/m3)_air/(kg/m3)_aer]
                !-------------------------------------------------------------------------
                v = mmr(i,k,so4_ndx) * rho_air/rho_sulf
                !-------------------------------------------------------------------------
                ! calculate the number density of aerosol (aerosols/cm3)
                ! assuming a lognormal distribution
                ! n  = (aerosols/cm3)
                ! dm = geometric mean diameter
                !
                ! because only the dry mass of the aerosols is known, we
                ! use the mean dry radius
                !-------------------------------------------------------------------------
                n  = v * (6._r8/pi)*(1._r8/(dm_sulf**3._r8))*n_exp
                !-------------------------------------------------------------------------
                ! find surface area of aerosols using dm_wet, log_sd 
                !  (increase of sd due to RH is negligible)
                ! and number density calculated above as distribution
                ! parameters
                ! sfc = surface area of wet aerosols (cm^2/cm^3)
                !-------------------------------------------------------------------------
                s_exp    = exp(2._r8*log_sd_sulf*log_sd_sulf)
                sfc_sulf = n * pi * (dm_sulf_wet**2._r8) * s_exp

             else
                !-------------------------------------------------------------------------
                !  if so4 not simulated, use off-line sulfate and calculate as above
                !  convert sulfate vmr to volume density of aerosol (cm^3_aerosol/cm^3_air)           
                !-------------------------------------------------------------------------
                v = sulfate(i,k) * m(i,k) * mw_so4 / (avo * rho_sulf) *1.e6_r8
                n  = v * (6._r8/pi)*(1._r8/(dm_sulf**3._r8))*n_exp
                s_exp    = exp(2._r8*log_sd_sulf*log_sd_sulf)
                sfc_sulf = n * pi * (dm_sulf_wet**2._r8) * s_exp

             end if
          end if

          !-------------------------------------------------------------------------
          ! ammonium nitrate (follow same procedure as sulfate, using size and density of sulfate)
          !-------------------------------------------------------------------------
          if( nit_ndx > 0 ) then
             v = mmr(i,k,nit_ndx) * rho_air/rho_sulf
             n  = v * (6._r8/pi)*(1._r8/(dm_sulf**3._r8))*n_exp
             s_exp   = exp(2._r8*log_sd_sulf*log_sd_sulf)
             sfc_nit = n * pi * (dm_sulf_wet**2._r8) * s_exp
          else
             sfc_nit = 0._r8
          end if

          !-------------------------------------------------------------------------
          ! hydrophylic organic carbon (follow same procedure as sulfate)
          !-------------------------------------------------------------------------
          if( oc2_ndx > 0 ) then
             v = mmr(i,k,oc2_ndx) * rho_air/rho_orgc
             n  = v * (6._r8/pi)*(1._r8/(dm_orgc**3))*n_exp
             s_exp    = exp(2._r8*log_sd_orgc*log_sd_orgc)
             sfc_oc   = n * pi * (dm_orgc_wet**2._r8) * s_exp
          else
             sfc_oc = 0._r8
          end if

          !-------------------------------------------------------------------------
          ! secondary organic carbon (follow same procedure as sulfate)
          !-------------------------------------------------------------------------
          if( soa_ndx > 0 ) then
             v = mmr(i,k,soa_ndx) * rho_air/rho_orgc
             n  = v * (6._r8/pi)*(1._r8/(dm_orgc**3._r8))*n_exp
             s_exp     = exp(2._r8*log_sd_orgc*log_sd_orgc)
             sfc_soa   = n * pi * (dm_orgc_wet**2._r8) * s_exp
          else
             sfc_soa = 0._r8
          end if
          if( soai_ndx > 0 ) then
             v = mmr(i,k,soai_ndx) * rho_air/rho_orgc
             n  = v * (6._r8/pi)*(1._r8/(dm_orgc**3._r8))*n_exp
             s_exp     = exp(2._r8*log_sd_orgc*log_sd_orgc)
             sfc_soai   = n * pi * (dm_orgc_wet**2._r8) * s_exp
          else
             sfc_soai = 0._r8
          end if
          if( soam_ndx > 0 ) then
             v = mmr(i,k,soam_ndx) * rho_air/rho_orgc
             n  = v * (6._r8/pi)*(1._r8/(dm_orgc**3._r8))*n_exp
             s_exp     = exp(2._r8*log_sd_orgc*log_sd_orgc)
             sfc_soam   = n * pi * (dm_orgc_wet**2._r8) * s_exp
          else
             sfc_soam = 0._r8
          end if
          if( soab_ndx > 0 ) then
             v = mmr(i,k,soab_ndx) * rho_air/rho_orgc
             n  = v * (6._r8/pi)*(1._r8/(dm_orgc**3._r8))*n_exp
             s_exp     = exp(2._r8*log_sd_orgc*log_sd_orgc)
             sfc_soab   = n * pi * (dm_orgc_wet**2._r8) * s_exp
          else
             sfc_soab = 0._r8
          end if
          if( soat_ndx > 0 ) then
             v = mmr(i,k,soat_ndx) * rho_air/rho_orgc
             n  = v * (6._r8/pi)*(1._r8/(dm_orgc**3._r8))*n_exp
             s_exp     = exp(2._r8*log_sd_orgc*log_sd_orgc)
             sfc_soat   = n * pi * (dm_orgc_wet**2._r8) * s_exp
          else
             sfc_soat = 0._r8
          end if
          if( soax_ndx > 0 ) then
             v = mmr(i,k,soax_ndx) * rho_air/rho_orgc
             n  = v * (6._r8/pi)*(1._r8/(dm_orgc**3._r8))*n_exp
             s_exp     = exp(2._r8*log_sd_orgc*log_sd_orgc)
             sfc_soax   = n * pi * (dm_orgc_wet**2._r8) * s_exp
          else
             sfc_soax = 0._r8
          end if
          sfc_soa = sfc_soa + sfc_soai + sfc_soam + sfc_soab + sfc_soat + sfc_soax 

          !-------------------------------------------------------------------------
          ! black carbon (follow same procedure as sulfate)
          !-------------------------------------------------------------------------
          if( cb2_ndx > 0 ) then
             v = mmr(i,k,cb2_ndx) * rho_air/rho_bc
             n  = v * (6._r8/pi)*(1._r8/(dm_bc**3._r8))*n_exp
             s_exp     = exp(2._r8*log_sd_bc*log_sd_bc)
             sfc_bc   = n * pi * (dm_bc_wet**2._r8) * s_exp
          else
             sfc_bc = 0._r8
          end if

          sfc(i,k,:) = (/ sfc_sulf, sfc_nit, sfc_oc, sfc_soa, sfc_bc /)
          dm_aer(i,k,:) = (/ dm_sulf_wet,dm_sulf_wet,dm_orgc_wet,dm_orgc_wet,dm_bc_wet /)

          !-------------------------------------------------------------------------
          !  	... add up total surface area density for output
          !-------------------------------------------------------------------------
          sad_total(i,k) = sfc_sulf + sfc_nit + sfc_oc + sfc_soa + sfc_bc

       enddo col_loop
    enddo ver_loop

  end subroutine aero_model_surfarea

  !=============================================================================
  !=============================================================================
  subroutine aero_model_gasaerexch( loffset, ncol, lchnk, delt, latndx, lonndx, reaction_rates, &
                                    tfld, pmid, pdel, mbar, relhum, &
                                    zm,  qh2o, cwat, cldfr, cldnum, &
                                    airdens, invariants, del_h2so4_gasprod,  &
                                    vmr0, vmr, pbuf )

    use chem_mods,   only : gas_pcnst
    use mo_aerosols, only : aerosols_formation, has_aerosols
    use mo_setsox,   only : setsox, has_sox
    use mo_setsoa,   only : setsoa, has_soa

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
    
    ! These are declared here so that arguments remain consistent with modal aerosol version.
    integer, intent(in)  ::  latndx(pcols)                         ! chunk lat indicies
    integer, intent(in)  ::  lonndx(pcols)                         ! chunk lon indicies

    ! local vars 

    real(r8) :: vmrcw(ncol,pver,gas_pcnst)            ! cloud-borne aerosol (vmr)

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

    if( has_soa ) then
       call setsoa( ncol, lchnk, delt, reaction_rates, tfld, airdens, vmr, pbuf)
    endif

    if( has_aerosols ) then
       call aerosols_formation( ncol, lchnk, tfld, relhum, vmr )
    endif

  end subroutine aero_model_gasaerexch

  !=============================================================================
  !=============================================================================
  subroutine aero_model_emissions( state, cam_in )
    use seasalt_model, only: seasalt_emis, seasalt_indices
    use dust_model,    only: dust_emis, dust_indices
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
       do m=1,dust_nbin
          mm = dust_indices(m)
          sflx(:ncol)=sflx(:ncol)+cam_in%cflx(:ncol,mm)
          call outfld(trim(dust_names(m))//'SF',cam_in%cflx(:,mm),pcols, lchnk)
       enddo
       call outfld('DSTSFMBL',sflx(:),pcols,lchnk)
       call outfld('LND_MBL',soil_erod_tmp(:),pcols, lchnk )
    endif

    if (sslt_active) then
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

end module aero_model
