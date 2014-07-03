
      module mo_waccm_hrates

      use shr_kind_mod,      only : r8 => shr_kind_r8
      use cam_logfile,       only : iulog

      implicit none

      save

      real(r8), parameter :: secpday       = 86400._r8
      real(r8), parameter :: daypsec       = 1._r8/secpday
      real(r8), parameter :: aur_therm     = 807._r8
      real(r8), parameter :: jkcal         = 4184._r8
      real(r8), parameter :: aur_heat_eff  = .05_r8
      real(r8), parameter :: aur_hconst    = 1.e3_r8*jkcal*aur_therm*aur_heat_eff

      real(r8) :: max_zen_angle

      private
      public :: waccm_hrates, init_hrates, has_hrates

      integer :: id_co2, id_o2, id_o3, id_o2_1d, id_o2_1s, id_h2o, id_o, id_h
      logical :: has_hrates

      contains
   
      subroutine init_hrates( )
        use mo_chem_utls, only : get_spc_ndx
        use cam_history,  only : addfld, phys_decomp
        use ppgrid,       only : pver
        use ref_pres,     only : ptop_ref, psurf_ref


        implicit none

        integer :: ids(8)
        character(len=128) :: attr  ! netcdf variable attribute

        id_co2   = get_spc_ndx( 'CO2' )
        id_o2    = get_spc_ndx( 'O2' )
        id_o3    = get_spc_ndx( 'O3' )
        id_o2_1d = get_spc_ndx( 'O2_1D' )
        id_o2_1s = get_spc_ndx( 'O2_1S' )
        id_h2o   = get_spc_ndx( 'H2O' )
        id_o     = get_spc_ndx( 'O' )
        id_h     = get_spc_ndx( 'H' )

        ids = (/ id_co2, id_o2, id_o3, id_o2_1d, id_o2_1s, id_h2o, id_o, id_h /)

        has_hrates = all( ids(:) > 0 ) .and. ptop_ref < 0.0004_r8 * psurf_ref

        if (.not. has_hrates) return

        call addfld( 'CPAIR', 'J/K/kg', pver, 'I', 'specific heat cap air', phys_decomp )
        call addfld( 'QRS_AUR', 'K/s', pver, 'I', 'total auroral heating rate', phys_decomp )
        call addfld( 'QRS_CO2NIR', 'K/s', pver, 'I', 'co2 nir heating rate', phys_decomp )
        call addfld( 'QTHERMAL', 'K/s', pver, 'I', 'non-euv photolysis heating rate', phys_decomp )
        call addfld( 'QRS_MLT', 'K/s', pver, 'I', 'Total heating rate (unmerged with tropospheric RT heating)', phys_decomp )

        attr = 'O2 + hv -> O1D + O3P solar heating rate < 200nm'
        call addfld( 'QRS_SO2A', 'K/s ', pver, 'I', trim(attr), phys_decomp )
        attr = 'O2 + hv -> O3P + O3P solar heating rate < 200nm'
        call addfld( 'QRS_SO2B', 'K/s ', pver, 'I', trim(attr), phys_decomp )
        attr = 'O3 + hv -> O1D + O2_1S solar heating rate < 200nm'
        call addfld( 'QRS_SO3A', 'K/s ', pver, 'I', trim(attr), phys_decomp )
        attr = 'O3 + hv -> O3P + O2 solar heating rate < 200nm'
        call addfld( 'QRS_SO3B', 'K/s ', pver, 'I', trim(attr), phys_decomp )
        attr = 'O2 + hv -> 2*O3P solar heating rate > 200nm'
        call addfld( 'QRS_LO2B', 'K/s ', pver, 'I', trim(attr), phys_decomp )
        attr = 'O3 + hv -> O1D + O2_1S solar heating rate > 200nm'
        call addfld( 'QRS_LO3A', 'K/s ', pver, 'I', trim(attr), phys_decomp )
        attr = 'O3 + hv -> O3P + O2 solar heating rate > 200nm'
        call addfld( 'QRS_LO3B', 'K/s ', pver, 'I', trim(attr), phys_decomp )
        attr = 'Total O3 solar heating > 200nm'
        call addfld( 'QRS_LO3',  'K/s ', pver, 'I', trim(attr), phys_decomp )
        attr = 'total euv heating rate'
        call addfld( 'QRS_EUV', 'K/s', pver, 'I', trim(attr), phys_decomp )
        attr = 'total jo2 euv photolysis rate'
        call addfld( 'JO2_EUV', '/s', pver, 'I', trim(attr), phys_decomp )

      end subroutine init_hrates

      subroutine waccm_hrates(ncol, state, asdir, bot_mlt_lev, qrs_tot, pbuf )
!-----------------------------------------------------------------------
!     ... computes the short wavelength heating rates
!-----------------------------------------------------------------------

      use chem_mods,         only : nabscol, nfs, gas_pcnst, rxntot, indexm
      use ppgrid,            only : pcols, pver
      use physconst,         only : rga, mbarv, cpairv
      use constituents,      only : pcnst
      use mo_gas_phase_chemdr,only: map2chm
      use mo_photo,          only : set_ub_col, setcol
      use mo_jlong,          only : jlong
      use mo_jshort,         only : jshort
      use mo_jeuv,           only : heuv
      use mo_cph,            only : cph
      use mo_heatnirco2,     only : heatnirco2
      use mo_airglow,        only : airglow
      use mo_aurora,         only : aurora
      use mo_setrxt,         only : setrxt_hrates
      use mo_adjrxt,         only : adjrxt
      use mo_usrrxt,         only : usrrxt_hrates
      use mo_setinv,         only : setinv
      use mo_mass_xforms,    only : mmr2vmr
      use physics_types,     only : physics_state
      use phys_grid,         only : get_rlat_all_p, get_rlon_all_p, &
                                    get_lat_all_p, get_lon_all_p
      use mo_mean_mass,      only : set_mean_mass
      use set_cp,            only : calc_cp
      use cam_history,       only : outfld
      use shr_orb_mod,       only : shr_orb_decl
      use time_manager,      only : get_curr_calday
      use cam_control_mod,   only : lambm0, eccen, mvelpp, obliqr
      use mo_constants,      only : r2d
      use short_lived_species,only: get_short_lived_species
      use physics_buffer,    only : physics_buffer_desc
      use phys_control,      only : waccmx_is

!-----------------------------------------------------------------------
!        ... dummy arguments
!-----------------------------------------------------------------------
      integer,             intent(in)  ::  ncol                  ! number columns in chunk
      type(physics_state), intent(in)  ::  state                 ! physics state structure
      real(r8),            intent(in)  ::  asdir(pcols)          ! shortwave, direct albedo
      integer,             intent(in)  ::  bot_mlt_lev           ! lowest model level where MLT heating is needed
      real(r8),            intent(out) ::  qrs_tot(pcols,pver)   ! total heating (K/s)
      type(physics_buffer_desc), pointer :: pbuf(:)

!-----------------------------------------------------------------------
!     	... local variables
!-----------------------------------------------------------------------
      integer             :: lchnk                 ! chunk index
      real(r8), parameter :: m2km  = 1.e-3_r8
      real(r8), parameter :: Pa2mb = 1.e-2_r8

      integer      ::  i, k, m, n
      integer      ::  kbot_hrates
      real(r8)     ::  esfact
      real(r8)     ::  sza                                           ! solar zenith angle (degrees)
      integer      ::  latndx(pcols)                                 ! chunk lat indicies
      integer      ::  lonndx(pcols)                                 ! chunk lon indicies
      real(r8)     ::  invariants(ncol,pver,nfs)
      real(r8)     ::  col_dens(ncol,pver,nabscol)                   ! column densities (molecules/cm^2)
      real(r8)     ::  col_delta(ncol,0:pver,nabscol)                ! layer column densities (molecules/cm^2)
      real(r8)     ::  vmr(ncol,pver,gas_pcnst)                      ! xported species (vmr)
      real(r8)     ::  reaction_rates(ncol,pver,rxntot)              ! reaction rates
      real(r8)     ::  mmr(pcols,pver,gas_pcnst)                     ! chem working concentrations (kg/kg)
      real(r8)     ::  h2ovmr(ncol,pver)                             ! water vapor concentration (mol/mol)
      real(r8)     ::  mbar(ncol,pver)                               ! mean wet atmospheric mass (kg/mole)
      real(r8)     ::  zmid(ncol,pver)                               ! midpoint geopotential (km)
      real(r8)     ::  cpair(ncol,pver)                              ! specific heat capacity (J/K/kg)
      real(r8)     ::  cphrate(ncol,pver)                            ! chemical pot heat rate (K/s)
      real(r8)     ::  aghrate(ncol,pver)                            ! airglow heat rate (K/s)
      real(r8)     ::  qrs_col(pver,4)                               ! column thermal heating < 200nm
      real(r8)     ::  qrl_col(pver,4)                               ! column thermal heating > 200nm
      real(r8)     ::  qrs(ncol,pver,4)                              ! chunk thermal heating < 200nm
      real(r8)     ::  qrl(ncol,pver,4)                              ! chunk thermal heating > 200nm
      real(r8)     ::  euv_hrate_col(pver)                           ! column euv thermal heating rate
      real(r8)     ::  co2_hrate_col(pver)                           ! column co2 nir heating rate
      real(r8)     ::  euv_hrate(ncol,pver)                          ! chunk euv thermal heating rate
      real(r8)     ::  aur_hrate(ncol,pver)                          ! chunk auroral heating rate
      real(r8)     ::  co2_hrate(ncol,pver)                          ! chunk co2 nir heating rate
      real(r8)     ::  o2mmr(ncol,pver)                              ! chunk o2 concentration (kg/kg)
      real(r8)     ::  ommr(ncol,pver)                               ! chunk o concentration (kg/kg)
      real(r8)     ::  fac1(pver)                                    ! work array
      real(r8)     ::  colo3(pver)                                   ! vertical o3 column density
      real(r8)     ::  zarg(pver)                                    ! vertical height array
      real(r8)     ::  parg(pver)                                    ! vertical pressure array (hPa)
      real(r8)     ::  tline(pver)                                   ! vertical temperature array
      real(r8)     ::  eff_alb(pver)                                 ! albedo
      real(r8)     ::  mw(pver)                                      ! atms molecular weight
      real(r8)     ::  n2_line(pver)                                 ! n2 density (mol/mol)
      real(r8)     ::  o_line(pver)                                  ! o density (mol/mol)
      real(r8)     ::  o2_line(pver)                                 ! o2 density (mol/mol)
      real(r8)     ::  o3_line(pver)                                 ! o3 density (mol/mol)
      real(r8)     ::  co2_line(pver)                                ! co2 density (mol/mol)
      real(r8)     ::  scco2(pver)                                   ! co2 slant column concentration (molec/cm^2)
      real(r8)     ::  scco2i(pver)                                  ! co2 slant column concentration (molec/cm^2)
      real(r8)     ::  occ(pver)                                     ! o density (molecules/cm^3)
      real(r8)     ::  o2cc(pver)                                    ! o2 density (molecules/cm^3)
      real(r8)     ::  co2cc(pver)                                   ! co2 density (molecules/cm^3)
      real(r8)     ::  n2cc(pver)                                    ! n2 density (molecules/cm^3)
      real(r8)     ::  o3cc(pver)                                    ! o3 density (molecules/cm^3)
      real(r8)     ::  cparg(pver)                                   ! specific heat capacity
      real(r8)     ::  zen_angle(ncol)                               ! solar zenith angles (radians)
      real(r8)     ::  zsurf(ncol)                                   ! surface height (m)
      real(r8)     ::  rlats(ncol)                                   ! chunk latitudes (radians)
      real(r8)     ::  rlons(ncol)                                   ! chunk longitudes (radians)
      real(r8)     ::  calday                                        ! day of year
      real(r8)     ::  delta                                         ! solar declination (radians)
      logical      ::  do_diag

      qrs_tot(:ncol,:) = 0._r8
      if (.not. has_hrates) return
      
!-------------------------------------------------------------------------      
!        ... set maximum zenith angle - higher value for higher top model
!-------------------------------------------------------------------------      
      if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then 
         max_zen_angle = 116._r8
      else
         max_zen_angle = 97.01_r8 ! degrees
      endif

!-----------------------------------------------------------------------      
!        ... get chunk latitudes and longitudes
!-----------------------------------------------------------------------      
      lchnk = state%lchnk

      call get_lat_all_p( lchnk, ncol, latndx )
      call get_lon_all_p( lchnk, ncol, lonndx )
      call get_rlat_all_p( lchnk, ncol, rlats )
      call get_rlon_all_p( lchnk, ncol, rlons )

!-----------------------------------------------------------------------      
!        ... set lower limit for heating rates which is now dictated by radheat module
!-----------------------------------------------------------------------      
      kbot_hrates = bot_mlt_lev
      kbot_hrates = min( kbot_hrates,pver )
!     write(iulog,*) 'hrates: kbot_hrates = ',kbot_hrates

!-----------------------------------------------------------------------      
!        ... calculate cosine of zenith angle then cast back to angle
!-----------------------------------------------------------------------      
      calday = get_curr_calday()
      call zenith( calday, rlats, rlons, zen_angle, ncol )
      zen_angle(:) = acos( zen_angle(:) )

!-----------------------------------------------------------------------      
!        ... map incoming concentrations to working array
!-----------------------------------------------------------------------      
      do m = 1,pcnst
         n = map2chm(m)
         if( n > 0 ) then
            do k = 1,pver
               mmr(:ncol,k,n) = state%q(:ncol,k,m)
            end do
         end if
      end do
      call get_short_lived_species( mmr, lchnk, ncol, pbuf )

!-----------------------------------------------------------------------      
!        ... set atmosphere mean mass
!-----------------------------------------------------------------------      
      if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then 
        do k = 1,pver
          mbar(:ncol,k) = mbarv(:ncol,k,lchnk)
        enddo
      else      
        call set_mean_mass( ncol, mmr, mbar )
      endif
!
!-----------------------------------------------------------------------      
!        ... xform from mmr to vmr
!-----------------------------------------------------------------------      
      call mmr2vmr( mmr, vmr, mbar, ncol )
!-----------------------------------------------------------------------      
!        ... xform water vapor from mmr to vmr
!-----------------------------------------------------------------------      
      do k = 1,pver
         h2ovmr(:ncol,k) = vmr(:ncol,k,id_h2o)
      end do
!-----------------------------------------------------------------------      
!        ... xform geopotential height from m to km 
!            and pressure from Pa to mb
!-----------------------------------------------------------------------      
      zsurf(:ncol) = rga * state%phis(:ncol)
      do k = 1,pver
         zmid(:ncol,k) = m2km * (state%zm(:ncol,k) + zsurf(:ncol))
      end do

!-----------------------------------------------------------------------      
!        ... set the "invariants"
!-----------------------------------------------------------------------      
      call setinv( invariants, state%t, h2ovmr, vmr, state%pmid, ncol, lchnk, pbuf )

!-----------------------------------------------------------------------      
!        ... set the column densities at the upper boundary
!-----------------------------------------------------------------------      
      call set_ub_col( col_delta, vmr, invariants, state%pint(:,1), state%pdel, ncol, lchnk )

!-----------------------------------------------------------------------      
!       ...  set rates for "tabular" and user specified reactions
!-----------------------------------------------------------------------      
      do m = 1,rxntot
         do k = 1,pver
            reaction_rates(:,k,m) = 0._r8
         end do
      end do
      call setrxt_hrates( reaction_rates, state%t, invariants(1,1,indexm), ncol, kbot_hrates )
      call usrrxt_hrates( reaction_rates, state%t, state%t, state%t, invariants, &
                          h2ovmr, state%pmid, invariants(:,:,indexm), ncol, kbot_hrates )
      call adjrxt( reaction_rates, invariants, invariants(1,1,indexm), ncol )
      
!-----------------------------------------------------------------------      
!     	... set cp array
!-----------------------------------------------------------------------      
      if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then 
        do k = 1, pver
           cpair(:ncol,k) = cpairv(:ncol,k,lchnk)
        enddo
      else      
        call calc_cp( ncol, vmr, cpair )
      endif

      call outfld( 'CPAIR', cpair, ncol, lchnk )
#ifdef HRATES_DEBUG
      write(iulog,*) ' '
      write(iulog,*) '---------------------------------------------'
      write(iulog,*) 'waccm_hrates: cp at lchnk = ',lchnk
      write(iulog,'(1p,5g15.7)') cpair(1,:)
      write(iulog,*) '---------------------------------------------'
      write(iulog,*) ' '
#endif

!-----------------------------------------------------------------------      
!     	... set the earth-sun distance factor
!-----------------------------------------------------------------------      
      call shr_orb_decl( calday, eccen, mvelpp, lambm0, obliqr  , &
                         delta, esfact )
!-----------------------------------------------------------------------      
!     	... set the column densities
!-----------------------------------------------------------------------      
      call setcol( col_delta, col_dens, vmr, state%pdel,  ncol )
!-----------------------------------------------------------------------
!        ... compute the thermal heating rates
!-----------------------------------------------------------------------      
      do m = 1,4
         do k = 1,pver
            qrs(:,k,m) = 0._r8
            qrl(:,k,m) = 0._r8
         end do
      end do
      do k = 1,pver
         euv_hrate(:,k) = 0._r8
         co2_hrate(:,k) = 0._r8
      end do
column_loop : &
      do i = 1,ncol
         sza = zen_angle(i)*r2d
         if( sza < max_zen_angle ) then
            zarg(:)     = zmid(i,:)
            parg(:)     = Pa2mb*state%pmid(i,:)
            colo3(:)    = col_dens(i,:,1)
            tline(:)    = state%t(i,:)
            eff_alb(:)  = asdir(i)
            o_line(:)   = vmr(i,:,id_o)
            o2_line(:)  = vmr(i,:,id_o2)
            co2_line(:) = vmr(i,:,id_co2)
            n2_line(:)  = 1._r8 - (o_line(:) + o2_line(:) + vmr(i,:,id_h))
            o3_line(:)  = vmr(i,:,id_o3)
            occ(:)      = o_line(:) * invariants(i,:,indexm)
            o2cc(:)     = o2_line(:) * invariants(i,:,indexm)
            co2cc(:)    = co2_line(:) * invariants(i,:,indexm)
            n2cc(:)     = n2_line(:) * invariants(i,:,indexm)
            o3cc(:)     = o3_line(:) * invariants(i,:,indexm)
            mw(:)       = mbar(i,:)
            cparg(:)    = cpair(i,:)
            do_diag     = .false.
            call jshort( pver, sza, o2_line, o3_line, o2cc, &
                         o3cc, tline, zarg, mw, qrs_col, &
                         cparg, lchnk, i, co2cc, scco2, do_diag )
            call jlong( pver, sza, eff_alb, parg, tline, &
                        mw, o2_line, o3_line, colo3, qrl_col, &
                        cparg, kbot_hrates )
            do m = 1,4
               qrs(i,pver:1:-1,m) = qrs_col(:,m) * esfact
            end do
            do m = 2,4
               qrl(i,:,m) = qrl_col(:,m) * esfact
            end do
            call heuv( pver, sza, occ, o2cc, n2cc, &
                       o_line, o2_line, n2_line, cparg, mw, &
                       zarg, euv_hrate_col, kbot_hrates )
            euv_hrate(i,:) = euv_hrate_col(:) * esfact
            scco2i(1:pver) = scco2(pver:1:-1)
            call heatnirco2( co2_line, scco2i, state%pmid(i,:kbot_hrates), co2_hrate_col, kbot_hrates, &
                             zarg, sza )
#ifdef HRATES_DEBUG
            write(iulog,*) '==================================='
            write(iulog,*) 'hrates: diagnostics for heatco2nir'
            write(iulog,*) 'hrates: co2_line'
            write(iulog,'(1p,5g15.7)') co2_line(:)
            write(iulog,*) 'hrates: scco2'
            write(iulog,'(1p,5g15.7)') scco2i(:)
            write(iulog,*) 'hrates: co2_hrate'
            write(iulog,'(1p,5g15.7)') co2_hrate_col(:)
            write(iulog,*) '==================================='
#endif
            co2_hrate(i,:kbot_hrates) = co2_hrate_col(:kbot_hrates) * esfact * daypsec
         end if
      end do column_loop


      call outfld( 'QRS_SO2A', qrs(:,:,1), ncol, lchnk )
      call outfld( 'QRS_SO2B', qrs(:,:,2), ncol, lchnk )
      call outfld( 'QRS_SO3A', qrs(:,:,3), ncol, lchnk )
      call outfld( 'QRS_SO3B', qrs(:,:,4), ncol, lchnk )
      call outfld( 'QRS_LO2B', qrl(:,:,2), ncol, lchnk )
      call outfld( 'QRS_LO3A', qrl(:,:,3), ncol, lchnk )
      call outfld( 'QRS_LO3B', qrl(:,:,4), ncol, lchnk )
      call outfld( 'QRS_LO3',  qrl(:,:,3)+qrl(:,:,4), ncol, lchnk )
      call outfld( 'QRS_EUV', euv_hrate(:,:), ncol, lchnk )
      call outfld( 'QRS_CO2NIR', co2_hrate(:,:), ncol, lchnk )

!-----------------------------------------------------------------------      
!     	... chemical pot heating rate
!-----------------------------------------------------------------------      
      call cph( cphrate, vmr, reaction_rates, cpair, mbar, &
                kbot_hrates, ncol, lchnk )

!-----------------------------------------------------------------------      
!     	... auroral ion production
!-----------------------------------------------------------------------      
      do k = 1,pver
         o2mmr(:ncol,k) = mmr(:ncol,k,id_o2)
         ommr(:ncol,k)  = mmr(:ncol,k,id_o)
      end do
      call aurora( state%t, o2mmr, ommr, mbar, rlats, &
                   aur_hrate, cpair, state%pmid, lchnk, calday, &
                   ncol, rlons )
      do k = 1,pver
         aur_hrate(:,k)  = aur_hrate(:,k)/invariants(:,k,indexm)
      end do
      call outfld( 'QRS_AUR', aur_hrate(:,:), ncol, lchnk )

!-----------------------------------------------------------------------      
!     	... airglow heating rate
!-----------------------------------------------------------------------      
      call airglow( aghrate, vmr(1,1,id_o2_1s), vmr(1,1,id_o2_1d), reaction_rates, cpair, &
                    ncol, lchnk )

!-----------------------------------------------------------------------      
!     	... form total heating rate
!-----------------------------------------------------------------------      
      do k = 1,kbot_hrates
         qrs_tot(:ncol,k) = qrs(:,k,1) + qrs(:,k,2) + qrs(:,k,3) + qrs(:,k,4) &
                          + qrl(:,k,1) + qrl(:,k,2) + qrl(:,k,3) + qrl(:,k,4)
      end do
      call outfld( 'QTHERMAL', qrs_tot, pcols, lchnk )
      do k = 1,kbot_hrates
         qrs_tot(:ncol,k) = qrs_tot(:ncol,k) &
                          + cphrate(:,k) + euv_hrate(:,k) + aur_hrate(:,k) + co2_hrate(:,k)
      end do
      call outfld( 'QRS_MLT', qrs_tot, pcols, lchnk )

      end subroutine waccm_hrates

      end module mo_waccm_hrates
