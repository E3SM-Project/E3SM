module ice_import_export

  use shr_kind_mod    , only: r8 => shr_kind_r8, cl=>shr_kind_cl
  use shr_sys_mod     , only: shr_sys_abort, shr_sys_flush
  use ice_kinds_mod   , only: int_kind, dbl_kind, char_len_long, log_kind
  use ice_constants   , only: c0, c1, puny, tffresh, spval_dbl
  use ice_constants   , only: field_loc_center, field_type_scalar
  use ice_constants   , only: field_type_vector, c100
  use ice_blocks      , only: block, get_block, nx_block, ny_block
  use ice_flux        , only: strairxt, strairyt, strocnxt, strocnyt           
  use ice_flux        , only: alvdr, alidr, alvdf, alidf, Tref, Qref, Uref
  use ice_flux        , only: flat, fsens, flwout, evap, fswabs, fhocn, fswthru                                          
  use ice_flux        , only: fresh, fsalt, zlvl, uatm, vatm, potT, Tair, Qa
  use ice_flux        , only: rhoa, swvdr, swvdf, swidr, swidf, flw, frain
  use ice_flux        , only: fsnow, uocn, vocn, sst, ss_tltx, ss_tlty, frzmlt
  use ice_flux        , only: sss, tf, wind, fsw, init_flux_atm, init_flux_ocn, faero
  use ice_state       , only: vice, vsno, aice, trcr, filename_aero, filename_iage
  use ice_state       , only: filename_volpn, filename_FY, filename_lvl
  use ice_state       , only: tr_aero, tr_iage, tr_FY, tr_pond, tr_lvl 
  use ice_domain      , only: nblocks, blocks_ice, halo_info, distrb_info, profile_barrier
  use ice_domain_size , only: nx_global, ny_global, block_size_x, block_size_y, max_blocks
  use ice_grid        , only: tlon, tlat, tarea, tmask, anglet, hm, ocn_gridcell_frac
  use ice_grid        , only: grid_type, t2ugrid_vector, gridcpl_file
  use ice_boundary    , only: ice_HaloUpdate 
  use ice_fileunits   , only: nu_diag
  use ice_prescribed_mod
  use ice_cpl_indices
  use perf_mod        , only: t_startf, t_stopf, t_barrierf
  implicit none
  public

  !==============================================================================
contains
  !==============================================================================

  subroutine ice_import( x2i )

    !-----------------------------------------------------
    ! Arguments
    real(r8), intent(inout) :: x2i(:,:)
    !
    ! Local variables
    integer     :: i, j, iblk, n
    integer     :: ilo, ihi, jlo, jhi !beginning and end of physical domain
    type(block) :: this_block         ! block information for current block
    integer,parameter                :: nflds=15,nfldv=6
    real (kind=dbl_kind),allocatable :: aflds(:,:,:,:)
    real (kind=dbl_kind)             :: workx, worky
    logical (kind=log_kind)          :: first_call = .true.
    !-----------------------------------------------------

    ! Note that the precipitation fluxes received  from the coupler
    ! are in units of kg/s/m^2 which is what CICE requires.
    ! Note also that the read in below includes only values needed
    ! by the thermodynamic component of CICE.  Variables uocn, vocn,
    ! ss_tltx, and ss_tlty are excluded. Also, because the SOM and
    ! DOM don't  compute SSS.   SSS is not read in and is left at
    ! the initilized value (see ice_flux.F init_coupler_flux) of
    ! 34 ppt

    ! Use aflds to gather the halo updates of multiple fields
    ! Need to separate the scalar from the vector halo updates

    allocate(aflds(nx_block,ny_block,nflds,nblocks))
    aflds = c0

    n=0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)         
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi

       do j = jlo, jhi
          do i = ilo, ihi

             n = n+1
             aflds(i,j, 1,iblk)   = x2i(index_x2i_So_t,n)
             aflds(i,j, 2,iblk)   = x2i(index_x2i_So_s,n)
             aflds(i,j, 3,iblk)   = x2i(index_x2i_Sa_z,n)
             aflds(i,j, 4,iblk)   = x2i(index_x2i_Sa_ptem,n)
             aflds(i,j, 5,iblk)   = x2i(index_x2i_Sa_tbot,n)
             aflds(i,j, 6,iblk)   = x2i(index_x2i_Sa_shum,n)
             aflds(i,j, 7,iblk)   = x2i(index_x2i_Sa_dens,n)
             aflds(i,j, 8,iblk)   = x2i(index_x2i_Fioo_q,n)
             aflds(i,j, 9,iblk)   = x2i(index_x2i_Faxa_swvdr,n)
             aflds(i,j,10,iblk)   = x2i(index_x2i_Faxa_swndr,n)
             aflds(i,j,11,iblk)   = x2i(index_x2i_Faxa_swvdf,n)
             aflds(i,j,12,iblk)   = x2i(index_x2i_Faxa_swndf,n)
             aflds(i,j,13,iblk)   = x2i(index_x2i_Faxa_lwdn,n)
             aflds(i,j,14,iblk)   = x2i(index_x2i_Faxa_rain,n)
             aflds(i,j,15,iblk)   = x2i(index_x2i_Faxa_snow,n)

          enddo    !i
       enddo    !j

    enddo        !iblk

    if (.not.prescribed_ice) then
       call t_startf ('cice_imp_halo')
       call ice_HaloUpdate(aflds, halo_info, field_loc_center, &
            field_type_scalar)
       call t_stopf ('cice_imp_halo')
    endif

    !$OMP PARALLEL DO PRIVATE(iblk,i,j)
    do iblk = 1, nblocks
       do j = 1,ny_block
          do i = 1,nx_block
             sst  (i,j,iblk)   = aflds(i,j, 1,iblk)
             sss  (i,j,iblk)   = aflds(i,j, 2,iblk)
             zlvl (i,j,iblk)   = aflds(i,j, 3,iblk)
             potT (i,j,iblk)   = aflds(i,j, 4,iblk)
             Tair (i,j,iblk)   = aflds(i,j, 5,iblk)
             Qa   (i,j,iblk)   = aflds(i,j, 6,iblk)
             rhoa (i,j,iblk)   = aflds(i,j, 7,iblk)
             frzmlt (i,j,iblk) = aflds(i,j, 8,iblk)
             swvdr(i,j,iblk)   = aflds(i,j, 9,iblk)
             swidr(i,j,iblk)   = aflds(i,j,10,iblk)
             swvdf(i,j,iblk)   = aflds(i,j,11,iblk)
             swidf(i,j,iblk)   = aflds(i,j,12,iblk)
             flw  (i,j,iblk)   = aflds(i,j,13,iblk)
             frain(i,j,iblk)   = aflds(i,j,14,iblk)
             fsnow(i,j,iblk)   = aflds(i,j,15,iblk)
          enddo    !i
       enddo    !j
    enddo        !iblk
    !$OMP END PARALLEL DO

    deallocate(aflds)
    allocate(aflds(nx_block,ny_block,nfldv,nblocks))
    aflds = c0

    n=0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi

       do j = jlo, jhi
          do i = ilo, ihi
             n = n+1
             aflds(i,j, 1,iblk)   = x2i(index_x2i_So_u,n)
             aflds(i,j, 2,iblk)   = x2i(index_x2i_So_v,n)
             aflds(i,j, 3,iblk)   = x2i(index_x2i_Sa_u,n)
             aflds(i,j, 4,iblk)   = x2i(index_x2i_Sa_v,n)
             aflds(i,j, 5,iblk)   = x2i(index_x2i_So_dhdx,n)
             aflds(i,j, 6,iblk)   = x2i(index_x2i_So_dhdy,n)
          enddo
       enddo
    enddo

    if (.not.prescribed_ice) then
       call t_startf ('cice_imp_halo')
       call ice_HaloUpdate(aflds, halo_info, field_loc_center, &
            field_type_vector)
       call t_stopf ('cice_imp_halo')
    endif

    !$OMP PARALLEL DO PRIVATE(iblk,i,j)
    do iblk = 1, nblocks
       do j = 1,ny_block
          do i = 1,nx_block
             uocn (i,j,iblk)   = aflds(i,j, 1,iblk)
             vocn (i,j,iblk)   = aflds(i,j, 2,iblk)
             uatm (i,j,iblk)   = aflds(i,j, 3,iblk)
             vatm (i,j,iblk)   = aflds(i,j, 4,iblk)
             ss_tltx(i,j,iblk) = aflds(i,j, 5,iblk)
             ss_tlty(i,j,iblk) = aflds(i,j, 6,iblk)
          enddo    !i
       enddo    !j
    enddo        !iblk
    !$OMP END PARALLEL DO

    deallocate(aflds)

    !-------------------------------------------------------
    ! Set aerosols from coupler 
    !-------------------------------------------------------

    n=0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)         
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi

       do j = jlo, jhi
          do i = ilo, ihi

             n = n+1

#ifdef MODAL_AER
            !HW++ modal treatment

            ! BC species 1 (=intersitial/external BC)
            faero(i,j,1,iblk) = x2i(index_x2i_Faxa_bcphodry,n) &
                              + x2i(index_x2i_Faxa_bcphidry,n)

            ! BC species 2 (=cloud_water/within-ice BC)
            faero(i,j,2,iblk) = x2i(index_x2i_Faxa_bcphiwet,n)

            ! Combine all of the dust into one category
            faero(i,j,3,iblk) = x2i(index_x2i_Faxa_dstwet1,n) &
                              + x2i(index_x2i_Faxa_dstdry1,n) &
                              + x2i(index_x2i_Faxa_dstwet2,n) &
                              + x2i(index_x2i_Faxa_dstdry2,n) &
                              + x2i(index_x2i_Faxa_dstwet3,n) &
                              + x2i(index_x2i_Faxa_dstdry3,n) &
                              + x2i(index_x2i_Faxa_dstwet4,n) &
                              + x2i(index_x2i_Faxa_dstdry4,n)
            !mgf--
#else
            ! bulk treatment

             faero(i,j,1,iblk) = x2i(index_x2i_Faxa_bcphodry,n)

             faero(i,j,2,iblk) = x2i(index_x2i_Faxa_bcphidry,n) &
                  + x2i(index_x2i_Faxa_bcphiwet,n)
             ! Combine all of the dust into one category
             faero(i,j,3,iblk) = x2i(index_x2i_Faxa_dstwet1,n) &
                  + x2i(index_x2i_Faxa_dstdry1,n) &
                  + x2i(index_x2i_Faxa_dstwet2,n) &
                  + x2i(index_x2i_Faxa_dstdry2,n) &
                  + x2i(index_x2i_Faxa_dstwet3,n) &
                  + x2i(index_x2i_Faxa_dstdry3,n) &
                  + x2i(index_x2i_Faxa_dstwet4,n) &
                  + x2i(index_x2i_Faxa_dstdry4,n)

#endif
!HW +++
          enddo    !i
       enddo    !j

    enddo        !iblk


    !-----------------------------------------------------------------
    ! rotate zonal/meridional vectors to local coordinates
    ! compute data derived quantities
    !-----------------------------------------------------------------

    ! Vector fields come in on T grid, but are oriented geographically
    ! need to rotate to pop-grid FIRST using ANGLET
    ! then interpolate to the U-cell centers  (otherwise we
    ! interpolate across the pole)
    ! use ANGLET which is on the T grid !

    call t_startf ('cice_imp_ocn')
    !$OMP PARALLEL DO PRIVATE(iblk,i,j,workx,worky)
    do iblk = 1, nblocks

       do j = 1,ny_block
          do i = 1,nx_block

             ! ocean
             workx      = uocn  (i,j,iblk) ! currents, m/s 
             worky      = vocn  (i,j,iblk)
             uocn(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) & ! convert to POP grid 
                  + worky*sin(ANGLET(i,j,iblk))
             vocn(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) &
                  - workx*sin(ANGLET(i,j,iblk))

             workx      = ss_tltx  (i,j,iblk)           ! sea sfc tilt, m/m
             worky      = ss_tlty  (i,j,iblk)
             ss_tltx(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) & ! convert to POP grid 
                  + worky*sin(ANGLET(i,j,iblk))
             ss_tlty(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) &
                  - workx*sin(ANGLET(i,j,iblk))

             sst(i,j,iblk) = sst(i,j,iblk) - Tffresh       ! sea sfc temp (C)
             Tf (i,j,iblk) = -1.8_dbl_kind                 ! hardwired for NCOM
             !         Tf (i,j,iblk) = -depressT*sss(i,j,iblk)       ! freezing temp (C)
             !         Tf (i,j,iblk) = -depressT*max(sss(i,j,iblk),ice_ref_salinity)

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    call t_stopf ('cice_imp_ocn')

    ! Interpolate ocean dynamics variables from T-cell centers to 
    ! U-cell centers.

    if (.not.prescribed_ice) then
       call t_startf ('cice_imp_t2u')
       call t2ugrid_vector(uocn)
       call t2ugrid_vector(vocn)
       call t2ugrid_vector(ss_tltx)
       call t2ugrid_vector(ss_tlty)
       call t_stopf ('cice_imp_t2u')
    end if

    ! Atmosphere variables are needed in T cell centers in
    ! subroutine stability and are interpolated to the U grid
    ! later as necessary.

    call t_startf ('cice_imp_atm')
    !$OMP PARALLEL DO PRIVATE(iblk,i,j,workx,worky)
    do iblk = 1, nblocks
       do j = 1, ny_block
          do i = 1, nx_block

             ! atmosphere
             workx      = uatm(i,j,iblk) ! wind velocity, m/s
             worky      = vatm(i,j,iblk) 
             uatm (i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) & ! convert to POP grid
                             + worky*sin(ANGLET(i,j,iblk))   ! note uatm, vatm, wind
             vatm (i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) & ! are on the T-grid here
                  - workx*sin(ANGLET(i,j,iblk))

             wind (i,j,iblk) = sqrt(uatm(i,j,iblk)**2 + vatm(i,j,iblk)**2)
             fsw  (i,j,iblk) = swvdr(i,j,iblk) + swvdf(i,j,iblk) &
                             + swidr(i,j,iblk) + swidf(i,j,iblk)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    call t_stopf ('cice_imp_atm')

  end subroutine ice_import

  !===============================================================================

  subroutine ice_export( i2x )   

    !-----------------------------------------------------
    !
    ! Arguments
    real(r8), intent(inout) :: i2x(:,:)
    !
    ! Local Variables
    integer :: i, j, iblk, n, ij 
    integer :: ilo, ihi, jlo, jhi !beginning and end of physical domain
    integer (kind=int_kind)                                :: icells ! number of ocean/ice cells
    integer (kind=int_kind), dimension (nx_block*ny_block) :: indxi  ! compressed indices in i
    integer (kind=int_kind), dimension (nx_block*ny_block) :: indxj  ! compressed indices in i

    real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: &
         Tsrf  &      ! surface temperature
         ,  tauxa &      ! atmo/ice stress
         ,  tauya &
         ,  tauxo &      ! ice/ocean stress
         ,  tauyo &
         ,  ailohi       ! fractional ice area

    real (kind=dbl_kind) :: &
         workx, worky           ! tmps for converting grid

    type(block)        :: this_block                           ! block information for current block
    logical :: flag
    !-----------------------------------------------------

    flag=.false.

    !calculate ice thickness from aice and vice. Also
    !create Tsrf from the first tracer (trcr) in ice_state.F

    !$OMP PARALLEL DO PRIVATE(iblk,i,j,workx,worky)
    do iblk = 1, nblocks
       do j = 1, ny_block
          do i = 1, nx_block

             ! ice fraction
             ailohi(i,j,iblk) = min(aice(i,j,iblk), c1)

             ! surface temperature
             Tsrf(i,j,iblk)  = Tffresh + trcr(i,j,1,iblk)             !Kelvin (original ???)

             ! wind stress  (on POP T-grid:  convert to lat-lon)
             workx = strairxT(i,j,iblk)                             ! N/m^2
             worky = strairyT(i,j,iblk)                             ! N/m^2
             tauxa(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) &
                             - worky*sin(ANGLET(i,j,iblk))
             tauya(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) &
                             + workx*sin(ANGLET(i,j,iblk))

             ! ice/ocean stress (on POP T-grid:  convert to lat-lon)
             workx = -strocnxT(i,j,iblk)                            ! N/m^2
             worky = -strocnyT(i,j,iblk)                            ! N/m^2
             tauxo(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) &
                             - worky*sin(ANGLET(i,j,iblk))
             tauyo(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) &
                             + workx*sin(ANGLET(i,j,iblk))

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    do iblk = 1, nblocks
       do j = 1, ny_block
          do i = 1, nx_block
             if (tmask(i,j,iblk) .and. ailohi(i,j,iblk) < c0 ) then
                flag = .true.
             endif
          end do
       end do
    end do
    if (flag) then
       do iblk = 1, nblocks
          do j = 1, ny_block
             do i = 1, nx_block
                if (tmask(i,j,iblk) .and. ailohi(i,j,iblk) < c0 ) then
                   write(nu_diag,*) &
                        ' (ice) send: ERROR ailohi < 0.0 ',i,j,ailohi(i,j,iblk)
                   call shr_sys_flush(nu_diag)
                endif
             end do
          end do
       end do
    endif

    ! Fill export state i2x_i

    i2x(:,:) = spval_dbl

    n=0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)         
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi

       do j = jlo, jhi
          do i = ilo, ihi

             n = n+1

             !-------states-------------------- 
             i2x(index_i2x_Si_ifrac ,n)    = ailohi(i,j,iblk)   

             if ( tmask(i,j,iblk) .and. ailohi(i,j,iblk) > c0 ) then
                !-------states-------------------- 
                i2x(index_i2x_Si_t     ,n)    = Tsrf(i,j,iblk)
                i2x(index_i2x_Si_avsdr ,n)    = alvdr(i,j,iblk)
                i2x(index_i2x_Si_anidr ,n)    = alidr(i,j,iblk)
                i2x(index_i2x_Si_avsdf ,n)    = alvdf(i,j,iblk)
                i2x(index_i2x_Si_anidf ,n)    = alidf(i,j,iblk)
                i2x(index_i2x_Si_u10  ,n)     = Uref(i,j,iblk)
                i2x(index_i2x_Si_tref  ,n)    = Tref(i,j,iblk)
                i2x(index_i2x_Si_qref  ,n)    = Qref(i,j,iblk)
                i2x(index_i2x_Si_snowh ,n)    = vsno(i,j,iblk) &
                     / ailohi(i,j,iblk)

                !--- a/i fluxes computed by ice
                i2x(index_i2x_Faii_taux ,n)   = tauxa(i,j,iblk)    
                i2x(index_i2x_Faii_tauy ,n)   = tauya(i,j,iblk)    
                i2x(index_i2x_Faii_lat  ,n)   = flat(i,j,iblk)     
                i2x(index_i2x_Faii_sen  ,n)   = fsens(i,j,iblk)    
                i2x(index_i2x_Faii_lwup ,n)   = flwout(i,j,iblk)   
                i2x(index_i2x_Faii_evap ,n)   = evap(i,j,iblk)     
                i2x(index_i2x_Faii_swnet,n)   = fswabs(i,j,iblk)

                !--- i/o fluxes computed by ice
                i2x(index_i2x_Fioi_melth,n)   = fhocn(i,j,iblk)
                i2x(index_i2x_Fioi_swpen,n)   = fswthru(i,j,iblk) ! hf from melting          
                i2x(index_i2x_Fioi_meltw,n)   = fresh(i,j,iblk)   ! h2o flux from melting    ???
                i2x(index_i2x_Fioi_salt ,n)   = fsalt(i,j,iblk)   ! salt flux from melting   ???
                i2x(index_i2x_Fioi_taux ,n)   = tauxo(i,j,iblk)   ! stress : i/o zonal       ???
                i2x(index_i2x_Fioi_tauy ,n)   = tauyo(i,j,iblk)   ! stress : i/o meridional  ???
             end if
          enddo    !i
       enddo    !j
    enddo        !iblk

  end subroutine ice_export

end module ice_import_export
