!  SVN:$Id: ice_tracers.F90 -1   $
!=======================================================================
! Indices and flags associated with the tracer infrastructure. 
! Grid-dependent and max_trcr-dependent arrays are declared in ice_state.F90.
!
! author Elizabeth C. Hunke, LANL

      module ice_colpkg_tracers

      use ice_kinds_mod

      implicit none
      save

      private
      public :: colpkg_compute_tracers

      integer (kind=int_kind), public :: &
         ntrcr     ! number of tracers in use

      integer (kind=int_kind), public :: &
         nbtrcr    ! number of bgc tracers in use
      
      integer (kind=int_kind), public :: &
         nt_Tsfc  , & ! ice/snow temperature
         nt_qice  , & ! volume-weighted ice enthalpy (in layers)
         nt_qsno  , & ! volume-weighted snow enthalpy (in layers)
         nt_sice  , & ! volume-weighted ice bulk salinity (CICE grid layers)
         nt_fbri  , & ! volume fraction of ice with dynamic salt (hinS/vicen*aicen)
         nt_iage  , & ! volume-weighted ice age
         nt_FY    , & ! area-weighted first-year ice area
         nt_alvl  , & ! level ice area fraction
         nt_vlvl  , & ! level ice volume fraction
         nt_apnd  , & ! melt pond area fraction
         nt_hpnd  , & ! melt pond depth
         nt_ipnd  , & ! melt pond refrozen lid thickness
         nt_aero  , & ! starting index for aerosols in ice
         nt_bgc_N_sk    , & ! algae (skeletal layer)
         nt_bgc_C_sk    , & ! 
         nt_bgc_chl_sk  , & ! 
         nt_bgc_Nit_sk  , & ! nutrients (skeletal layer) 
         nt_bgc_Am_sk   , & ! 
         nt_bgc_Sil_sk  , & !
         nt_bgc_DMSPp_sk, & ! trace gases (skeletal layer)
         nt_bgc_DMSPd_sk, & ! 
         nt_bgc_DMS_sk  , & ! 
         nt_bgc_Nit_ml  , & ! nutrients (ocean mixed layer) 
         nt_bgc_Am_ml   , & ! 
         nt_bgc_Sil_ml  , & !
         nt_bgc_DMSP_ml , & ! trace gases (ocean mixed layer)
         nt_bgc_DMS_ml

      logical (kind=log_kind), public :: &
         tr_iage     , & ! if .true., use age tracer
         tr_FY       , & ! if .true., use first-year area tracer
         tr_lvl      , & ! if .true., use level ice tracer
         tr_pond     , & ! if .true., use melt pond tracer
         tr_pond_cesm, & ! if .true., use cesm pond tracer
         tr_pond_lvl , & ! if .true., use level-ice pond tracer
         tr_pond_topo, & ! if .true., use explicit topography-based ponds
         tr_aero     , & ! if .true., use aerosol tracers
         tr_brine        ! if .true., brine height differs from ice thickness

!=======================================================================

      contains

!=======================================================================

! Compute tracer fields.
! Given atrcrn = aicen*trcrn (or vicen*trcrn, vsnon*trcrn), compute trcrn.
!
! author: William H. Lipscomb, LANL

      subroutine colpkg_compute_tracers (ntrcr,     trcr_depend,    &
                                         atrcrn,    aicen,          &
                                         vicen,     vsnon,          &
                                         trcr_base, n_trcr_strata,  &
                                         nt_strata, trcrn)

      use ice_constants_colpkg, only: c0, c1, puny, Tocnfrz

      integer (kind=int_kind), intent(in) :: &
         ntrcr                 ! number of tracers in use

      integer (kind=int_kind), dimension (ntrcr), intent(in) :: &
         trcr_depend, & ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon
         n_trcr_strata  ! number of underlying tracer layers

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         trcr_base      ! = 0 or 1 depending on tracer dependency
                        ! argument 2:  (1) aice, (2) vice, (3) vsno

      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nt_strata      ! indices of underlying tracer layers

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         atrcrn    ! aicen*trcrn or vicen*trcrn or vsnon*trcrn

      real (kind=dbl_kind), intent(in) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (ntrcr), intent(out) :: &
         trcrn     ! ice tracers

      ! local variables

      integer (kind=int_kind) :: &
         it,     & ! tracer index
         itl,    & ! tracer index
         ntr,    & ! tracer index
         k         ! loop index

      real (kind=dbl_kind), dimension(3) :: &
         divisor   ! base quantity on which tracers are carried

      real (kind=dbl_kind) :: &
         work      ! temporary scalar

      !-----------------------------------------------------------------
      ! Compute new tracers
      !-----------------------------------------------------------------

      do it = 1, ntrcr
         divisor(1) = trcr_base(it,1)*aicen
         divisor(2) = trcr_base(it,2)*vicen
         divisor(3) = trcr_base(it,3)*vsnon

         if (trcr_depend(it) == 0) then ! ice area tracers
            if (aicen > puny) then  
               trcrn(it) = atrcrn(it) / aicen
            else
               trcrn(it) = c0
               if (it == nt_Tsfc) trcrn(it) = Tocnfrz  ! surface temperature
            endif

         else

            work = c0
            do k = 1, 3
               if (divisor(k) > c0) then
                  work = atrcrn(it) / divisor(k)
               endif
               if (trcr_base(it,k) > c0) then
                  ! nonzero default values could be put in an array
                  if (it == nt_fbri) work = c1       ! brine fraction
               endif
            enddo
            trcrn(it) = work                ! save
            if (n_trcr_strata(it) > 0) then          ! additional tracer layers
               do itl = 1, n_trcr_strata(it)
                  ntr = nt_strata(it,itl)
                  if (trcrn(ntr) > c0) then
                      trcrn(it) = trcrn(it) / trcrn(ntr)
                  else
                      trcrn(it) = c0
                  endif
               enddo
            endif

         endif ! trcr_depend=0

      enddo

    end subroutine colpkg_compute_tracers

!=======================================================================

      end module ice_colpkg_tracers

!=======================================================================
