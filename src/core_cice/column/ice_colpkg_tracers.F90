!  SVN:$Id: ice_tracers.F90 -1   $
!=======================================================================
! Indices and flags associated with the tracer infrastructure. 
! Grid-dependent and max_trcr-dependent arrays are declared in ice_state.F90.
!
! author Elizabeth C. Hunke, LANL

      module ice_colpkg_tracers

      use ice_kinds_mod
      use ice_colpkg_shared, only: max_algae, max_dic, max_doc, max_don, &
          max_fe, max_aero, max_nbtrcr

      implicit none
      save

      private
      public :: colpkg_compute_tracers

      integer (kind=int_kind), public :: &
         ntrcr   , &  ! number of tracers in use
         ntrcr_o      ! number of non-bio tracers in use

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
         nt_bgc_Nit,   & ! nutrients  
         nt_bgc_Am,    & ! 
         nt_bgc_Sil,   & !
         nt_bgc_DMSPp, & ! trace gases (skeletal layer)
         nt_bgc_DMSPd, & ! 
         nt_bgc_DMS,   & ! 
         nt_bgc_PON,   & ! zooplankton and detritus 
         nt_bgc_hum,   & ! humic material 
         nt_zbgc_frac, & ! fraction of tracer in the mobile phase
         nt_bgc_S        ! Bulk salinity in fraction ice with dynamic salinity (Bio grid)

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

      !-----------------------------------------------------------------
      !  biogeochemistry
      !-----------------------------------------------------------------

      logical (kind=log_kind), public :: & 
         tr_bgc_S,       & ! if .true., use zsalinity
         tr_zaero,       & ! if .true., black carbon is tracers  (n_zaero)
         tr_bgc_Nit,     & ! if .true. Nitrate tracer in ice 
         tr_bgc_N,       & ! if .true., algal nitrogen tracers  (n_algae)
         tr_bgc_DON,     & ! if .true., DON pools are tracers  (n_don)
         tr_bgc_C,       & ! if .true., algal carbon tracers + DOC and DIC 
         tr_bgc_chl,     & ! if .true., algal chlorophyll tracers 
         tr_bgc_Am,      & ! if .true., ammonia/um as nutrient tracer 
         tr_bgc_Sil,     & ! if .true., silicon as nutrient tracer 
         tr_bgc_DMS,     & ! if .true., DMS as  tracer 
         tr_bgc_Fe,      & ! if .true., Fe as  tracer 
         tr_bgc_PON,     & ! if .true., PON as tracer 
         tr_bgc_hum        ! if .true., humic material as tracer 

      integer (kind=int_kind), public :: &
         nbtrcr,         & ! number of bgc tracers in use
         nbtrcr_sw,      & ! number of bgc tracers which impact shortwave
         nlt_chl_sw        ! points to total chla in trcrn_sw

      integer (kind=int_kind), dimension(max_aero), public :: &
         nlt_zaero_sw       ! points to aerosol in trcrn_sw
  
      integer (kind=int_kind), dimension(max_algae), public :: &
         nlt_bgc_N      , & ! algae 
         nlt_bgc_C      , & ! 
         nlt_bgc_chl   

      integer (kind=int_kind), dimension(max_doc), public :: &
         nlt_bgc_DOC        ! disolved organic carbon

      integer (kind=int_kind), dimension(max_don), public :: &
         nlt_bgc_DON        !

      integer (kind=int_kind), dimension(max_dic), public :: &
         nlt_bgc_DIC        ! disolved inorganic carbon

      integer (kind=int_kind), dimension(max_fe), public :: &
         nlt_bgc_Fed    , & !
         nlt_bgc_Fep        !

      integer (kind=int_kind), dimension(max_aero), public :: &
         nlt_zaero          ! non-reacting layer aerosols

      integer (kind=int_kind), public :: &
         nlt_bgc_Nit   ,   & ! nutrients  
         nlt_bgc_Am    ,   & ! 
         nlt_bgc_Sil   ,   & !
         nlt_bgc_DMSPp ,   & ! trace gases (skeletal layer)
         nlt_bgc_DMSPd ,   & ! 
         nlt_bgc_DMS   ,   & ! 
         nlt_bgc_PON   ,   & ! zooplankton and detritus
         nlt_bgc_hum         ! humic material

      integer (kind=int_kind), dimension(max_algae), public :: &  
         nt_bgc_N , & ! diatoms, phaeocystis, pico/small   
         nt_bgc_C , & ! diatoms, phaeocystis, pico/small   
         nt_bgc_chl   ! diatoms, phaeocystis, pico/small 

      integer (kind=int_kind), dimension(max_doc), public :: &  
         nt_bgc_DOC      !  dissolved organic carbon

      integer (kind=int_kind), dimension(max_don), public :: & 
         nt_bgc_DON         !  dissolved organic nitrogen

      integer (kind=int_kind), dimension(max_dic), public :: &  
         nt_bgc_DIC         !  dissolved inorganic carbon

      integer (kind=int_kind), dimension(max_fe), public :: & 
         nt_bgc_Fed,     & !  dissolved iron
         nt_bgc_Fep        !  particulate iron

      integer (kind=int_kind), dimension(max_aero), public :: &  
         nt_zaero       !  black carbon and other aerosols
      
      integer (kind=int_kind), dimension(max_nbtrcr), public :: &
         bio_index_o         ! relates nlt_bgc_NO to ocean concentration index
                             ! see ocean_bio_all

      integer (kind=int_kind), dimension(max_nbtrcr), public :: &
         bio_index           ! relates bio indices, ie.  nlt_bgc_N to nt_bgc_N 

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
            if (vicen <= c0 .and. it == nt_fbri) trcrn(it) = c1

         endif ! trcr_depend=0

      enddo

    end subroutine colpkg_compute_tracers

!=======================================================================

      end module ice_colpkg_tracers

!=======================================================================
