!=======================================================================
!BOP
!
! !MODULE: ice_ocean - ocean mixed layer internal to sea ice model
!
! !DESCRIPTION:
!
! Ocean mixed layer calculation (internal to sea ice model).
! Allows heat storage in ocean for uncoupled runs.
!
! !REVISION HISTORY:
!  SVN:$Id: ice_ocean.F90 131 2008-05-30 16:53:40Z eclare $
!
! authors:   John Weatherly, CRREL
!            C.M. Bitz, UW
!            Elizabeth C. Hunke, LANL
!            Bruce P. Briegleb, NCAR
!            William H. Lipscomb, LANL
!
! 2004: Block structure added by William Lipscomb
! 2005: Ocean-to-atmosphere fluxes added as 3D arrays, William Lipscomb
! 2006: Converted to free source form (F90) by Elizabeth Hunke
!
! !INTERFACE:
!
      module ice_ocean
!
! !USES:
!
      use ice_kinds_mod
      use ice_constants
!
!EOP
!
      implicit none
      save

      logical (kind=log_kind) :: &
         oceanmixed_ice           ! if true, use ocean mixed layer

      real (kind=dbl_kind), parameter :: &
         cprho = cp_ocn*rhow

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !ROUTINE: ocean_mixed_layer - compute SST and freeze/melt potential
!
! !DESCRIPTION:
!
! Compute the mixed layer heat balance and update the SST.
! Compute the energy available to freeze or melt ice.
! NOTE: SST changes due to fluxes through the ice are computed in
!       ice_therm_vertical.
!
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine ocean_mixed_layer (dt)
!
! !USES:
!
      use ice_blocks
      use ice_domain
      use ice_state, only: aice, uvel, vvel
      use ice_flux
      use ice_grid, only: tmask
      use ice_atmo, only: atmo_boundary_layer, atmbndy, atmo_boundary_const
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step
!
!EOP
!
      real (kind=dbl_kind) :: &
         TsfK , & ! surface temperature (K)
         swabs    ! surface absorbed shortwave heat flux (W/m^2)

      real (kind=dbl_kind), parameter :: &
         frzmlt_max = c1000   ! max magnitude of frzmlt (W/m^2)

      integer (kind=int_kind) :: &
         i, j           , & ! horizontal indices
         ij             , & ! combined ij index
         iblk           , & ! block index
         ilo,ihi,jlo,jhi    ! beginning and end of physical domain

      real (kind=dbl_kind), dimension(nx_block,ny_block) :: &
         delt  , & ! potential temperature difference   (K)
         delq  , & ! specific humidity difference   (kg/kg)
         shcoef, & ! transfer coefficient for sensible heat
         lhcoef    ! transfer coefficient for latent heat

      integer (kind=int_kind), save :: &
         icells    ! number of ocean cells

      integer (kind=int_kind), dimension(nx_block*ny_block), save :: &
         indxi, indxj    ! compressed indices for ocean cells

      type (block) :: &
         this_block           ! block information for current block

      do iblk = 1, nblocks

      !-----------------------------------------------------------------
      ! Identify ocean cells.
      ! Set fluxes to zero in land cells.
      !-----------------------------------------------------------------

         icells = 0
         do j = 1, ny_block
         do i = 1, nx_block
            if (tmask(i,j,iblk)) then
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
            else
               sst       (i,j,iblk) = c0
               frzmlt    (i,j,iblk) = c0
               flwout_ocn(i,j,iblk) = c0
               fsens_ocn (i,j,iblk) = c0
               flat_ocn  (i,j,iblk) = c0
               evap_ocn  (i,j,iblk) = c0
            endif
         enddo                  ! i
         enddo                  ! j

      !-----------------------------------------------------------------
      ! Compute boundary layer quantities
      !-----------------------------------------------------------------

         if (trim(atmbndy) == 'constant') then
            call atmo_boundary_const (nx_block,  ny_block,   &
                                      'ice',     icells,     &
                                      indxi,     indxj,      &
                                      uatm       (:,:,iblk), &   
                                      vatm       (:,:,iblk), &   
                                      wind       (:,:,iblk), &   
                                      rhoa       (:,:,iblk), &
                                      strairx_ocn(:,:,iblk), & 
                                      strairy_ocn(:,:,iblk), & 
                                      lhcoef     (:,:),      &
                                      shcoef     (:,:) )
         else ! default
            call atmo_boundary_layer (nx_block,  ny_block,   &
                                      'ocn',     icells,     &
                                      indxi,     indxj,      &
                                      sst        (:,:,iblk), &    
                                      potT       (:,:,iblk), &
                                      uatm       (:,:,iblk), &   
                                      vatm       (:,:,iblk), &   
                                      uvel       (:,:,iblk), &   
                                      vvel       (:,:,iblk), &   
                                      wind       (:,:,iblk), &   
                                      zlvl       (:,:,iblk), &   
                                      Qa         (:,:,iblk), &     
                                      rhoa       (:,:,iblk), &
                                      strairx_ocn(:,:,iblk), & 
                                      strairy_ocn(:,:,iblk), & 
                                      Uref_ocn   (:,:,iblk), &
                                      Tref_ocn   (:,:,iblk), & 
                                      Qref_ocn   (:,:,iblk), & 
                                      delt       (:,:),      &    
                                      delq       (:,:),      &
                                      lhcoef     (:,:),      &
                                      shcoef     (:,:) )
         endif

      !-----------------------------------------------------------------
      ! Ocean  albedo
      ! For now, assume albedo = albocn in each spectral band.
      !-----------------------------------------------------------------

         alvdr_ocn(:,:,iblk) = albocn
         alidr_ocn(:,:,iblk) = albocn
         alvdf_ocn(:,:,iblk) = albocn
         alidf_ocn(:,:,iblk) = albocn

      !-----------------------------------------------------------------
      ! Compute ocean fluxes and update SST
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         ! shortwave radiative flux
         swabs = (c1-alvdr_ocn(i,j,iblk)) * swvdr(i,j,iblk) &
               + (c1-alidr_ocn(i,j,iblk)) * swidr(i,j,iblk) &
               + (c1-alvdf_ocn(i,j,iblk)) * swvdf(i,j,iblk) &
               + (c1-alidf_ocn(i,j,iblk)) * swidf(i,j,iblk) 

         ! ocean surface temperature in Kelvin
         TsfK = sst(i,j,iblk) + Tffresh

         ! longwave radiative flux
         flwout_ocn(i,j,iblk) = -stefan_boltzmann * TsfK**4

         ! downward latent and sensible heat fluxes
         fsens_ocn(i,j,iblk) =  shcoef(i,j) * delt(i,j)
         flat_ocn (i,j,iblk) =  lhcoef(i,j) * delq(i,j)
         evap_ocn (i,j,iblk) = -flat_ocn(i,j,iblk) / Lvap

         ! Compute sst change due to exchange with atm/ice above
         sst(i,j,iblk) = sst(i,j,iblk) + dt * ( &
              (fsens_ocn(i,j,iblk) + flat_ocn(i,j,iblk) + flwout_ocn(i,j,iblk) &
             + flw(i,j,iblk) + swabs) * (c1-aice(i,j,iblk)) &
             + fhocn(i,j,iblk) + fswthru(i,j,iblk))         &  ! these are *aice
             / (cprho*hmix(i,j,iblk))

         ! adjust qdp if cooling of mixed layer would occur when sst <= Tf
         if (sst(i,j,iblk) <= Tf(i,j,iblk) .and. qdp(i,j,iblk) > c0) qdp(i,j,iblk) = c0

         ! computed T change due to exchange with deep layers:
         sst(i,j,iblk) = sst(i,j,iblk) - qdp(i,j,iblk)*dt/(cprho*hmix(i,j,iblk))

         ! compute potential to freeze or melt ice
         frzmlt(i,j,iblk) = (Tf(i,j,iblk)-sst(i,j,iblk))*cprho*hmix(i,j,iblk)/dt
         frzmlt(i,j,iblk) = min(max(frzmlt(i,j,iblk),-frzmlt_max),frzmlt_max)

         ! if sst is below freezing, reset sst to Tf
         if (sst(i,j,iblk) <= Tf(i,j,iblk)) sst(i,j,iblk) = Tf(i,j,iblk)

      enddo                     ! ij
      enddo                     ! iblk

      end subroutine ocean_mixed_layer

!=======================================================================

      end module ice_ocean

!=======================================================================
