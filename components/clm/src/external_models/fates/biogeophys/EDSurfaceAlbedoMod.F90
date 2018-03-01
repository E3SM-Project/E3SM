module EDSurfaceRadiationMod
   
   !-------------------------------------------------------------------------------------
   ! EDSurfaceRadiation
   !
   ! This module contains function and type definitions for all things related
   ! to radiative transfer in ED modules at the land surface.
   !
   !-------------------------------------------------------------------------------------

#include "shr_assert.h"
   
  use EDTypesMod        , only : ed_patch_type, ed_site_type
  use EDTypesMod        , only : maxPatchesPerSite
  use EDTypesMod        , only : maxpft
  use FatesConstantsMod , only : r8 => fates_r8
  use FatesInterfaceMod , only : bc_in_type
  use FatesInterfaceMod , only : bc_out_type
  use FatesInterfaceMod , only : hlm_numSWb
  use FatesInterfaceMod , only : numpft
  use EDTypesMod        , only : maxSWb
  use EDTypesMod        , only : nclmax
  use EDTypesMod        , only : nlevleaf
  use EDCanopyStructureMod, only: calc_areaindex
  use FatesGlobals      , only : fates_log
  
  ! CIME globals
  use shr_log_mod       , only : errMsg => shr_log_errMsg

  implicit none

  private
  public :: ED_Norman_Radiation  ! Surface albedo and two-stream fluxes
  public :: ED_SunShadeFracs
  
  logical :: DEBUG = .false.  ! for debugging this module

  
  real(r8), public  :: albice(maxSWb) = &       ! albedo land ice by waveband (1=vis, 2=nir)
        (/ 0.80_r8, 0.55_r8 /)

  ! INTERF-TODO: THIS NEEDS SOME CONSISTENCY AND SHOULD BE SET IN THE INTERFACE
  ! WITH OTHER DIMENSIONS
  integer, parameter :: ipar = 1          ! The band index for PAR

contains
   
  subroutine ED_Norman_Radiation (nsites, sites, bc_in, bc_out )
      !

      !
      ! !USES:
      use EDPftvarcon       , only : EDPftvarcon_inst
      use EDtypesMod        , only : ed_patch_type
      use EDTypesMod        , only : ed_site_type


      ! !ARGUMENTS:
      
      integer,            intent(in)            :: nsites 
      type(ed_site_type), intent(inout), target :: sites(nsites)      ! FATES site vector
      type(bc_in_type),   intent(in)            :: bc_in(nsites)
      type(bc_out_type),  intent(inout)         :: bc_out(nsites)


      ! !LOCAL VARIABLES:
      ! ============================================================================
      ! ED/NORMAN RADIATION DECS
      ! ============================================================================
      type (ed_patch_type) , pointer :: currentPatch
      integer  :: radtype, L, ft, j, ifp
      integer  :: iter                                          ! Iteration index
      integer  :: irep                                          ! Flag to exit iteration loop
      real(r8) :: sb
      real(r8) :: error                                         ! Error check
      real(r8) :: down_rad, up_rad                              ! Iterative solution do Dif_dn and Dif_up
      real(r8) :: ftweight(nclmax,maxpft,nlevleaf)
      real(r8) :: k_dir(maxpft)                              ! Direct beam extinction coefficient
      real(r8) :: tr_dir_z(nclmax,maxpft,nlevleaf)         ! Exponential transmittance of direct beam radiation through a single layer
      real(r8) :: tr_dif_z(nclmax,maxpft,nlevleaf)         ! Exponential transmittance of diffuse radiation through a single layer
      real(r8) :: forc_dir(maxPatchesPerSite,maxSWb)
      real(r8) :: forc_dif(maxPatchesPerSite,maxSWb)
      real(r8) :: weighted_dir_tr(nclmax)
      real(r8) :: weighted_fsun(nclmax)
      real(r8) :: weighted_dif_ratio(nclmax,maxSWb)
      real(r8) :: weighted_dif_down(nclmax)
      real(r8) :: weighted_dif_up(nclmax)
      real(r8) :: refl_dif(nclmax,maxpft,nlevleaf,maxSWb)  ! Term for diffuse radiation reflected by laye
      real(r8) :: tran_dif(nclmax,maxpft,nlevleaf,maxSWb)  ! Term for diffuse radiation transmitted by layer
      real(r8) :: dif_ratio(nclmax,maxpft,nlevleaf,maxSWb) ! Ratio of upward to forward diffuse fluxes
      real(r8) :: Dif_dn(nclmax,maxpft,nlevleaf)           ! Forward diffuse flux onto canopy layer J (W/m**2 ground area)
      real(r8) :: Dif_up(nclmax,maxpft,nlevleaf)           ! Upward diffuse flux above canopy layer J (W/m**2 ground area)
      real(r8) :: lai_change(nclmax,maxpft,nlevleaf)       ! Forward diffuse flux onto canopy layer J (W/m**2 ground area)
      real(r8) :: f_not_abs(maxpft,maxSWb)                   ! Fraction reflected + transmitted. 1-absorbtion.
      real(r8) :: Abs_dir_z(maxpft,nlevleaf)
      real(r8) :: Abs_dif_z(maxpft,nlevleaf)
      real(r8) :: abs_rad(maxSWb)                               !radiation absorbed by soil
      real(r8) :: tr_soili                                      ! Radiation transmitted to the soil surface.
      real(r8) :: tr_soild                                      ! Radiation transmitted to the soil surface.
      real(r8) :: phi1b(maxPatchesPerSite,maxpft)      ! Radiation transmitted to the soil surface.
      real(r8) :: phi2b(maxPatchesPerSite,maxpft)
      real(r8) :: laisum                                        ! cumulative lai+sai for canopy layer (at middle of layer)
      real(r8) :: angle

      real(r8),parameter :: tolerance = 0.000000001_r8
      real(r8), parameter :: pi   = 3.141592654                 ! PI

      
      integer, parameter :: max_diag_nlevleaf = 4
      integer, parameter :: diag_nlevleaf = min(nlevleaf,max_diag_nlevleaf)  ! for diagnostics, write a small number of leaf layers

      real(r8) :: denom
      real(r8) :: lai_reduction(2)

      integer  :: fp,iv,s      ! array indices
      integer  :: ib               ! waveband number
      real(r8) :: cosz             ! 0.001 <= coszen <= 1.000
      real(r8) :: chil(maxPatchesPerSite)     ! -0.4 <= xl <= 0.6
      real(r8) :: gdir(maxPatchesPerSite)    ! leaf projection in solar direction (0 to 1)

      !-----------------------------------------------------------------------

      associate(&
            rhol         =>    EDPftvarcon_inst%rhol                     , & ! Input:  [real(r8) (:)   ] leaf reflectance: 1=vis, 2=nir
            rhos         =>    EDPftvarcon_inst%rhos                     , & ! Input:  [real(r8) (:)   ] stem reflectance: 1=vis, 2=nir
            taul         =>    EDPftvarcon_inst%taul                     , & ! Input:  [real(r8) (:)   ] leaf transmittance: 1=vis, 2=nir
            taus         =>    EDPftvarcon_inst%taus                     , & ! Input:  [real(r8) (:)   ] stem transmittance: 1=vis, 2=nir
            xl           =>    EDPftvarcon_inst%xl)                          ! Input:  [real(r8) (:)   ] ecophys const - leaf/stem orientation index
            
!            albd         =>    surfalb_inst%albd_patch         , & ! Output: [real(r8) (:,:) ] surface albedo (direct) (USED IN LND2ATM,BALANCE_CHECK)
!            albi         =>    surfalb_inst%albi_patch         , & ! Output: [real(r8) (:,:) ] surface albedo (diffuse) (LND2ATM,BALANCE_CHECK)
!            fabd         =>    surfalb_inst%fabd_patch         , & ! Output: [real(r8) (:,:) ] flux absorbed by canopy per unit direct flux (BALANCE_CHECK)
!            fabi         =>    surfalb_inst%fabi_patch         , & ! Output: [real(r8) (:,:) ] flux absorbed by canopy per unit diffuse flux (BALANCE_CHECK)
!            ftdd         =>    surfalb_inst%ftdd_patch         , & ! Output: [real(r8) (:,:) ] down direct flux below canopy per unit direct flx (BALANCE_CHECK)
!            ftid         =>    surfalb_inst%ftid_patch         , & ! Output: [real(r8) (:,:) ] down diffuse flux below canopy per unit direct flx (BALANCE_CHECK)
!            ftii         =>    surfalb_inst%ftii_patch         , & ! Output: [real(r8) (:,:) ] down diffuse flux below canopy per unit diffuse flx (BALANCE_CHECK)

        ! -------------------------------------------------------------------------------
        ! TODO (mv, 2014-10-29) the filter here is different than below 
        ! this is needed to have the VOC's be bfb - this needs to be
        ! re-examined int he future
        ! RGK,2016-08-06: FATES is still incompatible with VOC emission module
        ! -------------------------------------------------------------------------------


        do s = 1, nsites

           ifp = 0
           currentpatch => sites(s)%oldest_patch
           do while (associated(currentpatch))  
              ifp = ifp+1
              
              currentPatch%f_sun      (:,:,:) = 0._r8
              currentPatch%fabd_sun_z (:,:,:) = 0._r8
              currentPatch%fabd_sha_z (:,:,:) = 0._r8
              currentPatch%fabi_sun_z (:,:,:) = 0._r8
              currentPatch%fabi_sha_z (:,:,:) = 0._r8
              currentPatch%fabd       (:)     = 0._r8
              currentPatch%fabi       (:)     = 0._r8

              if(bc_in(s)%filter_vegzen_pa(ifp))then

                 weighted_dir_tr(:)   = 0._r8
                 weighted_dif_down(:) = 0._r8
                 weighted_dif_up(:)   = 0._r8
                 bc_out(s)%albd_parb(ifp,:)            = 0._r8  ! output HLM
                 bc_out(s)%albi_parb(ifp,:)            = 0._r8  ! output HLM
                 bc_out(s)%fabi_parb(ifp,:)            = 0._r8  ! output HLM
                 bc_out(s)%fabd_parb(ifp,:)            = 0._r8  ! output HLM
                 tr_dir_z(:,:,:)      = 0._r8   
                 tr_dif_z(:,:,:)      = 0._r8
                 ftweight(:,:,:)      = 0._r8
                 lai_change(:,:,:)    = 0._r8
                 Dif_up(:,:,:)        = 0._r8
                 Dif_dn(:,:,:)        = 0._r8
                 refl_dif(:,:,:,:)    = 0.0_r8
                 tran_dif(:,:,:,:)    = 0.0_r8
                 dif_ratio(:,:,:,:)   = 0.0_r8
                 bc_out(s)%ftdd_parb(ifp,:)            = 1._r8 ! output HLM
                 bc_out(s)%ftid_parb(ifp,:)            = 1._r8 ! output HLM
                 bc_out(s)%ftii_parb(ifp,:)            = 1._r8 ! output HLM
                 
                 if (maxval(currentPatch%nrad(1,:))==0)then
                    !there are no leaf layers in this patch. it is effectively bare ground. 
                    ! no radiation is absorbed  
                    bc_out(s)%fabd_parb(ifp,:) = 0.0_r8
                    bc_out(s)%fabi_parb(ifp,:) = 0.0_r8
                    do ib = 1,hlm_numSWb
                       bc_out(s)%albd_parb(ifp,ib) = bc_in(s)%albgr_dir_rb(ib)
                       bc_out(s)%albd_parb(ifp,ib) = bc_in(s)%albgr_dif_rb(ib)
                       bc_out(s)%ftdd_parb(ifp,ib)= 1.0_r8
                       bc_out(s)%ftid_parb(ifp,ib)= 1.0_r8
                       bc_out(s)%ftii_parb(ifp,ib)= 1.0_r8
                    enddo
                 else
                    
                    ! Is this pft/canopy layer combination present in this patch?
                    do L = 1,nclmax
                       do ft = 1,numpft
                          currentPatch%present(L,ft) = 0
                          do  iv = 1, currentPatch%nrad(L,ft)
                             if (currentPatch%canopy_area_profile(L,ft,iv) > 0._r8)then
                                currentPatch%present(L,ft) = 1
                                !I think 'present' is only used here...
                             endif
                          end do !iv
                       end do !ft
                    end do !L

                    do radtype = 1,2 !do this once for one unit of diffuse, and once for one unit of direct radiation
                       do ib = 1,hlm_numSWb
                          if (radtype == 1) then
                             ! Set the hypothetical driving radiation. We do this once for a single unit of direct and
                             ! once for a single unit of diffuse radiation.
                             forc_dir(ifp,ib) = 1.00_r8
                             forc_dif(ifp,ib) = 0.00_r8
                          else !dif
                             forc_dir(ifp,ib) = 0.00_r8
                             forc_dif(ifp,ib) = 1.00_r8
                          end if
                       end do !ib

                       !Extract information that needs to be provided by ED into local array.
                       ftweight(:,:,:) = 0._r8
                       do L = 1,currentPatch%NCL_p
                          do ft = 1,numpft
                             do  iv = 1, currentPatch%nrad(L,ft)
                                !this is already corrected for area in CLAP
                                ftweight(L,ft,iv) = currentPatch%canopy_area_profile(L,ft,iv) 
                             end do  !iv
                          end do  !ft1
                       end do  !L
                       if (sum(ftweight(1,:,1))<0.999_r8)then
                          write(fates_log(),*) 'canopy not full',ftweight(1,:,1)
                       endif
                       if (sum(ftweight(1,:,1))>1.0001_r8)then
                          write(fates_log(),*) 'canopy too full',ftweight(1,:,1)
                       endif

                       !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
                       ! Direct beam extinction coefficient, k_dir. PFT specific.
                       !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
                       cosz = max(0.001_r8, bc_in(s)%coszen_pa(ifp)) !copied from previous radiation code...
                       do ft = 1,numpft
                          sb = (90._r8 - (acos(cosz)*180/pi)) * (pi / 180._r8)
                          chil(ifp) = xl(ft) !min(max(xl(ft), -0.4_r8), 0.6_r8 )
                          if (abs(chil(ifp)) <= 0.01_r8) then
                             chil(ifp) = 0.01_r8
                          end if
                          phi1b(ifp,ft) = 0.5_r8 - 0.633_r8*chil(ifp) - 0.330_r8*chil(ifp)*chil(ifp)
                          phi2b(ifp,ft) = 0.877_r8 * (1._r8 - 2._r8*phi1b(ifp,ft)) !0 = horiz leaves, 1 - vert leaves.
                          gdir(ifp) = phi1b(ifp,ft) + phi2b(ifp,ft) * sin(sb)
                          !how much direct light penetrates a singleunit of lai?
                          k_dir(ft) = gdir(ifp) / sin(sb)
                       end do !FT
                       
                       do L = 1,currentPatch%NCL_p !start at the top canopy layer (1 is the top layer.)
                          weighted_dir_tr(L) = 0.0_r8
                          weighted_fsun(L) = 0._r8
                          weighted_dif_ratio(L,1:hlm_numSWb) = 0._r8
                          !Each canopy layer (canopy, understorey) has multiple 'parallel' pft's
                          do ft =1,numpft
                             if (currentPatch%present(L,ft) == 1)then !only do calculation if there are the appropriate leaves.
                                !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
                                ! Diffuse transmittance, tr_dif, do each layer with thickness elai_z.
                                ! Estimated do nine sky angles in increments of 10 degrees
                                ! PFT specific...
                                !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
                                tr_dif_z(L,ft,:) = 0._r8
                                do iv = 1,currentPatch%nrad(L,ft)
                                   do j = 1,9
                                      angle = (5._r8 + (j - 1) * 10._r8) * 3.142 / 180._r8
                                      gdir(ifp) = phi1b(ifp,ft) + phi2b(ifp,ft) * sin(angle) !This line is redundant FIX(RF,032414). 
                                      tr_dif_z(L,ft,iv) = tr_dif_z(L,ft,iv) + exp(-gdir(ifp) / sin(angle) * &
                                           (currentPatch%elai_profile(L,ft,iv)+currentPatch%esai_profile(L,ft,iv))) * &
                                           sin(angle)*cos(angle)
                                   end do
                                   
                                   tr_dif_z(L,ft,iv) = tr_dif_z(L,ft,iv) * 2._r8 * (10.00*pi/180._r8)
                                   
                                end do

                                !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
                                ! Direct beam transmittance, tr_dir_z, uses cumulative LAI above layer J to give
                                ! unscattered direct beam onto layer J. do each PFT section.
                                ! This is just an  decay curve based on k_dir. (leaf & sun angle)
                                !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
                                if (L==1)then
                                   tr_dir_z(L,ft,1) = 1._r8
                                else   
                                   tr_dir_z(L,ft,1)  = weighted_dir_tr(L-1)
                                endif
                                laisum = 0.00_r8                         
                                !total direct beam getting to the bottom of the top canopy.
                                do iv = 1,currentPatch%nrad(L,ft)
                                   laisum = laisum + currentPatch%elai_profile(L,ft,iv)+currentPatch%esai_profile(L,ft,iv)
                                   lai_change(L,ft,iv) = 0.0_r8
                                   if (( ftweight(L,ft,iv+1)  >  0.0_r8 ) .and. ( ftweight(L,ft,iv+1)  <  ftweight(L,ft,iv) ))then
                                      !where there is a partly empty leaf layer, some fluxes go straight through.
                                      lai_change(L,ft,iv) = ftweight(L,ft,iv)-ftweight(L,ft,iv+1)
                                   endif
                                   if (ftweight(L,ft,iv+1) - ftweight(L,ft,iv) > 1.e-10_r8)then
                                      write(fates_log(),*) 'lower layer has more coverage. This is wrong' , &
                                           ftweight(L,ft,iv),ftweight(L,ft,iv+1),ftweight(L,ft,iv+1)-ftweight(L,ft,iv)
                                   endif
                                   
                                   !n.b. in theory lai_change could be calculated daily in the ED code.
                                   !This is light coming striaght through the canopy.
                                   if (L==1)then
                                      tr_dir_z(L,ft,iv+1) = exp(-k_dir(ft) * laisum)* &
                                           (ftweight(L,ft,iv)/ftweight(L,ft,1))
                                   else   
                                      tr_dir_z(L,ft,iv+1) = weighted_dir_tr(L-1)*exp(-k_dir(ft) * laisum)* &
                                           (ftweight(L,ft,iv)/ftweight(L,ft,1))
                                   endif
                                   
                                   if (iv == 1)then
                                      !this is the top layer.
                                      tr_dir_z(L,ft,iv+1) = tr_dir_z(L,ft,iv+1) + tr_dir_z(L,ft,iv) * &
                                           ((ftweight(L,ft,1)-ftweight(L,ft,iv))/ftweight(L,ft,1))                               
                                   else
                                      !the lai_change(iv) affects the light incident on layer iv+2 not iv+1
                                      ! light coming from the layer above (iv-1) goes through iv and onto iv+1. 
                                      if (lai_change(L,ft,iv-1) > 0.0_r8)then
                                         tr_dir_z(L,ft,iv+1) = tr_dir_z(L,ft,iv+1) + tr_dir_z(L,ft,iv)* &
                                              lai_change(L,ft,iv-1) / ftweight(L,ft,1)
                                         tr_dir_z(L,ft,iv+1) = tr_dir_z(L,ft,iv+1) + tr_dir_z(L,ft,iv-1)* &
                                              (ftweight(L,ft,1)-ftweight(L,ft,iv-1))/ftweight(L,ft,1)
                                      else
                                         !account fot the light that comes striaght down from unfilled layers above.
                                         tr_dir_z(L,ft,iv+1) = tr_dir_z(L,ft,iv+1) + tr_dir_z(L,ft,iv) * &
                                              ((ftweight(L,ft,1)-ftweight(L,ft,iv))/ftweight(L,ft,1))
                                      endif
                                   endif
                                end do
                                
                                !add up all the weighted contributions from the different PFT columns.
                                weighted_dir_tr(L) = weighted_dir_tr(L) + tr_dir_z(L,ft,currentPatch%nrad(L,ft)+1)*ftweight(L,ft,1)

                                !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
                                ! Sunlit and shaded fraction of leaf layer
                                !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
                                
                                !laisum = 0._r8
                                do iv = 1,currentPatch%nrad(L,ft)
                                   ! Cumulative leaf area. Original code uses cumulative lai do layer.
                                   ! Now use cumulative lai at center of layer.
                                   ! Same as tr_dir_z calcualtions, but in the middle of the layer? FIX(RF,032414)-WHY?
                                   if (iv  ==  1) then
                                      laisum = 0.5_r8 * (currentPatch%elai_profile(L,ft,iv)+currentPatch%esai_profile(L,ft,iv))
                                   else
                                      laisum = laisum + currentPatch%elai_profile(L,ft,iv)+currentPatch%esai_profile(L,ft,iv)
                                   end if
                                   
                                   
                                   if (L == 1)then !top canopy layer
                                      currentPatch%f_sun(L,ft,iv) = exp(-k_dir(ft) * laisum)* &
                                           (ftweight(L,ft,iv)/ftweight(L,ft,1))
                                   else
                                      currentPatch%f_sun(L,ft,iv) = weighted_fsun(L-1)* exp(-k_dir(ft) * laisum)* &
                                           (ftweight(L,ft,iv)/ftweight(L,ft,1))
                                   endif
                                   
                                   if ( iv > 1 ) then  ! becasue we are looking at this layer (not the next)
                                      ! we only ever add fluxes if iv>1
                                      if (lai_change(L,ft,iv-1) > 0.0_r8)then
                                         currentPatch%f_sun(L,ft,iv) = currentPatch%f_sun(L,ft,iv) + &
                                              currentPatch%f_sun(L,ft,iv) * &
                                              lai_change(L,ft,iv-1)/ftweight(L,ft,1)
                                         currentPatch%f_sun(L,ft,iv) = currentPatch%f_sun(L,ft,iv) + &
                                              currentPatch%f_sun(L,ft,iv-1) * &
                                              (ftweight(L,ft,1)-ftweight(L,ft,iv-1))/ftweight(L,ft,1)
                                      else
                                         currentPatch%f_sun(L,ft,iv) = currentPatch%f_sun(L,ft,iv) + &
                                              currentPatch%f_sun(L,ft,iv-1) * &
                                              (ftweight(L,ft,1)-ftweight(L,ft,iv))/ftweight(L,ft,1)
                                      endif
                                   endif
                                   
                                end do !iv
                                weighted_fsun(L) = weighted_fsun(L) + currentPatch%f_sun(L,ft,currentPatch%nrad(L,ft))* &
                                     ftweight(L,ft,1)
                                
                                ! instance where the first layer ftweight is used a proxy for the whole column. FTWA
                                ! this is possibly a source of slight error. If we use the ftweight at the top of the PFT column,
                                ! then we willl underestimate fsun, but if we use ftweight at the bottom of the column, we will
                                ! underestimate it. Really, we should be tracking the release of direct light from the column as it tapers
                                ! towards the ground. Is that necessary to get energy closure? It would be quite hard...
                             endif !present.
                          end do!pft loop
                       end do !L
                       
                       do L = currentPatch%NCL_p,1, -1 !start at the bottom and work up.
                          do ft = 1,numpft
                             if (currentPatch%present(L,ft) == 1)then
                                !==============================================================================!
                                ! Iterative solution do scattering
                                !==============================================================================!
                                
                                do ib = 1,hlm_numSWb !vis, nir
                                   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
                                   ! Leaf scattering coefficient and terms do diffuse radiation reflected
                                   ! and transmitted by a layer
                                   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
                                   f_not_abs(ft,ib) = rhol(ft,ib) + taul(ft,ib) !leaf level fraction NOT absorbed.
                                   !tr_dif_z is a term that uses the LAI in each layer, whereas rhol and taul do not,
                                   !because they are properties of leaf surfaces and not of the leaf matrix.
                                   do iv = 1,currentPatch%nrad(L,ft)
                                      !How much diffuse light is intercepted and then reflected?
                                      refl_dif(L,ft,iv,ib) = (1._r8 - tr_dif_z(L,ft,iv)) * rhol(ft,ib)
                                      !How much diffuse light in this layer is transmitted?
                                      tran_dif(L,ft,iv,ib) = (1._r8 - tr_dif_z(L,ft,iv)) * &
                                            taul(ft,ib) + tr_dif_z(L,ft,iv)
                                   end do
                                   
                                   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
                                   ! Ratio of upward to forward diffuse fluxes, dif_ratio
                                   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
                                   ! Soil diffuse reflectance (ratio of down to up radiation).
                                   iv = currentPatch%nrad(L,ft) + 1
                                   if (L  == currentPatch%NCL_p)then !nearest the soil
                                      dif_ratio(L,ft,iv,ib) = bc_in(s)%albgr_dif_rb(ib)
                                   else
                                      dif_ratio(L,ft,iv,ib) = weighted_dif_ratio(L+1,ib)
                                   end if
                                   ! Canopy layers, working upwardfrom soil with dif_ratio(iv+1) known
                                   ! FIX(RF,032414) ray tracing eqution - need to find derivation of this...
                                   ! for each unit going down, there are x units going up.
                                   do iv = currentPatch%nrad(L,ft),1, -1
                                      dif_ratio(L,ft,iv,ib) = dif_ratio(L,ft,iv+1,ib) * &
                                            tran_dif(L,ft,iv,ib)*tran_dif(L,ft,iv,ib) / &
                                            (1._r8 - dif_ratio(L,ft,iv+1,ib) * refl_dif(L,ft,iv,ib)) &
                                            + refl_dif(L,ft,iv,ib)
                                      dif_ratio(L,ft,iv,ib) = dif_ratio(L,ft,iv,ib) * &
                                            ftweight(L,ft,iv)/ftweight(L,ft,1)
                                      dif_ratio(L,ft,iv,ib) = dif_ratio(L,ft,iv,ib) + dif_ratio(L,ft,iv+1,ib) * &
                                           (ftweight(L,ft,1)-ftweight(L,ft,iv))/ftweight(L,ft,1)
                                   end do
                                   weighted_dif_ratio(L,ib) = weighted_dif_ratio(L,ib) + &
                                         dif_ratio(L,ft,1,ib) * ftweight(L,ft,1)
                                   !instance where the first layer ftweight is used a proxy for the whole column. FTWA
                                end do!hlm_numSWb
                             endif ! currentPatch%present
                          end do!ft
                       end do!L
                       
                       do ib = 1,hlm_numSWb
                          Dif_dn(:,:,:) = 0.00_r8
                          Dif_up(:,:,:) = 0.00_r8
                          do L = 1, currentPatch%NCL_p !work down from the top of the canopy.
                             weighted_dif_down(L) = 0._r8
                             do ft = 1, numpft
                                if (currentPatch%present(L,ft) == 1)then
                                   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
                                   ! First estimates do downward and upward diffuse flux
                                   !
                                   ! Dif_dn =  forward diffuse flux onto layer J
                                   ! Dif_up =  Upward diffuse flux above layer J
                                   !
                                   ! Solved here without direct beam radiation and using dif_ratio = Dif_up / Dif_dn
                                   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
                                   ! downward diffuse flux onto the top surface of the canopy
                                   
                                   if (L == 1)then
                                      Dif_dn(L,ft,1) = forc_dif(ifp,ib)
                                   else
                                      Dif_dn(L,ft,1) = weighted_dif_down(L-1)
                                   end if
                                   ! forward diffuse flux within the canopy and at soil, working forward through canopy
                                   do iv = 1,currentPatch%nrad(L,ft)
                                      denom = refl_dif(L,ft,iv,ib) *  dif_ratio(L,ft,iv,ib)
                                      denom = 1._r8 - denom
                                      Dif_dn(L,ft,iv+1) = Dif_dn(L,ft,iv) * tran_dif(L,ft,iv,ib) / &
                                           denom *ftweight(L,ft,iv)/ftweight(L,ft,1)
                                      if (iv > 1)then
                                         if (lai_change(L,ft,iv-1) > 0.0_r8)then
                                            !here we are thinking about whether the layer above had an laichange,
                                            !but calculating the flux onto the layer below.
                                            Dif_dn(L,ft,iv+1) = Dif_dn(L,ft,iv+1)+ Dif_dn(L,ft,iv)* &
                                                 lai_change(L,ft,iv-1)/ftweight(L,ft,1)
                                            Dif_dn(L,ft,iv+1) = Dif_dn(L,ft,iv+1)+ Dif_dn(L,ft,iv-1)* &
                                                 (ftweight(L,ft,1)-ftweight(L,ft,iv-1)/ftweight(L,ft,1))
                                         else
                                            Dif_dn(L,ft,iv+1) = Dif_dn(L,ft,iv+1) + Dif_dn(L,ft,iv) * &
                                                 (ftweight(L,ft,1)-ftweight(L,ft,iv))/ftweight(L,ft,1)
                                         endif
                                      else
                                         Dif_dn(L,ft,iv+1)    = Dif_dn(L,ft,iv+1) + Dif_dn(L,ft,iv) * &
                                              (ftweight(L,ft,1)-ftweight(L,ft,iv))/ftweight(L,ft,1)
                                      endif
                                   end do
                                   
                                   weighted_dif_down(L) = weighted_dif_down(L) + Dif_dn(L,ft,currentPatch%nrad(L,ft)+1) * &
                                        ftweight(L,ft,1)
                                   
                                   !instance where the first layer ftweight is used a proxy for the whole column. FTWA
                                endif !present
                             end do !ft
                             if (L == currentPatch%NCL_p.and.currentPatch%NCL_p > 1)then !is the the (incomplete) understorey?
                                !Add on the radiation going through the canopy gaps.
                                weighted_dif_down(L) = weighted_dif_down(L) + weighted_dif_down(L-1)*(1.0-sum(ftweight(L,:,1)))
                                !instance where the first layer ftweight is used a proxy for the whole column. FTWA
                             endif
                          end do !L
                          
                          do L = currentPatch%NCL_p,1 ,-1 !work up from the bottom.
                             weighted_dif_up(L) = 0._r8
                             do ft = 1, numpft
                                if (currentPatch%present(L,ft) == 1)then
                                   !Bounce diffuse radiation off soil surface.
                                   iv = currentPatch%nrad(L,ft) + 1
                                   if (L==currentPatch%NCL_p)then !is this the bottom layer ?
                                      Dif_up(L,ft,iv) =bc_in(s)%albgr_dif_rb(ib) * Dif_dn(L,ft,iv)
                                   else
                                      Dif_up(L,ft,iv) = weighted_dif_up(L+1)
                                   end if
                                   ! Upward diffuse flux within the canopy and above the canopy, working upward through canopy
                                   
                                   do iv = currentPatch%nrad(L,ft), 1, -1
                                      if (lai_change(L,ft,iv) > 0.0_r8)then
                                         Dif_up(L,ft,iv) = dif_ratio(L,ft,iv,ib) * Dif_dn(L,ft,iv) * &
                                              ftweight(L,ft,iv) / ftweight(L,ft,1)
                                         Dif_up(L,ft,iv) = Dif_up(L,ft,iv) + Dif_up(L,ft,iv+1) * &
                                              tran_dif(L,ft,iv,ib) * lai_change(L,ft,iv)/ftweight(L,ft,1)
                                         Dif_up(L,ft,iv) = Dif_up(L,ft,iv) + Dif_up(L,ft,iv+1) * &
                                              (ftweight(L,ft,1)-ftweight(L,ft,iv))/ftweight(L,ft,1)
                                         !nb is this the right constuction?
                                         ! the radiation that hits the empty space is not reflected.
                                      else
                                         Dif_up(L,ft,iv) = dif_ratio(L,ft,iv,ib) * Dif_dn(L,ft,iv) * ftweight(L,ft,iv)
                                         Dif_up(L,ft,iv) = Dif_up(L,ft,iv) + Dif_up(L,ft,iv+1) * (1.0_r8-ftweight(L,ft,iv))
                                      endif
                                   end do
                                   
                                   weighted_dif_up(L) = weighted_dif_up(L) + Dif_up(L,ft,1) * ftweight(L,ft,1)
                                   !instance where the first layer ftweight is used a proxy for the whole column. FTWA
                                endif !present
                             end do !ft
                             if (L == currentPatch%NCL_p.and.currentPatch%NCL_p > 1)then !is this the (incomplete) understorey?
                                !Add on the radiation coming up through the canopy gaps.
                                !diffuse to diffuse
                                weighted_dif_up(L) = weighted_dif_up(L) +(1.0-sum(ftweight(L,1:numpft,1))) * &
                                     weighted_dif_down(L-1) * bc_in(s)%albgr_dif_rb(ib)
                                !direct to diffuse
                                weighted_dif_up(L) = weighted_dif_up(L) + forc_dir(ifp,ib) * &
                                     weighted_dir_tr(L-1) * (1.0-sum(ftweight(L,1:numpft,1)))*bc_in(s)%albgr_dir_rb(ib)
                             endif
                          end do !L
                          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
                          ! 3. Iterative calculation of forward and upward diffuse fluxes, iNCL_puding
                          ! scattered direct beam
                          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

                          ! Flag to exit iteration loop: 0 = exit and 1 = iterate
                          irep = 1
                          ! Iteration loop
                          iter = 0
                          do while(irep ==1 .and. iter<50)
                             
                             iter = iter + 1
                             irep = 0
                             do L = 1,currentPatch%NCL_p !working from the top down
                                weighted_dif_down(L) = 0._r8
                                do ft =1,numpft
                                   if (currentPatch%present(L,ft) == 1)then
                                      ! forward diffuse flux within the canopy and at soil, working forward through canopy
                                      ! with Dif_up -from previous iteration-. Dif_dn(1) is the forward diffuse flux onto the canopy.
                                      ! Note: down = forward flux onto next layer
                                      if (L == 1)then !is this the top layer?
                                         Dif_dn(L,ft,1) = forc_dif(ifp,ib)
                                      else
                                         Dif_dn(L,ft,1) = weighted_dif_down(L-1)
                                      end if
                                      down_rad = 0._r8
                                      
                                      do iv = 1, currentPatch%nrad(L,ft)
                                         
                                         down_rad = Dif_dn(L,ft,iv) * tran_dif(L,ft,iv,ib) + &
                                              Dif_up(L,ft,iv+1) * refl_dif(L,ft,iv,ib) + &
                                              forc_dir(ifp,ib) * tr_dir_z(L,ft,iv) * (1.00_r8 - &
                                              exp(-k_dir(ft) * (currentPatch%elai_profile(L,ft,iv)+ &
                                              currentPatch%esai_profile(L,ft,iv)))) * taul(ft,ib)
                                         down_rad = down_rad *(ftweight(L,ft,iv)/ftweight(L,ft,1))
                                         
                                         if (iv > 1)then
                                            if (lai_change(L,ft,iv-1) > 0.0_r8)then
                                               down_rad = down_rad + Dif_dn(L,ft,iv)   * lai_change(L,ft,iv-1)/ftweight(L,ft,1)
                                               down_rad = down_rad + Dif_dn(L,ft,iv-1) * (ftweight(L,ft,1)-ftweight(L,ft,iv-1))/ &
                                                    ftweight(L,ft,1)
                                            else
                                               down_rad = down_rad + Dif_dn(L,ft,iv)   * (ftweight(L,ft,1)-ftweight(L,ft,iv))/ &
                                                    ftweight(L,ft,1)
                                            endif
                                         else
                                            down_rad = down_rad + Dif_dn(L,ft,iv)   * (ftweight(L,ft,1)-ftweight(L,ft,iv))/ &
                                                 ftweight(L,ft,1)
                                         endif
                                         
                                         !this is just Dif down, plus refl up, plus dir intercepted and turned into dif... ,
                                         if (abs(down_rad - Dif_dn(L,ft,iv+1)) > tolerance)then
                                            irep = 1
                                         end if
                                         Dif_dn(L,ft,iv+1) = down_rad
                                         
                                      end do !iv
                                      
                                      weighted_dif_down(L) = weighted_dif_down(L) + Dif_dn(L,ft,currentPatch%nrad(L,ft)+1) * &
                                           ftweight(L,ft,1)

                                   endif !present
                                end do!ft
                                if (L == currentPatch%NCL_p.and.currentPatch%NCL_p > 1)then !is this the (incomplete) understorey?
                                   weighted_dif_down(L) = weighted_dif_down(L) + weighted_dif_down(L-1) * &
                                         (1.0-sum(ftweight(L,1:numpft,1)))
                                end if
                             end do ! do L loop

                             do L = 1, currentPatch%NCL_p ! working from the top down.
                                weighted_dif_up(L) = 0._r8
                                do ft =1,numpft
                                   if (currentPatch%present(L,ft) == 1)then
                                      ! Upward diffuse flux at soil or from lower canopy (forward diffuse and unscattered direct beam)
                                      iv = currentPatch%nrad(L,ft) + 1
                                      if (L==currentPatch%NCL_p)then  !In the bottom canopy layer, reflect off the soil
                                         Dif_up(L,ft,iv) = Dif_dn(L,ft,iv) *bc_in(s)%albgr_dif_rb(ib) + &
                                              forc_dir(ifp,ib) * tr_dir_z(L,ft,iv) *bc_in(s)%albgr_dir_rb(ib)
                                      else      !In the other canopy layers, reflect off the underlying vegetation.
                                         Dif_up(L,ft,iv) =  weighted_dif_up(L+1)
                                      end if
                                      
                                      ! Upward diffuse flux within and above the canopy, working upward through canopy
                                      ! with Dif_dn from previous interation.  Note: up = upward flux above current layer
                                      do iv = currentPatch%nrad(L,ft),1,-1
                                         !this is radiation up, by layer transmittance, by
                                         
                                         !reflection of the lower layer,
                                         up_rad = Dif_dn(L,ft,iv) * refl_dif(L,ft,iv,ib)
                                         up_rad = up_rad + forc_dir(ifp,ib) * tr_dir_z(L,ft,iv) * (1.00_r8 - exp(-k_dir(ft) * &
                                              (currentPatch%elai_profile(L,ft,iv) + currentPatch%esai_profile(L,ft,iv)))) * &
                                              rhol(ft,ib)
                                         up_rad = up_rad + Dif_up(L,ft,iv+1) * tran_dif(L,ft,iv,ib)
                                         up_rad = up_rad * ftweight(L,ft,iv)/ftweight(L,ft,1)
                                         up_rad = up_rad + Dif_up(L,ft,iv+1) *(ftweight(L,ft,1)-ftweight(L,ft,iv))/ftweight(L,ft,1)
                                         ! THE LOWER LAYER FLUX IS HOMOGENIZED, SO WE DON"T CONSIDER THE LAI_CHANGE HERE...
                                         
                                         if (abs(up_rad - Dif_up(L,ft,iv)) > tolerance) then !are we close to the tolerance level?
                                            irep = 1
                                         end if
                                         Dif_up(L,ft,iv) = up_rad
                                         
                                      end do  !iv
                                      weighted_dif_up(L) = weighted_dif_up(L) + Dif_up(L,ft,1) * ftweight(L,ft,1)
                                   end if !present
                                end do!ft

                                if (L == currentPatch%NCL_p.and.currentPatch%NCL_p > 1)then  !is this the (incomplete) understorey?
                                   !Add on the radiation coming up through the canopy gaps.
                                   weighted_dif_up(L) = weighted_dif_up(L) +(1.0_r8-sum(ftweight(L,1:numpft,1))) * &
                                        weighted_dif_down(L-1) * bc_in(s)%albgr_dif_rb(ib)
                                   weighted_dif_up(L) = weighted_dif_up(L) + forc_dir(ifp,ib) * &
                                        weighted_dir_tr(L-1) * (1.0_r8-sum(ftweight(L,1:numpft,1)))*bc_in(s)%albgr_dir_rb(ib)
                                end if
                             end do!L
                          end do ! do while over iter
                          
                          abs_rad(ib) = 0._r8
                          tr_soili = 0._r8
                          tr_soild = 0._r8
                          do L = 1, currentPatch%NCL_p !working from the top down.
                             abs_dir_z(:,:) = 0._r8
                             abs_dif_z(:,:) = 0._r8
                             do ft =1,numpft
                                if (currentPatch%present(L,ft) == 1)then
                                   !==============================================================================!
                                   ! Compute absorbed flux densities
                                   !==============================================================================!
                                   
                                   ! Absorbed direct beam and diffuse do leaf layers
                                   do iv = 1, currentPatch%nrad(L,ft)
                                      Abs_dir_z(ft,iv) = ftweight(L,ft,iv)* forc_dir(ifp,ib) * tr_dir_z(L,ft,iv) * &
                                           (1.00_r8 - exp(-k_dir(ft) * (currentPatch%elai_profile(L,ft,iv)+ &
                                           currentPatch%esai_profile(L,ft,iv)))) * (1.00_r8 - f_not_abs(ft,ib))
                                      Abs_dif_z(ft,iv) = ftweight(L,ft,iv)* ((Dif_dn(L,ft,iv) + &
                                           Dif_up(L,ft,iv+1)) * (1.00_r8 - tr_dif_z(L,ft,iv)) * &
                                           (1.00_r8 - f_not_abs(ft,ib)))
                                   end do
                                   
                                   ! Absorbed direct beam and diffuse do soil
                                   if (L == currentPatch%NCL_p)then
                                      iv = currentPatch%nrad(L,ft) + 1
                                      Abs_dif_z(ft,iv) = ftweight(L,ft,1)*Dif_dn(L,ft,iv) * (1.0_r8 -bc_in(s)%albgr_dif_rb(ib))
                                      Abs_dir_z(ft,iv) = ftweight(L,ft,1)*forc_dir(ifp,ib) * &
                                           tr_dir_z(L,ft,iv) * (1.0_r8 -bc_in(s)%albgr_dir_rb(ib))
                                      tr_soild = tr_soild + ftweight(L,ft,1)*forc_dir(ifp,ib) * tr_dir_z(L,ft,iv)
                                      tr_soili = tr_soili + ftweight(L,ft,1)*Dif_dn(L,ft,iv)
                                   end if
                                   ! Absorbed radiation, shaded and sunlit portions of leaf layers
                                   !here we get one unit of diffuse radiation... how much of
                                   !it is absorbed?
                                   if (ib == 1) then ! only set the absorbed PAR for the visible light band. 
                                      do iv = 1, currentPatch%nrad(L,ft)
                                         if (radtype==1) then
                                            if ( DEBUG ) then
                                               write(fates_log(),*) 'EDsurfAlb 730 ',Abs_dif_z(ft,iv),currentPatch%f_sun(L,ft,iv)
                                               write(fates_log(),*) 'EDsurfAlb 731 ', currentPatch%fabd_sha_z(L,ft,iv), &
                                                    currentPatch%fabd_sun_z(L,ft,iv)
                                            endif
                                            currentPatch%fabd_sha_z(L,ft,iv) = Abs_dif_z(ft,iv) * &
                                                 (1._r8 - currentPatch%f_sun(L,ft,iv))
                                            currentPatch%fabd_sun_z(L,ft,iv) = Abs_dif_z(ft,iv) * &
                                                 currentPatch%f_sun(L,ft,iv) + &
                                                 Abs_dir_z(ft,iv)
                                         else
                                            currentPatch%fabi_sha_z(L,ft,iv) = Abs_dif_z(ft,iv) * &
                                                 (1._r8 - currentPatch%f_sun(L,ft,iv))
                                            currentPatch%fabi_sun_z(L,ft,iv) = Abs_dif_z(ft,iv) * &
                                                 currentPatch%f_sun(L,ft,iv)
                                         endif
                                         if ( DEBUG ) then
                                            write(fates_log(),*) 'EDsurfAlb 740 ', currentPatch%fabd_sha_z(L,ft,iv), &
                                                 currentPatch%fabd_sun_z(L,ft,iv)
                                         endif
                                      end do
                                   endif ! ib 
                                   
                                   !==============================================================================!
                                   ! Sum fluxes
                                   !==============================================================================!
                                   ! Solar radiation absorbed by ground
                                   iv = currentPatch%nrad(L,ft) + 1
                                   if (L==currentPatch%NCL_p)then
                                      abs_rad(ib) = abs_rad(ib) +  (Abs_dir_z(ft,iv) + Abs_dif_z(ft,iv))
                                   end if
                                   ! Solar radiation absorbed by vegetation and sunlit/shaded leaves
                                   do iv = 1,currentPatch%nrad(L,ft)
                                      if (radtype == 1)then
                                         currentPatch%fabd(ib) = currentPatch%fabd(ib) + Abs_dir_z(ft,iv)+Abs_dif_z(ft,iv)
                                         ! bc_out(s)%fabd_parb(ifp,ib) = currentPatch%fabd(ib)
                                      else
                                         currentPatch%fabi(ib) = currentPatch%fabi(ib) + Abs_dif_z(ft,iv)
                                         ! bc_out(s)%fabi_parb(ifp,ib) = currentPatch%fabi(ib)
                                      endif
                                   end do
                                   ! Albefor
                                   if (L==1)then !top canopy layer.
                                      if (radtype == 1)then
                                         bc_out(s)%albd_parb(ifp,ib) = bc_out(s)%albd_parb(ifp,ib) + &
                                              Dif_up(L,ft,1) * ftweight(L,ft,1)
                                      else
                                         bc_out(s)%albi_parb(ifp,ib) = bc_out(s)%albi_parb(ifp,ib) + &
                                              Dif_up(L,ft,1) * ftweight(L,ft,1)
                                      end if
                                   end if
                                end if ! present
                             end do !ft
                             if (radtype == 1)then
                                bc_out(s)%fabd_parb(ifp,ib) = currentPatch%fabd(ib)
                             else
                                bc_out(s)%fabi_parb(ifp,ib) = currentPatch%fabi(ib)
                             endif
                             
                             
                             !radiation absorbed from fluxes through unfilled part of lower canopy.
                             if (currentPatch%NCL_p > 1.and.L == currentPatch%NCL_p)then 
                                abs_rad(ib) = abs_rad(ib) + weighted_dif_down(L-1) * &
                                     (1.0_r8-sum(ftweight(L,1:numpft,1)))*(1.0_r8-bc_in(s)%albgr_dif_rb(ib))
                                abs_rad(ib) = abs_rad(ib) + forc_dir(ifp,ib) * weighted_dir_tr(L-1) * &
                                     (1.0_r8-sum(ftweight(L,1:numpft,1)))*(1.0_r8-bc_in(s)%albgr_dir_rb(ib))
                                tr_soili = tr_soili + weighted_dif_down(L-1) * (1.0_r8-sum(ftweight(L,1:numpft,1)))
                                tr_soild = tr_soild + forc_dir(ifp,ib) * weighted_dir_tr(L-1) * (1.0_r8-sum(ftweight(L,1:numpft,1)))
                             endif
                             
                             if (radtype == 1)then
                                currentPatch%tr_soil_dir(ib) = tr_soild
                                currentPatch%tr_soil_dir_dif(ib) = tr_soili
                                currentPatch%sabs_dir(ib)     = abs_rad(ib)
                                bc_out(s)%ftdd_parb(ifp,ib)  = tr_soild
                                bc_out(s)%ftid_parb(ifp,ib) =  tr_soili
                             else
                                currentPatch%tr_soil_dif(ib) = tr_soili
                                currentPatch%sabs_dif(ib)     = abs_rad(ib)
                                bc_out(s)%ftii_parb(ifp,ib) =  tr_soili
                             end if
                             
                          end do!l
                          
                          
                          !==============================================================================!
                          ! Conservation check
                          !==============================================================================!
                          ! Total radiation balance: absorbed = incoming - outgoing
                          
                          if (radtype == 1)then
                             error = abs(currentPatch%sabs_dir(ib) - (currentPatch%tr_soil_dir(ib) * &
                                  (1.0_r8-bc_in(s)%albgr_dir_rb(ib)) + &
                                  currentPatch%tr_soil_dir_dif(ib) * (1.0_r8-bc_in(s)%albgr_dif_rb(ib))))
                             if ( abs(error) > 0.0001)then
                                write(fates_log(),*)'dir ground absorption error',ifp,s,error,currentPatch%sabs_dir(ib), &
                                     currentPatch%tr_soil_dir(ib)* &
                                     (1.0_r8-bc_in(s)%albgr_dir_rb(ib)),currentPatch%NCL_p,ib,sum(ftweight(1,1:numpft,1))
                                write(fates_log(),*) 'albedos',currentPatch%sabs_dir(ib) ,currentPatch%tr_soil_dir(ib), &
                                     (1.0_r8-bc_in(s)%albgr_dir_rb(ib)),currentPatch%lai
                                
                                do ft =1,3
                                   iv = currentPatch%nrad(1,ft) + 1
                                   write(fates_log(),*) 'abs soil fluxes', Abs_dir_z(ft,iv),Abs_dif_z(ft,iv)
                                end do
                                
                             end if
                          else
                             if ( abs(currentPatch%sabs_dif(ib)-(currentPatch%tr_soil_dif(ib) * &
                                  (1.0_r8-bc_in(s)%albgr_dif_rb(ib)))) > 0.0001)then
                                write(fates_log(),*)'dif ground absorption error',ifp,s,currentPatch%sabs_dif(ib) , &
                                     (currentPatch%tr_soil_dif(ib)* &
                                     (1.0_r8-bc_in(s)%albgr_dif_rb(ib))),currentPatch%NCL_p,ib,sum(ftweight(1,1:numpft,1))
                             endif
                          endif
                          
                          if (radtype == 1)then
                             error = (forc_dir(ifp,ib) + forc_dif(ifp,ib)) - &
                                  (bc_out(s)%fabd_parb(ifp,ib)  + bc_out(s)%albd_parb(ifp,ib) + currentPatch%sabs_dir(ib))
                          else
                             error = (forc_dir(ifp,ib) + forc_dif(ifp,ib)) - &
                                  (bc_out(s)%fabi_parb(ifp,ib)  + bc_out(s)%albi_parb(ifp,ib) + currentPatch%sabs_dif(ib))
                          endif
                          lai_reduction(:) = 0.0_r8
                          do L = 1, currentPatch%NCL_p
                             do ft =1,numpft
                                if (currentPatch%present(L,ft) == 1)then
                                   do iv = 1, currentPatch%nrad(L,ft)
                                      if (lai_change(L,ft,iv) > 0.0_r8)then
                                         lai_reduction(L) = max(lai_reduction(L),lai_change(L,ft,iv))
                                      endif
                                   enddo
                                endif
                             enddo
                          enddo
                          
                          if (radtype == 1)then
                             !here we are adding a within-ED radiation scheme tolerance, and then adding the diffrence onto the albedo
                             !it is important that the lower boundary for this is ~1000 times smaller than the tolerance in surface albedo. 
                             if (abs(error)  >  1.e-9_r8 .and. abs(error) < 0.15_r8)then
                                bc_out(s)%albd_parb(ifp,ib) = bc_out(s)%albd_parb(ifp,ib) + error
                                !this terms adds the error back on to the albedo. While this is partly inexcusable, it is 
                                ! in the medium term a solution that
                                ! prevents the model from crashing with small and occasional energy balances issues.
                                ! These are extremely difficult to debug, many have been solved already, leading
                                ! to the complexity of this code, but where the system generates occasional errors, we
                                ! will deal with them for now.
                             end if
                             if (abs(error)  >  0.15_r8)then
                                write(fates_log(),*) 'Large Dir Radn consvn error',error ,ifp,ib
                                write(fates_log(),*) 'diags', bc_out(s)%albd_parb(ifp,ib), bc_out(s)%ftdd_parb(ifp,ib), &
                                     bc_out(s)%ftid_parb(ifp,ib), bc_out(s)%fabd_parb(ifp,ib)
                                write(fates_log(),*) 'lai_change',lai_change(currentpatch%ncl_p,1:numpft,1:diag_nlevleaf)
                                write(fates_log(),*) 'elai',currentpatch%elai_profile(currentpatch%ncl_p,1:numpft,1:diag_nlevleaf)
                                write(fates_log(),*) 'esai',currentpatch%esai_profile(currentpatch%ncl_p,1:numpft,1:diag_nlevleaf)
                                write(fates_log(),*) 'ftweight',ftweight(1,1:numpft,1:diag_nlevleaf)
                                write(fates_log(),*) 'cp',currentPatch%area, currentPatch%patchno
                                write(fates_log(),*) 'bc_in(s)%albgr_dir_rb(ib)',bc_in(s)%albgr_dir_rb(ib)
                                
                                bc_out(s)%albd_parb(ifp,ib) = bc_out(s)%albd_parb(ifp,ib) + error
                             end if
                          else
                             
                             if (abs(error)  >  1.e-9_r8 .and. abs(error) < 0.15_r8)then
                                bc_out(s)%albi_parb(ifp,ib) = bc_out(s)%albi_parb(ifp,ib) + error
                             end if
                             
                             if (abs(error)  >  0.15_r8)then
                                write(fates_log(),*)  '>5% Dif Radn consvn error',error ,ifp,ib
                                write(fates_log(),*) 'diags', bc_out(s)%albi_parb(ifp,ib), bc_out(s)%ftii_parb(ifp,ib), &
                                     bc_out(s)%fabi_parb(ifp,ib)
                                write(fates_log(),*) 'lai_change',lai_change(currentpatch%ncl_p,1:numpft,1:diag_nlevleaf)
                                write(fates_log(),*) 'elai',currentpatch%elai_profile(currentpatch%ncl_p,1:numpft,1:diag_nlevleaf)
                                write(fates_log(),*) 'esai',currentpatch%esai_profile(currentpatch%ncl_p,1:numpft,1:diag_nlevleaf)
                                write(fates_log(),*) 'ftweight',ftweight(currentpatch%ncl_p,1:numpft,1:diag_nlevleaf)
                                write(fates_log(),*) 'cp',currentPatch%area, currentPatch%patchno
                                write(fates_log(),*) 'bc_in(s)%albgr_dif_rb(ib)',bc_in(s)%albgr_dif_rb(ib)
                                write(fates_log(),*) 'rhol',rhol(1:numpft,:)
                                write(fates_log(),*) 'ftw',sum(ftweight(1,1:numpft,1)),ftweight(1,1:numpft,1)
                                write(fates_log(),*) 'present',currentPatch%present(1,1:numpft)
                                write(fates_log(),*) 'CAP',currentPatch%canopy_area_profile(1,1:numpft,1)
                                
                                bc_out(s)%albi_parb(ifp,ib) = bc_out(s)%albi_parb(ifp,ib) + error
                             end if
                             
                             if (radtype == 1)then
                                error = (forc_dir(ifp,ib) + forc_dif(ifp,ib)) - &
                                     (bc_out(s)%fabd_parb(ifp,ib)  + bc_out(s)%albd_parb(ifp,ib) + currentPatch%sabs_dir(ib))
                             else
                                error = (forc_dir(ifp,ib) + forc_dif(ifp,ib)) - &
                                     (bc_out(s)%fabi_parb(ifp,ib)  + bc_out(s)%albi_parb(ifp,ib) + currentPatch%sabs_dif(ib))
                             endif
                             
                             if (abs(error)  >  0.00000001_r8)then
                                write(fates_log(),*)  'there is still error after correction',error ,ifp,ib
                             end if
                             
                          end if
                          
                       end do !hlm_numSWb
                       
                    enddo ! rad-type
                 endif ! is there vegetation? 

              end if    ! if the vegetation and zenith filter is active
              currentPatch => currentPatch%younger
           end do       ! Loop linked-list patches
        enddo           ! Loop Sites
        
      end associate
      return
    end subroutine ED_Norman_Radiation
    
 ! ======================================================================================

 subroutine ED_SunShadeFracs(nsites, sites,bc_in,bc_out)
    
    implicit none

    ! Arguments
    integer,intent(in)                      :: nsites
    type(ed_site_type),intent(inout),target :: sites(nsites)
    type(bc_in_type),intent(in)             :: bc_in(nsites)
    type(bc_out_type),intent(inout)         :: bc_out(nsites)
    

    ! locals
    type (ed_patch_type),pointer :: cpatch   ! c"urrent" patch
    real(r8)          :: sunlai
    real(r8)          :: shalai
    real(r8)          :: elai
    integer           :: CL
    integer           :: FT
    integer           :: iv
    integer           :: s
    integer           :: ifp
    

    do s = 1,nsites

       ifp = 0
       cpatch => sites(s)%oldest_patch

       do while (associated(cpatch))                 
          
          ifp=ifp+1
          
          if( DEBUG ) write(fates_log(),*) 'edsurfRad_5600',ifp,s,cpatch%NCL_p,numpft
          
          ! zero out various datas
          cpatch%ed_parsun_z(:,:,:) = 0._r8
          cpatch%ed_parsha_z(:,:,:) = 0._r8
          cpatch%ed_laisun_z(:,:,:) = 0._r8     
          cpatch%ed_laisha_z(:,:,:) = 0._r8

          bc_out(s)%fsun_pa(ifp) = 0._r8

          sunlai  = 0._r8
          shalai  = 0._r8

          ! Loop over patches to calculate laisun_z and laisha_z for each layer.
          ! Derive canopy laisun, laisha, and fsun from layer sums.
          ! If sun/shade big leaf code, nrad=1 and fsun_z(p,1) and tlai_z(p,1) from
          ! SurfaceAlbedo is canopy integrated so that layer value equals canopy value.
          
          ! cpatch%f_sun is calculated in the surface_albedo routine...
          
          do CL = 1, cpatch%NCL_p
             do FT = 1,numpft

                if( DEBUG ) write(fates_log(),*) 'edsurfRad_5601',CL,FT,cpatch%nrad(CL,ft)
                
                do iv = 1, cpatch%nrad(CL,ft) !NORMAL CASE. 
                   
                   ! FIX(SPM,040114) - existing comment
                   ! ** Should this be elai or tlai? Surely we only do radiation for elai? 
                   
                   cpatch%ed_laisun_z(CL,ft,iv) = cpatch%elai_profile(CL,ft,iv) * &
                         cpatch%f_sun(CL,ft,iv)
                   
                   if ( DEBUG ) write(fates_log(),*) 'edsurfRad 570 ',cpatch%elai_profile(CL,ft,iv)
                   if ( DEBUG ) write(fates_log(),*) 'edsurfRad 571 ',cpatch%f_sun(CL,ft,iv)
                   
                   cpatch%ed_laisha_z(CL,ft,iv) = cpatch%elai_profile(CL,ft,iv) * &
                         (1._r8 - cpatch%f_sun(CL,ft,iv))
                   
                end do
                
                !needed for the VOC emissions, etc. 
                sunlai = sunlai + sum(cpatch%ed_laisun_z(CL,ft,1:cpatch%nrad(CL,ft)))
                shalai = shalai + sum(cpatch%ed_laisha_z(CL,ft,1:cpatch%nrad(CL,ft)))
                
             end do
          end do
          
          if(sunlai+shalai > 0._r8)then
             bc_out(s)%fsun_pa(ifp) = sunlai / (sunlai+shalai) 
          else
             bc_out(s)%fsun_pa(ifp) = 0._r8
          endif
          
          if(bc_out(s)%fsun_pa(ifp) > 1._r8)then
             write(fates_log(),*) 'too much leaf area in profile',  bc_out(s)%fsun_pa(ifp), &
                   cpatch%lai,sunlai,shalai
          endif

          elai = calc_areaindex(cpatch,'elai')

          bc_out(s)%laisun_pa(ifp) = elai*bc_out(s)%fsun_pa(ifp)
          bc_out(s)%laisha_pa(ifp) = elai*(1.0_r8-bc_out(s)%fsun_pa(ifp))

         ! Absorbed PAR profile through canopy
         ! If sun/shade big leaf code, nrad=1 and fluxes from SurfaceAlbedo
         ! are canopy integrated so that layer values equal big leaf values.
         
         if ( DEBUG ) write(fates_log(),*) 'edsurfRad 645 ',cpatch%NCL_p,numpft
         
         do CL = 1, cpatch%NCL_p
            do FT = 1,numpft
               
               if ( DEBUG ) write(fates_log(),*) 'edsurfRad 649 ',cpatch%nrad(CL,ft)
               
               do iv = 1, cpatch%nrad(CL,ft)
                  
                  if ( DEBUG ) then
                     write(fates_log(),*) 'edsurfRad 653 ', cpatch%ed_parsun_z(CL,ft,iv)
                     write(fates_log(),*) 'edsurfRad 654 ', bc_in(s)%solad_parb(ifp,ipar)
                     write(fates_log(),*) 'edsurfRad 655 ', bc_in(s)%solai_parb(ifp,ipar)
                     write(fates_log(),*) 'edsurfRad 656 ', cpatch%fabd_sun_z(CL,ft,iv)
                     write(fates_log(),*) 'edsurfRad 657 ', cpatch%fabi_sun_z(CL,ft,iv)
                  endif
                  
                  cpatch%ed_parsun_z(CL,ft,iv) = &
                        bc_in(s)%solad_parb(ifp,ipar)*cpatch%fabd_sun_z(CL,ft,iv) + &
                        bc_in(s)%solai_parb(ifp,ipar)*cpatch%fabi_sun_z(CL,ft,iv) 
                  
                  if ( DEBUG )write(fates_log(),*) 'edsurfRad 663 ', cpatch%ed_parsun_z(CL,ft,iv)
                  
                  cpatch%ed_parsha_z(CL,ft,iv) = &
                        bc_in(s)%solad_parb(ifp,ipar)*cpatch%fabd_sha_z(CL,ft,iv) + &
                        bc_in(s)%solai_parb(ifp,ipar)*cpatch%fabi_sha_z(CL,ft,iv)          
                  
                  if ( DEBUG ) write(fates_log(),*) 'edsurfRad 669 ', cpatch%ed_parsha_z(CL,ft,iv)
                  
               end do !iv
            end do !FT
         end do !CL
         
         cpatch => cpatch%younger
      enddo
      
      
   enddo
   return
   
end subroutine ED_SunShadeFracs


!      ! MOVE TO THE INTERFACE
!      subroutine ED_CheckSolarBalance(g,filter_nourbanp,num_nourbanp,fsa,fsr,forc_solad,forc_solai)


!         implicit none
!         integer,intent(in),dimension(:)    :: gridcell     ! =>    gridcell index
!         integer,intent(in),dimension(:)    :: filter_nourbanp ! => patch filter for non-urban points
!         integer, intent(in)                :: num_nourbanp !       number of patches in non-urban points in patch  filter
!         real(r8),intent(in),dimension(:,:) :: forc_solad   ! =>    atm2lnd_inst%forc_solad_grc, direct radiation (W/m**2
!         real(r8),intent(in),dimension(:,:) :: forc_solai   ! =>    atm2lnd_inst%forc_solai_grc, diffuse radiation (W/m**2)
!         real(r8),intent(in),dimension(:,:) :: fsa          ! =>    solarabs_inst%fsa_patch, solar radiation absorbed (total) (W/m**2)
!         real(r8),intent(in),dimension(:,:) :: fsr          ! =>    solarabs_inst%fsr_patch, solar radiation reflected (W/m**2)      

!         integer :: p
!         integer :: fp
!         integer :: g
!         real(r8) :: errsol

!         do fp = 1,num_nourbanp
!            p = filter_nourbanp(fp)
!            g = gridcell(p)
!            errsol = (fsa(p) + fsr(p)  - (forc_solad(g,1) + forc_solad(g,2) + forc_solai(g,1) + forc_solai(g,2)))
!            if(abs(errsol) > 0.1_r8)then
!               write(fates_log(),*) 'sol error in surf rad',p,g, errsol
!            endif
!         end do
!         return
!      end subroutine ED_CheckSolarBalance
   

end module EDSurfaceRadiationMod
