module EDSurfaceAlbedoMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Performs surface albedo calculations
  !
  ! !PUBLIC TYPES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_assert_mod , only : shr_assert
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use clm_varpar     , only : numrad, nclmax
  use decompMod      , only : bounds_type

  implicit none
  save
  private

  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: ED_Norman_Radiation  ! Surface albedo and two-stream fluxes
  !
  ! !PUBLIC DATA MEMBERS:
  ! The CLM default albice values are too high.
  ! Full-spectral albedo for land ice is ~0.5 (Paterson, Physics of Glaciers, 1994, p. 59)
  ! This is the value used in CAM3 by Pritchard et al., GRL, 35, 2008.

  real(r8), public  :: albice(numrad) = &       ! albedo land ice by waveband (1=vis, 2=nir)
       (/ 0.80_r8, 0.55_r8 /)
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine ED_Norman_Radiation (bounds, filter_vegsol, num_vegsol, coszen, &
       surfalb_vars)
    !
    ! !DESCRIPTION:
    ! Two-stream fluxes for canopy radiative transfer
    ! Use two-stream approximation of Dickinson (1983) Adv Geophysics
    ! 25:305-353 and Sellers (1985) Int J Remote Sensing 6:1335-1372
    ! to calculate fluxes absorbed by vegetation, reflected by vegetation,
    ! and transmitted through vegetation for unit incoming direct or diffuse
    ! flux given an underlying surface with known albedo.
    ! Calculate sunlit and shaded fluxes as described by
    ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 and extended to
    ! a multi-layer canopy to calculate APAR profile
    !
    ! !USES:
    use clm_varctl        , only : iulog
    use EcophysconType    , only : ecophyscon
    use EDtypesMod        , only : patch, numpft_ed, nlevcan_ed, gridCellEdState
    use PatchType         , only : pft  
    use EDVecPatchType    , only : EDpft
    use SurfaceAlbedoType , only : surfalb_type
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type)  , intent(in)    :: bounds                 ! bounds
    integer            , intent(in)    :: filter_vegsol(:)       ! filter for vegetated pfts with coszen>0
    integer            , intent(in)    :: num_vegsol             ! number of vegetated pfts where coszen>0
    real(r8)           , intent(in)    :: coszen( bounds%begp: ) ! cosine solar zenith angle for next time step [pft]
    type(surfalb_type) , intent(inout) :: surfalb_vars 
    !
    ! !LOCAL VARIABLES:
    ! ============================================================================
    ! ED/NORMAN RADIATION DECS
    ! ============================================================================
    type (patch) , pointer :: currentPatch

    integer  :: radtype, L, ft, g ,j
    integer  :: iter                                          ! Iteration index
    integer  :: irep                                          ! Flag to exit iteration loop
    real(r8) :: sb
    real(r8) :: error                                         ! Error check
    real(r8) :: down_rad, up_rad                              ! Iterative solution do Dif_dn and Dif_up
    real(r8) :: ftweight(nclmax,numpft_ed,nlevcan_ed)
    real(r8) :: k_dir(numpft_ed)                              ! Direct beam extinction coefficient
    real(r8) :: tr_dir_z(nclmax,numpft_ed,nlevcan_ed)         ! Exponential transmittance of direct beam radiation through a single layer
    real(r8) :: tr_dif_z(nclmax,numpft_ed,nlevcan_ed)         ! Exponential transmittance of diffuse radiation through a single layer
    real(r8) :: forc_dir(bounds%begp:bounds%endp,numrad)
    real(r8) :: forc_dif(bounds%begp:bounds%endp,numrad)
    real(r8) :: weighted_dir_tr(nclmax)
    real(r8) :: weighted_fsun(nclmax)
    real(r8) :: weighted_dif_ratio(nclmax,numrad)
    real(r8) :: weighted_dif_down(nclmax)
    real(r8) :: weighted_dif_up(nclmax)
    real(r8) :: refl_dif(nclmax,numpft_ed,nlevcan_ed,numrad)  ! Term for diffuse radiation reflected by laye
    real(r8) :: tran_dif(nclmax,numpft_ed,nlevcan_ed,numrad)  ! Term for diffuse radiation transmitted by layer
    real(r8) :: dif_ratio(nclmax,numpft_ed,nlevcan_ed,numrad) ! Ratio of upward to forward diffuse fluxes
    real(r8) :: Dif_dn(nclmax,numpft_ed,nlevcan_ed)           ! Forward diffuse flux onto canopy layer J (W/m**2 ground area)
    real(r8) :: Dif_up(nclmax,numpft_ed,nlevcan_ed)           ! Upward diffuse flux above canopy layer J (W/m**2 ground area)
    real(r8) :: lai_change(nclmax,numpft_ed,nlevcan_ed)       ! Forward diffuse flux onto canopy layer J (W/m**2 ground area)

    real(r8) :: f_not_abs(numpft_ed,numrad)                   ! Fraction reflected + transmitted. 1-absorbtion.
    real(r8) :: tolerance
    real(r8) :: Abs_dir_z(numpft_ed,nlevcan_ed)
    real(r8) :: Abs_dif_z(numpft_ed,nlevcan_ed)
    real(r8) :: abs_rad(numrad)                               !radiation absorbed by soil
    real(r8) :: tr_soili                                      ! Radiation transmitted to the soil surface.
    real(r8) :: tr_soild                                      ! Radiation transmitted to the soil surface.
    real(r8) :: phi1b(bounds%begp:bounds%endp,numpft_ed)      ! Radiation transmitted to the soil surface.
    real(r8) :: phi2b(bounds%begp:bounds%endp,numpft_ed)
    real(r8) :: laisum                                        ! cumulative lai+sai for canopy layer (at middle of layer)

    real(r8) :: angle
    real(r8), parameter :: pi   = 3.141592654                 ! PI
    real(r8) :: denom
    real(r8) :: lai_reduction(2)

    integer  :: fp,p,c,iv                                     ! array indices
    integer  :: ib                                            ! waveband number
    real(r8) :: cosz                                          ! 0.001 <= coszen <= 1.000
    real(r8) :: chil(bounds%begp:bounds%endp)                 ! -0.4 <= xl <= 0.6
    real(r8) :: gdir(bounds%begp:bounds%endp)                 ! leaf projection in solar direction (0 to 1)
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    ! What is this about? (FIX(RF,032414))
    SHR_ASSERT_ALL((ubound(coszen) == (/bounds%endp/)),         errMsg(__FILE__, __LINE__))

    associate(                                                &
         rhol         =>    ecophyscon%rhol                 , & ! Input:  [real(r8) (:)   ] leaf reflectance: 1=vis, 2=nir
         rhos         =>    ecophyscon%rhos                 , & ! Input:  [real(r8) (:)   ] stem reflectance: 1=vis, 2=nir
         taul         =>    ecophyscon%taul                 , & ! Input:  [real(r8) (:)   ] leaf transmittance: 1=vis, 2=nir
         taus         =>    ecophyscon%taus                 , & ! Input:  [real(r8) (:)   ] stem transmittance: 1=vis, 2=nir
         xl           =>    ecophyscon%xl                   , & ! Input:  [real(r8) (:)   ] ecophys const - leaf/stem orientation index

         albgrd       =>    surfalb_vars%albgrd_col         , & ! Input:  [real(r8) (:,:) ] ground albedo (direct) (column-level)
         albgri       =>    surfalb_vars%albgri_col         , & ! Input:  [real(r8) (:,:) ] ground albedo (diffuse)(column-level)
         albd         =>    surfalb_vars%albd_patch         , & ! Output: [real(r8) (:,:) ] surface albedo (direct)
         albi         =>    surfalb_vars%albi_patch         , & ! Output: [real(r8) (:,:) ] surface albedo (diffuse)
         fabd         =>    surfalb_vars%fabd_patch         , & ! Output: [real(r8) (:,:) ] flux absorbed by canopy per unit direct flux
         fabd_sun     =>    surfalb_vars%fabd_sun_patch     , & ! Output: [real(r8) (:,:) ] flux absorbed by sunlit canopy per unit direct flux
         fabd_sha     =>    surfalb_vars%fabd_sha_patch     , & ! Output: [real(r8) (:,:) ] flux absorbed by shaded canopy per unit direct flux
         fabi         =>    surfalb_vars%fabi_patch         , & ! Output: [real(r8) (:,:) ] flux absorbed by canopy per unit diffuse flux
         fabi_sun     =>    surfalb_vars%fabi_sun_patch     , & ! Output: [real(r8) (:,:) ] flux absorbed by sunlit canopy per unit diffuse flux
         fabi_sha     =>    surfalb_vars%fabi_sha_patch     , & ! Output: [real(r8) (:,:) ] flux absorbed by shaded canopy per unit diffuse flux
         ftdd         =>    surfalb_vars%ftdd_patch         , & ! Output: [real(r8) (:,:) ] down direct flux below canopy per unit direct flx
         ftid         =>    surfalb_vars%ftid_patch         , & ! Output: [real(r8) (:,:) ] down diffuse flux below canopy per unit direct flx
         ftii         =>    surfalb_vars%ftii_patch         , & ! Output: [real(r8) (:,:) ] down diffuse flux below canopy per unit diffuse flx
         nrad         =>    surfalb_vars%nrad_patch         , & ! Input:  [integer  (:)   ] number of canopy layers, above snow for radiative transfer
         fabd_sun_z   =>    surfalb_vars%fabd_sun_z_patch   , & ! Output: [real(r8) (:,:) ] absorbed sunlit leaf direct  PAR (per unit lai+sai) for each canopy layer
         fabd_sha_z   =>    surfalb_vars%fabd_sha_z_patch   , & ! Output: [real(r8) (:,:) ] absorbed shaded leaf direct  PAR (per unit lai+sai) for each canopy layer
         fabi_sun_z   =>    surfalb_vars%fabi_sun_z_patch   , & ! Output: [real(r8) (:,:) ] absorbed sunlit leaf diffuse PAR (per unit lai+sai) for each canopy layer
         fabi_sha_z   =>    surfalb_vars%fabi_sha_z_patch   , & ! Output: [real(r8) (:,:) ] absorbed shaded leaf diffuse PAR (per unit lai+sai) for each canopy layer
         fsun_z       =>    surfalb_vars%fsun_z_patch       , & ! Output: [real(r8) (:,:) ] sunlit fraction of canopy layer
         vcmaxcintsun =>    surfalb_vars%vcmaxcintsun_patch , & ! Output: [real(r8) (:)   ] leaf to canopy scaling coefficient, sunlit leaf vcmax
         vcmaxcintsha =>    surfalb_vars%vcmaxcintsha_patch , & ! Output: [real(r8) (:)   ] leaf to canopy scaling coefficient, shaded leaf vcmax
         
         ED_patch     =>    EDpft%ED_patch                    &
         )

      !================================================================
      ! NORMAN RADIATION CODE
      ! ============================================================================
      ! FIX(SPM,032414) refactor this...too long for one routine.

      tolerance = 0.000000001_r8 ! FIX(SPM,032414) make this a param

      do fp = 1,num_vegsol
         p = filter_vegsol(fp)
         c = pft%column(p)
         if(ED_patch(p) == 1)then ! We have vegetation...
            g = pft%gridcell(p)
            currentPatch => gridCellEdState(g)%spnt%oldest_patch
            do while(p /= currentPatch%clm_pno)
               currentPatch => currentPatch%younger
            enddo

            if(associated(currentPatch))then
               !zero all of the matrices used here to reduce potential for errors.
               weighted_dir_tr(:) = 0._r8
               weighted_dif_down(:) = 0._r8
               weighted_dif_up(:) = 0._r8
               albd(p,:)  = 0._r8
               albi(p,:)  = 0._r8
               currentPatch%f_sun(:,:,:) = 0._r8
               tr_dir_z(:,:,:) = 0._r8
               tr_dif_z(:,:,:) = 0._r8
               ftweight(:,:,:) = 0._r8
               lai_change(:,:,:) = 0._r8
               Dif_up(:,:,:) = 0._r8
               Dif_dn(:,:,:) = 0._r8
               refl_dif(:,:,:,:) = 0.0_r8
               tran_dif(:,:,:,:) = 0.0_r8
               dif_ratio(:,:,:,:) = 0.0_r8
               currentPatch%fabd_sun_z(:,:,:) = 0._r8
               currentPatch%fabd_sha_z(:,:,:) = 0._r8
               currentPatch%fabi_sun_z(:,:,:) = 0._r8
               currentPatch%fabi_sha_z(:,:,:) = 0._r8
               currentPatch%fabd(:) = 0._r8
               currentPatch%fabi(:) = 0._r8

               ! Is this pft/canopy layer combination present in this patch?
               do L = 1,nclmax
                  do ft = 1,numpft_ed
                     currentPatch%present(L,ft) = 0
                     do  iv = 1, currentPatch%nrad(L,ft)
                        if(currentPatch%canopy_area_profile(L,ft,iv) > 0._r8)then
                           currentPatch%present(L,ft) = 1
                           !I think 'present' is only used here...
                        endif
                     end do !iv
                  end do !ft
               end do !L
               g = currentPatch%siteptr%clmgcell

               do radtype = 1,2 !do this once for one unit of diffuse, and once for one unit of direct radiation
                  do ib = 1,numrad
                     if(radtype == 1) then
                        ! Set the hypothetical driving radiation. We do this once for a single unit of direct and
                        ! once for a single unit of diffuse radiation.
                        forc_dir(p,ib) = 1.00_r8
                        forc_dif(p,ib) = 0.00_r8
                     else !dif
                        forc_dir(p,ib) = 0.00_r8
                        forc_dif(p,ib) = 1.00_r8
                     end if
                  end do !ib

                  !Extract information that needs to be provided by ED into local array.
                  ftweight(:,:,:) = 0._r8
                  do L = 1,currentPatch%NCL_p
                     do ft = 1,numpft_ed
                        do  iv = 1, currentPatch%nrad(L,ft)
                           !this is already corrected for area in CLAP
                           ftweight(L,ft,iv) = currentPatch%canopy_area_profile(L,ft,iv) 
                        end do  !iv
                     end do  !ft
                  end do  !L

                  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
                  ! Direct beam extinction coefficient, k_dir. PFT specific.
                  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
                  cosz = max(0.001_r8, coszen(p)) !copied from previous radiation code...
                  do ft = 1,numpft_ed
                     sb = (90._r8 - (acos(cosz)*180/pi)) * (pi / 180._r8)
                     chil(p) = xl(ft) !min(max(xl(ft), -0.4_r8), 0.6_r8 )
                     if (abs(chil(p)) <= 0.01_r8) then
                        chil = 0.01_r8
                     end if
                     phi1b(p,ft) = 0.5_r8 - 0.633_r8*chil(p) - 0.330_r8*chil(p)*chil(p)
                     phi2b(p,ft) = 0.877_r8 * (1._r8 - 2._r8*phi1b(p,ft)) !0 = horiz leaves, 1 - vert leaves.
                     gdir(p) = phi1b(p,ft) + phi2b(p,ft) * sin(sb)
                     !how much direct light penetrates a singleunit of lai?
                     k_dir(ft) = gdir(p) / sin(sb)
                  end do !FT

                  do L = 1,currentPatch%NCL_p !start at the top canopy layer (1 is the top layer.)
                     weighted_dir_tr(L) = 0.0_r8
                     weighted_fsun(L) = 0._r8
                     weighted_dif_ratio(L,1:numrad) = 0._r8
                     !Each canopy layer (canopy, understorey) has multiple 'parallel' pft's
                     do ft =1,numpft_ed
                        if(currentPatch%present(L,ft) == 1)then !only do calculation if there are the appropriate leaves.
                           !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
                           ! Diffuse transmittance, tr_dif, do each layer with thickness elai_z.
                           ! Estimated do nine sky angles in increments of 10 degrees
                           ! PFT specific...
                           !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
                           tr_dif_z(L,ft,:) = 0._r8
                           do iv = 1,currentPatch%nrad(L,ft)
                              do j = 1,9
                                 angle = (5._r8 + (j - 1) * 10._r8) * 3.142 / 180._r8
                                 gdir(p) = phi1b(p,ft) + phi2b(p,ft) * sin(angle) !This line is redundant FIX(RF,032414). 
                                 tr_dif_z(L,ft,iv) = tr_dif_z(L,ft,iv) + exp(-gdir(p) / sin(angle) * &
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
                           !THIS CAN BE COLLAPSED DOWN THE SAME WAY AS F_SUN
                           if(L==1)then
                              laisum = 0.00_r8
                              tr_dir_z(L,ft,1) = 1._r8
                              !total direct beam getting to the bottom of the top canopy.
                              do iv = 1,currentPatch%nrad(L,ft)
                                 laisum = laisum + currentPatch%elai_profile(L,ft,iv)+currentPatch%esai_profile(L,ft,iv)
                                 lai_change(L,ft,iv) = 0.0_r8
                                 if (( ftweight(L,ft,iv+1)  >  0.0_r8 ) .and. ( ftweight(L,ft,iv+1)  <  ftweight(L,ft,iv) ))then
                                    !where there is a partly empty leaf layer, some fluxes go straight through.
                                    lai_change(L,ft,iv) = ftweight(L,ft,iv)-ftweight(L,ft,iv+1)
                                 endif
                                 !n.b. in theory lai_change could be calculated daily in the ED code.
                                 tr_dir_z(L,ft,iv+1) = exp(-k_dir(ft) * laisum)* &
                                      (ftweight(L,ft,iv)/ftweight(L,ft,1))
                                 if(iv > 1)then
                                    !the lai_change(iv) affects the light incident on layer iv+2 not iv+1
                                    if(lai_change(L,ft,iv-1) > 0.0_r8)then
                                       tr_dir_z(L,ft,iv+1) = tr_dir_z(L,ft,iv+1) + tr_dir_z(L,ft,iv)* &
                                            lai_change(L,ft,iv-1) / ftweight(L,ft,1)
                                       tr_dir_z(L,ft,iv+1) = tr_dir_z(L,ft,iv+1) + tr_dir_z(L,ft,iv-1)* &
                                            (ftweight(L,ft,1)-ftweight(L,ft,iv-1))/ftweight(L,ft,1)
                                    else
                                       !account fot the light that comes striaght down from unfilled layers above.
                                       tr_dir_z(L,ft,iv+1) = tr_dir_z(L,ft,iv+1) + tr_dir_z(L,ft,iv) * &
                                            ((ftweight(L,ft,1)-ftweight(L,ft,iv))/ftweight(L,ft,1))
                                    endif
                                 else
                                    !account fot the light that comes striaght down from unfilled layers above.
                                    tr_dir_z(L,ft,iv+1) = tr_dir_z(L,ft,iv+1) + tr_dir_z(L,ft,iv) * &
                                         ((ftweight(L,ft,1)-ftweight(L,ft,iv))/ftweight(L,ft,1))
                                 endif
                              end do
                              !add up all the weighted contributions from the different PFT columns.
                              weighted_dir_tr(L) = weighted_dir_tr(L) + tr_dir_z(L,ft,currentPatch%nrad(L,ft)+1)* &
                                   ftweight(L,ft,1)
                           else
                              laisum = 0._r8

                              tr_dir_z(L,ft,1)  = weighted_dir_tr(L-1)
                              !total direct beam getting to the bottom of the top canopy.
                              !+1 Changed this Nov 13. FIX(RF,032414). nrad+1+1 seems to make no sense?
                              do iv = 1,currentPatch%nrad(L,ft) 
                                 lai_change(L,ft,iv) = 0.0_r8
                                 if(ftweight(L,ft,iv+1) > 0.0_r8.and.ftweight(L,ft,iv+1) < ftweight(L,ft,iv))then
                                    !where there is a partly empty leaf layer, some fluxes go straight through.
                                    lai_change(L,ft,iv) = ftweight(L,ft,iv)-ftweight(L,ft,iv+1)
                                 endif
                                 if(ftweight(L,ft,iv+1) > ftweight(L,ft,iv))then
                                    write(iulog,*) 'lower layer has more coverage. This is wrong' , &
                                         ftweight(L,ft,iv),ftweight(L,ft,iv+1)
                                 endif
                                 laisum = laisum + currentPatch%elai_profile(L,ft,iv)+currentPatch%esai_profile(L,ft,iv)
                                 tr_dir_z(L,ft,iv+1) = weighted_dir_tr(L-1)*exp(-k_dir(ft) * laisum)* &
                                      (ftweight(L,ft,iv)/ftweight(L,ft,1))
                                 if(iv > 1)then !the lai_change(iv) affects the light incident on layer iv+2 not iv+1
                                    if(lai_change(L,ft,iv-1) > 0.0_r8)then
                                       tr_dir_z(L,ft,iv+1) = tr_dir_z(L,ft,iv+1) + tr_dir_z(L,ft,iv)* &
                                            lai_change(L,ft,iv-1)/ftweight(L,ft,1)
                                       tr_dir_z(L,ft,iv+1) = tr_dir_z(L,ft,iv+1) + tr_dir_z(L,ft,iv-1)* &
                                            ((ftweight(L,ft,1)-ftweight(L,ft,iv-1))/ftweight(L,ft,1))
                                    else
                                       !account fot the light that comes striaght down from unfilled layers above.
                                       tr_dir_z(L,ft,iv+1) = tr_dir_z(L,ft,iv+1) + tr_dir_z(L,ft,iv) * &
                                            ((ftweight(L,ft,1)-ftweight(L,ft,iv))/ftweight(L,ft,1))
                                    endif
                                 else
                                    !account fot the light that comes striaght down from unfilled layers above.
                                    tr_dir_z(L,ft,iv+1) = tr_dir_z(L,ft,iv+1) + tr_dir_z(L,ft,iv) * &
                                         ((ftweight(L,ft,1)-ftweight(L,ft,iv))/ftweight(L,ft,1))
                                 endif
                              end do

                              weighted_dir_tr(L) = weighted_dir_tr(L) + tr_dir_z(L,ft,currentPatch%nrad(L,ft)+1)*ftweight(L,ft,1)

                           end if

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


                              if(L == 1)then !top canopy layer
                                 currentPatch%f_sun(L,ft,iv) = exp(-k_dir(ft) * laisum)* &
                                      (ftweight(L,ft,iv)/ftweight(L,ft,1))
                              else
                                 currentPatch%f_sun(L,ft,iv) = weighted_fsun(L-1)* exp(-k_dir(ft) * laisum)* &
                                      (ftweight(L,ft,iv)/ftweight(L,ft,1))
                              endif

                              if(iv > 1)then  !becasue we are looking at this layer (not the next) we only ever add fluxes if iv>1
                                 if(lai_change(L,ft,iv-1) > 0.0_r8)then
                                    currentPatch%f_sun(L,ft,iv) = currentPatch%f_sun(L,ft,iv) + currentPatch%f_sun(L,ft,iv) * &
                                         lai_change(L,ft,iv-1)/ftweight(L,ft,1)
                                    currentPatch%f_sun(L,ft,iv) = currentPatch%f_sun(L,ft,iv) + currentPatch%f_sun(L,ft,iv-1) * &
                                         (ftweight(L,ft,1)-ftweight(L,ft,iv-1))/ftweight(L,ft,1)
                                 else
                                    currentPatch%f_sun(L,ft,iv) = currentPatch%f_sun(L,ft,iv) + currentPatch%f_sun(L,ft,iv-1) * &
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
                     do ft = 1,numpft_ed
                        if(currentPatch%present(L,ft) == 1)then
                           !==============================================================================!
                           ! Iterative solution do scattering
                           !==============================================================================!

                           do ib = 1,numrad !vis, nir
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
                                 tran_dif(L,ft,iv,ib) = (1._r8 - tr_dif_z(L,ft,iv)) * taul(ft,ib) + tr_dif_z(L,ft,iv)
                              end do

                              !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
                              ! Ratio of upward to forward diffuse fluxes, dif_ratio
                              !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
                              ! Soil diffuse reflectance (ratio of down to up radiation).
                              iv = currentPatch%nrad(L,ft) + 1
                              if(L  == currentPatch%NCL_p)then !nearest the soil
                                 dif_ratio(L,ft,iv,ib) = albgri(c,ib)
                              else
                                 dif_ratio(L,ft,iv,ib) = weighted_dif_ratio(L+1,ib)
                              end if
                              ! Canopy layers, working upwardfrom soil with dif_ratio(iv+1) known
                              ! FIX(RF,032414) ray tracing eqution - need to find derivation of this...
                              ! for each unit going down, there are x units going up.
                              do iv = currentPatch%nrad(L,ft),1, -1
                                 dif_ratio(L,ft,iv,ib) = dif_ratio(L,ft,iv+1,ib) * tran_dif(L,ft,iv,ib)*tran_dif(L,ft,iv,ib) / &
                                      (1._r8 - dif_ratio(L,ft,iv+1,ib) * refl_dif(L,ft,iv,ib)) + refl_dif(L,ft,iv,ib)
                                 dif_ratio(L,ft,iv,ib) = dif_ratio(L,ft,iv,ib) * ftweight(L,ft,iv)/ftweight(L,ft,1)
                                 dif_ratio(L,ft,iv,ib) = dif_ratio(L,ft,iv,ib) + dif_ratio(L,ft,iv+1,ib)* &
                                      (ftweight(L,ft,1)-ftweight(L,ft,iv))/ftweight(L,ft,1)
                              end do
                              weighted_dif_ratio(L,ib) = weighted_dif_ratio(L,ib) + dif_ratio(L,ft,1,ib) * ftweight(L,ft,1)
                              !instance where the first layer ftweight is used a proxy for the whole column. FTWA
                           end do!numrad
                        endif ! currentPatch%present
                     end do!ft
                  end do!L

                  do ib = 1,numrad
                     Dif_dn(:,:,:) = 0.00_r8
                     Dif_up(:,:,:) = 0.00_r8
                     do L = 1, currentPatch%NCL_p !work down from the top of the canopy.
                        weighted_dif_down(L) = 0._r8
                        do ft = 1, numpft_ed
                           if(currentPatch%present(L,ft) == 1)then
                              !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
                              ! First estimates do downward and upward diffuse flux
                              !
                              ! Dif_dn =  forward diffuse flux onto layer J
                              ! Dif_up =  Upward diffuse flux above layer J
                              !
                              ! Solved here without direct beam radiation and using dif_ratio = Dif_up / Dif_dn
                              !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
                              ! downward diffuse flux onto the top surface of the canopy

                              if(L == 1)then
                                 Dif_dn(L,ft,1) = forc_dif(p,ib)
                              else
                                 Dif_dn(L,ft,1) = weighted_dif_down(L-1)
                              end if
                              ! forward diffuse flux within the canopy and at soil, working forward through canopy
                              do iv = 1,currentPatch%nrad(L,ft)
                                 denom = refl_dif(L,ft,iv,ib) *  dif_ratio(L,ft,iv,ib)
                                 denom = 1._r8 - denom
                                 Dif_dn(L,ft,iv+1) = Dif_dn(L,ft,iv) * tran_dif(L,ft,iv,ib) / &
                                      denom *ftweight(L,ft,iv)/ftweight(L,ft,1)
                                 if(iv > 1)then
                                    if(lai_change(L,ft,iv-1) > 0.0_r8)then
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
                        if(L == currentPatch%NCL_p.and.currentPatch%NCL_p > 1)then !is the the (incomplete) understorey?
                           !Add on the radiation going through the canopy gaps.
                           weighted_dif_down(L) = weighted_dif_down(L) + weighted_dif_down(L-1)*(1.0-sum(ftweight(L,:,1)))
                           !instance where the first layer ftweight is used a proxy for the whole column. FTWA
                        endif
                     end do !L

                     do L = currentPatch%NCL_p,1 ,-1 !work up from the bottom.
                        weighted_dif_up(L) = 0._r8
                        do ft = 1, numpft_ed
                           if(currentPatch%present(L,ft) == 1)then
                              !Bounce diffuse radiation off soil surface.
                              iv = currentPatch%nrad(L,ft) + 1
                              if(L==currentPatch%NCL_p)then !is this the bottom layer ?
                                 Dif_up(L,ft,iv) =albgri(c,ib) * Dif_dn(L,ft,iv)
                              else
                                 Dif_up(L,ft,iv) = weighted_dif_up(L+1)
                              end if
                              ! Upward diffuse flux within the canopy and above the canopy, working upward through canopy

                              do iv = currentPatch%nrad(L,ft), 1, -1
                                 if(lai_change(L,ft,iv) > 0.0_r8)then
                                    Dif_up(L,ft,iv) = dif_ratio(L,ft,iv,ib) * Dif_dn(L,ft,iv)*ftweight(L,ft,iv)/ftweight(L,ft,1)
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
                        if(L == currentPatch%NCL_p.and.currentPatch%NCL_p > 1)then !is this the (incomplete) understorey?
                           !Add on the radiation coming up through the canopy gaps.
                           !diffuse to diffuse
                           weighted_dif_up(L) = weighted_dif_up(L) +(1.0-sum(ftweight(L,:,1))) * &
                                weighted_dif_down(L-1) * albgri(c,ib)
                           !direct to diffuse
                           weighted_dif_up(L) = weighted_dif_up(L) + forc_dir(p,ib) * &
                                weighted_dir_tr(L-1) * (1.0-sum(ftweight(L,:,1)))*albgrd(c,ib)
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
                           do ft =1,numpft_ed
                              if(currentPatch%present(L,ft) == 1)then
                                 ! forward diffuse flux within the canopy and at soil, working forward through canopy
                                 ! with Dif_up -from previous iteration-. Dif_dn(1) is the forward diffuse flux onto the canopy.
                                 ! Note: down = forward flux onto next layer
                                 if(L == 1)then !is this the top layer?
                                    Dif_dn(L,ft,1) = forc_dif(p,ib)
                                 else
                                    Dif_dn(L,ft,1) = weighted_dif_down(L-1)
                                 end if
                                 down_rad = 0._r8

                                 !IS THE ORDERING OF THIS THE PROBLEM!!!? (just before thanksgiving)
                                 do iv = 1, currentPatch%nrad(L,ft)

                                    down_rad = Dif_dn(L,ft,iv) * tran_dif(L,ft,iv,ib) + &
                                         Dif_up(L,ft,iv+1) * refl_dif(L,ft,iv,ib) + &
                                         forc_dir(p,ib) * tr_dir_z(L,ft,iv) * (1.00_r8 - &
                                         exp(-k_dir(ft) * (currentPatch%elai_profile(L,ft,iv)+ &
                                         currentPatch%esai_profile(L,ft,iv)))) * taul(ft,ib)
                                    down_rad = down_rad *(ftweight(L,ft,iv)/ftweight(L,ft,1))

                                    if(iv > 1)then
                                       if(lai_change(L,ft,iv-1) > 0.0_r8)then
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
                           if(L == currentPatch%NCL_p.and.currentPatch%NCL_p > 1)then !is this the (incomplete) understorey?
                              weighted_dif_down(L) = weighted_dif_down(L) + weighted_dif_down(L-1)*(1.0-sum(ftweight(L,:,1)))
                           end if
                        end do ! do L loop

                        do L = 1, currentPatch%NCL_p ! working from the top down.
                           weighted_dif_up(L) = 0._r8
                           do ft =1,numpft_ed
                              if(currentPatch%present(L,ft) == 1)then
                                 ! Upward diffuse flux at soil or from lower canopy (forward diffuse and unscattered direct beam)
                                 iv = currentPatch%nrad(L,ft) + 1
                                 if(L==currentPatch%NCL_p)then  !In the bottom canopy layer, reflect off the soil
                                    Dif_up(L,ft,iv) = Dif_dn(L,ft,iv) *albgri(c,ib) + &
                                         forc_dir(p,ib) * tr_dir_z(L,ft,iv) *albgrd(c,ib)
                                 else      !In the other canopy layers, reflect off the underlying vegetation.
                                    Dif_up(L,ft,iv) =  weighted_dif_up(L+1)
                                 end if

                                 ! Upward diffuse flux within and above the canopy, working upward through canopy
                                 ! with Dif_dn from previous interation.  Note: up = upward flux above current layer
                                 do iv = currentPatch%nrad(L,ft),1,-1
                                    !this is radiation up, by layer transmittance, by

                                    !reflection of the lower layer,
                                    up_rad = Dif_dn(L,ft,iv) * refl_dif(L,ft,iv,ib)
                                    up_rad = up_rad + forc_dir(p,ib) * tr_dir_z(L,ft,iv) * (1.00_r8 - exp(-k_dir(ft) * &
                                         (currentPatch%elai_profile(L,ft,iv)+currentPatch%esai_profile(L,ft,iv)))) * rhol(ft,ib)
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

                           if(L == currentPatch%NCL_p.and.currentPatch%NCL_p > 1)then  !is this the (incomplete) understorey?
                              !Add on the radiation coming up through the canopy gaps.
                              weighted_dif_up(L) = weighted_dif_up(L) +(1.0_r8-sum(ftweight(L,:,1))) * &
                                   weighted_dif_down(L-1) * albgri(c,ib)
                              weighted_dif_up(L) = weighted_dif_up(L) + forc_dir(p,ib) * &
                                   weighted_dir_tr(L-1) * (1.0_r8-sum(ftweight(L,:,1)))*albgrd(c,ib)
                           end if
                        end do!L
                     end do ! do while over iter

                     abs_rad(ib) = 0._r8
                     tr_soili = 0._r8
                     tr_soild = 0._r8
                     do L = 1, currentPatch%NCL_p !working from the top down.
                        abs_dir_z(:,:) = 0._r8
                        abs_dif_z(:,:) = 0._r8
                        do ft =1,numpft_ed
                           if(currentPatch%present(L,ft) == 1)then
                              !==============================================================================!
                              ! Compute absorbed flux densities
                              !==============================================================================!

                              ! Absorbed direct beam and diffuse do leaf layers
                              do iv = 1, currentPatch%nrad(L,ft)
                                 Abs_dir_z(ft,iv) = ftweight(L,ft,iv)* forc_dir(p,ib) * tr_dir_z(L,ft,iv) * &
                                      (1.00_r8 - exp(-k_dir(ft) * (currentPatch%elai_profile(L,ft,iv)+ &
                                      currentPatch%esai_profile(L,ft,iv)))) * (1.00_r8 - f_not_abs(ft,ib))
                                 Abs_dif_z(ft,iv) = ftweight(L,ft,iv)* ((Dif_dn(L,ft,iv) + &
                                      Dif_up(L,ft,iv+1)) * (1.00_r8 - tr_dif_z(L,ft,iv)) * &
                                      (1.00_r8 - f_not_abs(ft,ib)))
                              end do

                              ! Absorbed direct beam and diffuse do soil
                              if(L == currentPatch%NCL_p)then
                                 iv = currentPatch%nrad(L,ft) + 1
                                 Abs_dif_z(ft,iv) = ftweight(L,ft,1)*Dif_dn(L,ft,iv) * (1.0_r8 -albgri(c,ib))
                                 Abs_dir_z(ft,iv) = ftweight(L,ft,1)*forc_dir(p,ib) * &
                                      tr_dir_z(L,ft,iv) * (1.0_r8 -albgrd(c,ib))
                                 tr_soild = tr_soild + ftweight(L,ft,1)*forc_dir(p,ib) * tr_dir_z(L,ft,iv)
                                 tr_soili = tr_soili + ftweight(L,ft,1)*Dif_dn(L,ft,iv)
                              end if
                              ! Absorbed radiation, shaded and sunlit portions of leaf layers
                              !here we get one unit of diffuse radiation... how much of
                              !it is absorbed?
                              do iv = 1, currentPatch%nrad(L,ft)
                                 if(radtype==1)then
                                    currentPatch%fabd_sha_z(L,ft,iv) = Abs_dif_z(ft,iv) * (1._r8 - currentPatch%f_sun(L,ft,iv))
                                    currentPatch%fabd_sun_z(L,ft,iv) = Abs_dif_z(ft,iv) * currentPatch%f_sun(L,ft,iv) + &
                                         Abs_dir_z(ft,iv)
                                 else
                                    currentPatch%fabi_sha_z(L,ft,iv) = Abs_dif_z(ft,iv) * (1._r8 - currentPatch%f_sun(L,ft,iv))
                                    currentPatch%fabi_sun_z(L,ft,iv) = Abs_dif_z(ft,iv) * currentPatch%f_sun(L,ft,iv)
                                 end if
                              end do

                              !==============================================================================!
                              ! Sum fluxes
                              !==============================================================================!
                              ! Solar radiation absorbed by ground
                              iv = currentPatch%nrad(L,ft) + 1
                              if(L==currentPatch%NCL_p)then
                                 abs_rad(ib) = abs_rad(ib) +  (Abs_dir_z(ft,iv) + Abs_dif_z(ft,iv))
                              end if
                              ! Solar radiation absorbed by vegetation and sunlit/shaded leaves
                              do iv = 1,currentPatch%nrad(L,ft)
                                 if(radtype == 1)then
                                    currentPatch%fabd(ib) = currentPatch%fabd(ib) + Abs_dir_z(ft,iv)+Abs_dif_z(ft,iv)
                                    fabd(p,ib) = currentPatch%fabd(ib)
                                 else
                                    currentPatch%fabi(ib) = currentPatch%fabi(ib) + Abs_dif_z(ft,iv)
                                    fabi(p,ib) = currentPatch%fabi(ib)
                                 endif
                              end do
                              ! Albefor
                              if(L==1)then !top canopy layer.
                                 if(radtype == 1)then
                                    albd(p,ib) = albd(p,ib) + Dif_up(L,ft,1) * ftweight(L,ft,1)
                                 else
                                    albi(p,ib) = albi(p,ib) + Dif_up(L,ft,1) * ftweight(L,ft,1)
                                 end if
                              end if
                           end if ! present
                        end do !ft

                        !radiation absorbed from fluxes through unfilled part of lower canopy.
                        if(currentPatch%NCL_p > 1.and.L == currentPatch%NCL_p)then 
                           abs_rad(ib) = abs_rad(ib) + weighted_dif_down(L-1) * &
                                (1.0_r8-sum(ftweight(L,:,1)))*(1.0_r8-albgri(c,ib))
                           abs_rad(ib) = abs_rad(ib) + forc_dir(p,ib) * weighted_dir_tr(L-1) * &
                                (1.0_r8-sum(ftweight(L,:,1)))*(1.0_r8-albgrd(c,ib))
                           tr_soili = tr_soili + weighted_dif_down(L-1) * (1.0_r8-sum(ftweight(L,:,1)))
                           tr_soild = tr_soild + forc_dir(p,ib) * weighted_dir_tr(L-1) * (1.0_r8-sum(ftweight(L,:,1)))
                        endif

                        if(radtype == 1)then
                           currentPatch%tr_soil_dir(ib) = tr_soild
                           currentPatch%tr_soil_dir_dif(ib) = tr_soili
                           currentPatch%sabs_dir(ib)     = abs_rad(ib)
                           ftdd(p,ib)  = tr_soild
                           ftid(p,ib) =  tr_soili
                        else
                           currentPatch%tr_soil_dif(ib) = tr_soili
                           currentPatch%sabs_dif(ib)     = abs_rad(ib)
                           ftii(p,ib) =  tr_soili
                        end if

                     end do!l


                     !==============================================================================!
                     ! Conservation check
                     !==============================================================================!
                     ! Total radiation balance: absorbed = incoming - outgoing

                     if(radtype == 1)then
                        error = abs(currentPatch%sabs_dir(ib)-(currentPatch%tr_soil_dir(ib)*(1.0_r8-albgrd(c,ib))+ &
                             currentPatch%tr_soil_dir_dif(ib)*(1.0_r8-albgri(c,ib))))
                        if( abs(error) > 0.0001)then
                           write(iulog,*)'dir ground absorption error',p,g,error,currentPatch%sabs_dir(ib), &
                                currentPatch%tr_soil_dir(ib)* &
                                (1.0_r8-albgrd(c,ib)),currentPatch%NCL_p,ib,sum(ftweight(1,:,1))
                           write(iulog,*) 'albedos',currentPatch%sabs_dir(ib) ,currentPatch%tr_soil_dir(ib), &
                                (1.0_r8-albgrd(c,ib)),currentPatch%lai

                           do ft =1,3
                              iv = currentPatch%nrad(1,ft) + 1
                              write(iulog,*) 'abs soil fluxes', Abs_dir_z(ft,iv),Abs_dif_z(ft,iv)
                           end do

                        end if
                     else
                        if( abs(currentPatch%sabs_dif(ib)-(currentPatch%tr_soil_dif(ib)*(1.0_r8-albgri(c,ib)))) > 0.0001)then
                           write(iulog,*)'dif ground absorption error',p,g,&
                                currentPatch%sabs_dif(ib),&
                                (currentPatch%tr_soil_dif(ib)*(1.0_r8-albgri(c,ib))),currentPatch%NCL_p,ib,sum(ftweight(1,:,1))
                        endif
                     endif

                     if(radtype == 1)then
                        error = (forc_dir(p,ib) + forc_dif(p,ib)) - (fabd(p,ib)  + albd(p,ib) + currentPatch%sabs_dir(ib))
                     else
                        error = (forc_dir(p,ib) + forc_dif(p,ib)) - (fabi(p,ib)  + albi(p,ib) + currentPatch%sabs_dif(ib))
                     endif
                     lai_reduction(:) = 0.0_r8
                     do L = 1, currentPatch%NCL_p
                        do ft =1,numpft_ed
                           if(currentPatch%present(L,ft) == 1)then
                              do iv = 1, currentPatch%nrad(L,ft)
                                 if(lai_change(L,ft,iv) > 0.0_r8)then
                                    lai_reduction(L) = max(lai_reduction(L),lai_change(L,ft,iv))
                                 endif
                              enddo
                           endif
                        enddo
                     enddo


                     if(radtype == 1)then

                        if (abs(error)  >  0.0001_r8)then
                           write(iulog,*) 'Dir Radn consvn error',error ,p,ib*10+radtype,lai_reduction
                           albd(p,ib) = albd(p,ib) + error 
                           !this terms adds the error back on to the albedo. While this is partly inexcusable, it is 
                           ! in the medium term a solution that
                           ! prevents the model from crashing with small and occasional energy balances issues.
                           ! These are extremely difficult to debug, many have been solved already, leading
                           ! to the complexity of this code, but where the system generates occasional errors, we
                           ! will deal with them for now.
                        end if
                        if (abs(error)  >  0.1_r8)then
                           write(iulog,*) 'Large Dir Radn consvn error',error ,p,ib
                           albd(p,ib) = albd(p,ib) + error
                        end if
                     else

                        if (abs(error)  >  0.0001_r8)then
                           write(iulog,*)  'Dif Radn consvn error',error ,p,ib
                           albi(p,ib) = albi(p,ib) + error
                        end if
                        if (abs(error)  >  0.1_r8)then
                           write(iulog,*)  '>10% Dif Radn consvn error',error ,p,ib
                           albi(p,ib) = albi(p,ib) + error
                        end if

                     end if

                  end do !numrad

               enddo ! rad-type
            endif !associated
         endif ! EDPATCH
      enddo ! loop over fp and indirection to p

    end associate

  end subroutine ED_Norman_Radiation

end module EDSurfaceAlbedoMod
