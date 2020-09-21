module ErosionMod 

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate erosion induced soil particulate C, N and P fluxes 
  !
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use decompMod         , only : bounds_type
  use clm_time_manager  , only : get_step_size
  use elm_varcon        , only : dzsoi_decomp
  use elm_varpar        , only : ndecomp_pools, nlevdecomp
  use CNCarbonFluxType  , only : carbonflux_type
  use CNCarbonStateType , only : carbonstate_type
  use CNDecompCascadeConType , only : decomp_cascade_con
  use CNNitrogenFluxType  , only : nitrogenflux_type
  use CNNitrogenStateType , only : nitrogenstate_type
  use PhosphorusFluxType  , only : phosphorusflux_type
  use PhosphorusStateType , only : phosphorusstate_type
  use SedFluxType       , only : sedflux_type
  use SoilStateType     , only : soilstate_type
  use ColumnDataType    , only : col_cs, col_ns, col_ps
  use ColumnDataType    , only : col_cf, col_nf, col_pf
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: ErosionFluxes           ! Calculate erosion introduced CNP fluxes 
  ! !MODULE CONSTANTS
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine ErosionFluxes (bounds, num_soilc, filter_soilc, &
    soilstate_vars, sedflux_vars, carbonstate_vars, nitrogenstate_vars, & 
    phosphorusstate_vars, carbonflux_vars, nitrogenflux_vars, &
    phosphorusflux_vars)
    !
    ! !DESCRIPTION:
    ! Calculate erosion introduced soil C, N, P fluxes 
    !
    ! !USES:
    use clm_varctl      , only : iulog
    use spmdMod         , only : masterproc
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc       ! number of column soil points in column filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    type(sedflux_type)       , intent(in)    :: sedflux_vars
    type(carbonstate_type)   , intent(in)    :: carbonstate_vars
    type(nitrogenstate_type) , intent(in)    :: nitrogenstate_vars
    type(phosphorusstate_type), intent(in)   :: phosphorusstate_vars
    type(carbonflux_type)    , intent(inout) :: carbonflux_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
    type(phosphorusflux_type), intent(inout) :: phosphorusflux_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: c, fc, l, j                              ! indices
    integer  :: kc, kn, kp                               ! new layer indices
    integer  :: k1, k2                                   ! new layer indices     
    real(r8) :: dh, dhn, dhp                             ! erosion height (m)
    real(r8) :: dm_ero, dm_yld, dm                       ! erosion mass (kg/m2)
    real(r8) :: soc_yld                                  ! SOC erosion (gC/m2)
    real(r8) :: zsoi_tot, zsoi_top                       ! soil layer height (m)
    real(r8) :: soc, son                                 ! soil OC and ON content (%)
    real(r8) :: ctot, ntot, ptot                         ! total particulate C, N and P density (g/m3)
    real(r8) :: pn2poc, pp2poc                           ! particulate particle N/C and N/P weight ratios
    real(r8) :: decomp_cpools_vr_new(1:ndecomp_pools)    ! new SOM C pools (gC/m3)
    real(r8) :: decomp_npools_vr_new(1:ndecomp_pools)    ! new SOM N pools (gN/m3)
    real(r8) :: decomp_ppools_vr_new(1:ndecomp_pools)    ! new SOM P pools (gP/m3)
    real(r8) :: labilep_vr_new                           ! new labile mineral P (gP/m3)
    real(r8) :: secondp_vr_new                           ! new secondary mineral P (gP/m3)
    real(r8) :: occlp_vr_new                             ! new occluded mineral P (gP/m3)
    real(r8) :: primp_vr_new                             ! new primary mineral P (gP/m3)
    real(r8) :: dt                                       ! radiation time step (seconds)
    character(len=32) :: subname = 'ErosionFluxes'       ! subroutine name
    !-----------------------------------------------------------------------

    associate(                                                        &
         decomp_cpools_vr =>    col_cs%decomp_cpools_vr       , & ! Input: [real(r8) (:,:,:) ] SOM C pools [gC/m3]
         decomp_npools_vr =>    col_ns%decomp_npools_vr       , & ! Input: [real(r8) (:,:,:) ] SOM N pools [gN/m3]
         decomp_ppools_vr =>    col_ps%decomp_ppools_vr       , & ! Input: [real(r8) (:,:,:) ] SOM P pools [gP/m3]
         labilep_vr       =>    col_ps%labilep_vr             , & ! Input: [real(r8) (:,:) ] soil labile mineral P [gP/m3]
         secondp_vr       =>    col_ps%secondp_vr             , & ! Input: [real(r8) (:,:) ] soil secondary mineral P [gP/m3]
         occlp_vr         =>    col_ps%occlp_vr               , & ! Input: [real(r8) (:,:) ] soil occluded mineral P [gP/m3]
         primp_vr         =>    col_ps%primp_vr               , & ! Input: [real(r8) (:,:) ] soil primary mineral P [gP/m3]

         flx_sed_ero      =>    sedflux_vars%sed_ero_col      , & ! Input: [real(r8) (:) ] sed total detachment (kg/m2/s)
         flx_sed_yld      =>    sedflux_vars%sed_yld_col      , & ! Input: [real(r8) (:) ] sed total yield (kg/m2/s)
         
         bd               =>    soilstate_vars%bd_col         , & ! Input: [real(r8) (:,:) ] soil bulk density (kg/m3)

         labilep_erode    =>    col_pf%labilep_erode          , & ! Output: [real(r8) (:) ] labile mineral P detachment (gP/m2/s)
         labilep_deposit  =>    col_pf%labilep_deposit        , & ! Output: [real(r8) (:) ] labile mineral P hillslope redeposition (gP/m2/s)
         labilep_yield_vr =>    col_pf%labilep_yield_vr       , & ! Output: [real(r8) (:,:) ] vertically-resolved labile mineral P loss (gP/m3/s)
         secondp_erode    =>    col_pf%secondp_erode          , & ! Output: [real(r8) (:) ] secondary mineral P detachment (gP/m2/s)
         secondp_deposit  =>    col_pf%secondp_deposit        , & ! Output: [real(r8) (:) ] secondary mineral P hillslope redeposition (gP/m2/s)
         secondp_yield_vr =>    col_pf%secondp_yield_vr       , & ! Output: [real(r8) (:,:) ] vertically-resolved secondary mineral P loss (gP/m3/s)
         occlp_erode      =>    col_pf%occlp_erode            , & ! Output: [real(r8) (:) ] occluded mineral P detachment (gP/m2/s)
         occlp_deposit    =>    col_pf%occlp_deposit          , & ! Output: [real(r8) (:) ] occluded mineral P hillslope redeposition (gP/m2/s)
         occlp_yield_vr   =>    col_pf%occlp_yield_vr         , & ! Output: [real(r8) (:,:) ] vertically-resolved occluded mineral P loss (gP/m3/s)
         primp_erode      =>    col_pf%primp_erode            , & ! Output: [real(r8) (:) ] primary mineral P detachment (gP/m2/s)
         primp_deposit    =>    col_pf%primp_deposit          , & ! Output: [real(r8) (:) ] primary mineral P hillslope redeposition (gP/m2/s)
         primp_yield_vr   =>    col_pf%primp_yield_vr         , & ! Output: [real(r8) (:,:) ] vertically-resolved primary mineral P loss (gP/m3/s)

         cpools_erode     =>    col_cf%decomp_cpools_erode    , & ! Output: [real(r8) (:,:) ] vertically-integrated decomposing C detachment (gC/m2/s)
         cpools_deposit   =>    col_cf%decomp_cpools_deposit  , & ! Output: [real(r8) (:,:) ] vertically-integrated decomposing C redeposition (gC/m2/s)
         cpools_yield_vr  =>    col_cf%decomp_cpools_yield_vr , & ! Output: [real(r8) (:,:,:) ] vertically-resolved decomposing C loss (gC/m3/s)
                  
         npools_erode     =>    col_nf%decomp_npools_erode    , & ! Output: [real(r8) (:,:) ] vertically-integrated decomposing N detachment (gN/m2/s)
         npools_deposit   =>    col_nf%decomp_npools_deposit  , & ! Output: [real(r8) (:,:) ] vertically-integrated decomposing N redeposition (gN/m2/s)
         npools_yield_vr  =>    col_nf%decomp_npools_yield_vr , & ! Output: [real(r8) (:,:,:) ] vertically-resolved decomposing N loss (gN/m3/s)
         
         ppools_erode     =>    col_pf%decomp_ppools_erode    , & ! Output: [real(r8) (:,:) ] vertically-integrated decomposing P detachment (gP/m2/s)
         ppools_deposit   =>    col_pf%decomp_ppools_deposit  , & ! Output: [real(r8) (:,:) ] vertically-integrated decomposing P redeposition (gP/m2/s)
         ppools_yield_vr  =>    col_pf%decomp_ppools_yield_vr   & ! Output: [real(r8) (:,:,:) ] vertically-resolved decomposing P loss (gP/m3/s)
         )

         dt = real( get_step_size(), r8 )

         ! soil col filters
         do fc = 1, num_soilc
            c = filter_soilc(fc)

            ! initialization
            labilep_erode(c)        = 0._r8
            labilep_deposit(c)      = 0._r8
            labilep_yield_vr(c,:)   = 0._r8
            secondp_erode(c)        = 0._r8
            secondp_deposit(c)      = 0._r8
            secondp_yield_vr(c,:)   = 0._r8
            occlp_erode(c)          = 0._r8
            occlp_deposit(c)        = 0._r8
            occlp_yield_vr(c,:)     = 0._r8
            primp_erode(c)          = 0._r8
            primp_deposit(c)        = 0._r8
            primp_yield_vr(c,:)     = 0._r8

            cpools_erode(c,:)       = 0._r8
            cpools_deposit(c,:)     = 0._r8
            cpools_yield_vr(c,:,:)  = 0._r8
            
            npools_erode(c,:)       = 0._r8
            npools_deposit(c,:)     = 0._r8
            npools_yield_vr(c,:,:)  = 0._r8
            
            ppools_erode(c,:)       = 0._r8
            ppools_deposit(c,:)     = 0._r8
            ppools_yield_vr(c,:,:)  = 0._r8    

            dm_ero = max(0._r8, flx_sed_ero(c)) * dt
            dm_yld = max(0._r8, flx_sed_yld(c)) * dt
            ! calculate erosion heights 
            dh = 0._r8
            dm = dm_yld
            soc_yld = 0._r8
            kc = 1
            do j = 1, nlevdecomp
               ctot = 0._r8
               do l = 1, ndecomp_pools
                  if ( decomp_cascade_con%is_soil(l) .and. all(decomp_cpools_vr(c,:,l)>=0) ) then
                     ctot = ctot + decomp_cpools_vr(c,j,l)
                  end if
               end do
               if (dm<=bd(c,j)*dzsoi_decomp(j)) then
                  dh = dh + dm/bd(c,j) 
                  soc_yld = soc_yld + ctot*dm/bd(c,j)
                  exit
               end if
               dh = dh + dzsoi_decomp(j) 
               dm = dm - bd(c,j)*dzsoi_decomp(j)
               soc_yld = soc_yld + ctot*dzsoi_decomp(j)
               kc = kc + 1
            end do

            ! SOC content in eroded soil (only mud in surface)
            if (dh>0._r8 .and. dh<=dzsoi_decomp(1)) then
               soc = 0.1_r8 * soc_yld / 1460._r8 / dh
            else if (dh>dzsoi_decomp(1)) then
               soc = 0.1_r8 * soc_yld / dm_yld
            else
               soc = 0._r8
            end if
            son = 0.116_r8 * soc - 0.019_r8     ! Beusen et al. (2005)
            if (soc>0._r8 .and. son>0._r8) then
               pn2poc = son / soc 
            else
               pn2poc = 1._r8 / 10.2_r8         ! Beusen et al. (2005)
            end if
            pp2poc = 1._r8 / 22_r8              ! Meybeck (1982)

            ! equivalent erosion heights for N and P
            dm = pn2poc * soc_yld
            dhn = 0._r8
            kn = 1
            do j = 1, nlevdecomp
               ntot = 0._r8
               do l = 1, ndecomp_pools
                  if ( decomp_cascade_con%is_soil(l) .and. all(decomp_npools_vr(c,:,l)>=0._r8) ) then
                     ntot = ntot + decomp_npools_vr(c,j,l)   
                  end if
               end do
               if (dm<=ntot*dzsoi_decomp(j)) then
                  dhn = dhn + dm/ntot
                  exit
               end if
               dhn = dhn + dzsoi_decomp(j)
               dm = dm - ntot*dzsoi_decomp(j)
               kn = kn + 1
            end do

            dm = pp2poc * soc_yld
            dhp = 0._r8
            kp = 1
            do j = 1, nlevdecomp
               ptot = 0._r8
               do l = 1, ndecomp_pools
                  if ( decomp_cascade_con%is_soil(l) .and. all(decomp_ppools_vr(c,:,l)>=0._r8) ) then
                     ptot = ptot + decomp_ppools_vr(c,j,l)
                  end if
               end do
               if ( all(labilep_vr(c,:)>=0._r8) ) then
                  ptot = ptot + labilep_vr(c,j)
               end if
               if ( all(secondp_vr(c,:)>=0._r8) ) then
                  ptot = ptot + secondp_vr(c,j)
               end if
               if ( all(occlp_vr(c,:)>=0._r8) ) then
                  ptot = ptot + occlp_vr(c,j)
               end if
               if ( all(primp_vr(c,:)>=0._r8) ) then
                  ptot = ptot + primp_vr(c,j)
               end if
               if (dm<=ptot*dzsoi_decomp(j)) then
                  dhp = dhp + dm/ptot
                  exit
               end if
               dhp = dhp + dzsoi_decomp(j)
               dm = dm - ptot*dzsoi_decomp(j)
               kp = kp + 1
            end do

            ! particulate C loss
            if ( soc_yld > 0._r8 .and. kc <= nlevdecomp) then
               k1 = kc
               do j = 1, nlevdecomp
                  ! update bottom index
                  if (k1<=nlevdecomp) then
                     zsoi_tot = sum(dzsoi_decomp(1:k1))
                     zsoi_top = sum(dzsoi_decomp(1:j-1)) + dh
                     if (zsoi_top+dzsoi_decomp(j)<zsoi_tot) then
                        k2 = k1
                     else
                        k2 = k1 + 1
                     end if
                  end if
                  if (k1==k2) then
                     do l = 1, ndecomp_pools
                        if ( all(decomp_cpools_vr(c,:,l)>=0._r8) ) then
                           if (k1<=nlevdecomp) then
                              decomp_cpools_vr_new(l) = decomp_cpools_vr(c,k1,l)
                           else
                              decomp_cpools_vr_new(l) = 0._r8
                           end if
                        else
                           decomp_cpools_vr_new(l) = decomp_cpools_vr(c,j,l)
                        end if
                     end do
                  else
                     do l = 1, ndecomp_pools
                        if ( all(decomp_cpools_vr(c,:,l)>=0._r8) ) then
                           if (k2<=nlevdecomp) then
                              decomp_cpools_vr_new(l) = decomp_cpools_vr(c,k1,l)*(zsoi_tot-zsoi_top)/dzsoi_decomp(j) + &
                                 decomp_cpools_vr(c,k2,l)*(1._r8-(zsoi_tot-zsoi_top)/dzsoi_decomp(j)) 
                           else
                              decomp_cpools_vr_new(l) = decomp_cpools_vr(c,k1,l)*(zsoi_tot-zsoi_top)/dzsoi_decomp(j)
                           end if
                        else
                           decomp_cpools_vr_new(l) = decomp_cpools_vr(c,j,l)
                        end if
                     end do
                  end if
                  do l = 1, ndecomp_pools
                     if ( decomp_cascade_con%is_soil(l) ) then
                        cpools_yield_vr(c,j,l) = (decomp_cpools_vr(c,j,l)-decomp_cpools_vr_new(l))/dt
                        cpools_erode(c,l) = cpools_erode(c,l) + cpools_yield_vr(c,j,l)* &
                           dzsoi_decomp(j)*dm_ero/dm_yld
                        cpools_deposit(c,l) = cpools_deposit(c,l) + cpools_erode(c,l) - &
                           cpools_yield_vr(c,j,l)*dzsoi_decomp(j)
                     end if
                  end do
                  k1 = k2     ! update top index
               end do
            end if

            ! particulate N loss
            if ( soc_yld > 0._r8 .and. kn <= nlevdecomp ) then
               k1 = kn
               do j = 1, nlevdecomp 
                  ! update bottom index
                  if (k1<=nlevdecomp) then
                     zsoi_tot = sum(dzsoi_decomp(1:k1))
                     zsoi_top = sum(dzsoi_decomp(1:j-1)) + dhn
                     if (zsoi_top+dzsoi_decomp(j)<zsoi_tot) then
                        k2 = k1
                     else
                        k2 = k1 + 1
                     end if
                  end if
                  if (k1==k2) then
                     do l = 1, ndecomp_pools
                        if ( all(decomp_npools_vr(c,:,l)>=0._r8) ) then
                           if (k1<=nlevdecomp) then
                              decomp_npools_vr_new(l) = decomp_npools_vr(c,k1,l)
                           else
                              decomp_npools_vr_new(l) = 0._r8
                           end if
                        else
                           decomp_npools_vr_new(l) = decomp_npools_vr(c,j,l)
                        end if
                     end do
                  else
                     do l = 1, ndecomp_pools
                        if ( all(decomp_npools_vr(c,:,l)>=0._r8) ) then
                           if (k2<=nlevdecomp) then
                              decomp_npools_vr_new(l) = decomp_npools_vr(c,k1,l)*(zsoi_tot-zsoi_top)/dzsoi_decomp(j) + & 
                                 decomp_npools_vr(c,k2,l)*(1._r8-(zsoi_tot-zsoi_top)/dzsoi_decomp(j))
                           else
                              decomp_npools_vr_new(l) = decomp_npools_vr(c,k1,l)*(zsoi_tot-zsoi_top)/dzsoi_decomp(j)
                           end if
                        else
                           decomp_npools_vr_new(l) = decomp_npools_vr(c,j,l)
                        end if
                     end do
                  end if
                  do l = 1, ndecomp_pools
                     if ( decomp_cascade_con%is_soil(l) ) then
                        npools_yield_vr(c,j,l) = (decomp_npools_vr(c,j,l)-decomp_npools_vr_new(l))/dt 
                        npools_erode(c,l) = npools_erode(c,l) + npools_yield_vr(c,j,l)* &
                           dzsoi_decomp(j)*dm_ero/dm_yld 
                        npools_deposit(c,l) = npools_deposit(c,l) + npools_erode(c,l) - &
                           npools_yield_vr(c,j,l)*dzsoi_decomp(j)
                     end if
                  end do
                  k1 = k2     ! update top index
               end do
            end if

            ! particulate P loss
            if ( soc_yld > 0._r8 .and. kp <= nlevdecomp ) then
               k1 = kp
               do j = 1, nlevdecomp
                  ! update bottom index
                  if (k1<=nlevdecomp) then
                     zsoi_tot = sum(dzsoi_decomp(1:k1))
                     zsoi_top = sum(dzsoi_decomp(1:j-1)) + dhp
                     if (zsoi_top+dzsoi_decomp(j)<zsoi_tot) then
                        k2 = k1
                     else
                        k2 = k1 + 1
                     end if
                  end if
                  if (k1==k2) then
                     do l = 1, ndecomp_pools
                        if ( all(decomp_ppools_vr(c,:,l)>=0._r8) ) then
                           if (k1<=nlevdecomp) then
                              decomp_ppools_vr_new(l) = decomp_ppools_vr(c,k1,l)
                           else
                              decomp_ppools_vr_new(l) = 0._r8
                           end if
                        else
                           decomp_ppools_vr_new(l) = decomp_ppools_vr(c,j,l)
                        end if
                     end do
                     if ( all(labilep_vr(c,:)>=0._r8) ) then
                        if (k1<=nlevdecomp) then
                           labilep_vr_new = labilep_vr(c,k1)
                        else
                           labilep_vr_new = 0._r8
                        end if
                     else
                        labilep_vr_new = labilep_vr(c,j) 
                     end if
                     if ( all(secondp_vr(c,:)>=0._r8) ) then
                        if (k1<=nlevdecomp) then
                           secondp_vr_new = secondp_vr(c,k1)
                        else
                           secondp_vr_new = 0._r8
                        end if
                     else
                        secondp_vr_new = secondp_vr(c,j)
                     end if
                     if ( all(occlp_vr(c,:)>=0._r8) ) then
                        if (k1<=nlevdecomp) then
                           occlp_vr_new = occlp_vr(c,k1)
                        else
                           occlp_vr_new = 0._r8
                        end if
                     else
                        occlp_vr_new = occlp_vr(c,j)
                     end if
                     if ( all(primp_vr(c,:)>=0._r8) ) then
                        if (k1<=nlevdecomp) then
                           primp_vr_new = primp_vr(c,k1)
                        else
                           primp_vr_new = primp_vr(c,nlevdecomp)
                        end if
                     else
                        primp_vr_new = primp_vr(c,j)
                     end if
                  else
                     do l = 1, ndecomp_pools
                        if ( all(decomp_ppools_vr(c,:,l)>=0._r8) ) then
                           if (k2<=nlevdecomp) then
                              decomp_ppools_vr_new(l) = decomp_ppools_vr(c,k1,l)*(zsoi_tot-zsoi_top)/dzsoi_decomp(j) + &
                                 decomp_ppools_vr(c,k2,l)*(1._r8-(zsoi_tot-zsoi_top)/dzsoi_decomp(j)) 
                           else
                              decomp_ppools_vr_new(l) = decomp_ppools_vr(c,k1,l)*(zsoi_tot-zsoi_top)/dzsoi_decomp(j) 
                           end if
                        else
                           decomp_ppools_vr_new(l) = decomp_ppools_vr(c,j,l)
                        end if
                     end do
                     if ( all(labilep_vr(c,:)>=0._r8) ) then
                        if (k2<=nlevdecomp) then
                           labilep_vr_new = labilep_vr(c,k1)*(zsoi_tot-zsoi_top)/dzsoi_decomp(j) + &
                              labilep_vr(c,k2)*(1._r8-(zsoi_tot-zsoi_top)/dzsoi_decomp(j)) 
                        else
                           labilep_vr_new = labilep_vr(c,k1)*(zsoi_tot-zsoi_top)/dzsoi_decomp(j)
                        end if
                     else
                        labilep_vr_new = labilep_vr(c,j)
                     end if
                     if ( all(secondp_vr(c,:)>=0._r8) ) then
                        if (k2<=nlevdecomp) then
                           secondp_vr_new = secondp_vr(c,k1)*(zsoi_tot-zsoi_top)/dzsoi_decomp(j) + &
                              secondp_vr(c,k2)*(1._r8-(zsoi_tot-zsoi_top)/dzsoi_decomp(j))
                        else
                           secondp_vr_new = secondp_vr(c,k1)*(zsoi_tot-zsoi_top)/dzsoi_decomp(j)
                        end if
                     else
                        secondp_vr_new = secondp_vr(c,j)
                     end if
                     if ( all(occlp_vr(c,:)>=0._r8) ) then
                        if (k2<=nlevdecomp) then
                           occlp_vr_new = occlp_vr(c,k1)*(zsoi_tot-zsoi_top)/dzsoi_decomp(j) + &
                              occlp_vr(c,k2)*(1._r8-(zsoi_tot-zsoi_top)/dzsoi_decomp(j))
                        else
                           occlp_vr_new = occlp_vr(c,k1)*(zsoi_tot-zsoi_top)/dzsoi_decomp(j)
                        end if
                     else
                        occlp_vr_new = occlp_vr(c,j)
                     end if
                     if ( all(primp_vr(c,:)>=0._r8) ) then
                        if (k2<=nlevdecomp) then
                           primp_vr_new = primp_vr(c,k1)*(zsoi_tot-zsoi_top)/dzsoi_decomp(j) + &
                              primp_vr(c,k2)*(1._r8-(zsoi_tot-zsoi_top)/dzsoi_decomp(j))
                        else
                           primp_vr_new = primp_vr(c,k1)*(zsoi_tot-zsoi_top)/dzsoi_decomp(j) + &
                              primp_vr(c,nlevdecomp)*(1._r8-(zsoi_tot-zsoi_top)/dzsoi_decomp(j))
                        end if
                     else
                        primp_vr_new = primp_vr(c,j)
                     end if
                  end if
                  do l = 1, ndecomp_pools
                     if ( decomp_cascade_con%is_soil(l) ) then
                        ppools_yield_vr(c,j,l) = (decomp_ppools_vr(c,j,l)-decomp_ppools_vr_new(l))/dt
                        ppools_erode(c,l) = ppools_erode(c,l) + ppools_yield_vr(c,j,l)* &
                           dzsoi_decomp(j)*dm_ero/dm_yld
                        ppools_deposit(c,l) = ppools_deposit(c,l) + ppools_erode(c,l) - &
                           ppools_yield_vr(c,j,l)*dzsoi_decomp(j)
                     end if
                  end do
                  labilep_yield_vr(c,j) = (labilep_vr(c,j)-labilep_vr_new)/dt
                  labilep_erode(c) = labilep_erode(c) + labilep_yield_vr(c,j)* &
                     dzsoi_decomp(j)*dm_ero/dm_yld
                  labilep_deposit(c) = labilep_deposit(c) + labilep_erode(c) - &
                     labilep_yield_vr(c,j)*dzsoi_decomp(j)

                  secondp_yield_vr(c,j) = (secondp_vr(c,j)-secondp_vr_new)/dt
                  secondp_erode(c) = secondp_erode(c) + secondp_yield_vr(c,j)* &
                     dzsoi_decomp(j)*dm_ero/dm_yld
                  secondp_deposit(c) = secondp_deposit(c) + secondp_erode(c) - &
                     secondp_yield_vr(c,j)*dzsoi_decomp(j)

                  occlp_yield_vr(c,j) = (occlp_vr(c,j)-occlp_vr_new)/dt
                  occlp_erode(c) = occlp_erode(c) + occlp_yield_vr(c,j)* &
                     dzsoi_decomp(j)*dm_ero/dm_yld
                  occlp_deposit(c) = occlp_deposit(c) + occlp_erode(c) - &
                     occlp_yield_vr(c,j)*dzsoi_decomp(j)

                  primp_yield_vr(c,j) = (primp_vr(c,j)-primp_vr_new)/dt
                  primp_erode(c) = primp_erode(c) + primp_yield_vr(c,j)* &
                     dzsoi_decomp(j)*dm_ero/dm_yld
                  primp_deposit(c) = primp_deposit(c) + primp_erode(c) - &
                     primp_yield_vr(c,j)*dzsoi_decomp(j)

                  k1 = k2     ! update top index
               end do
            end if

         end do
    
    end associate

  end subroutine ErosionFluxes

end module ErosionMod 
