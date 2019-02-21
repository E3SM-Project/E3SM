module CNErosionMod 

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate the fluxes of SOMC, SOMN and SOMP induced by sediment flux 
  !
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use decompMod         , only : bounds_type
  use clm_time_manager  , only : get_step_size
  use clm_varcon        , only : dzsoi_decomp
  use clm_varpar        , only : ndecomp_pools, nlevdecomp
  use CNCarbonFluxType  , only : carbonflux_type
  use CNCarbonStateType , only : carbonstate_type
  use CNDecompCascadeConType , only : decomp_cascade_con
  use CNNitrogenFluxType  , only : nitrogenflux_type
  use CNNitrogenStateType , only : nitrogenstate_type
  use PhosphorusFluxType  , only : phosphorusflux_type
  use PhosphorusStateType , only : phosphorusstate_type
  use SedFluxType       , only : sedflux_type
  use SoilStateType     , only : soilstate_type
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CNErosionFluxes         ! Calculate erosion introduced CN fluxes 
  ! !MODULE CONSTANTS
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CNErosionFluxes (bounds, num_soilc, filter_soilc, &
    soilstate_vars, sedflux_vars, carbonstate_vars, nitrogenstate_vars, & 
    phosphorusstate_vars, carbonflux_vars, nitrogenflux_vars, &
    phosphorusflux_vars)
    !
    ! !DESCRIPTION:
    ! Calculate soil erosion introduced soil CN fluxes 
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
    integer  :: k, kn, kp                                ! new layer indices
    integer  :: k1, k2                                   ! new layer indices     
    real(r8) :: dh, dhn, dhp                             ! erosion height (m)
    real(r8) :: dm_ero, dm_yld, dm                       ! erosion mass (kg/m2)
    real(r8) :: soc_yld                                  ! SOC erosion (gC/m2)
    real(r8) :: erode_tmp, deposit_tmp                   ! temporal variable
    real(r8) :: zsoi_tot, zsoi_top                       ! soil layer height (m)
    real(r8) :: soc, son                                 ! soil OC and ON content (%)
    real(r8) :: decomp_ctot                              ! total C density (gC/m3)
    real(r8) :: decomp_ntot                              ! total N density (gN/m3)
    real(r8) :: decomp_ptot                              ! total P density (gP/m3)
    real(r8) :: som_n2c, som_p2c                         ! eroded SOM N/C and N/P weight ratios
    real(r8) :: decomp_cpool_vr_new                      ! new SOM C pool (gC/m3)
    real(r8) :: decomp_npool_vr_new                      ! new SOM N pool (gN/m3)
    real(r8) :: decomp_ppool_vr_new                      ! new SOM P pool (gP/m3)
    real(r8) :: dt                                       ! radiation time step (seconds)
    character(len=32) :: subname = 'CNErosionFluxes'     ! subroutine name
    !-----------------------------------------------------------------------

    associate(                                                        &
         decomp_cpools_vr =>    carbonstate_vars%decomp_cpools_vr_col       , & ! Input: [real(r8) (:) ] SOM C pools [gC/m3]
         decomp_npools_vr =>    nitrogenstate_vars%decomp_npools_vr_col     , & ! Input: [real(r8) (:) ] SOM N pools [gN/m3]
         decomp_ppools_vr =>    phosphorusstate_vars%decomp_ppools_vr_col   , & ! Input: [real(r8) (:) ] SOM P pools [gP/m3]

         flx_sed_ero      =>    sedflux_vars%flx_sed_ero_col        , & ! Input: [real(r8) (:) ] sed total detachment (kg/m2/s)
         flx_sed_yld      =>    sedflux_vars%flx_sed_yld_col        , & ! Input: [real(r8) (:) ] sed total yield (kg/m2/s)
         
         bd               =>    soilstate_vars%bd_col               , & ! Input: [real(r8) (:,:) ] soil bulk density (kg/m3)

         cpools_erode     =>    carbonflux_vars%decomp_cpools_erode_col         , & ! Output: [real(r8) (:,:) ] vertically-integrated decomposing C detachment (gC/m2/s)
         cpools_deposit   =>    carbonflux_vars%decomp_cpools_deposit_col       , & ! Output: [real(r8) (:,:) ] vertically-integrated decomposing C redeposition (gC/m2/s)
         cpools_yield_vr  =>    carbonflux_vars%decomp_cpools_yield_vr_col      , & ! Output: [real(r8) (:,:,:) ] vertically-resolved decomposing C loss (gC/m3/s)
                  
         npools_erode     =>    nitrogenflux_vars%decomp_npools_erode_col       , & ! Output: [real(r8) (:,:) ] vertically-integrated decomposing N detachment (gN/m2/s)
         npools_deposit   =>    nitrogenflux_vars%decomp_npools_deposit_col     , & ! Output: [real(r8) (:,:) ] vertically-integrated decomposing N redeposition (gN/m2/s)
         npools_yield_vr  =>    nitrogenflux_vars%decomp_npools_yield_vr_col    , & ! Output: [real(r8) (:,:,:) ] vertically-resolved decomposing N loss (gN/m3/s)
         
         ppools_erode     =>    phosphorusflux_vars%decomp_ppools_erode_col     , & ! Output: [real(r8) (:,:) ] vertically-integrated decomposing P detachment (gP/m2/s)
         ppools_deposit   =>    phosphorusflux_vars%decomp_ppools_deposit_col   , & ! Output: [real(r8) (:,:) ] vertically-integrated decomposing P redeposition (gP/m2/s)
         ppools_yield_vr  =>    phosphorusflux_vars%decomp_ppools_yield_vr_col    & ! Output: [real(r8) (:,:,:) ] vertically-resolved decomposing P loss (gP/m3/s)
         )

         dt = real( get_step_size(), r8 )

         ! soil col filters
         do fc = 1, num_soilc
            c = filter_soilc(fc)

            ! initialization
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
            k = 1
            do j = 1, nlevdecomp
               decomp_ctot = 0._r8
               do l = 1, ndecomp_pools
                  if ( decomp_cascade_con%is_soil(l) ) then
                     decomp_ctot = decomp_ctot + decomp_cpools_vr(c,j,l)
                  end if
               end do
               if (dm<=bd(c,j)*dzsoi_decomp(j)) then
                  dh = dh + dm/bd(c,j) 
                  soc_yld = soc_yld + decomp_ctot*dm/bd(c,j)
                  exit
               end if
               dh = dh + dzsoi_decomp(j) 
               dm = dm - bd(c,j)*dzsoi_decomp(j)
               soc_yld = soc_yld + decomp_ctot*dzsoi_decomp(j)
               k = k + 1
            end do

            ! SOC content in eroded soil (only mud in surface)
            if (dh>0._r8 .and. dh<=dzsoi_decomp(1)) then
               soc = 0.1 * soc_yld / 1460._r8 / dh
            else if (dh>dzsoi_decomp(1)) then
               soc = 0.1 * soc_yld / dm_yld
            else
               soc = 0._r8
            end if
            son = 0.116 * soc - 0.019           ! Beusen et al. (2005)
            if (soc>0._r8 .and. son>0._r8) then
               som_n2c = son / soc 
            else
               som_n2c = 1._r8 / 10.2           ! Beusen et al. (2005)
            end if
            som_p2c = 1._r8 / 22                ! Meybeck (1982)

            ! equivalent erosion heights for N and P
            dm = som_n2c * soc_yld
            dhn = 0._r8
            kn = 1
            do j = 1, nlevdecomp
               decomp_ntot = 0._r8
               do l = 1, ndecomp_pools
                  if ( decomp_cascade_con%is_soil(l) ) then
                     decomp_ntot = decomp_ntot + decomp_npools_vr(c,j,l)   
                  end if
               end do
               if (dm<=decomp_ntot*dzsoi_decomp(j)) then
                  dhn = dhn + dm/decomp_ntot
                  exit
               end if
               dhn = dhn + dzsoi_decomp(j)
               dm = dm - decomp_ntot*dzsoi_decomp(j)
               kn = kn + 1
            end do

            dm = som_p2c * soc_yld
            dhp = 0._r8
            kp = 1
            do j = 1, nlevdecomp
               decomp_ptot = 0._r8
               do l = 1, ndecomp_pools
                  if ( decomp_cascade_con%is_soil(l) ) then
                     decomp_ptot = decomp_ptot + decomp_ppools_vr(c,j,l)
                  end if
               end do
               if (dm<=decomp_ptot*dzsoi_decomp(j)) then
                  dhp = dhp + dm/decomp_ptot
                  exit
               end if
               dhp = dhp + dzsoi_decomp(j)
               dm = dm - decomp_ptot*dzsoi_decomp(j)
               kp = kp + 1
            end do

            do l = 1, ndecomp_pools
               if ( decomp_cascade_con%is_soil(l) ) then
                  ! C, N, P loss from soil column
                  if ( soc_yld > 0._r8 .and. k <= nlevdecomp ) then
                     cpools_erode(c,l) = 0._r8 
                     cpools_deposit(c,l) = 0._r8 
                     k1 = k
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
                           if (k1<=nlevdecomp) then
                              decomp_cpool_vr_new = decomp_cpools_vr(c,k1,l)
                           else
                              decomp_cpool_vr_new = 0._r8
                           end if
                        else
                           if (k2<=nlevdecomp) then
                              decomp_cpool_vr_new = decomp_cpools_vr(c,k1,l)*(zsoi_tot-zsoi_top)/dzsoi_decomp(j) + &
                                 decomp_cpools_vr(c,k2,l)*(1._r8-(zsoi_tot-zsoi_top)/dzsoi_decomp(j))
                           else
                              decomp_cpool_vr_new = decomp_cpools_vr(c,k1,l)*(zsoi_tot-zsoi_top)/dzsoi_decomp(j)
                           end if
                        end if
                        cpools_yield_vr(c,j,l) = (decomp_cpools_vr(c,j,l) - decomp_cpool_vr_new) / dt
                        erode_tmp = cpools_yield_vr(c,j,l)*dzsoi_decomp(j)*dm_ero/dm_yld
                        deposit_tmp = cpools_yield_vr(c,j,l)*dzsoi_decomp(j)*(dm_ero-dm_yld)/dm_yld
                        cpools_erode(c,l) = cpools_erode(c,l) + erode_tmp 
                        cpools_deposit(c,l) = cpools_deposit(c,l) + deposit_tmp
                        k1 = k2     ! update top index
                     end do
                  end if
                  if ( soc_yld > 0._r8 .and. kn <= nlevdecomp ) then
                     npools_erode(c,l) = 0._r8 
                     npools_deposit(c,l) = 0._r8 
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
                           if (k1<=nlevdecomp) then
                              decomp_npool_vr_new = decomp_npools_vr(c,k1,l)
                           else
                              decomp_npool_vr_new = 0._r8
                           end if
                        else
                           if (k2<=nlevdecomp) then
                              decomp_npool_vr_new = decomp_npools_vr(c,k1,l)*(zsoi_tot-zsoi_top)/dzsoi_decomp(j) + &
                                 decomp_npools_vr(c,k2,l)*(1._r8-(zsoi_tot-zsoi_top)/dzsoi_decomp(j))
                           else
                              decomp_npool_vr_new = decomp_npools_vr(c,k1,l)*(zsoi_tot-zsoi_top)/dzsoi_decomp(j)
                           end if
                        end if
                        npools_yield_vr(c,j,l) = (decomp_npools_vr(c,j,l) - decomp_npool_vr_new) / dt
                        erode_tmp = npools_yield_vr(c,j,l)*dzsoi_decomp(j)*dm_ero/dm_yld
                        deposit_tmp = npools_yield_vr(c,j,l)*dzsoi_decomp(j)*(dm_ero-dm_yld)/dm_yld
                        npools_erode(c,l) = npools_erode(c,l) + erode_tmp
                        npools_deposit(c,l) = npools_deposit(c,l) + deposit_tmp
                        k1 = k2     ! update top index
                     end do
                  end if 
                  if ( soc_yld > 0._r8 .and. kp <= nlevdecomp ) then
                     ppools_erode(c,l) = 0._r8 
                     ppools_deposit(c,l) = 0._r8 
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
                           if (k1<=nlevdecomp) then
                              decomp_ppool_vr_new = decomp_ppools_vr(c,k1,l)
                           else
                              decomp_ppool_vr_new = 0._r8
                           end if
                        else
                           if (k2<=nlevdecomp) then
                              decomp_ppool_vr_new = decomp_ppools_vr(c,k1,l)*(zsoi_tot-zsoi_top)/dzsoi_decomp(j) + &
                                 decomp_ppools_vr(c,k2,l)*(1._r8-(zsoi_tot-zsoi_top)/dzsoi_decomp(j))
                           else
                              decomp_ppool_vr_new = decomp_ppools_vr(c,k1,l)*(zsoi_tot-zsoi_top)/dzsoi_decomp(j)
                           end if
                        end if
                        ppools_yield_vr(c,j,l) = (decomp_ppools_vr(c,j,l) - decomp_ppool_vr_new) / dt
                        erode_tmp = ppools_yield_vr(c,j,l)*dzsoi_decomp(j)*dm_ero/dm_yld
                        deposit_tmp = ppools_yield_vr(c,j,l)*dzsoi_decomp(j)*(dm_ero-dm_yld)/dm_yld
                        ppools_erode(c,l) = ppools_erode(c,l) + erode_tmp
                        ppools_deposit(c,l) = ppools_deposit(c,l) + deposit_tmp
                        k1 = k2     ! update top index
                     end do
                  end if
               end if
            end do

         end do
    
    end associate

  end subroutine CNErosionFluxes

end module CNErosionMod 
