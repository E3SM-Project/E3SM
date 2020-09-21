module atm2lndMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle atm2lnd forcing
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use shr_megan_mod  , only : shr_megan_mechcomps_n
  use elm_varpar     , only : numrad, ndst, nlevgrnd !ndst = number of dust bins.
  use elm_varcon     , only : rair, grav, cpair, hfus, tfrz, spval
  use clm_varctl     , only : iulog, use_c13, use_cn, use_lch4, iulog
  use seq_drydep_mod , only : n_drydep, drydep_method, DD_XLND
  use abortutils     , only : endrun
  use decompMod      , only : bounds_type
  use atm2lndType    , only : atm2lnd_type
  use LandunitType   , only : lun_pp                
  use ColumnType     , only : col_pp                
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: downscale_forcings           ! Downscale atm forcing fields from gridcell to column
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: downscale_longwave          ! Downscale longwave radiation from gridcell to column
  private :: build_normalization         ! Compute normalization factors so that downscaled fields are conservative
  private :: check_downscale_consistency ! Check consistency of downscaling
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine downscale_forcings(bounds,num_do_smb_c,filter_do_smb_c, atm2lnd_vars)
    !
    ! !DESCRIPTION:
    ! Downscale atmospheric forcing fields from gridcell to column
    !
    ! Downscaling is done over columns defined by filter_do_smb_c. But we also do direct copies
    ! of gridcell-level forcings into column-level forcings over all other active columns.
    !
    ! !USES:
    use elm_varcon      , only : rair, cpair, grav, lapse_glcmec
    use elm_varcon      , only : glcmec_rain_snow_threshold
    use landunit_varcon , only : istice_mec 
    use clm_varctl      , only : glcmec_downscale_rain_snow_convert
    use domainMod       , only : ldomain
    use QsatMod         , only : Qsat
    !
    ! !ARGUMENTS:
    type(bounds_type)  , intent(in)    :: bounds  
    integer            , intent(in)    :: num_do_smb_c       ! number of columns in filter_do_smb_c
    integer            , intent(in)    :: filter_do_smb_c(:) ! filter_do_smb_c giving columns over which downscaling should be done   
    type(atm2lnd_type) , intent(inout) :: atm2lnd_vars
    !
    ! !LOCAL VARIABLES:
    integer :: g, l, c, fc         ! indices
    integer :: clo, cc

    ! temporaries for topo downscaling
    real(r8) :: hsurf_g,hsurf_c,Hbot
    real(r8) :: zbot_g, tbot_g, pbot_g, thbot_g, qbot_g, qs_g, es_g
    real(r8) :: zbot_c, tbot_c, pbot_c, thbot_c, qbot_c, qs_c, es_c
    real(r8) :: egcm_c, rhos_c
    real(r8) :: dum1,   dum2

    character(len=*), parameter :: subname = 'downscale_forcings'
    !-----------------------------------------------------------------------

    associate(&
         ! Gridcell-level non-downscaled fields:
         forc_t_g     => atm2lnd_vars%forc_t_not_downscaled_grc    , & ! Input:  [real(r8) (:)]  atmospheric temperature (Kelvin)        
         forc_th_g    => atm2lnd_vars%forc_th_not_downscaled_grc   , & ! Input:  [real(r8) (:)]  atmospheric potential temperature (Kelvin)
         forc_q_g     => atm2lnd_vars%forc_q_not_downscaled_grc    , & ! Input:  [real(r8) (:)]  atmospheric specific humidity (kg/kg)   
         forc_pbot_g  => atm2lnd_vars%forc_pbot_not_downscaled_grc , & ! Input:  [real(r8) (:)]  atmospheric pressure (Pa)               
         forc_rho_g   => atm2lnd_vars%forc_rho_not_downscaled_grc  , & ! Input:  [real(r8) (:)]  atmospheric density (kg/m**3)           
         forc_rain_g  => atm2lnd_vars%forc_rain_not_downscaled_grc , & ! Input:  [real(r8) (:)]  rain rate [mm/s]
         forc_snow_g  => atm2lnd_vars%forc_snow_not_downscaled_grc , & ! Input:  [real(r8) (:)]  snow rate [mm/s]
         
         ! Column-level downscaled fields:
         forc_t_c     => atm2lnd_vars%forc_t_downscaled_col        , & ! Output: [real(r8) (:)]  atmospheric temperature (Kelvin)        
         forc_th_c    => atm2lnd_vars%forc_th_downscaled_col       , & ! Output: [real(r8) (:)]  atmospheric potential temperature (Kelvin)
         forc_q_c     => atm2lnd_vars%forc_q_downscaled_col        , & ! Output: [real(r8) (:)]  atmospheric specific humidity (kg/kg)   
         forc_pbot_c  => atm2lnd_vars%forc_pbot_downscaled_col     , & ! Output: [real(r8) (:)]  atmospheric pressure (Pa)               
         forc_rho_c   => atm2lnd_vars%forc_rho_downscaled_col      , & ! Output: [real(r8) (:)]  atmospheric density (kg/m**3)           
         forc_rain_c  => atm2lnd_vars%forc_rain_downscaled_col     , & ! Output: [real(r8) (:)]  rain rate [mm/s]
         forc_snow_c  => atm2lnd_vars%forc_snow_downscaled_col       & ! Output: [real(r8) (:)]  snow rate [mm/s]
         )
      
      ! Initialize column forcing (needs to be done for ALL active columns)
      do c = bounds%begc,bounds%endc
         if (col_pp%active(c)) then
            g = col_pp%gridcell(c)

            forc_t_c(c)     = forc_t_g(g)
            forc_th_c(c)    = forc_th_g(g)
            forc_q_c(c)     = forc_q_g(g)
            forc_pbot_c(c)  = forc_pbot_g(g)
            forc_rho_c(c)   = forc_rho_g(g)
            forc_rain_c(c)  = forc_rain_g(g)
            forc_snow_c(c)  = forc_snow_g(g)
         end if
      end do

      ! Downscale forc_t, forc_th, forc_q, forc_pbot, and forc_rho to columns.
      ! For glacier_mec columns the downscaling is based on surface elevation.
      ! For other columns the downscaling is a simple copy (above).
      do fc = 1, num_do_smb_c
         c = filter_do_smb_c(fc)
         l = col_pp%landunit(c)
         g = col_pp%gridcell(c)

         ! This is a simple downscaling procedure 
         ! Note that forc_hgt, forc_u, and forc_v are not downscaled.

         hsurf_g = ldomain%topo(g)                       ! gridcell sfc elevation
         hsurf_c = col_pp%glc_topo(c)                       ! column sfc elevation
         tbot_g  = forc_t_g(g)                           ! atm sfc temp
         thbot_g = forc_th_g(g)                          ! atm sfc pot temp
         qbot_g  = forc_q_g(g)                           ! atm sfc spec humid
         pbot_g  = forc_pbot_g(g)                        ! atm sfc pressure
         zbot_g  = atm2lnd_vars%forc_hgt_grc(g)          ! atm ref height
         zbot_c  = zbot_g
         tbot_c  = tbot_g-lapse_glcmec*(hsurf_c-hsurf_g) ! sfc temp for column
         Hbot    = rair*0.5_r8*(tbot_g+tbot_c)/grav      ! scale ht at avg temp
         pbot_c  = pbot_g*exp(-(hsurf_c-hsurf_g)/Hbot)   ! column sfc press

         ! Derivation of potential temperature calculation:
         ! 
         ! The textbook definition would be:
         ! thbot_c = tbot_c * (p0/pbot_c)^(rair/cpair)
         ! 
         ! Note that pressure is related to scale height as:
         ! pbot_c = p0 * exp(-zbot_c/H)
         !
         ! Using Hbot in place of H, we get:
         ! pbot_c = p0 * exp(-zbot_c/Hbot)
         !
         ! Plugging this in to the textbook definition, then manipulating, we get:
         ! thbot_c = tbot_c * (p0/(p0*exp(-zbot_c/Hbot)))^(rair/cpair)
         !         = tbot_c * (1/exp(-zbot_c/Hbot))^(rair/cpair)
         !         = tbot_c * (exp(zbot_c/Hbot))^(rair/cpair)
         !         = tbot_c * exp((zbot_c/Hbot) * (rair/cpair))

         thbot_c= tbot_c*exp((zbot_c/Hbot)*(rair/cpair))  ! pot temp calc

         call Qsat(tbot_g,pbot_g,es_g,dum1,qs_g,dum2)
         call Qsat(tbot_c,pbot_c,es_c,dum1,qs_c,dum2)

         qbot_c = qbot_g*(qs_c/qs_g)
         egcm_c = qbot_c*pbot_c/(0.622+0.378*qbot_c)
         rhos_c = (pbot_c-0.378*egcm_c) / (rair*tbot_c)

         forc_t_c(c)    = tbot_c
         forc_th_c(c)   = thbot_c
         forc_q_c(c)    = qbot_c
         forc_pbot_c(c) = pbot_c
         forc_rho_c(c)  = rhos_c

         ! Optionally, convert rain to snow or vice versa based on tbot_c
         ! Note: This conversion does not conserve energy.
         !       It would be better to compute the net latent energy associated
         !        with the conversion and to apply it as a pseudo-flux.
         if (glcmec_downscale_rain_snow_convert) then
            if (tbot_c > glcmec_rain_snow_threshold) then  ! too warm for snow
               forc_rain_c(c) = forc_rain_c(c) + forc_snow_c(c)
               forc_snow_c(c) = 0._r8
            else                                           ! too cold for rain
               forc_snow_c(c) = forc_rain_c(c) + forc_snow_c(c)
               forc_rain_c(c) = 0._r8
            endif
         endif   ! glcmec_downscale_rain_snow_convert

      end do

      call downscale_longwave(bounds, num_do_smb_c, filter_do_smb_c, atm2lnd_vars)

      call check_downscale_consistency(bounds, atm2lnd_vars)

    end associate

  end subroutine downscale_forcings

  !-----------------------------------------------------------------------
  subroutine downscale_longwave(bounds, num_do_smb_c, filter_do_smb_c, atm2lnd_vars)
    !
    ! !DESCRIPTION:
    ! Downscale longwave radiation from gridcell to column
    ! Must be done AFTER temperature downscaling
    !
    ! !USES:
    use domainMod       , only : ldomain
    use landunit_varcon , only : istice_mec 
    use elm_varcon      , only : lapse_glcmec
    use clm_varctl      , only : glcmec_downscale_longwave
    !
    ! !ARGUMENTS:
    type(bounds_type)  , intent(in)    :: bounds  
    integer            , intent(in)    :: num_do_smb_c       ! number of columns in filter_do_smb_c
    integer            , intent(in)    :: filter_do_smb_c(:) ! filter_do_smb_c giving columns over which downscaling should be done (currently glcmec columns)
    type(atm2lnd_type) , intent(inout) :: atm2lnd_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: c,l,g,fc     ! indices
    real(r8) :: hsurf_c      ! column-level elevation (m)
    real(r8) :: hsurf_g      ! gridcell-level elevation (m)

    real(r8), dimension(bounds%begg : bounds%endg) :: sum_lwrad_g    ! weighted sum of column-level lwrad
    real(r8), dimension(bounds%begg : bounds%endg) :: sum_wts_g      ! sum of weights that contribute to sum_lwrad_g
    real(r8), dimension(bounds%begg : bounds%endg) :: lwrad_norm_g   ! normalization factors
    real(r8), dimension(bounds%begg : bounds%endg) :: newsum_lwrad_g ! weighted sum of column-level lwrad after normalization

    character(len=*), parameter :: subname = 'downscale_longwave'
    !-----------------------------------------------------------------------

    associate(&
         ! Gridcell-level fields:
         forc_t_g     => atm2lnd_vars%forc_t_not_downscaled_grc    , & ! Input:  [real(r8) (:)]  atmospheric temperature (Kelvin)        
         forc_lwrad_g => atm2lnd_vars%forc_lwrad_not_downscaled_grc, & ! Input:  [real(r8) (:)]  downward longwave (W/m**2)
         
         ! Column-level (downscaled) fields:
         forc_t_c     => atm2lnd_vars%forc_t_downscaled_col        , & ! Input:  [real(r8) (:)]  atmospheric temperature (Kelvin)        
         forc_lwrad_c => atm2lnd_vars%forc_lwrad_downscaled_col      & ! Output: [real(r8) (:)]  downward longwave (W/m**2)
         )
    
      ! Initialize column forcing (needs to be done for ALL active columns)
      do c = bounds%begc, bounds%endc
         if (col_pp%active(c)) then
            g = col_pp%gridcell(c)
            forc_lwrad_c(c) = forc_lwrad_g(g)
         end if
      end do

      ! Optionally, downscale the longwave radiation, conserving energy
      if (glcmec_downscale_longwave) then

         ! Initialize variables related to normalization
         do g = bounds%begg, bounds%endg
            sum_lwrad_g(g) = 0._r8
            sum_wts_g(g) = 0._r8
            newsum_lwrad_g(g) = 0._r8
         end do

         ! Do the downscaling
         do fc = 1, num_do_smb_c
            c = filter_do_smb_c(fc)
            l = col_pp%landunit(c)
            g = col_pp%gridcell(c)

            hsurf_g = ldomain%topo(g)
            hsurf_c = col_pp%glc_topo(c)

            ! Here we assume that deltaLW = (dLW/dT)*(dT/dz)*deltaz
            ! We get dLW/dT = 4*eps*sigma*T^3 = 4*LW/T from the Stefan-Boltzmann law,
            ! evaluated at the mean temp.
            ! We assume the same temperature lapse rate as above.

            forc_lwrad_c(c) = forc_lwrad_g(g) - &
                 4.0_r8 * forc_lwrad_g(g)/(0.5_r8*(forc_t_c(c)+forc_t_g(g))) * &
                 lapse_glcmec * (hsurf_c - hsurf_g)

            ! Keep track of the gridcell-level weighted sum for later normalization.
            !
            ! This gridcell-level weighted sum just includes points for which we do the
            ! downscaling (e.g., glc_mec points). Thus the contributing weights
            ! generally do not add to 1. So to do the normalization properly, we also
            ! need to keep track of the weights that have contributed to this sum.
            sum_lwrad_g(g) = sum_lwrad_g(g) + col_pp%wtgcell(c)*forc_lwrad_c(c)
            sum_wts_g(g) = sum_wts_g(g) + col_pp%wtgcell(c)
         end do


         ! Normalize forc_lwrad_c(c) to conserve energy

         call build_normalization(orig_field=forc_lwrad_g(bounds%begg:bounds%endg), &
              sum_field=sum_lwrad_g(bounds%begg:bounds%endg), &
              sum_wts=sum_wts_g(bounds%begg:bounds%endg), &
              norms=lwrad_norm_g(bounds%begg:bounds%endg))

         do fc = 1, num_do_smb_c
            c = filter_do_smb_c(fc)
            l = col_pp%landunit(c)
            g = col_pp%gridcell(c)

            forc_lwrad_c(c) = forc_lwrad_c(c) * lwrad_norm_g(g)
            newsum_lwrad_g(g) = newsum_lwrad_g(g) + col_pp%wtgcell(c)*forc_lwrad_c(c)
         end do


         ! Make sure that, after normalization, the grid cell mean is conserved

         do g = bounds%begg, bounds%endg
            if (sum_wts_g(g) > 0._r8) then
               if (abs((newsum_lwrad_g(g) / sum_wts_g(g)) - forc_lwrad_g(g)) > 1.e-8_r8) then
                  write(iulog,*) 'g, newsum_lwrad_g, sum_wts_g, forc_lwrad_g: ', &
                       g, newsum_lwrad_g(g), sum_wts_g(g), forc_lwrad_g(g)
                  call endrun(msg=' ERROR: Energy conservation error downscaling longwave'//&
                       errMsg(__FILE__, __LINE__))
               end if
            end if
         end do

      end if    ! glcmec_downscale_longwave

    end associate

  end subroutine downscale_longwave

  !-----------------------------------------------------------------------
  subroutine build_normalization(orig_field, sum_field, sum_wts, norms)
    !
    ! !DESCRIPTION:
    ! Build an array of normalization factors that can be applied to a downscaled forcing
    ! field, in order to force the mean of the new field to be the same as the mean of
    ! the old field (for conservation).
    !
    ! This allows for the possibility that only a subset of columns are downscaled. Only
    ! the columns that are adjusted should be included in the weighted sum, sum_field;
    ! sum_wts gives the sum of contributing weights on the grid cell level. 

    ! For example, if a grid cell has an original forcing value of 1.0, and contains 4
    ! columns with the following weights on the gridcell, and the following values after
    ! normalization:
    !
    !       col #:    1     2     3     4
    !      weight:  0.1   0.2   0.3   0.4
    ! downscaled?:  yes   yes    no    no
    !       value:  0.9   1.1   1.0   1.0
    !
    ! Then we would have:
    ! orig_field(g) = 1.0
    ! sum_field(g) = 0.1*0.9 + 0.2*1.1 = 0.31
    ! sum_wts(g) = 0.1 + 0.2 = 0.3
    ! norms(g) = 1.0 / (0.31 / 0.3) = 0.9677
    !
    ! The field can then be normalized as:
    !              forc_lwrad_c(c) = forc_lwrad_c(c) * lwrad_norm_g(g)
    !   where lwrad_norm_g is the array of norms computed by this routine

    !
    ! !ARGUMENTS:
    real(r8), intent(in)  :: orig_field(:)  ! the original field, at the grid cell level
    real(r8), intent(in)  :: sum_field(:)   ! the new weighted sum across columns (dimensioned by grid cell)
    real(r8), intent(in)  :: sum_wts(:)     ! sum of the weights used to create sum_field (dimensioned by grid cell)
    real(r8), intent(out) :: norms(:)       ! computed normalization factors
    !-----------------------------------------------------------------------

    SHR_ASSERT((size(orig_field) == size(norms)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT((size(sum_field) == size(norms)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT((size(sum_wts) == size(norms)), errMsg(__FILE__, __LINE__))

    where (sum_wts == 0._r8)
       ! Avoid divide by zero; if sum_wts is 0, then the normalization doesn't matter,
       ! because the adjusted values won't affect the grid cell mean.
       norms = 1.0_r8

    elsewhere (sum_field == 0._r8)
       ! Avoid divide by zero; this should only happen if the gridcell-level value is 0,
       ! in which case the normalization doesn't matter
       norms = 1.0_r8

    elsewhere
       ! The standard case
       norms = orig_field / (sum_field / sum_wts)

    end where

  end subroutine build_normalization


  !-----------------------------------------------------------------------
  subroutine check_downscale_consistency(bounds, atm2lnd_vars)
    !
    ! !DESCRIPTION:
    ! Check consistency of downscaling
    !
    ! Note that this operates over more than just the filter used for the downscaling,
    ! because it checks some things outside that filter.
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type) , intent(in) :: bounds  
    type(atm2lnd_type), intent(in) :: atm2lnd_vars
    !
    ! !LOCAL VARIABLES:
    integer :: g, l, c    ! indices
    character(len=*), parameter :: subname = 'check_downscale_consistency'
    !-----------------------------------------------------------------------

    associate(&
         ! Gridcell-level fields:
         forc_t_g     => atm2lnd_vars%forc_t_not_downscaled_grc     , & ! Input:  [real(r8) (:)]  atmospheric temperature (Kelvin)        
         forc_th_g    => atm2lnd_vars%forc_th_not_downscaled_grc    , & ! Input:  [real(r8) (:)]  atmospheric potential temperature (Kelvin)
         forc_q_g     => atm2lnd_vars%forc_q_not_downscaled_grc     , & ! Input:  [real(r8) (:)]  atmospheric specific humidity (kg/kg)   
         forc_pbot_g  => atm2lnd_vars%forc_pbot_not_downscaled_grc  , & ! Input:  [real(r8) (:)]  atmospheric pressure (Pa)               
         forc_rho_g   => atm2lnd_vars%forc_rho_not_downscaled_grc   , & ! Input:  [real(r8) (:)]  atmospheric density (kg/m**3)           
         forc_rain_g  => atm2lnd_vars%forc_rain_not_downscaled_grc  , & ! Input:  [real(r8) (:)]  rain rate [mm/s]
         forc_snow_g  => atm2lnd_vars%forc_snow_not_downscaled_grc  , & ! Input:  [real(r8) (:)]  snow rate [mm/s]
         forc_lwrad_g => atm2lnd_vars%forc_lwrad_not_downscaled_grc , & ! Input:  [real(r8) (:)]  downward longwave (W/m**2)
         
         ! Column-level (downscaled) fields:
         forc_t_c     => atm2lnd_vars%forc_t_downscaled_col         , & ! Input:  [real(r8) (:)]  atmospheric temperature (Kelvin)        
         forc_th_c    => atm2lnd_vars%forc_th_downscaled_col        , & ! Input:  [real(r8) (:)]  atmospheric potential temperature (Kelvin)
         forc_q_c     => atm2lnd_vars%forc_q_downscaled_col         , & ! Input:  [real(r8) (:)]  atmospheric specific humidity (kg/kg)   
         forc_pbot_c  => atm2lnd_vars%forc_pbot_downscaled_col      , & ! Input:  [real(r8) (:)]  atmospheric pressure (Pa)               
         forc_rho_c   => atm2lnd_vars%forc_rho_downscaled_col       , & ! Input:  [real(r8) (:)]  atmospheric density (kg/m**3)           
         forc_rain_c  => atm2lnd_vars%forc_rain_downscaled_col      , & ! Input:  [real(r8) (:)]  rain rate [mm/s]
         forc_snow_c  => atm2lnd_vars%forc_snow_downscaled_col      , & ! Input:  [real(r8) (:)]  snow rate [mm/s]
         forc_lwrad_c => atm2lnd_vars%forc_lwrad_downscaled_col       & ! Input:  [real(r8) (:)]  downward longwave (W/m**2)
         )
    
      ! Make sure that, for urban points, the column-level forcing fields are identical to
      ! the gridcell-level forcing fields. This is needed because the urban-specific code
      ! sometimes uses the gridcell-level forcing fields (and it would take a large
      ! refactor to change this to use column-level fields).
      
      do c = bounds%begc, bounds%endc
         if (col_pp%active(c)) then
            l = col_pp%landunit(c)
            g = col_pp%gridcell(c)

            if (lun_pp%urbpoi(l)) then
               if (forc_t_c(c)     /= forc_t_g(g)    .or. &
                    forc_th_c(c)    /= forc_th_g(g)   .or. &
                    forc_q_c(c)     /= forc_q_g(g)    .or. &
                    forc_pbot_c(c)  /= forc_pbot_g(g) .or. &
                    forc_rho_c(c)   /= forc_rho_g(g)  .or. &
                    forc_rain_c(c)  /= forc_rain_g(g) .or. &
                    forc_snow_c(c)  /= forc_snow_g(g) .or. &
                    forc_lwrad_c(c) /= forc_lwrad_g(g)) then
                  write(iulog,*) subname//' ERROR: column-level forcing differs from gridcell-level forcing for urban point'
                  write(iulog,*) 'c, g = ', c, g
                  call endrun(msg=errMsg(__FILE__, __LINE__))
               end if  ! inequal
            end if  ! urbpoi
         end if  ! active
      end do

    end associate

  end subroutine check_downscale_consistency

end module atm2lndMod
