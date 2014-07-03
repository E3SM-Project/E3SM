  !! Determine the stratifrom cloud fractions using the CAM routines. This will return the
  !! ice and liquid cloud fractions as well as the minimum relative humidity for the onset
  !! of liquid clouds.
  !!
  !! NOTE: This is just a stub for models that don't use cloud fraction. It should be replaced
  !! be a new routine in a file of the same name in the model directory if the model needs
  !! cloud fraction. This routine needs to be in its own file to avoid circular references when
  !! using the CAM cloud fraction routines (see cirrus model).
  !!
  !!  @version Aug-2010 
  !!  @author  Chuck Bardeen 
  subroutine CARMA_CloudFraction(carma, cstate, cam_in, state, icol, cldfrc, rhcrit, rc)
    use carma_precision_mod
    use carma_enums_mod
    use carma_constants_mod
    use carma_types_mod
    use carma_model_mod
    use carma_flags_mod
    use carmaelement_mod
    use carmagas_mod
    use carmagroup_mod
    use carmasolute_mod
    use carmastate_mod
    use carma_mod
  
    use shr_kind_mod,   only: r8 => shr_kind_r8
    use physics_types,  only: physics_state, physics_ptend, set_wet_to_dry, &
                              set_dry_to_wet
    use constituents,   only: cnst_get_ind
    use abortutils,     only: endrun

    use camsrfexch,       only : cam_in_t
    use ppgrid,           only : pcols, pver, pverp
    use cldwat2m_macro,   only : astG_RHU_single, astG_PDF_single, aist_single, CAMstfrac
  
    type(carma_type)        :: carma            !! the carma object
    type(carmastate_type)   :: cstate           !! the carma state object
    type(cam_in_t)          :: cam_in
    type(physics_state)     :: state            !! physics state variables
    integer                 :: icol             !! column index
    real(r8)                :: cldfrc(pver)     !! total cloud fraction [fraction]
    real(r8)                :: rhcrit(pver)     !! realtive humidity for onset of liquid clouds [fraction]
    integer                 :: rc               !! return code, negative indicates failure

    real(r8)     :: liqcldf(pver)         ! liquid cloud fraction [fraction]
    real(r8)     :: icecldf(pver)         ! ice cloud fraction [fraction]
    real(r8)     :: Ga                    ! dU/da
    real(r8)     :: ssl               
    real(r8)     :: ssi               
    real(r8)     :: qi(pver)              ! ice mass mixing ratio (kg/kg)
    real(r8)     :: mmr(pver)             ! ice mass mixing ratio (kg/kg)
    integer      :: ixcldice
    integer      :: ielem
    integer      :: igroup
    integer      :: ibin
    integer      :: iz
    integer      :: ienconc
    logical      :: is_ice
    logical      :: is_cloud
    logical      :: do_detrain
    
    rc = RC_OK
    
    call CARMA_Get(carma, rc, do_detrain=do_detrain) 
    if (rc < RC_OK) call endrun('CARMA_CloudFraction::CARMA_Get failed.')
 
    ! Get the cloud ice mmr. For the cloud fraction, we only want to include the
    ! ice that is considered in-cloud.
    qi = 0._f
    
    do ielem = 1, NELEM
       
       call CARMAELEMENT_Get(carma, ielem, rc, igroup=igroup)
       if (rc < RC_OK) call endrun('CARMA_CloudFraction::CARMAELEMENT_Get failed.')
       
       call CARMAGROUP_Get(carma, igroup, rc, ienconc=ienconc, is_ice=is_ice, is_cloud=is_cloud)
       if (rc < RC_OK) call endrun('CARMA_CloudFraction::CARMAGROUP_Get failed.')
       
       ! Is this an ice cloud?
       if ((ielem == ienconc) .and. (is_ice) .and. (is_cloud)) then
          
          do ibin = 1, NBIN
             
             ! Get the mass mixing ration for this bin.
             call CARMASTATE_GetBin(cstate, ielem, ibin, mmr, rc)
             if (rc < RC_OK) call endrun('CARMA_CloudFraction::CARMASTATE_GetBin failed.')
             
             ! Add it to the existing ice.
             qi = qi + mmr
             
             ! Add in the detrained ice for this bin.
             if (do_detrain) then
                call CARMASTATE_GetDetrain(cstate, ielem, ibin, mmr, rc)
                if (rc < RC_OK) call endrun('CARMA_CloudFraction::CARMASTATE_GetDetrain failed.')
                
                ! Add it to the existing ice.
                qi = qi + mmr
             end if
          end do
       end if
    end do

    
    ! Calculate the cloud fractions.
    do iz = 1, pver
    
      ! Get a supersaturation that has not been scaled based upon the cloud
      ! fraction.
      call supersat_nocldf(carma, cstate, iz, I_GAS_H2O, ssi, ssl, rc)
    
      ! Get the liquid cloud fraction and the onset humidity for liquid clouds.
      !
      ! NOTE: There is also a PDF based routine, but for now it isn't being used. If
      ! it starts to be used, then a general routine astG_single should be written.
      if (CAMstfrac) then
        call astG_RHU_single(ssl + 1._f, state%pmid(icol, iz), state%q(icol, iz, 1), &
                cam_in%landfrac(icol), cam_in%snowhland(icol), liqcldf(iz), Ga, rhcrit(iz))
      else
        call astG_PDF_single(ssl + 1._f, state%pmid(icol, iz), state%q(icol, iz, 1), &
                cam_in%landfrac(icol), cam_in%snowhland(icol), liqcldf(iz), Ga, rhcrit(iz))
      end if

      ! Now get the ice cloud fraction.
      call aist_single(state%q(icol, iz, 1), state%t(icol, iz), state%pmid(icol, iz), &
              qi(iz), cam_in%landfrac(icol), cam_in%snowhland(icol), icecldf(iz))
    end do
              
    ! Calculate an overall cloud fraction. This may vary depending upon the model,
    ! but defaults to minimum overlap (Wilson and Ballard, 1999). This may not be
    ! the same as the assumptions made by the CAM cloud scheme.
    cldfrc(:) = min(1.0_f, icecldf(:) + liqcldf(:))

    ! For the cirrus model, we just want to use the ice cloud fraction.
!    cldfrc(:) = icecldf(:)

    ! Don't let the cloud fraction get too small.
    cldfrc(:) = max(CLDFRC_MIN, cldfrc(:))
    
    return
  end subroutine CARMA_CloudFraction


