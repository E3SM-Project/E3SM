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
  
    use physics_types,    only : physics_state
    use camsrfexch,       only : cam_in_t
    use ppgrid,           only : pver
  
    type(carma_type)        :: carma            !! the carma object
    type(carmastate_type)   :: cstate           !! the carma state object
    type(cam_in_t)          :: cam_in
    type(physics_state)     :: state            !! physics state variables
    integer                 :: icol             !! column index
    real(kind=f)            :: cldfrc(pver)     !! total cloud fraction [fraction]
    real(kind=f)            :: rhcrit(pver)     !! realtive humidity for onset of liquid clouds [fraction]
    integer                 :: rc               !! return code, negative indicates failure

    rc = RC_OK

    cldfrc(:) = 1._f
    rhcrit(:) = 1._f
    
    return
  end subroutine CARMA_CloudFraction


