!! The CARMAGROUP module contains configuration information about a CARMA partcile.
!!
!! NOTE: Because of the way Fortran handles pointers and allocations, it is much
!! simpiler to have these methods directly access the group array that is in the
!! CARMA object rather than having this as its own objects. Some compilers (like
!! IBM on AIX do not by default automatically deallocate automatically created
!! derived types that contain allocations. This can result in memory leaks that
!! are difficult to find.
!!
!! These calls are written like they are part of CARMA, but they are called
!! CARMAGROUP and kept by themselves in their own file to make it easier to keep
!! track of what is required when adding an attribute to a group.
!!
!!  @version July-2009 
!!  @author  Chuck Bardeen 
module carmagroup_mod

  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod

  ! CARMA explicitly declares all variables. 
  implicit none

  ! All CARMA variables and procedures are private except those explicitly declared to be public.
  private

  ! Declare the public methods.
  public CARMAGROUP_Create
  public CARMAGROUP_Destroy
  public CARMAGROUP_Get
  public CARMAGROUP_Print

contains

  subroutine CARMAGROUP_Create(carma, igroup, name, rmin, rmrat, ishape, eshape, is_ice, &
      rc, irhswell, irhswcomp, refidx, do_mie, do_wetdep, do_drydep, do_vtran, solfac, scavcoef, shortname, &
      cnsttype, maxbin, ifallrtn, is_cloud, rmassmin, imiertn, is_sulfate, dpc_threshold)
    type(carma_type), intent(inout)             :: carma               !! the carma object
    integer, intent(in)                         :: igroup              !! the group index
    character(*), intent(in)                    :: name                !! the group name, maximum of 255 characters
    real(kind=f), intent(in)                    :: rmin                !! the minimum radius, can be specified [cm]
    real(kind=f), intent(in)                    :: rmrat               !! the volume ratio between bins
    integer, intent(in)                         :: ishape              !! the type of the particle shape [I_SPHERE | I_HEXAGON | I_CYLINDER]
    real(kind=f), intent(in)                    :: eshape              !! the aspect ratio of the particle shape (length/diameter)
    logical, intent(in)                         :: is_ice              !! is this an ice particle?
    integer, intent(out)                        :: rc                  !! return code, negative indicates failure
    integer, optional, intent(in)               :: irhswell            !! the parameterization for particle swelling from relative humidity [I_FITZGERALD | I_GERBER]
    integer, optional, intent(in)               :: irhswcomp           !! the composition for particle swelling from relative humidity [I_FITZGERALD | I_GERBER]
    complex(kind=f), optional, intent(in)       :: refidx(carma%f_NWAVE) !! refractive index for the particle
    logical, optional, intent(in)               :: do_mie              !! do mie calculations?
    logical, optional, intent(in)               :: do_wetdep           !! do wet deposition for this particle?
    logical, optional, intent(in)               :: do_drydep           !! do dry deposition for this particle?
    logical, optional, intent(in)               :: do_vtran            !! do sedimentation for this particle?
    real(kind=f), intent(in), optional          :: solfac              !! the solubility factor for wet deposition
    real(kind=f), intent(in), optional          :: scavcoef            !! the scavenging coefficient for wet deposition
    character(*), optional, intent(in)          :: shortname           !! the group shortname, maximum of 6 characters
    integer, optional, intent(in)               :: cnsttype            !! constituent type in parent model [I_CNSTTYPE_PROGNOSTIC | I_CNSTTYPE_DIAGNOSTIC]
    integer, optional, intent(in)               :: maxbin              !! bin number of the last prognostic bin, the remaining bins are diagnostic
    integer, optional, intent(in)               :: ifallrtn            !! fall velocity routine [I_FALLRTN_STD | I_FALLRTN_STD_SHAPE | I_FALLRTN_HEYMSFIELD2010 | I_FALLRTN_ACKERMAN_DROP | I_FALLRTN_ACKERMAN_ICE]
    logical, optional, intent(in)               :: is_cloud            !! is this a cloud particle?
    real(kind=f), optional, intent(in)          :: rmassmin            !! the minimum mass, when used overrides rmin[g]
    integer, optional, intent(in)               :: imiertn             !! mie routine [I_MIERTN_TOON1981 | I_MIERTN_BOHREN1983]
    logical, optional, intent(in)               :: is_sulfate          !! is this a sulfate particle?
    real(kind=f), optional, intent(in)          :: dpc_threshold       !! convergence criteria for particle concentration [fraction]

    ! Local variables
    integer                               :: ier
    
    ! Assume success.
    rc = RC_OK
    
    ! Make sure there are enough groups allocated.
    if (igroup > carma%f_NGROUP) then
      if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMAGROUP_Add:: ERROR - The specifed group (", &
        igroup, ") is larger than the number of groups (", carma%f_NGROUP, ")."
      rc = RC_ERROR
      return
    end if
    
    allocate( &
      carma%f_group(igroup)%f_r(carma%f_NBIN), &
      carma%f_group(igroup)%f_rmass(carma%f_NBIN), &
      carma%f_group(igroup)%f_vol(carma%f_NBIN), &
      carma%f_group(igroup)%f_dr(carma%f_NBIN), &
      carma%f_group(igroup)%f_dm(carma%f_NBIN), &
      carma%f_group(igroup)%f_rmassup(carma%f_NBIN), &
      carma%f_group(igroup)%f_rup(carma%f_NBIN), &
      carma%f_group(igroup)%f_rlow(carma%f_NBIN), &
      carma%f_group(igroup)%f_icorelem(carma%f_NELEM), &
      carma%f_group(igroup)%f_arat(carma%f_NBIN), &
      carma%f_group(igroup)%f_rrat(carma%f_NBIN), &
      stat=ier) 
    if(ier /= 0) then
        if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMAGROUP_Add: ERROR allocating, status=", ier
      rc = RC_ERROR
      return
    end if
    
    ! Initialize
    carma%f_group(igroup)%f_r(:)        = 0._f
    carma%f_group(igroup)%f_rmass(:)    = 0._f
    carma%f_group(igroup)%f_vol(:)      = 0._f
    carma%f_group(igroup)%f_dr(:)       = 0._f
    carma%f_group(igroup)%f_dm(:)       = 0._f
    carma%f_group(igroup)%f_rmassup(:)  = 0._f
    carma%f_group(igroup)%f_rup(:)      = 0._f
    carma%f_group(igroup)%f_rlow(:)     = 0._f
    carma%f_group(igroup)%f_icorelem(:) = 0
    carma%f_group(igroup)%f_ifallrtn    = I_FALLRTN_STD
    carma%f_group(igroup)%f_imiertn     = I_MIERTN_TOON1981
    carma%f_group(igroup)%f_is_cloud    = .false.
    carma%f_group(igroup)%f_is_sulfate  = .false.
    carma%f_group(igroup)%f_dpc_threshold = 0._f
    

    ! Any optical properties?
    if (carma%f_NWAVE > 0) then
      allocate( &
        carma%f_group(igroup)%f_refidx(carma%f_NWAVE), &
        carma%f_group(igroup)%f_qext(carma%f_NWAVE,carma%f_NBIN), &
        carma%f_group(igroup)%f_ssa(carma%f_NWAVE,carma%f_NBIN), &
        carma%f_group(igroup)%f_asym(carma%f_NWAVE,carma%f_NBIN), &
        stat=ier) 
      if(ier /= 0) then
        if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMAGROUP_Add: ERROR allocating, status=", ier
        rc = RC_ERROR
        return
      endif

      ! Initialize
      carma%f_group(igroup)%f_refidx(:) = (0._f, 0._f)
      carma%f_group(igroup)%f_qext(:,:) = 0._f
      carma%f_group(igroup)%f_ssa(:,:)  = 0._f
      carma%f_group(igroup)%f_asym(:,:) = 0._f
    end if
    

    ! Save off the settings.
    carma%f_group(igroup)%f_name        = name
    carma%f_group(igroup)%f_rmin        = rmin
    carma%f_group(igroup)%f_rmrat       = rmrat
    carma%f_group(igroup)%f_ishape      = ishape
    carma%f_group(igroup)%f_eshape      = eshape
    carma%f_group(igroup)%f_is_ice      = is_ice
    
    
    ! Defaults for optional parameters
    carma%f_group(igroup)%f_irhswell    = 0
    carma%f_group(igroup)%f_do_mie      = .false.
    carma%f_group(igroup)%f_do_wetdep   = .false.
    carma%f_group(igroup)%f_grp_do_drydep   = .false.
    carma%f_group(igroup)%f_grp_do_vtran  = .true.
    carma%f_group(igroup)%f_solfac      = 0.3_f
    carma%f_group(igroup)%f_scavcoef    = 0.1_f
    carma%f_group(igroup)%f_shortname   = ""
    carma%f_group(igroup)%f_cnsttype    = I_CNSTTYPE_PROGNOSTIC
    carma%f_group(igroup)%f_maxbin      = carma%f_NBIN
    carma%f_group(igroup)%f_rmassmin    = 0.0_f
    
    ! Set optional parameters.
    if (present(irhswell))   carma%f_group(igroup)%f_irhswell     = irhswell
    if (present(irhswcomp))  carma%f_group(igroup)%f_irhswcomp    = irhswcomp
    if (present(refidx))     carma%f_group(igroup)%f_refidx(:)    = refidx(:)
    if (present(do_mie))     carma%f_group(igroup)%f_do_mie       = do_mie
    if (present(do_wetdep))  carma%f_group(igroup)%f_do_wetdep    = do_wetdep
    if (present(do_drydep))  carma%f_group(igroup)%f_grp_do_drydep  = do_drydep
    if (present(do_vtran))   carma%f_group(igroup)%f_grp_do_vtran = do_vtran
    if (present(solfac))     carma%f_group(igroup)%f_solfac       = solfac
    if (present(scavcoef))   carma%f_group(igroup)%f_scavcoef     = scavcoef
    if (present(shortname))  carma%f_group(igroup)%f_shortname    = shortname
    if (present(cnsttype))   carma%f_group(igroup)%f_cnsttype     = cnsttype
    if (present(maxbin))     carma%f_group(igroup)%f_maxbin       = maxbin
    if (present(ifallrtn))   carma%f_group(igroup)%f_ifallrtn     = ifallrtn
    if (present(is_cloud))   carma%f_group(igroup)%f_is_cloud     = is_cloud
    if (present(rmassmin))   carma%f_group(igroup)%f_rmassmin     = rmassmin
    if (present(imiertn))    carma%f_group(igroup)%f_imiertn      = imiertn
    if (present(is_sulfate)) carma%f_group(igroup)%f_is_sulfate   = is_sulfate
    if (present(dpc_threshold)) carma%f_group(igroup)%f_dpc_threshold = dpc_threshold

    
    ! Initialize other properties.
    carma%f_group(igroup)%f_nelem         = 0
    carma%f_group(igroup)%f_if_sec_mom    = .FALSE.
    carma%f_group(igroup)%f_ncore         = 0
    carma%f_group(igroup)%f_ienconc       = 0
    carma%f_group(igroup)%f_imomelem      = 0
    
    
    ! The area ratio is the ratio of the area of the shape to the area of the
    ! circumscribing circle. The radius ratio is the ratio between the radius
    ! of the longest dimension and the radius of the enclosing sphere.
    if (ishape .eq. I_HEXAGON) then
      carma%f_group(igroup)%f_arat(:) = 3._f * sqrt(3._f) / 2._f / PI 
      carma%f_group(igroup)%f_rrat(:) = ((4._f * PI / 9._f / sqrt(3._f)) ** (1._f / 3._f)) * eshape**(-1._f / 3._f)
    else if (ishape .eq. I_CYLINDER) then
      carma%f_group(igroup)%f_arat(:) = 1.0_f
      carma%f_group(igroup)%f_rrat(:) = ((2._f / 3._f) ** (1._f / 3._f)) * eshape**(-1._f / 3._f)
    else

      ! Default to a sphere.
      !
      ! NOTE: Should add code here to handle oblate and prolate spheroids.
      carma%f_group(igroup)%f_arat(:) = 1.0_f
      carma%f_group(igroup)%f_rrat(:) = 1.0_f
    end if
    
    return
  end subroutine CARMAGROUP_Create
    

  !! Deallocates the memory associated with a CARMAGROUP object.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  !!
  !! @see CARMAGROUP_Create
  subroutine CARMAGROUP_Destroy(carma, igroup, rc)
    type(carma_type), intent(inout)      :: carma         !! the carma object
    integer, intent(in)                  :: igroup        !! the group index
    integer, intent(out)                 :: rc            !! return code, negative indicates failure

    ! Local variables
    integer                              :: ier
    
    ! Assume success.
    rc = RC_OK
    
    ! Make sure there are enough groups allocated.
    if (igroup > carma%f_NGROUP) then
      if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMAGROUP_Destroy:: ERROR - The specifed group (", &
        igroup, ") is larger than the number of groups (", carma%f_NGROUP, ")."
      rc = RC_ERROR
      return
    end if
    
    if (allocated(carma%f_group(igroup)%f_refidx)) then
      deallocate( &
        carma%f_group(igroup)%f_refidx, &
        carma%f_group(igroup)%f_qext, &
        carma%f_group(igroup)%f_ssa, &
        carma%f_group(igroup)%f_asym, &
        stat=ier) 
      if(ier /= 0) then
        if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMAGROUP_Destroy: ERROR deallocating, status=", ier
        rc = RC_ERROR
        return
      endif
    endif
    
    ! Allocate dynamic data.
    if (allocated(carma%f_group(igroup)%f_r)) then
      deallocate( &
        carma%f_group(igroup)%f_r, &
        carma%f_group(igroup)%f_rmass, &
        carma%f_group(igroup)%f_vol, &
        carma%f_group(igroup)%f_dr, &
        carma%f_group(igroup)%f_dm, &
        carma%f_group(igroup)%f_rmassup, &
        carma%f_group(igroup)%f_rup, &
        carma%f_group(igroup)%f_rlow, &
        carma%f_group(igroup)%f_icorelem, &
        carma%f_group(igroup)%f_arat, &
        carma%f_group(igroup)%f_rrat, &
        stat=ier) 
      if(ier /= 0) then
        if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMAGROUP_Destroy: ERROR deallocating, status=", ier
        rc = RC_ERROR
        return
      endif
    endif

    return
  end subroutine CARMAGROUP_Destroy


  !! Gets information about a group.
  !!
  !! The group name and most other properties are available after a call to
  !! CARMAGROUP_Create(). After a call to CARMA_Initialize(), the bin
  !! dimensions and optical properties can be retrieved.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  !!
  !! @see CARMAGROUP_Create
  !! @see CARMA_GetGroup
  !! @see CARMA_Initialize 
  subroutine CARMAGROUP_Get(carma, igroup, rc, name, shortname, rmin, rmrat, ishape, eshape, is_ice, &
      irhswell, irhswcomp, cnsttype, r, rlow, rup, dr, rmass, dm, vol, qext, ssa, asym, do_mie, &
      do_wetdep, do_drydep, do_vtran, solfac, scavcoef, ienconc, refidx, ncore, icorelem, maxbin, &
      ifallrtn, is_cloud, rmassmin, arat, rrat, imiertn, is_sulfate, dpc_threshold)
      
    type(carma_type), intent(in)              :: carma                        !! the carma object
    integer, intent(in)                       :: igroup                       !! the group index
    integer, intent(out)                      :: rc                           !! return code, negative indicates failure
    character(len=*), optional, intent(out)   :: name                         !! the group name
    character(len=*), optional, intent(out)   :: shortname                    !! the group short name
    real(kind=f), optional, intent(out)       :: rmin                         !! the minimum radius [cm]
    real(kind=f), optional, intent(out)       :: rmrat                        !! the volume ratio between bins
    integer, optional, intent(out)            :: ishape                       !! the type of the particle shape
    real(kind=f), optional, intent(out)       :: eshape                       !! the aspect ratio of the particle shape
    logical, optional, intent(out)            :: is_ice                       !! is this an ice particle?
    integer, optional, intent(out)            :: irhswell                     !! the parameterization for particle swelling from relative humidity
    integer, optional, intent(out)            :: irhswcomp                    !! the composition for particle swelling from relative humidity
    integer, optional, intent(out)            :: cnsttype                     !! constituent type in the parent model
    real(kind=f), intent(out), optional       :: r(carma%f_NBIN)                !! the bin radius [cm]
    real(kind=f), intent(out), optional       :: rlow(carma%f_NBIN)             !! the bin radius lower bound [cm]
    real(kind=f), intent(out), optional       :: rup(carma%f_NBIN)              !! the bin radius upper bound [cm]
    real(kind=f), intent(out), optional       :: dr(carma%f_NBIN)               !! the bin width in radius space [cm]
    real(kind=f), intent(out), optional       :: rmass(carma%f_NBIN)            !! the bin mass [g]
    real(kind=f), intent(out), optional       :: dm(carma%f_NBIN)               !! the bin width in mass space [g]
    real(kind=f), intent(out), optional       :: vol(carma%f_NBIN)              !! the bin volume [cm<sup>3</sup>]
    real(kind=f), intent(out), optional       :: rrat(carma%f_NBIN)             !! the radius ratio (maximum dimension / radius of enclosing sphere)
    real(kind=f), intent(out), optional       :: arat(carma%f_NBIN)             !! the projected area ratio (area / area enclosing sphere)
    complex(kind=f), intent(out), optional    :: refidx(carma%f_NWAVE)          !! the refractive index at each wavelength
    real(kind=f), intent(out), optional       :: qext(carma%f_NWAVE,carma%f_NBIN) !! extinction efficiency
    real(kind=f), intent(out), optional       :: ssa(carma%f_NWAVE,carma%f_NBIN)  !! single scattering albedo
    real(kind=f), intent(out), optional       :: asym(carma%f_NWAVE,carma%f_NBIN) !! asymmetry factor
    logical, optional, intent(out)            :: do_mie                       !! do mie calculations?
    logical, optional, intent(out)            :: do_wetdep                    !! do wet deposition for this particle?
    logical, optional, intent(out)            :: do_drydep                    !! do dry deposition for this particle?
    logical, optional, intent(out)            :: do_vtran                     !! do sedimentation for this particle?
    real(kind=f), intent(out), optional       :: solfac                       !! the solubility factor for wet deposition
    real(kind=f), intent(out), optional       :: scavcoef                     !! the scavenging coefficient for wet deposition
    integer, intent(out), optional            :: ienconc                      !! Particle number conc. element for group
    integer, intent(out), optional            :: ncore                        !! Number of core mass elements for group
    integer, intent(out), optional            :: icorelem(carma%f_NELEM)        !! Element index of core mass elements for group
    integer, optional, intent(out)            :: maxbin                       !! the last prognostic bin in the group
    integer, optional, intent(out)            :: ifallrtn                     !! fall velocity routine [I_FALLRTN_STD | I_FALLRTN_STD_SHAPE | I_FALLRTN_HEYMSFIELD2010 | I_FALLRTN_ACKERMAN_DROP | I_FALLRTN_ACKERMAN_ICE]
    logical, optional, intent(out)            :: is_cloud                     !! is this a cloud particle?
    real(kind=f), optional, intent(out)       :: rmassmin                     !! the minimum mass [g]
    integer, optional, intent(out)            :: imiertn                      !! mie routine [I_MIERTN_TOON1981 | I_MIERTN_BOHREN1983]
    logical, optional, intent(out)            :: is_sulfate                   !! is this a sulfate particle?
    real(kind=f), optional, intent(out)       :: dpc_threshold                !! convergence criteria for particle concentration [fraction]

    ! Assume success.
    rc = RC_OK

    ! Make sure there are enough groups allocated.
    if (igroup > carma%f_NGROUP) then
      if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMAGROUP_Get:: ERROR - The specifed group (", &
        igroup, ") is larger than the number of groups (", carma%f_NGROUP, ")."
      rc = RC_ERROR
      return
    end if
      
    ! Return any requested properties of the group.
    if (present(name))         name         = carma%f_group(igroup)%f_name
    if (present(shortname))    shortname    = carma%f_group(igroup)%f_shortname
    if (present(rmin))         rmin         = carma%f_group(igroup)%f_rmin
    if (present(rmrat))        rmrat        = carma%f_group(igroup)%f_rmrat
    if (present(ishape))       ishape       = carma%f_group(igroup)%f_ishape
    if (present(eshape))       eshape       = carma%f_group(igroup)%f_eshape
    if (present(is_ice))       is_ice       = carma%f_group(igroup)%f_is_ice
    if (present(irhswell))     irhswell     = carma%f_group(igroup)%f_irhswell
    if (present(irhswcomp))    irhswcomp    = carma%f_group(igroup)%f_irhswcomp
    if (present(cnsttype))     cnsttype     = carma%f_group(igroup)%f_cnsttype
    if (present(r))            r(:)         = carma%f_group(igroup)%f_r(:)
    if (present(rlow))         rlow(:)      = carma%f_group(igroup)%f_rlow(:)
    if (present(rup))          rup(:)       = carma%f_group(igroup)%f_rup(:)
    if (present(dr))           dr(:)        = carma%f_group(igroup)%f_dr(:)
    if (present(rmass))        rmass(:)     = carma%f_group(igroup)%f_rmass(:)
    if (present(rrat))         rrat(:)      = carma%f_group(igroup)%f_rrat(:)
    if (present(arat))         arat(:)      = carma%f_group(igroup)%f_arat(:)
    if (present(dm))           dm(:)        = carma%f_group(igroup)%f_dm(:)
    if (present(vol))          vol(:)       = carma%f_group(igroup)%f_vol(:)
    if (present(do_mie))       do_mie       = carma%f_group(igroup)%f_do_mie
    if (present(do_wetdep))    do_wetdep    = carma%f_group(igroup)%f_do_wetdep
    if (present(do_drydep))    do_drydep    = carma%f_group(igroup)%f_grp_do_drydep
    if (present(do_vtran))     do_vtran     = carma%f_group(igroup)%f_grp_do_vtran
    if (present(solfac))       solfac       = carma%f_group(igroup)%f_solfac
    if (present(scavcoef))     scavcoef     = carma%f_group(igroup)%f_scavcoef
    if (present(ienconc))      ienconc      = carma%f_group(igroup)%f_ienconc
    if (present(ncore))        ncore        = carma%f_group(igroup)%f_ncore
    if (present(icorelem))     icorelem     = carma%f_group(igroup)%f_icorelem(:)
    if (present(maxbin))       maxbin       = carma%f_group(igroup)%f_maxbin
    if (present(ifallrtn))     ifallrtn     = carma%f_group(igroup)%f_ifallrtn
    if (present(is_cloud))     is_cloud     = carma%f_group(igroup)%f_is_cloud
    if (present(rmassmin))     rmassmin     = carma%f_group(igroup)%f_rmassmin
    if (present(imiertn))      imiertn      = carma%f_group(igroup)%f_imiertn
    if (present(is_sulfate))   is_sulfate   = carma%f_group(igroup)%f_is_sulfate
    if (present(dpc_threshold)) dpc_threshold = carma%f_group(igroup)%f_dpc_threshold
    
    if (carma%f_NWAVE == 0) then
      if (present(refidx) .or. present(qext) .or. present(ssa) .or. present(asym)) then
        if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMAGROUP_Get: ERROR no optical properties defined."
        rc = RC_ERROR
        return
      end if
    else
      if (present(refidx))     refidx(:)    = carma%f_group(igroup)%f_refidx(:)
      if (present(qext))       qext(:,:)    = carma%f_group(igroup)%f_qext(:,:)
      if (present(ssa))        ssa(:,:)     = carma%f_group(igroup)%f_ssa(:,:)
      if (present(asym))       asym(:,:)    = carma%f_group(igroup)%f_asym(:,:)
    end if
    
    return
  end subroutine CARMAGROUP_Get
  
  
  
  !! Prints information about a group.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  !!
  !! @see CARMAGROUP_Get
  subroutine CARMAGROUP_Print(carma, igroup, rc)
    type(carma_type), intent(in)              :: carma              !! the carma object
    integer, intent(in)                       :: igroup             !! the group index
    integer, intent(out)                      :: rc                 !! return code, negative indicates failure
    
    ! Local variables
    integer                                   :: i
    character(len=CARMA_NAME_LEN)             :: name               ! name
    character(len=CARMA_SHORT_NAME_LEN)       :: shortname          ! shortname
    real(kind=f)                              :: rmin               ! the minimum radius [cm]
    real(kind=f)                              :: rmrat              ! the volume ratio between bins
    integer                                   :: ishape             ! the type of the particle shape
    real(kind=f)                              :: eshape             ! the aspect ratio of the particle shape
    logical                                   :: is_ice             ! is this an ice particle?
    integer                                   :: irhswell           ! the parameterization for particle swelling from relative humidity
    integer                                   :: irhswcomp          ! the composition for particle swelling from relative humidity
    integer                                   :: cnsttype           ! constituent type in the parent model
    real(kind=f)                              :: r(carma%f_NBIN)      ! the bin radius [m]
    real(kind=f)                              :: dr(carma%f_NBIN)     ! the bin width in radius space [m]
    real(kind=f)                              :: rmass(carma%f_NBIN)  ! the bin mass [kg]
    real(kind=f)                              :: dm(carma%f_NBIN)     ! the bin width in mass space [kg]
    real(kind=f)                              :: vol(carma%f_NBIN)    ! the bin volume [m<sup>3</sup>]
    integer                                   :: ifallrtn           ! fall velocity routine [I_FALLRTN_STD | I_FALLRTN_STD_SHAPE | I_FALLRTN_HEYMSFIELD2010 | I_FALLRTN_ACKERMAN_DROP | I_FALLRTN_ACKERMAN_ICE]
    logical                                   :: is_cloud           ! is this a cloud particle?
    real(kind=f)                              :: rmassmin           ! the minimum mass [g]
    logical                                   :: do_mie             ! do mie calculations?
    logical                                   :: do_wetdep          ! do wet deposition for this particle?
    logical                                   :: do_drydep          ! do dry deposition for this particle?
    logical                                   :: do_vtran           ! do sedimentation for this particle?
    integer                                   :: imiertn            ! mie velocity routine
    logical                                   :: is_sulfate         ! is this a sulfate particle?
    real(kind=f)                              :: dpc_threshold      ! convergence criteria for particle concentration [fraction]

    ! Assume success.
    rc = RC_OK

    ! Test out the Get method.
    if (carma%f_do_print) then
      call CARMAGROUP_Get(carma, igroup, rc, name=name, shortname=shortname, rmin=rmin, rmrat=rmrat, ishape=ishape, eshape=eshape, &
                        is_ice=is_ice, is_cloud=is_cloud, irhswell=irhswell, irhswcomp=irhswcomp, cnsttype=cnsttype, r=r, dr=dr, &
                        rmass=rmass, dm=dm, vol=vol, ifallrtn=ifallrtn, rmassmin=rmassmin, do_mie=do_mie, do_wetdep=do_wetdep, &
                        do_drydep=do_drydep, do_vtran=do_vtran, imiertn=imiertn)
      if (rc < 0) return

    
      write(carma%f_LUNOPRT,*) "    name          : ", trim(name)
      write(carma%f_LUNOPRT,*) "    shortname     : ", trim(shortname)
      write(carma%f_LUNOPRT,*) "    rmin          : ", rmin, " (cm)"
      write(carma%f_LUNOPRT,*) "    rmassmin      : ", rmassmin, " (g)"
      write(carma%f_LUNOPRT,*) "    rmrat         : ", rmrat
      write(carma%f_LUNOPRT,*) "    dpc_threshold : ", dpc_threshold

      select case(ishape)
        case (I_SPHERE)
          write(carma%f_LUNOPRT,*) "    ishape        :    spherical"
        case (I_HEXAGON)
          write(carma%f_LUNOPRT,*) "    ishape        :    hexagonal"
        case (I_CYLINDER)
          write(carma%f_LUNOPRT,*) "    ishape        :    cylindrical"
        case default
          write(carma%f_LUNOPRT,*) "    ishape        :    unknown, ", ishape
      end select

      write(carma%f_LUNOPRT,*) "    eshape        : ", eshape
      write(carma%f_LUNOPRT,*) "    is_ice        : ", is_ice
      write(carma%f_LUNOPRT,*) "    is_cloud      : ", is_cloud
      write(carma%f_LUNOPRT,*) "    is_sulfate    : ", is_sulfate
      
      write(carma%f_LUNOPRT,*) "    do_drydep     : ", do_drydep
      write(carma%f_LUNOPRT,*) "    do_mie        : ", do_mie
      write(carma%f_LUNOPRT,*) "    do_vtran      : ", do_vtran
      write(carma%f_LUNOPRT,*) "    do_wetdep     : ", do_wetdep
      
      select case(irhswell)
        case (0)
          write(carma%f_LUNOPRT,*) "    irhswell      :    none"
        case (I_FITZGERALD)
          write(carma%f_LUNOPRT,*) "    irhswell      :    Fitzgerald"
        case (I_GERBER)
          write(carma%f_LUNOPRT,*) "    irhswell      :    Gerber"
        case default
          write(carma%f_LUNOPRT,*) "    irhswell      :    unknown, ", irhswell
      end select

      select case(irhswcomp)
        case (0)
          write(carma%f_LUNOPRT,*) "    irhswcomp     :    none"

        case (I_SWF_NH42SO4)
          write(carma%f_LUNOPRT,*) "    irhswcomp     :    (NH4)2SO4 (Fitzgerald)"
        case (I_SWF_NH4NO3)
          write(carma%f_LUNOPRT,*) "    irhswcomp     :    NH4NO3 (Fitzgerald)"
        case (I_SWF_NANO3)
          write(carma%f_LUNOPRT,*) "    irhswcomp     :    NaNO3 (Fitzgerald)"
        case (I_SWF_NH4CL)
          write(carma%f_LUNOPRT,*) "    irhswcomp     :    NH4Cl (Fitzgerald)"
        case (I_SWF_CACL2)
          write(carma%f_LUNOPRT,*) "    irhswcomp     :    CaCl2 (Fitzgerald)"
        case (I_SWF_NABR)
          write(carma%f_LUNOPRT,*) "    irhswcomp     :    NaBr (Fitzgerald)"
        case (I_SWF_NACL)
          write(carma%f_LUNOPRT,*) "    irhswcomp     :    NaCl (Fitzgerald)"
        case (I_SWF_MGCL2)
          write(carma%f_LUNOPRT,*) "    irhswcomp     :    MgCl2 (Fitzgerald)"
        case (I_SWF_LICL)
          write(carma%f_LUNOPRT,*) "    irhswcomp     :    LiCl (Fitzgerald)"

        case (I_SWG_NH42SO4)
          write(carma%f_LUNOPRT,*) "    irhswcomp     :    (NH4)2SO4 (Gerber)"
        case (I_SWG_RURAL)
          write(carma%f_LUNOPRT,*) "    irhswcomp     :    Rural (Gerber)"
        case (I_SWG_SEA_SALT)
          write(carma%f_LUNOPRT,*) "    irhswcomp     :    Sea Salt (Gerber)"
        case (I_SWG_URBAN)
          write(carma%f_LUNOPRT,*) "    irhswcomp     :    Urban (Gerber)"

        case default
          write(carma%f_LUNOPRT,*) "    irhswell      :    unknown, ", irhswcomp
      end select
      
      select case(cnsttype)
        case (0)
          write(carma%f_LUNOPRT,*) "    cnsttype      :    none"
        case (I_CNSTTYPE_PROGNOSTIC)
          write(carma%f_LUNOPRT,*) "    cnsttype      :    prognostic"
         case (I_CNSTTYPE_DIAGNOSTIC)
          write(carma%f_LUNOPRT,*) "    cnsttype      :    diagnostic"
        case default
          write(carma%f_LUNOPRT,*) "    cnsttype      :    unknown, ", cnsttype
      end select

      select case(ifallrtn)
        case (I_FALLRTN_STD)
          write(carma%f_LUNOPRT,*) "    ifallrtn      :    standard"
        case (I_FALLRTN_STD_SHAPE)
          write(carma%f_LUNOPRT,*) "    ifallrtn      :    standard (shape)"
        case (I_FALLRTN_HEYMSFIELD2010)
          write(carma%f_LUNOPRT,*) "    ifallrtn      :    Heymsfield & Westbrook, 2010"
        case default
          write(carma%f_LUNOPRT,*) "    ifallrtn      :    unknown, ", ifallrtn
      end select

      select case(imiertn)
        case (I_MIERTN_TOON1981)
          write(carma%f_LUNOPRT,*) "    imiertn       :    Toon & Ackerman, 1981"
        case (I_MIERTN_BOHREN1983)
          write(carma%f_LUNOPRT,*) "    imiertn       :    Bohren & Huffman, 1983"
        case default
          write(carma%f_LUNOPRT,*) "    imiertn       :    unknown, ", imiertn
      end select
  
      write(carma%f_LUNOPRT,*)   
      write(carma%f_LUNOPRT,"('    ', a4, 5a12)") "bin",  "r",  "dr",  "rmass",  "dm",  "vol"
      write(carma%f_LUNOPRT,"('    ', a4, 5a12)") "",  "(cm)",  "(cm)",  "(g)",  "(g)",  "(cm3)"
     
      do i = 1, carma%f_NBIN
        write(carma%f_LUNOPRT, "('    ', i4,  5g12.3)") i, r(i), dr(i), rmass(i), dm(i), vol(i)
      end do
    end if
    
    return
  end subroutine CARMAGROUP_Print
end module
