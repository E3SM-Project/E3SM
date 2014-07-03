!! The CARMA state module contains the atmospheric data for use with the CARMA
!! module. This implementation has been customized to work within other model 
!! frameworks. CARMA adds a lot of extra state information (atmospheric
!! properties, fall velocities, coagulation kernels, growth kernels, ...) and
!! thus has a large memory footprint. Because only one column will be operated
!! upon at a time per thread, only one cstate object needs to be instantiated
!! at a time and each cstate object only represents one column. This keeps
!! the memory requirements of CARMA to a minimum.
!!
!! @version Feb-2009 
!! @author  Chuck Bardeen, Pete Colarco, Jamie Smith 
!
! NOTE: Documentation for this code can be generated automatically using f90doc,
! which is freely available from:
!   http://erikdemaine.org/software/f90doc/
! Comment lines with double comment characters are processed by f90doc, and there are
! some special characters added to the comments to control the documentation process.
! In addition to the special characters mentioned in the f990doc documentation, html
! formatting tags (e.g. <i></i>, <sup></sup>, ...) can also be added to the f90doc
! comments.
module carmastate_mod

  ! This module maps the parents models constants into the constants need by CARMA.
  ! NOTE: CARMA constants are in CGS units, while the parent models are typically in
  ! MKS units.
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod

  ! cstate explicitly declares all variables. 
  implicit none

  ! All cstate variables and procedures are private except those explicitly
  ! declared to be public.
  private
  
  ! Declare the public methods.
  public CARMASTATE_Create
  public CARMASTATE_CreateFromReference
  public CARMASTATE_Destroy
  public CARMASTATE_Get
  public CARMASTATE_GetBin
  public CARMASTATE_GetDetrain
  public CARMASTATE_GetGas
  public CARMASTATE_GetState
  public CARMASTATE_SetBin
  public CARMASTATE_SetDetrain
  public CARMASTATE_SetGas
  public CARMASTATE_SetState
  public CARMASTATE_Step
  
contains
  
  ! These are the methods that provide the interface between the parent model and
  ! the atmospheric state data of the CARMA microphysical model. There are many other
  ! methods that are not in this file that are used to implement the microphysical
  ! calculations needed by the CARMA model. These other methods are in effect private
  ! methods of the CARMA module, but are in individual files since that is the way that
  ! CARMA has traditionally been structured and where users may want to extend or
  ! replace code to affect the microphysics.

  !! Create the CARMASTATE object, which contains information about the
  !! atmospheric state. Internally, CARMA uses CGS units, but this interface uses
  !! MKS units which are more commonly used in parent models. The units and grid
  !! orientation depend on the grid type:
  !!
  !!  - igridh
  !!    - I_CART   : Cartesian coordinates, units in [m]
  !!    - I_LL     : Lat/Lon coordinates, units in [degrees]
  !!
  !!  - igridv
  !!    - I_CART   : Cartesian coordinates, units in [m], bottom at NZ=1
  !!    - I_SIG    : Sigma coordinates, unitless [P/P0], top at NZ=1
  !!    - I_HYBRID : Hybrid coordinates, unitless [~P/P0], top at NZ=1
  !!
  !! NOTE: The supplied CARMA object should already have been created, configured,
  !! and initialized.
  !!
  !! NOTE: The relative humidity is optional, but needs to be supplied if particles
  !! are subject to swelling based upon relative humidity. The specific humdity can
  !! can be specified instead. If both are specified, then the realtive humidity is
  !! used.
  !!
  !! @author Chuck Bardeen
  !! @version Feb-2009
  !! @see CARMA_Create
  !! @see CARMA_Initialize
  !! @see CARMASTATE_Destroy
  subroutine CARMASTATE_Create(cstate, carma_ptr, time, dtime, NZ, igridv, igridh,  &
      lat, lon, xc, dx, yc, dy, zc, zl, p, pl, t, rc, qh2o, relhum, told, radint)
    type(carmastate_type), intent(inout)    :: cstate      !! the carma state object
    type(carma_type), pointer               :: carma_ptr   !! (in) the carma object
    real(kind=f), intent(in)                :: time        !! the model time [s]
    real(kind=f), intent(in)                :: dtime       !! the timestep size [s]
    integer, intent(in)                     :: NZ          !! the number of vertical grid points
    integer, intent(in)                     :: igridv      !! vertical grid type
    integer, intent(in)                     :: igridh      !! horizontal grid type
    real(kind=f), intent(in)                :: lat         !! latitude at center [degrees north]
    real(kind=f), intent(in)                :: lon         !! longitude at center [degrees east]
    real(kind=f), intent(in)                :: xc(NZ)      !! x at center
    real(kind=f), intent(in)                :: dx(NZ)      !! ix width
    real(kind=f), intent(in)                :: yc(NZ)      !! y at center
    real(kind=f), intent(in)                :: dy(NZ)      !! y width
    real(kind=f), intent(in)                :: zc(NZ)      !! z at center
    real(kind=f), intent(in)                :: zl(NZ+1)    !! z at edge
    real(kind=f), intent(in)                :: p(NZ)       !! pressure at center [Pa]
    real(kind=f), intent(in)                :: pl(NZ+1)    !! presssure at edge [Pa]
    real(kind=f), intent(in)                :: t(NZ)       !! temperature at center [K]
    integer, intent(out)                    :: rc          !! return code, negative indicates failure
    real(kind=f), intent(in) , optional     :: qh2o(NZ)    !! specific humidity at center [mmr]
    real(kind=f), intent(in) , optional     :: relhum(NZ)  !! relative humidity at center [fraction]
    real(kind=f), intent(in) , optional     :: told(NZ)    !! previous temperature at center [K]
    real(kind=f), intent(in) , optional     :: radint(NZ,carma_ptr%f_NWAVE)  !! radiative intensity [W/m2/sr/cm]

    integer                                 :: iz
    real(kind=f)                            :: rvap
    real(kind=f)                            :: pvap_liq
    real(kind=f)                            :: pvap_ice
    real(kind=f)                            :: gc_cgs  

    ! Assume success.
    rc = RC_OK

    ! Save the defintion of the number of comonents involved in the microphysics.
    cstate%f_carma => carma_ptr

    ! Save the model timing.
    cstate%f_time       = time
    cstate%f_dtime_orig = dtime
    cstate%f_dtime      = dtime
    cstate%f_nretries   = 0
    
    ! Save the grid dimensions.
    cstate%f_NZ   = NZ
    cstate%f_NZP1 = NZ+1
    
    ! Save the grid definition.
    cstate%f_igridv = igridv
    cstate%f_igridh = igridh
    
    ! Store away the grid location information.
    cstate%f_lat  = lat
    cstate%f_lon  = lon
    
    ! Allocate all the dynamic variables related to state.
    call CARMASTATE_Allocate(cstate, rc)
    if (rc < 0) return
    
    cstate%f_xc(:)  = xc(:)
    cstate%f_dx(:)  = dx(:)
    cstate%f_yc(:)  = yc(:)
    cstate%f_dy(:)  = dy(:)        
    cstate%f_zc(:)  = zc(:)
    cstate%f_zl(:)  = zl(:)

    ! Store away the grid state, doing any necessary unit conversions from MKS to CGS.
    cstate%f_p(:)  = p(:)  * RPA2CGS    
    cstate%f_pl(:) = pl(:) * RPA2CGS    
    cstate%f_t(:)  = t(:)
    
    cstate%f_pcd(:,:,:)     = 0._f
    
    if (carma_ptr%f_do_substep) then
      if (present(told)) then
        cstate%f_told(:) = told
      else
        if (carma_ptr%f_do_print) write(carma_ptr%f_LUNOPRT,*) "CARMASTATE_Create: Error - Need to specify told when substepping."
        rc = RC_ERROR
        
        return
      end if
    end if 
    
    ! Calculate the metrics, ...
    ! if Cartesian coordinates were specifed, then the units need to be converted
    ! from MKS to CGS.
    if (cstate%f_igridh == I_CART) then
      cstate%f_xc = cstate%f_xc * RM2CGS
      cstate%f_dx = cstate%f_dx * RM2CGS
      cstate%f_yc = cstate%f_yc * RM2CGS
      cstate%f_dy = cstate%f_dy * RM2CGS
    end if
    
    if (cstate%f_igridv == I_CART) then
      cstate%f_zc = cstate%f_zc * RM2CGS
      cstate%f_zl = cstate%f_zl * RM2CGS
    end if
    
    ! Initialize the state of the atmosphere.
    call setupatm(carma_ptr, cstate, carma_ptr%f_do_fixedinit, rc)
    if (rc < 0) return
    
    ! Set the realtive humidity. If necessary, it will be calculated from
    ! the specific humidity.
    if (present(relhum)) then
      cstate%f_relhum(:) = relhum(:)
    else if (present(qh2o)) then
    
      ! Define gas constant for this gas
      rvap = RGAS/WTMOL_H2O

      ! Calculate relative humidity
      do iz = 1, NZ
        call vaporp_h2o_murphy2005(carma_ptr, cstate, iz, rc, pvap_liq, pvap_ice)
        if (rc < 0) return

        gc_cgs = qh2o(iz)*cstate%f_rhoa_wet(iz) / (cstate%f_zmet(iz)*cstate%f_xmet(iz)*cstate%f_ymet(iz))
        cstate%f_relhum(iz) = ( gc_cgs * rvap * t(iz)) / pvap_liq
      enddo
    end if
    
    ! Need for vertical transport.
    !
    ! NOTE: How should these be set? Optional parameters?
    if (carma_ptr%f_do_vtran) then
      cstate%f_ftoppart(:,:)  = 0._f
      cstate%f_fbotpart(:,:)  = 0._f
      cstate%f_pc_topbnd(:,:) = 0._f
      cstate%f_pc_botbnd(:,:) = 0._f
    end if
        
    ! Radiative intensity for particle heating.
    !
    ! W/m2/sr/cm -> erg/s/cm2/sr/cm
    if (carma_ptr%f_do_grow) then
      if (present(radint)) cstate%f_radint(:,:) = radint(:,:) * 1e7_f / 1e4_f
    end if
    
    return
  end subroutine CARMASTATE_Create


  !! Create the CARMASTATE object, which contains information about the
  !! atmospheric state.
  !! 
  !! This call is similar to CARMASTATE_Create, but differs in that all the
  !! initialization happens here based on the the fixed state information provided rather
  !! than occurring in CARMASTATE_Step.
  !!
  !! This call should be done before CARMASTATE_Create when do_fixedinit has been
  !! specified. The temperatures and pressures specified here should be the reference
  !! state used for all columns, not an actual column from the model.
  !!
  !! A water vapor profile is optional, but is used whenever  either qh2o (preferred)
  !! or relhum have been provided. If this is not provided, then initialization will
  !! be done on a dry profile. If particle swelling occurs, initialization will be
  !! done on the wet radius; however, most of the initialized values will not get
  !! recalculated as the wet radius changes.
  !!
  !! CARMASTATE_Create should still be called again after this call with the actual
  !! column of state information from the model. The initialization will be done once 
  !! from the reference state, but the microphysical calculations will be done on the
  !! model state. Multiple CARMASTATE_Create ... CARMASTATE_Step calls can be done
  !! before a CARMASTATE_Destroy. This reduces the amount of memory allocations and
  !! when used with do_fixedinit, reduces the amount of time spent initializing.
  !!
  !! @author Chuck Bardeen
  !! @version June-2010
  !! @see CARMA_Create
  !! @see CARMA_Initialize
  !! @see CARMASTATE_Destroy
  subroutine CARMASTATE_CreateFromReference(cstate, carma_ptr, time, dtime, NZ, igridv, igridh,  &
      lat, lon, xc, dx, yc, dy, zc, zl, p, pl, t, rc, qh2o, relhum, qh2so4)
    type(carmastate_type), intent(inout)    :: cstate      !! the carma state object
    type(carma_type), pointer               :: carma_ptr   !! (in) the carma object
    real(kind=f), intent(in)                :: time        !! the model time [s]
    real(kind=f), intent(in)                :: dtime       !! the timestep size [s]
    integer, intent(in)                     :: NZ          !! the number of vertical grid points
    integer, intent(in)                     :: igridv      !! vertical grid type
    integer, intent(in)                     :: igridh      !! horizontal grid type
    real(kind=f), intent(in)                :: lat         !! latitude at center [degrees north]
    real(kind=f), intent(in)                :: lon         !! longitude at center [degrees east]
    real(kind=f), intent(in)                :: xc(NZ)      !! x at center
    real(kind=f), intent(in)                :: dx(NZ)      !! ix width
    real(kind=f), intent(in)                :: yc(NZ)      !! y at center
    real(kind=f), intent(in)                :: dy(NZ)      !! y width
    real(kind=f), intent(in)                :: zc(NZ)      !! z at center
    real(kind=f), intent(in)                :: zl(NZ+1)    !! z at edge
    real(kind=f), intent(in)                :: p(NZ)       !! pressure at center [Pa]
    real(kind=f), intent(in)                :: pl(NZ+1)    !! presssure at edge [Pa]
    real(kind=f), intent(in)                :: t(NZ)       !! temperature at center [K]
    integer, intent(out)                    :: rc          !! return code, negative indicates failure
    real(kind=f), intent(in) , optional     :: qh2o(NZ)    !! specific humidity at center [mmr]
    real(kind=f), intent(in) , optional     :: relhum(NZ)  !! relative humidity at center [fraction]
    real(kind=f), intent(in) , optional     :: qh2so4(NZ)  !! H2SO4 mass mixing ratio at center [mmr]
    
    integer                                 :: iz
    integer                                 :: igas
    real(kind=f)                            :: rvap
    real(kind=f)                            :: pvap_liq
    real(kind=f)                            :: pvap_ice
    real(kind=f)                            :: gc_cgs  

    ! Assume success.
    rc = RC_OK

    ! Save the defintion of the number of comonents involved in the microphysics.
    cstate%f_carma => carma_ptr

    ! Save the model timing.
    cstate%f_time       = time
    cstate%f_dtime_orig = dtime
    cstate%f_dtime      = dtime
    cstate%f_nretries   = 0
    
    ! Save the grid dimensions.
    cstate%f_NZ   = NZ
    cstate%f_NZP1 = NZ+1
    
    ! Save the grid definition.
    cstate%f_igridv = igridv
    cstate%f_igridh = igridh
    
    ! Store away the grid location information.
    cstate%f_lat  = lat
    cstate%f_lon  = lon
    
    ! Allocate all the dynamic variables related to state.
    call CARMASTATE_Allocate(cstate, rc)
    if (rc < 0) return
    
    cstate%f_xc(:)  = xc(:)
    cstate%f_dx(:)  = dx(:)
    cstate%f_yc(:)  = yc(:)
    cstate%f_dy(:)  = dy(:)        
    cstate%f_zc(:)  = zc(:)
    cstate%f_zl(:)  = zl(:)

    ! Store away the grid state, doing any necessary unit conversions from MKS to CGS.
    cstate%f_p(:)  = p(:)  * RPA2CGS    
    cstate%f_pl(:) = pl(:) * RPA2CGS    
    cstate%f_t(:)  = t(:)
    
    cstate%f_pcd(:,:,:)     = 0._f
    
    ! Calculate the metrics, ...
    ! if Cartesian coordinates were specifed, then the units need to be converted
    ! from MKS to CGS.
    if (cstate%f_igridh == I_CART) then
      cstate%f_xc = cstate%f_xc * RM2CGS
      cstate%f_dx = cstate%f_dx * RM2CGS
      cstate%f_yc = cstate%f_yc * RM2CGS
      cstate%f_dy = cstate%f_dy * RM2CGS
    end if
    
    if (cstate%f_igridv == I_CART) then
      cstate%f_zc = cstate%f_zc * RM2CGS
      cstate%f_zl = cstate%f_zl * RM2CGS
    end if
    
    ! Initialize the state of the atmosphere.
    call setupatm(carma_ptr, cstate, .false., rc)
    if (rc < 0) return

    ! If the model uses a gas, then set the relative and
    ! specific humidities.
    if (carma_ptr%f_igash2o /= 0) then
    
      if (present(qh2o)) then
        cstate%f_gc(:, carma_ptr%f_igash2o) = qh2o(:) * cstate%f_rhoa_wet(:)
      
        ! Define gas constant for this gas
        rvap = RGAS/WTMOL_H2O
  
        ! Calculate relative humidity
        do iz = 1, NZ
          call vaporp_h2o_murphy2005(carma_ptr, cstate, iz, rc, pvap_liq, pvap_ice)
          if (rc < 0) return
  
          gc_cgs = qh2o(iz) * cstate%f_rhoa_wet(iz) / (cstate%f_zmet(iz)*cstate%f_xmet(iz)*cstate%f_ymet(iz))
          cstate%f_relhum(iz) = (gc_cgs * rvap * t(iz)) / pvap_liq
        enddo
        
      else if (present(relhum)) then
        cstate%f_relhum(:) = relhum
        
        ! Define gas constant for this gas
        rvap = RGAS/WTMOL_H2O

        ! Calculate specific humidity
        do iz = 1, NZ
          call vaporp_h2o_murphy2005(carma_ptr, cstate, iz, rc, pvap_liq, pvap_ice)
          if (rc < 0) return
  
          gc_cgs = (rvap * t(iz)) / (pvap_liq * relhum(iz))
          cstate%f_gc(iz, carma_ptr%f_igash2o) = gc_cgs * (cstate%f_zmet(iz)*cstate%f_xmet(iz)*cstate%f_ymet(iz)) / cstate%f_rhoa_wet(iz) 
        enddo
      end if      
    end if

    ! If the model uses sulfuric acid, then set that gas concentration.
    if (carma_ptr%f_igash2so4 /= 0) then
      if (present(qh2so4)) then
        cstate%f_gc(:, carma_ptr%f_igash2so4) = qh2so4(:) * cstate%f_rhoa_wet(:)
      end if
    end if

    ! Determine the gas supersaturations.
    do iz = 1, cstate%f_NZ
      do igas = 1, cstate%f_carma%f_NGAS
        call supersat(cstate%f_carma, cstate, iz, igas, rc)
        if (rc < 0) return
      end do
    end do

    ! Need for vertical transport.
    !
    ! NOTE: How should these be set? Optional parameters?
    if (carma_ptr%f_do_vtran) then
      cstate%f_ftoppart(:,:) = 0._f
      cstate%f_fbotpart(:,:) = 0._f
      cstate%f_pc_topbnd(:,:) = 0._f
      cstate%f_pc_botbnd(:,:) = 0._f
    end if
    
    
    ! Now do the initialization that is normally done in CARMASTATE_Step. However
    ! here it is done using the reference atmosphere.
    
    ! Determine the particle densities.
    call rhopart(cstate%f_carma, cstate, rc)
    if (rc < 0) return
    
    ! Save off the wet radius and wet density as reference values to be used
    ! later to scale process rates based upon changes to the wet radius and
    ! wet density when particle swelling is used.
    cstate%f_r_ref(:,:,:)    = cstate%f_r_wet(:,:,:)
    cstate%f_rhop_ref(:,:,:) = cstate%f_rhop_wet(:,:,:)

    ! If configured for fixed initialization, then we will lose some accuracy
    ! in the calculation of the fall velocities, growth kernels, ... and in return
    ! will gain a significant performance by not having to initialize as often.
  
    ! Initialize the vertical transport.
    if (cstate%f_carma%f_do_vtran .or. cstate%f_carma%f_do_coag .or. cstate%f_carma%f_do_grow) then
      call setupvf(cstate%f_carma, cstate, rc)
      
      if (cstate%f_carma%f_do_vdiff) then
        call setupbdif(cstate%f_carma, cstate, rc)
      end if
    end if

    ! Intialize the nucleation, growth and evaporation.      
    if (cstate%f_carma%f_do_grow)  then
      call setupgrow(cstate%f_carma, cstate, rc)
      if (rc < 0) return

      call setupgkern(cstate%f_carma, cstate, rc)
      if (rc < 0) return
      
       call setupnuc(cstate%f_carma, cstate, rc)
      if (rc < 0) return
    end if
    
    ! Initialize the coagulation.
    if (cstate%f_carma%f_do_coag) then
      call setupckern(cstate%f_carma, cstate, rc)
      if (rc < 0) return
    end if
    
    return
  end subroutine CARMASTATE_CreateFromReference


  subroutine CARMASTATE_Allocate(cstate, rc)
    type(carmastate_type), intent(inout)  :: cstate
    integer, intent(out)                  :: rc
    
    ! Local Variables
    integer                               :: ier
    integer                               :: NZ
    integer                               :: NZP1
    integer                               :: NGROUP
    integer                               :: NELEM
    integer                               :: NBIN
    integer                               :: NGAS
    integer                               :: NWAVE
    
    ! Assume success.
    rc = RC_OK

    ! Check to see if the arrays are already allocated. If so, just reuse the
    ! existing allocations.
    
    ! Allocate the variables needed for setupatm.
    if (.not. (allocated(cstate%f_xmet))) then
    
      NZ      = cstate%f_NZ
      NZP1    = cstate%f_NZP1
      NGROUP  = cstate%f_carma%f_NGROUP
      NELEM   = cstate%f_carma%f_NELEM
      NBIN    = cstate%f_carma%f_NBIN
      NGAS    = cstate%f_carma%f_NGAS
      NWAVE   = cstate%f_carma%f_NWAVE
    
      allocate( &
        cstate%f_xmet(NZ), &
        cstate%f_ymet(NZ), &
        cstate%f_zmet(NZ), &
        cstate%f_zmetl(NZP1), &
        cstate%f_xc(NZ), &
        cstate%f_yc(NZ), &
        cstate%f_zc(NZ), &
        cstate%f_dx(NZ), &
        cstate%f_dy(NZ), &
        cstate%f_dz(NZ), &
        cstate%f_zl(NZP1), &
        cstate%f_pc(NZ,NBIN,NELEM), &
        cstate%f_pcd(NZ,NBIN,NELEM), &
        cstate%f_pc_surf(NBIN,NELEM), &
        cstate%f_sedimentationflux(NBIN,NELEM), &
        cstate%f_gc(NZ,NGAS), &
        cstate%f_cldfrc(NZ), &
        cstate%f_rhcrit(NZ), &
        cstate%f_rhop(NZ,NBIN,NGROUP), &
        cstate%f_r_wet(NZ,NBIN,NGROUP), &
        cstate%f_rlow_wet(NZ,NBIN,NGROUP), &
        cstate%f_rup_wet(NZ,NBIN,NGROUP), &
        cstate%f_rhop_wet(NZ,NBIN,NGROUP), &
        cstate%f_r_ref(NZ,NBIN,NGROUP), &
        cstate%f_rhop_ref(NZ,NBIN,NGROUP), &
        cstate%f_rhoa(NZ), &
        cstate%f_rhoa_wet(NZ), &
        cstate%f_t(NZ), &
        cstate%f_p(NZ), &
        cstate%f_pl(NZP1), &
        cstate%f_relhum(NZ), &
        cstate%f_wtpct(NZ), &
        cstate%f_rmu(NZ), &
        cstate%f_thcond(NZ), &
        cstate%f_thcondnc(NZ,NBIN,NGROUP), &
        cstate%f_dpc_sed(NBIN,NELEM), &
        cstate%f_pconmax(NZ,NGROUP), &
        cstate%f_pcl(NZ,NBIN,NELEM), &
        stat=ier)
      if (ier /= 0) then
        if (cstate%f_carma%f_do_print) write(cstate%f_carma%f_LUNOPRT, *) "CARMASTATE_Allocate::ERROR allocating atmosphere arrays, status=", ier
        rc = RC_ERROR
        return
      end if
      
      cstate%f_relhum(:)      = 0._f
      cstate%f_pc(:,:,:)      = 0._f
      cstate%f_pcd(:,:,:)     = 0._f
      cstate%f_pc_surf(:,:)   = 0._f
      cstate%f_sedimentationflux(:,:)   = 0._f
      cstate%f_cldfrc(:)      = 1._f
      cstate%f_rhcrit(:)      = 1._f
      
      ! Allocate the last fields if they are needed for substepping.
      if (cstate%f_carma%f_do_substep) then
        allocate( &
          cstate%f_gcl(NZ,NGAS), &
          cstate%f_d_gc(NZ,NGAS), &
          cstate%f_told(NZ), &
          cstate%f_d_t(NZ), &
          cstate%f_zsubsteps(NZ), &
          stat=ier)
        if (ier /= 0) then
          if (cstate%f_carma%f_do_print) write(cstate%f_carma%f_LUNOPRT, *) "CARMASTATE_Allocate::ERROR allocating stepping arrays, status=", ier
          rc = RC_ERROR
          return
        endif
      
        ! Initialize
        cstate%f_gcl(:,:)     = 0._f
        cstate%f_d_gc(:,:)    = 0._f
        cstate%f_told(:)      = 0._f
        cstate%f_d_t(:)       = 0._f
        cstate%f_zsubsteps(:) = 0._f

        ! When substepping is enabled, we want to initialize these statistics once for
        ! the life of the object.
        cstate%f_max_nsubstep = 0
        cstate%f_max_nretry   = 0._f
        cstate%f_nstep        = 0._f
        cstate%f_nsubstep     = 0
        cstate%f_nretry       = 0._f
      endif

      
      ! Allocate the variables needed for setupvf.
      !
      ! NOTE: Coagulation and dry deposition also need bpm, vf and re.
      if (cstate%f_carma%f_do_vtran .or. cstate%f_carma%f_do_coag .or. cstate%f_carma%f_do_grow .or. cstate%f_carma%f_do_drydep) then
        allocate( &
          cstate%f_bpm(NZ,NBIN,NGROUP), &
          cstate%f_vf(NZP1,NBIN,NGROUP), &
          cstate%f_re(NZ,NBIN,NGROUP), &
          cstate%f_dkz(NZP1,NBIN,NGROUP), &
          cstate%f_ftoppart(NBIN,NELEM), &
          cstate%f_fbotpart(NBIN,NELEM), &
          cstate%f_pc_topbnd(NBIN,NELEM), &
          cstate%f_pc_botbnd(NBIN,NELEM), &
          cstate%f_vd(NBIN, NGROUP), &
          stat=ier)
        if (ier /= 0) then
          if (cstate%f_carma%f_do_print) write(cstate%f_carma%f_LUNOPRT, *) "CARMASTATE_Allocate::ERROR allocating vertical transport arrays, status=", ier
          rc = RC_ERROR
          return
        endif

        ! Initialize
        cstate%f_bpm(:,:,:) = 0._f
        cstate%f_vf(:,:,:) = 0._f
        cstate%f_re(:,:,:) = 0._f
        cstate%f_dkz(:,:,:) = 0._f
        cstate%f_ftoppart(:,:) = 0._f
        cstate%f_fbotpart(:,:) = 0._f
        cstate%f_pc_topbnd(:,:) = 0._f
        cstate%f_pc_botbnd(:,:) = 0._f
        cstate%f_vd(:, :) = 0._f
      end if
      
      
      
      if (cstate%f_carma%f_NGAS > 0) then
        allocate( &
          cstate%f_pvapl(NZ,NGAS), &
          cstate%f_pvapi(NZ,NGAS), &
          cstate%f_supsatl(NZ,NGAS), &
          cstate%f_supsati(NZ,NGAS), &
          cstate%f_supsatlold(NZ,NGAS), &
          cstate%f_supsatiold(NZ,NGAS), &
          stat=ier)
        if (ier /= 0) then
          if (cstate%f_carma%f_do_print) write(cstate%f_carma%f_LUNOPRT, *) "CARMASTATE_Allocate::ERROR allocating gas arrays, status=", ier
          rc = RC_ERROR
          return
        endif
      end if
 
      
      if (cstate%f_carma%f_do_grow) then
        allocate( &
          cstate%f_diffus(NZ,NGAS), &
          cstate%f_rlhe(NZ,NGAS), &
          cstate%f_rlhm(NZ,NGAS), &
          cstate%f_surfctwa(NZ), &
          cstate%f_surfctiw(NZ), &
          cstate%f_surfctia(NZ), &
          cstate%f_akelvin(NZ,NGAS), &
          cstate%f_akelvini(NZ,NGAS), &
          cstate%f_ft(NZ,NBIN,NGROUP), &
          cstate%f_gro(NZ,NBIN,NGROUP),  &
          cstate%f_gro1(NZ,NBIN,NGROUP),  &
          cstate%f_gro2(NZ,NGROUP),  &
          cstate%f_scrit(NZ,NBIN,NGROUP), &
          cstate%f_rnuclg(NBIN,NGROUP,NGROUP),&
          cstate%f_rhompe(NBIN,NELEM), &
          cstate%f_rnucpe(NBIN,NELEM), &
          cstate%f_pc_nucl(NZ,NBIN,NELEM), &
          cstate%f_growpe(NBIN,NELEM), &
          cstate%f_evappe(NBIN,NELEM), &
          cstate%f_evcore(NELEM), &
          cstate%f_growlg(NBIN,NGROUP), &
          cstate%f_evaplg(NBIN,NGROUP), &
          cstate%f_gasprod(NGAS), &
          cstate%f_rlheat(NZ), &
          cstate%f_radint(NZ,NWAVE), &
          cstate%f_partheat(NZ), &
          cstate%f_dtpart(NZ,NBIN,NGROUP), &
          cstate%f_cmf(NBIN,NGROUP), &
          cstate%f_totevap(NBIN,NGROUP), &
          stat=ier)
        if (ier /= 0) then
          if (cstate%f_carma%f_do_print) write(cstate%f_carma%f_LUNOPRT, *) "CARMASTATE_Allocate::ERROR allocating growth arrays, status=", ier
          rc = RC_ERROR
          return
        endif
        
        cstate%f_radint(:,:) = 0._f
      end if
      
      if (cstate%f_carma%f_do_coag) then
        allocate( &
          cstate%f_coaglg(NZ,NBIN,NGROUP), &
          cstate%f_coagpe(NZ,NBIN,NELEM), &
          cstate%f_ckernel(NZ,NBIN,NBIN,NGROUP,NGROUP), &
          stat = ier)
        if (ier /= 0) then
          if (cstate%f_carma%f_do_print) write(cstate%f_carma%f_LUNOPRT, *) "CARMASTATE_Allocate::ERROR allocating coag arrays, status=", ier
          rc = RC_ERROR
          return
        end if

        ! Initialize
        cstate%f_coaglg(:,:,:) = 0._f
        cstate%f_coagpe(:,:,:) = 0._f
        cstate%f_ckernel(:,:,:,:,:) = 0._f
      end if
    end if
    
    return
  end subroutine CARMASTATE_Allocate
    

  !! The routine should be called when the carma state object is no longer needed.
  !! It deallocates any memory allocations made by CARMA during CARMASTATE_Create(), 
  !! and failure to call this routine could result in memory leaks.
  !!
  !! @author Chuck Bardeen
  !! @version Feb-2009
  !! @see CARMASTATE_Create
  subroutine CARMASTATE_Destroy(cstate, rc)
    type(carmastate_type), intent(inout)    :: cstate
    integer, intent(out)                    :: rc
    
    ! Local variables
    integer   :: ier
    
    ! Assume success.
    rc = RC_OK

    ! Check to see if the arrays are already allocated. If so, deallocate them.

    ! Allocate the variables needed for setupatm.
    if (allocated(cstate%f_xmet)) then
    
      deallocate( &
        cstate%f_xmet, &
        cstate%f_ymet, &
        cstate%f_zmet, &
        cstate%f_zmetl, &
        cstate%f_xc, &
        cstate%f_yc, &
        cstate%f_zc, &
        cstate%f_dx, &
        cstate%f_dy, &
        cstate%f_dz, &
        cstate%f_zl, &
        cstate%f_pc, &
        cstate%f_pcd, &
        cstate%f_pc_surf, &
        cstate%f_sedimentationflux, &
        cstate%f_gc, &
        cstate%f_cldfrc, &
        cstate%f_rhcrit, &
        cstate%f_rhop, &
        cstate%f_r_wet, &
        cstate%f_rlow_wet, &
        cstate%f_rup_wet, &
        cstate%f_rhop_wet, &
        cstate%f_r_ref, &
        cstate%f_rhop_ref, &
        cstate%f_rhoa, &
        cstate%f_rhoa_wet, &
        cstate%f_t, &
        cstate%f_p, &
        cstate%f_pl, &
        cstate%f_relhum, &
        cstate%f_wtpct, &
        cstate%f_rmu, &
        cstate%f_thcond, &
        cstate%f_thcondnc, &
        cstate%f_dpc_sed, &
        cstate%f_pconmax, &
        cstate%f_pcl, &
        stat=ier)
      if (ier /= 0) then
        if (cstate%f_carma%f_do_print) write(cstate%f_carma%f_LUNOPRT, *) "CARMASTATE_Destroy::ERROR deallocating atmosphere arrays, status=", ier
        rc = RC_ERROR
        return
      end if
      
      ! Allocate the last fields if they are needed for substepping stepping.
      if (allocated(cstate%f_gcl)) then
        deallocate( &
          cstate%f_gcl, &
          cstate%f_d_gc, &
          cstate%f_told, &
          cstate%f_d_t, &
          cstate%f_zsubsteps, &
          stat=ier)
        if (ier /= 0) then
          if (cstate%f_carma%f_do_print) write(cstate%f_carma%f_LUNOPRT, *) "CARMASTATE_Destroy::ERROR deallocating stepping arrays, status=", ier
          rc = RC_ERROR
          return
        endif
      endif
      
      ! Allocate the variables needed for setupvf.
      !
      ! NOTE: Coagulation also needs bpm, vf and re.
      if (allocated(cstate%f_bpm)) then
        deallocate( &
          cstate%f_bpm, &
          cstate%f_vf, &
          cstate%f_re, &
          cstate%f_dkz, &
          cstate%f_ftoppart, &
          cstate%f_fbotpart, &
          cstate%f_pc_topbnd, &
          cstate%f_pc_botbnd, &
          cstate%f_vd, &
          stat=ier)
        if (ier /= 0) then
          if (cstate%f_carma%f_do_print) write(cstate%f_carma%f_LUNOPRT, *) "CARMASTATE_Destroy::ERROR deallocating vertical transport arrays, status=", ier
          rc = RC_ERROR
          return
        endif
      end if
      
      if (allocated(cstate%f_diffus)) then
        deallocate( &
          cstate%f_diffus, &
          cstate%f_rlhe, &
          cstate%f_rlhm, &
          cstate%f_surfctwa, &
          cstate%f_surfctiw, &
          cstate%f_surfctia, &
          cstate%f_akelvin, &
          cstate%f_akelvini, &
          cstate%f_ft, &
          cstate%f_gro, &
          cstate%f_gro1, &
          cstate%f_gro2, &
          cstate%f_scrit, &
          cstate%f_rnuclg,&
          cstate%f_rnucpe, &
          cstate%f_rhompe, &
          cstate%f_pc_nucl, &
          cstate%f_growpe, &
          cstate%f_evappe, &
          cstate%f_evcore, &
          cstate%f_growlg, &
          cstate%f_evaplg, &
          cstate%f_gasprod, &
          cstate%f_rlheat, &
          cstate%f_radint, &
          cstate%f_partheat, &
          cstate%f_dtpart, &
          cstate%f_cmf, &
          cstate%f_totevap, &
          stat=ier)
        if (ier /= 0) then
          if (cstate%f_carma%f_do_print) write(cstate%f_carma%f_LUNOPRT, *) "CARMASTATE_Destroy::ERROR deallocating growth arrays, status=", ier
          rc = RC_ERROR
          return
        endif
      end if
      
      if (allocated(cstate%f_pvapl)) then
        deallocate( &
          cstate%f_pvapl, &
          cstate%f_pvapi, &
          cstate%f_supsatl, &
          cstate%f_supsati, &
          cstate%f_supsatlold, &
          cstate%f_supsatiold, &
          stat=ier)
        if (ier /= 0) then
          if (cstate%f_carma%f_do_print) write(cstate%f_carma%f_LUNOPRT, *) "CARMASTATE_Destroy::ERROR deallocating gas arrays, status=", ier
          rc = RC_ERROR
          return
        endif
      end if
      
      if (allocated(cstate%f_coaglg)) then
        deallocate( &
          cstate%f_coaglg, &
          cstate%f_coagpe, &
          cstate%f_ckernel, &
          stat = ier)
        if (ier /= 0) then
          if (cstate%f_carma%f_do_print) write(cstate%f_carma%f_LUNOPRT, *) "CARMASTATE_Destroy::ERROR deallocating coag arrays, status=", ier
          rc = RC_ERROR
          return
        end if
      end if
    end if
    
    return
  end subroutine CARMASTATE_Destroy


  !! The routine performs the main CARMA processing for one timestep of
  !! the parent model. The state variables should have all been set before
  !! calling CARMASTATE_Step(). When this routine returns, the state will
  !! have been updated to reflect the changes from the CARMA microphysics.
  !! If tendencies are desired, then the difference between the final and
  !! initial state will need to be computed by the caller.
  !!
  !! NIOTE: xxxfv, xxxram and xxxfrac need to be specified for dry deposition.
  !!
  !! @author Chuck Bardeen
  !! @version Feb-2009
  subroutine CARMASTATE_Step(cstate, rc, cldfrc, rhcrit, lndfv, ocnfv, icefv, lndram, ocnram, iceram, lndfrac, ocnfrac, icefrac)
    type(carmastate_type), intent(inout)  :: cstate
    integer, intent(out)                  :: rc
    real(kind=f), intent(in), optional    :: cldfrc(cstate%f_NZ)  !! cloud fraction [fraction]
    real(kind=f), intent(in), optional    :: rhcrit(cstate%f_NZ)  !! relative humidity for onset of liquid clouds [fraction]
    real(kind=f), intent(in), optional    :: lndfv                !! the surface friction velocity over land  [m/s]
    real(kind=f), intent(in), optional    :: ocnfv                !! the surface friction velocity over ocean  [m/s]
    real(kind=f), intent(in), optional    :: icefv                !! the surface friction velocity over ice  [m/s]
    real(kind=f), intent(in), optional    :: lndram               !! the aerodynamic resistance over land [s/m]
    real(kind=f), intent(in), optional    :: ocnram               !! the aerodynamic resistance over ocean [s/m]
    real(kind=f), intent(in), optional    :: iceram               !! the aerodynamic resistance over ice [s/m]
    real(kind=f), intent(in), optional    :: lndfrac              !! land fraction
    real(kind=f), intent(in), optional    :: ocnfrac              !! ocn fraction
    real(kind=f), intent(in), optional    :: icefrac              !! ice fraction

    
    integer                               :: iz     ! vertical index
    integer                               :: igas   ! gas index
    integer                               :: ielem
    integer                               :: ibin
    integer                               :: igroup
    logical                               :: swelling   ! Do any groups undergo partcile swelling?
    integer                               :: i1, i2, j1, j2
  
    ! Assume success.
    rc = RC_OK
    
    ! Store the cloud fraction if specified
    cstate%f_cldfrc(:) = 1._f
    cstate%f_rhcrit(:) = 1._f
    
    if (present(cldfrc)) cstate%f_cldfrc(:) = cldfrc(:)
    if (present(rhcrit)) cstate%f_rhcrit(:) = rhcrit(:)
    
    ! Determine the gas supersaturations.
    do iz = 1, cstate%f_NZ
      do igas = 1, cstate%f_carma%f_NGAS
        call supersat(cstate%f_carma, cstate, iz, igas, rc)
        if (rc < 0) return
      end do
    end do

    ! Determine the particle densities.
    call rhopart(cstate%f_carma, cstate, rc)
    if (rc < 0) return
    

    ! We have to hold off initialization until now, because the particle density
    ! (rhop) can not be determined until the particle masses are known (i.e. after
    ! CARMASTATE_SetBin), because rhop is used to determine the fall velocity.
    !
    ! NOTE: If configured for fixed initialization, then we will lose some accuracy
    ! in the calculation of the fall velocities, growth kernels, ... and in return
    ! will gain a significant performance by not having to initialize as often.
    if (.not. cstate%f_carma%f_do_fixedinit) then
    
      ! Initialize the vertical transport.
      if (cstate%f_carma%f_do_vtran .or. cstate%f_carma%f_do_coag .or. cstate%f_carma%f_do_grow) then
        call setupvf(cstate%f_carma, cstate, rc)

        if (cstate%f_carma%f_do_vdiff) then
          call setupbdif(cstate%f_carma, cstate, rc)
        end if
      end if
      
      ! intialize the dry deposition
      if (cstate%f_carma%f_do_drydep) then
        if (present(lndfv) .and. present(lndram) .and. present(lndfrac) .and. &
            present(ocnfv) .and. present(ocnram) .and. present(ocnfrac) .and. &
            present(icefv) .and. present(iceram) .and. present(icefrac)) then
        
          ! NOTE: Need to convert surfric and ram from mks to cgs units.
          call setupvdry(cstate%f_carma, cstate, &
            lndfv * 100._f, ocnfv * 100._f, icefv * 100._f, &
            lndram / 100._f, ocnram / 100._f, iceram / 100._f, &
            lndfrac, ocnfrac, icefrac, rc)
          if (rc < RC_OK) return
        else
          write(cstate%f_carma%f_LUNOPRT, *) "CARMASTATE_Step: do_drydep requires that the optional inputs xxxfv, xxxram and xxxfrac be provided."
          rc = RC_ERROR
          return
        end if
      end if
       
      ! Intialize the nucleation, growth and evaporation.      
      if (cstate%f_carma%f_do_grow)  then
        call setupgrow(cstate%f_carma, cstate, rc)
        if (rc < RC_OK) return
  
        call setupgkern(cstate%f_carma, cstate, rc)
        if (rc < RC_OK) return
        
         call setupnuc(cstate%f_carma, cstate, rc)
        if (rc < RC_OK) return
      end if
      
      ! Initialize the coagulation.
      if (cstate%f_carma%f_do_coag) then
        call setupckern(cstate%f_carma, cstate, rc)
        if (rc < RC_OK) return
      end if
    end if
    
    ! Calculate the impact of microphysics upon the state.
    call step(cstate%f_carma, cstate, rc)

    return
  end subroutine CARMASTATE_Step


  ! Query, Control and State I/O

  !! Gets the mass mixing ratio for the gas (igas). After a call to CARMA_Step(),
  !! the new mass mixing ratio of the gas can be retrieved.
  !!
  !! @author Chuck Bardeen
  !! @version Feb-2009
  !! @see CARMA_AddGas
  !! @see CARMA_GetGas
  !! @see CARMA_Step 
  !! @see CARMASTATE_SetGas
  subroutine CARMASTATE_Get(cstate, rc, max_nsubstep, max_nretry, nstep, nsubstep, nretry, zsubsteps, lat, lon)
    type(carmastate_type), intent(in)     :: cstate            !! the carma state object
    integer, intent(out)                  :: rc                !! return code, negative indicates failure
    integer, optional, intent(out)        :: max_nsubstep      !! maximum number of substeps in a step
    real(kind=f), optional, intent(out)   :: max_nretry        !! maximum number of retries in a step
    real(kind=f), optional, intent(out)   :: nstep             !! total number of steps taken
    integer, optional, intent(out)        :: nsubstep          !! total number of substeps taken
    real(kind=f), optional, intent(out)   :: nretry            !! total number of retries taken
    real(kind=f), optional, intent(out)   :: zsubsteps(cstate%f_NZ) !! number of substeps taken per vertical grid point
    real(kind=f), optional, intent(out)   :: lat               !! grid center latitude [deg]
    real(kind=f), optional, intent(out)   :: lon               !! grid center longitude [deg]
    
    ! Assume success.
    rc = RC_OK

    if (present(max_nsubstep)) max_nsubstep = cstate%f_max_nsubstep
    if (present(max_nretry))   max_nretry   = cstate%f_max_nretry
    if (present(nstep))        nstep        = cstate%f_nstep
    if (present(nsubstep))     nsubstep     = cstate%f_nsubstep
    if (present(nretry))       nretry       = cstate%f_nretry
    if (present(zsubsteps))    zsubsteps    = cstate%f_zsubsteps
    if (present(lat))          lat          = cstate%f_lat
    if (present(lon))          lon          = cstate%f_lon
    
    return
  end subroutine CARMASTATE_Get
  

  !! Gets the mass of the bins (ibin) for each particle element (ielem). After the
  !! CARMA_Step() call, new particle concentrations are determined. The number density
  !! and the nucleation rate are only calculated if the element is the number density
  !! element for the group.
  !!
  !! @author Chuck Bardeen
  !! @version Feb-2009
  !! @see CARMA_AddElement
  !! @see CARMA_AddGroup
  !! @see CARMA_Step 
  !! @see CARMASTATE_SetBin
  subroutine CARMASTATE_GetBin(cstate, ielem, ibin, mmr, rc, &
                               nmr, numberDensity, nucleationRate, r_wet, rhop_wet, &
                               surface, sedimentationflux, vf, vd, dtpart)
    type(carmastate_type), intent(in)     :: cstate         !! the carma state object
    integer, intent(in)                   :: ielem          !! the element index
    integer, intent(in)                   :: ibin           !! the bin index
    real(kind=f), intent(out)             :: mmr(cstate%f_NZ) !! the bin mass mixing ratio [kg/kg]
    integer, intent(out)                  :: rc             !! return code negative indicates failure
    real(kind=f), optional, intent(out)   :: nmr(cstate%f_NZ) !! number mixing ratio [#/kg]
    real(kind=f), optional, intent(out)   :: numberDensity(cstate%f_NZ)  !! number density [#/cm3]
    real(kind=f), optional, intent(out)   :: nucleationRate(cstate%f_NZ) !! nucleation rate [1/cm3/s]
    real(kind=f), optional, intent(out)   :: r_wet(cstate%f_NZ)          !! wet particle radius [cm]
    real(kind=f), optional, intent(out)   :: rhop_wet(cstate%f_NZ)       !! wet particle density [g/cm3]
    real(kind=f), optional, intent(out)   :: surface        !! particle mass on the surface [kg/m2]
    real(kind=f), optional, intent(out)   :: sedimentationflux         !! particle sedimentation mass flux to surface [kg/m2/s]
    real(kind=f), optional, intent(out)   :: vf(cstate%f_NZ+1) !! fall velocity [cm/s]
    real(kind=f), optional, intent(out)   :: vd             !! deposition velocity [cm/s]
    real(kind=f), optional, intent(out)   :: dtpart(cstate%f_NZ) !! delta particle temperature [K]
    
    integer                               :: ienconc        !! index of element that is the particle concentration for the group
    integer                               :: igroup         ! Group containing this bin

    ! Assume success.
    rc = RC_OK
    
    ! Determine the particle group for the bin.    
    igroup = cstate%f_carma%f_element(ielem)%f_igroup

    ! Make sure there are enough elements allocated.
    if (ielem > cstate%f_carma%f_NELEM) then
      if (cstate%f_carma%f_do_print) write(cstate%f_carma%f_LUNOPRT, *) "CARMASTATE_SetBin:: ERROR - The specifed element (", &
        ielem, ") is larger than the number of elements (", cstate%f_carma%f_NELEM, ")."
      rc = RC_ERROR
      return
    end if
    
    ! Make sure there are enough bins allocated.
    if (ibin > cstate%f_carma%f_NBIN) then
      if (cstate%f_carma%f_do_print) write(cstate%f_carma%f_LUNOPRT, *) "CARMA_SetBin:: ERROR - The specifed bin (", &
        ibin, ") is larger than the number of bins (", cstate%f_carma%f_NBIN, ")."
      rc = RC_ERROR
      return
    end if

    
    ! Use the specified mass mixing ratio and the air density to determine the mass
    ! of the particles in g/x/y/z.
    mmr(:) = cstate%f_pc(:, ibin, ielem) / cstate%f_rhoa_wet(:)


    ! Handle the special cases for different types of elements ...
    if ((cstate%f_carma%f_element(ielem)%f_itype == I_INVOLATILE) .or. (cstate%f_carma%f_element(ielem)%f_itype == I_VOLATILE)) then
      mmr(:) = mmr(:) * cstate%f_carma%f_group(igroup)%f_rmass(ibin)
    else if (cstate%f_carma%f_element(ielem)%f_itype == I_CORE2MOM) then
      mmr(:) = mmr(:) / cstate%f_carma%f_group(igroup)%f_rmass(ibin)
    end if
    
    ! If the number of particles in the group is less than the minimum value represented
    ! by CARMA, then return and mmr of 0.0 for all elements.
    ienconc = cstate%f_carma%f_group(igroup)%f_ienconc
!    where (cstate%f_pc(:, ibin, ienconc) <= SMALL_PC) mmr(:) = 0.0_f


    ! Do they also want the mass concentration of particles at the surface?
    if (present(surface)) then
      
      ! Convert from g/cm2 to kg/m2
      surface = cstate%f_pc_surf(ibin, ielem) * 1e4_f / 1e3_f

      ! Handle the special cases for different types of elements ...
      if ((cstate%f_carma%f_element(ielem)%f_itype == I_INVOLATILE) .or. (cstate%f_carma%f_element(ielem)%f_itype == I_VOLATILE)) then
        surface = surface * cstate%f_carma%f_group(igroup)%f_rmass(ibin)
      else if (cstate%f_carma%f_element(ielem)%f_itype == I_CORE2MOM) then
        surface = surface / cstate%f_carma%f_group(igroup)%f_rmass(ibin)
      end if
    end if
    
    ! Do they also want the mass flux of particles that sediment to the surface?
    if (present(sedimentationflux)) then
      
      ! Convert from g/cm2 to kg/m2
      sedimentationflux = cstate%f_sedimentationflux(ibin, ielem) * 1e4_f / 1e3_f

      ! Handle the special cases for different types of elements ...
      if ((cstate%f_carma%f_element(ielem)%f_itype == I_INVOLATILE) .or. (cstate%f_carma%f_element(ielem)%f_itype == I_VOLATILE)) then
        sedimentationflux = sedimentationflux * cstate%f_carma%f_group(igroup)%f_rmass(ibin)
      else if (cstate%f_carma%f_element(ielem)%f_itype == I_CORE2MOM) then
        sedimentationflux = sedimentationflux / cstate%f_carma%f_group(igroup)%f_rmass(ibin)
      end if
    end if
    
    ! If this is the partcile # element, then determine some other statistics.
    if (ienconc == ielem) then
      if (present(nmr))           nmr(:)             = (cstate%f_pc(:, ibin, ielem) / cstate%f_rhoa_wet(:)) * 1000._f
      if (present(numberDensity)) numberDensity(:)   = cstate%f_pc(:, ibin, ielem) / (cstate%f_xmet(:)*cstate%f_ymet(:)*cstate%f_zmet(:))
      if (present(r_wet))         r_wet(:)           = cstate%f_r_wet(:, ibin, igroup)
      if (present(rhop_wet))      rhop_wet(:)        = cstate%f_rhop_wet(:, ibin, igroup)

      if (cstate%f_carma%f_do_vtran) then
        if (present(vf))            vf(:)              = cstate%f_vf(:, ibin, igroup) / cstate%f_zmetl(:)
      else
        if (present(vf))            vf(:)              = CAM_FILL
      end if
      
      if (cstate%f_carma%f_do_drydep) then
        if (present(vd)) then
          if (cstate%f_igridv .eq. I_CART) then
            vd                 = cstate%f_vd(ibin, igroup) / cstate%f_zmetl(1)
          else
            vd                 = cstate%f_vd(ibin, igroup) / cstate%f_zmetl(cstate%f_NZP1)
          end if
        end if
      else 
        if (present(vd))          vd                 = CAM_FILL
      end if

      if (cstate%f_carma%f_do_grow) then
        if (present(nucleationRate)) nucleationRate(:) = cstate%f_pc_nucl(:, ibin, ielem) / (cstate%f_xmet(:)*cstate%f_ymet(:)*cstate%f_zmet(:)) / cstate%f_dtime
      else
        if (present(nucleationRate)) nucleationRate(:) = CAM_FILL
      end if

      if (cstate%f_carma%f_do_pheat) then
        if (present(dtpart))        dtpart(:)          = cstate%f_dtpart(:, ibin, igroup)
      else
        if (present(dtpart))        dtpart(:)          = CAM_FILL
      end if
    else
      if (present(nmr))            nmr(:)             = CAM_FILL
      if (present(numberDensity))  numberDensity(:)   = CAM_FILL
      if (present(nucleationRate)) nucleationRate(:)  = CAM_FILL
      if (present(r_wet))          r_wet(:)           = CAM_FILL
      if (present(rhop_wet))       rhop_wet(:)        = CAM_FILL
      if (present(dtpart))         dtpart(:)          = CAM_FILL
      if (present(vf))             vf(:)              = CAM_FILL
      if (present(vd))             vd                 = CAM_FILL
    end if
   
    return
  end subroutine CARMASTATE_GetBin
  

  !! Gets the mass of the detrained condensate for the bins (ibin) for each particle
  !! element (ielem) in the grid.
  !!
  !!
  !! @author Chuck Bardeen
  !! @version Feb-2009
  !! @see CARMA_AddElement
  !! @see CARMA_AddGroup
  !! @see CARMA_Step 
  !! @see CARMASTATE_SetDetrain
  subroutine CARMASTATE_GetDetrain(cstate, ielem, ibin, mmr, rc, nmr, numberDensity, r_wet, rhop_wet)
    type(carmastate_type), intent(in)     :: cstate         !! the carma state object
    integer, intent(in)                   :: ielem          !! the element index
    integer, intent(in)                   :: ibin           !! the bin index
    real(kind=f), intent(out)             :: mmr(cstate%f_NZ) !! the bin mass mixing ratio [kg/kg]
    integer, intent(out)                  :: rc             !! return code negative indicates failure
    real(kind=f), optional, intent(out)   :: nmr(cstate%f_NZ) !! number mixing ratio [#/kg]
    real(kind=f), optional, intent(out)   :: numberDensity(cstate%f_NZ)  !! number density [#/cm3]
    real(kind=f), optional, intent(out)   :: r_wet(cstate%f_NZ)          !! wet particle radius [cm]
    real(kind=f), optional, intent(out)   :: rhop_wet(cstate%f_NZ)       !! wet particle density [g/cm3]
    
    integer                               :: ienconc        !! index of element that is the particle concentration for the group
    integer                               :: igroup         ! Group containing this bin

    ! Assume success.
    rc = RC_OK
    
    ! Determine the particle group for the bin.    
    igroup = cstate%f_carma%f_element(ielem)%f_igroup

    ! Make sure there are enough elements allocated.
    if (ielem > cstate%f_carma%f_NELEM) then
      if (cstate%f_carma%f_do_print) write(cstate%f_carma%f_LUNOPRT, *) "CARMASTATE_SetDetrain:: ERROR - The specifed element (", &
        ielem, ") is larger than the number of elements (", cstate%f_carma%f_NELEM, ")."
      rc = RC_ERROR
      return
    end if
    
    ! Make sure there are enough bins allocated.
    if (ibin > cstate%f_carma%f_NBIN) then
      if (cstate%f_carma%f_do_print) write(cstate%f_carma%f_LUNOPRT, *) "CARMA_SetDetrainin:: ERROR - The specifed bin (", &
        ibin, ") is larger than the number of bins (", cstate%f_carma%f_NBIN, ")."
      rc = RC_ERROR
      return
    end if

    
    ! Use the specified mass mixing ratio and the air density to determine the mass
    ! of the particles in g/x/y/z.
    mmr(:) = cstate%f_pcd(:, ibin, ielem) / cstate%f_rhoa_wet(:)


    ! Handle the special cases for different types of elements ...
    if ((cstate%f_carma%f_element(ielem)%f_itype == I_INVOLATILE) .or. (cstate%f_carma%f_element(ielem)%f_itype == I_VOLATILE)) then
      mmr(:) = mmr(:) * cstate%f_carma%f_group(igroup)%f_rmass(ibin)
    else if (cstate%f_carma%f_element(ielem)%f_itype == I_CORE2MOM) then
      mmr(:) = mmr(:) / cstate%f_carma%f_group(igroup)%f_rmass(ibin)
    end if
       
    ! If this is the partcile # element, then determine some other statistics.
    ienconc = cstate%f_carma%f_group(igroup)%f_ienconc
    if (ienconc == ielem) then
      if (present(nmr))           nmr(:)             = (cstate%f_pcd(:, ibin, ielem) / cstate%f_rhoa_wet(:)) * 1000._f
      if (present(numberDensity)) numberDensity(:)   = cstate%f_pcd(:, ibin, ielem) / (cstate%f_xmet(:)*cstate%f_ymet(:)*cstate%f_zmet(:))
      if (present(r_wet))         r_wet(:)           = cstate%f_r_wet(:, ibin, igroup)
      if (present(rhop_wet))      rhop_wet(:)        = cstate%f_rhop_wet(:, ibin, igroup)
    else
      if (present(nmr))            nmr(:)             = CAM_FILL
      if (present(numberDensity))  numberDensity(:)   = CAM_FILL
    end if
    
   return
  end subroutine CARMASTATE_GetDetrain
  

  !! Gets the mass mixing ratio for the gas (igas). After a call to CARMA_Step(),
  !! the new mass mixing ratio of the gas can be retrieved.
  !!
  !! @author Chuck Bardeen
  !! @version Feb-2009
  !! @see CARMA_AddGas
  !! @see CARMA_GetGas
  !! @see CARMA_Step 
  !! @see CARMASTATE_SetGas
  subroutine CARMASTATE_GetGas(cstate, igas, mmr, rc, satice, satliq, eqice, eqliq, wtpct)
    type(carmastate_type), intent(in)     :: cstate            !! the carma state object
    integer, intent(in)                   :: igas              !! the gas index
    real(kind=f), intent(out)             :: mmr(cstate%f_NZ)    !! the gas mass mixing ratio [kg/kg]
    integer, intent(out)                  :: rc                !! return code, negative indicates failure
    real(kind=f), optional, intent(out)   :: satice(cstate%f_NZ) !! the gas saturation wrt ice
    real(kind=f), optional, intent(out)   :: satliq(cstate%f_NZ) !! the gas saturation wrt liquid
    real(kind=f), optional, intent(out)   :: eqice(cstate%f_NZ)  !! the gas vapor pressure wrt ice
    real(kind=f), optional, intent(out)   :: eqliq(cstate%f_NZ)  !! the gas vapor pressure wrt liquid
    real(kind=f), optional, intent(out)   :: wtpct(cstate%f_NZ)  !! weight percent aerosol composition

    ! Assume success.
    rc = RC_OK

    ! Make sure there are enough gases allocated.
    if (igas > cstate%f_carma%f_NGAS) then
      if (cstate%f_carma%f_do_print) write(cstate%f_carma%f_LUNOPRT, *) "CARMASTATE_GetGas:: ERROR - The specifed gas (", &
        igas, ") is larger than the number of gases (", cstate%f_carma%f_NGAS, ")."
      rc = RC_ERROR
      return
    end if
    
    ! Use the specified mass mixing ratio and the air density to determine the mass
    ! of the gas in g/x/y/z.
    mmr(:) = cstate%f_gc(:, igas) / cstate%f_rhoa_wet(:)

    if (present(satice)) satice(:) = cstate%f_supsati(:, igas) + 1._f
    if (present(satliq)) satliq(:) = cstate%f_supsatl(:, igas) + 1._f
    if (present(eqice))  eqice(:)  = cstate%f_pvapi(:, igas) / cstate%f_p(:)
    if (present(eqliq))  eqliq(:)  = cstate%f_pvapl(:, igas) / cstate%f_p(:)
    if (present(wtpct))  wtpct(:)  = cstate%f_wtpct(:)
    
    return
  end subroutine CARMASTATE_GetGas
  

  !! Gets information about the state of the atmosphere. After the CARMA_Step() call,
  !! a new atmospheric state is determined.
  !!
  !! @author Chuck Bardeen
  !! @version Feb-2009
  !! @see CARMA_Step 
  !! @see CARMASTATE_Create
  subroutine CARMASTATE_GetState(cstate, rc, t, p, rhoa_wet, rlheat)
    type(carmastate_type), intent(in)     :: cstate                !! the carma state object
    integer, intent(out)                  :: rc                    !! return code, negative indicates failure
    real(kind=f), optional, intent(out)   :: t(cstate%f_NZ)        !! the air temperature [K]
    real(kind=f), optional, intent(out)   :: p(cstate%f_NZ)        !! the air pressure [Pa]
    real(kind=f), optional, intent(out)   :: rhoa_wet(cstate%f_NZ) !! air density [kg m-3]
    real(kind=f), optional, intent(out)   :: rlheat(cstate%f_NZ)   !! latent heat [K/s]
    
    ! Assume success.
    rc = RC_OK

    ! Return the temperature, pressure, and/or density.
    if (present(t))         t(:) = cstate%f_t(:)
    
    ! DYNE -> Pa
    if (present(p))         p(:) = cstate%f_p(:) / RPA2CGS
    
    ! Convert rhoa from the scaled units to mks.
    if (present(rhoa_wet))  rhoa_wet(:) = (cstate%f_rhoa_wet(:) / (cstate%f_zmet(:)*cstate%f_xmet(:)*cstate%f_ymet(:))) * 1e6_f / 1e3_f
    
    if (present(rlheat))    rlheat(:) = cstate%f_rlheat(:)

    return
  end subroutine CARMASTATE_GetState
  

  !! Sets the mass of the bins (ibin) for each particle element (ielem) in the grid.
  !! This call should be made after CARMASTATE_Create() and before CARMA_Step().
  !!
  !! @author Chuck Bardeen
  !! @version Feb-2009
  !! @see CARMA_AddBin
  !! @see CARMA_Step 
  !! @see CARMASTATE_GetBin
  subroutine CARMASTATE_SetBin(cstate, ielem, ibin, mmr, rc, surface)
    type(carmastate_type), intent(inout)  :: cstate         !! the carma state object
    integer, intent(in)                   :: ielem          !! the element index
    integer, intent(in)                   :: ibin           !! the bin index
    real(kind=f), intent(in)              :: mmr(cstate%f_NZ) !! the bin mass mixing ratio [kg/kg]
    integer, intent(out)                  :: rc             !! return code, negative indicates failure
    real(kind=f), optional, intent(in)    :: surface        !! particles mass on the surface [kg/m2]
    
    integer                               :: igroup         ! Group containing this bin

    ! Assume success.
    rc = RC_OK
    
    ! Determine the particle group for the bin.    
    igroup = cstate%f_carma%f_element(ielem)%f_igroup

    ! Make sure there are enough elements allocated.
    if (ielem > cstate%f_carma%f_NELEM) then
      if (cstate%f_carma%f_do_print) write(cstate%f_carma%f_LUNOPRT, *) "CARMASTATE_SetBin:: ERROR - The specifed element (", &
        ielem, ") is larger than the number of elements (", cstate%f_carma%f_NELEM, ")."
      rc = RC_ERROR
      return
    end if
    
    ! Make sure there are enough bins allocated.
    if (ibin > cstate%f_carma%f_NBIN) then
      if (cstate%f_carma%f_do_print) write(cstate%f_carma%f_LUNOPRT, *) "CARMASTATE_SetBin:: ERROR - The specifed bin (", &
        ibin, ") is larger than the number of bins (", cstate%f_carma%f_NBIN, ")."
      rc = RC_ERROR
      return
    end if
    
    ! Use the specified mass mixing ratio and the air density to determine the mass
    ! of the particles in g/x/y/z.
    cstate%f_pc(:, ibin, ielem) = mmr(:) * cstate%f_rhoa_wet(:)
    
    ! Handle the special cases for different types of elements ...
    if ((cstate%f_carma%f_element(ielem)%f_itype == I_INVOLATILE) .or. (cstate%f_carma%f_element(ielem)%f_itype == I_VOLATILE)) then
      cstate%f_pc(:, ibin, ielem) = cstate%f_pc(:, ibin, ielem) / cstate%f_carma%f_group(igroup)%f_rmass(ibin)
    else if (cstate%f_carma%f_element(ielem)%f_itype == I_CORE2MOM) then
      cstate%f_pc(:, ibin, ielem) = cstate%f_pc(:, ibin, ielem) * cstate%f_carma%f_group(igroup)%f_rmass(ibin)
    end if
    
    ! If they specified an initial mass of particles on the surface, then use that
    ! value.
    if (present(surface)) then
      
      ! Convert from g/cm2 to kg/m2
      cstate%f_pc_surf(ibin, ielem) = surface / 1e4_f * 1e3_f

      ! Handle the special cases for different types of elements ...
      if ((cstate%f_carma%f_element(ielem)%f_itype == I_INVOLATILE) .or. (cstate%f_carma%f_element(ielem)%f_itype == I_VOLATILE)) then
        cstate%f_pc_surf(ibin, ielem) = cstate%f_pc_surf(ibin, ielem) / cstate%f_carma%f_group(igroup)%f_rmass(ibin)
      else if (cstate%f_carma%f_element(ielem)%f_itype == I_CORE2MOM) then
        cstate%f_pc_surf(ibin, ielem) = cstate%f_pc_surf(ibin, ielem) * cstate%f_carma%f_group(igroup)%f_rmass(ibin)
      end if
    else
      cstate%f_pc_surf(ibin, ielem) = 0.0_f
    end if
        
    return
  end subroutine CARMASTATE_SetBin
  

  !! Sets the mass of the detrained condensate for the bins (ibin) for each particle
  !! element (ielem) in the grid. This call should be made after CARMASTATE_Create()
  !! and before CARMA_Step().
  !!
  !! @author Chuck Bardeen
  !! @version May-2010
  !! @see CARMA_AddBin
  !! @see CARMA_Step 
  !! @see CARMASTATE_GetDetrain
  subroutine CARMASTATE_SetDetrain(cstate, ielem, ibin, mmr, rc)
    type(carmastate_type), intent(inout)  :: cstate         !! the carma state object
    integer, intent(in)                   :: ielem          !! the element index
    integer, intent(in)                   :: ibin           !! the bin index
    real(kind=f), intent(in)              :: mmr(cstate%f_NZ) !! the bin mass mixing ratio [kg/kg]
    integer, intent(out)                  :: rc             !! return code, negative indicates failure
    
    integer                               :: igroup         ! Group containing this bin

    ! Assume success.
    rc = RC_OK
    
    ! Determine the particle group for the bin.    
    igroup = cstate%f_carma%f_element(ielem)%f_igroup

    ! Make sure there are enough elements allocated.
    if (ielem > cstate%f_carma%f_NELEM) then
      if (cstate%f_carma%f_do_print) write(cstate%f_carma%f_LUNOPRT, *) "CARMASTATE_SetDetrain:: ERROR - The specifed element (", &
        ielem, ") is larger than the number of elements (", cstate%f_carma%f_NELEM, ")."
      rc = RC_ERROR
      return
    end if
    
    ! Make sure there are enough bins allocated.
    if (ibin > cstate%f_carma%f_NBIN) then
      if (cstate%f_carma%f_do_print) write(cstate%f_carma%f_LUNOPRT, *) "CARMASTATE_SetDetrain:: ERROR - The specifed bin (", &
        ibin, ") is larger than the number of bins (", cstate%f_carma%f_NBIN, ")."
      rc = RC_ERROR
      return
    end if
    
    ! Use the specified mass mixing ratio and the air density to determine the mass
    ! of the particles in g/x/y/z.
    cstate%f_pcd(:, ibin, ielem) = mmr(:) * cstate%f_rhoa_wet(:)
    
    ! Handle the special cases for different types of elements ...
    if ((cstate%f_carma%f_element(ielem)%f_itype == I_INVOLATILE) .or. (cstate%f_carma%f_element(ielem)%f_itype == I_VOLATILE)) then
      cstate%f_pcd(:, ibin, ielem) = cstate%f_pcd(:, ibin, ielem) / cstate%f_carma%f_group(igroup)%f_rmass(ibin)
    else if (cstate%f_carma%f_element(ielem)%f_itype == I_CORE2MOM) then
      cstate%f_pcd(:, ibin, ielem) = cstate%f_pcd(:, ibin, ielem) * cstate%f_carma%f_group(igroup)%f_rmass(ibin)
    end if
        
    return
  end subroutine CARMASTATE_SetDetrain
  


  !! Sets the mass of the gas (igas) in the grid. This call should be made after
  !! CARMASTATE_Create() and before CARMA_Step().
  !!
  !! @author Chuck Bardeen
  !! @version Feb-2009
  !! @see CARMA_AddGas
  !! @see CARMA_GetGas
  !! @see CARMA_InitializeStep 
  !! @see CARMA_Step 
  subroutine CARMASTATE_SetGas(cstate, igas, mmr, rc, mmr_old, satice_old, satliq_old)
    type(carmastate_type), intent(inout)  :: cstate         !! the carma object
    integer, intent(in)                   :: igas           !! the gas index
    real(kind=f), intent(in)              :: mmr(cstate%f_NZ) !! the gas mass mixing ratio [kg/kg]
    integer, intent(out)                  :: rc             !! return code, negative indicates failure
    real(kind=f), intent(in), optional    :: mmr_old(cstate%f_NZ) !! the previous gas mass mixing ratio [kg/kg]
    real(kind=f), intent(inout), optional :: satice_old(cstate%f_NZ) !! the previous gas saturation wrt ice, calculates if -1
    real(kind=f), intent(inout), optional :: satliq_old(cstate%f_NZ) !! the previous gas saturation wrt liquid, calculates if -1
    
    real(kind=f)                          :: tnew(cstate%f_NZ)
    integer                               :: iz
    logical                               :: calculateOld
    
    ! Assume success.
    rc = RC_OK

    ! Make sure there are enough gases allocated.
    if (igas > cstate%f_carma%f_NGAS) then
      if (cstate%f_carma%f_do_print) write(cstate%f_carma%f_LUNOPRT, *) "CARMASTATE_SetGas:: ERROR - The specifed gas (", &
        igas, ") is larger than the number of gases (", cstate%f_carma%f_NGAS, ")."
      rc = RC_ERROR
      return
    end if
    
    if (cstate%f_carma%f_do_substep) then
      if (.not. present(mmr_old)) then
        if (cstate%f_carma%f_do_print) write(cstate%f_carma%f_LUNOPRT,*) "CARMASTATE_SetGas: Error - Need to specify mmr_old, satic_old, satliq_old when substepping."
        rc = RC_ERROR
        
        return
        
      else
        cstate%f_gcl(:, igas) = mmr_old(:) * cstate%f_rhoa_wet(:) * cstate%f_t(:) / cstate%f_told(:)
      
        ! A value of -1 for the saturation ratio means that it needs to be calculated from the old temperature
        ! and the old gc.
        !
        ! NOTE: This is typically just a problem for the first step, so we just need to get close.
        calculateOld = .false.
        if (present(satice_old) .and. present(satliq_old)) then
          if (any(satice_old(:) == -1._f) .or. any(satliq_old(:) == -1._f)) calculateOld = .true.
        else 
          calculateOld = .true.
        end if
        
        if (calculateOld) then
          
          ! This is a bit of a hack, because of the way CARMA has the vapor pressure and saturation
          ! routines implemented.
          
          ! Temporarily set the temperature and gc of to the old state
          
          tnew(:)      = cstate%f_t(:)
          cstate%f_t(:)  = cstate%f_told(:)
       
          cstate%f_gc(:, igas) = mmr_old(:) * cstate%f_rhoa_wet(:)
          
          do iz = 1, cstate%f_NZ
            call supersat(cstate%f_carma, cstate, iz, igas, rc)
            if (rc < RC_OK) return
          
            if (present(satice_old)) then
              if (satice_old(iz) == -1._f) then
                cstate%f_supsatiold(iz, igas) = cstate%f_supsati(iz, igas)
              else
                cstate%f_supsatiold(iz, igas) = satice_old(iz) - 1._f
              endif
            else
              cstate%f_supsatiold(iz, igas) = cstate%f_supsati(iz, igas)
            end if
            
            if (present(satliq_old)) then
              if (satliq_old(iz) == -1._f) then
                cstate%f_supsatlold(iz, igas) = cstate%f_supsatl(iz, igas)
              else
                cstate%f_supsatlold(iz, igas) = satliq_old(iz) - 1._f
              endif
            else
              cstate%f_supsatlold(iz, igas) = cstate%f_supsatl(iz, igas)
            end if
          end do
          
          cstate%f_t(:) = tnew(:)
        
        else
          cstate%f_supsatiold(:, igas) = satice_old(:) - 1._f
          cstate%f_supsatlold(:, igas) = satliq_old(:) - 1._f
        end if
      end if
    end if

    ! Use the specified mass mixing ratio and the air density to determine the mass
    ! of the gas in g/x/y/z.
    cstate%f_gc(:, igas)  = mmr(:) * cstate%f_rhoa_wet(:)
    
    return
  end subroutine CARMASTATE_SetGas
  
  
  !! Sets information about the state of the atmosphere.
  !!
  !! @author Chuck Bardeen
  !! @version Feb-2009
  !! @see CARMA_Step 
  !! @see CARMASTATE_Create
  subroutine CARMASTATE_SetState(cstate, rc, t, rhoa_wet)
    type(carmastate_type), intent(inout)  :: cstate              !! the carma state object
    integer, intent(out)                  :: rc                  !! return code, negative indicates failure
    real(kind=f), optional, intent(in)    :: t(cstate%f_NZ)        !! the air temperature [K]
    real(kind=f), optional, intent(in)    :: rhoa_wet(cstate%f_NZ) !! air density [kg m-3]
    
    ! Assume success.
    rc = RC_OK

    ! Return the temperature or density.
    if (present(t))         cstate%f_t(:) = t(:)
    
    ! Convert rhoa from mks to the scaled units.
    if (present(rhoa_wet))  cstate%f_rhoa_wet(:) = (rhoa_wet(:) * (cstate%f_zmet(:)*cstate%f_xmet(:)*cstate%f_ymet(:))) / 1e6_f * 1e3_f
    
    return
  end subroutine CARMASTATE_SetState
end module
