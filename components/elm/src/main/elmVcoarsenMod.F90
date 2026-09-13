module elmVcoarsenMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Collapse ELM's multi-level history fields into a handful of 2-D fields,
  ! so an FME (Full Model Emulation) tape carries nothing but lat-lon images.
  !
  ! This is the ELM analogue of components/eam/src/control/eam_vcoarsen.F90
  ! and of the MPAS-Ocean fmeDepthCoarsening analysis member. All the
  ! vertical math is delegated to share/util/shr_vcoarsen_mod.F90; only the
  ! ELM-specific concerns live here: which data structure holds which field,
  ! the soil interface depths, and registering the outputs with ELM's history
  ! machinery so they average, remap and restart like any other field.
  !
  ! Two modes, because ELM's multi-level fields are of two different kinds:
  !
  !   1. SOIL PROFILE (type2d='levgrnd', nlevgrnd levels on the exponential
  !      soil grid). Collapsed by overlap-weighted averaging between depth
  !      boundaries given in metres -- the same operation MPAS-O's depth
  !      coarsening performs. n_out = size(vcoarsen_depth_bounds) - 1 output
  !      layers named <FLD>_0 .. <FLD>_{n_out-1} (0-based, matching EAM's
  !      vcoarsen naming). Registered on columns with the same
  !      l2g_scale_type='veg' the native field uses, so the column->gridcell
  !      average is identical.
  !
  !   2. RADIATION BANDS (type2d='numrad', 2 levels: visible, near-IR).
  !      There is nothing to integrate -- the two bands are distinct
  !      quantities -- so each is emitted under its own name,
  !      <FLD>_vis and <FLD>_nir. Registered on patches with
  !      c2l_scale_type='urbanf' to match the native fields.
  !
  ! Configuration (elm_inparm):
  !   vcoarsen_depth_bounds = 0.0, 0.07, 0.28, 1.00   ! m, strictly increasing
  !   vcoarsen_soil_flds    = 'TSOI', 'H2OSOI'
  !   vcoarsen_band_flds    = 'ALBD', 'ALBI'
  !
  ! The default bounds above are ERA5-Land's first three soil layers
  ! (L1 0-7 cm, L2 7-28 cm, L3 28-100 cm), so <FLD>_0/_1/_2 line up with the
  ! SM1/SM2/SM3 predictors an emulator is likely to be conditioned on.
  ! ERA5-Land's fourth layer is 100-289 cm; append 2.89 for a matching
  ! four-layer split. A Noah-MP-style 0/7/21/72 cm scheme is the other
  ! convention in common use. None of these align with ELM interfaces --
  ! the exponential soil grid has interfaces at 1.75, 4.51, 9.06, 16.55,
  ! 28.91, 49.29, 82.89 cm ... -- which is exactly why this collapses by
  ! DEPTH rather than by level-index range: the overlap weighting splits a
  ! partially covered ELM layer by the fraction that falls inside the bound.
  !
  ! Leave the field lists empty (the default) and the module is inert.
  !
  ! Adding a field means adding a case to get_soil_source / get_band_source
  ! below: ELM has no name->array registry the way EAM has constituents and
  ! the physics buffer, so the mapping is an explicit table.
  !
  ! Averaging caveat: the soil collapse is a THICKNESS-weighted mean, which
  ! is right for intensive quantities (a temperature, a volumetric water
  ! fraction). A per-layer MASS (kg/m2, e.g. SOILLIQ) would need a sum, not a
  ! mean -- that mode is deliberately not implemented rather than silently
  ! giving the wrong answer, so validate_soil_field rejects the fields it
  ! does not know.
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use abortutils     , only : endrun
  use spmdMod        , only : masterproc
  use elm_varctl     , only : iulog
  use elm_varcon     , only : spval, zisoi
  use elm_varpar     , only : nlevgrnd, numrad
  use decompMod      , only : bounds_type
  use SurfaceAlbedoType, only : surfalb_type
  !
  implicit none
  private
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: elm_vcoarsen_init    ! validate config, allocate, register fields
  public :: elm_vcoarsen_update  ! recompute the collapsed fields each step
  public :: elm_vcoarsen_active  ! is anything configured?
  !
  ! !PUBLIC DATA (namelist, read by controlMod into elm_inparm):
  integer, parameter, public :: vcoarsen_max_bounds = 21   ! => at most 20 layers
  integer, parameter, public :: vcoarsen_max_flds   = 20

  ! Depth boundaries in metres, strictly increasing, starting at or above 0.
  ! Entries at or below the "unset" sentinel terminate the list.
  real(r8), public :: vcoarsen_depth_bounds(vcoarsen_max_bounds) = -1.0_r8
  character(len=64), public :: vcoarsen_soil_flds(vcoarsen_max_flds) = ' '
  character(len=64), public :: vcoarsen_band_flds(vcoarsen_max_flds) = ' '
  !
  ! !PRIVATE DATA:
  integer :: n_out        = 0   ! number of coarsened soil layers
  integer :: n_soil_flds  = 0
  integer :: n_band_flds  = 0
  logical :: is_active    = .false.

  ! Output storage. One contiguous column per output field so a rank-1
  ! pointer into it can be handed to hist_addfld1d.
  !   soil_out(begc:endc, (i-1)*n_out + k)   field i, layer k
  !   band_out(begp:endp, (i-1)*numrad + k)  field i, band k
  real(r8), pointer :: soil_out(:,:) => null()
  real(r8), pointer :: band_out(:,:) => null()

  ! Soil interface depths replicated across columns, the layout
  ! shr_vcoarsen_avg_cols expects. zisoi is a single global profile in ELM,
  ! so this is built once at init rather than per call.
  real(r8), allocatable :: zi_cols(:,:)   ! (ncol, nlevgrnd+1)

  character(len=*), parameter :: unset_msg = &
       'vcoarsen_depth_bounds must be strictly increasing and non-negative'
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  logical function elm_vcoarsen_active()
    elm_vcoarsen_active = is_active
  end function elm_vcoarsen_active

  !-----------------------------------------------------------------------
  subroutine elm_vcoarsen_init(bounds)
    !
    ! !DESCRIPTION:
    ! Validate the namelist, allocate the output arrays and register the
    ! collapsed fields. Must run during the addfld phase -- after every
    ! component InitHistory and before htapes_fieldlist.
    !
    ! !USES:
    use histFileMod, only : hist_addfld1d
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: i, k, c, slot, nslots
    real(r8), pointer :: p1d(:)
    character(len=64)  :: fname
    character(len=256) :: lname
    character(len=32)  :: units
    character(len=*), parameter :: subname = 'elm_vcoarsen_init'
    !-----------------------------------------------------------------------

    ! --- count and validate the depth boundaries -------------------------
    n_out = 0
    do i = 1, vcoarsen_max_bounds
       if (vcoarsen_depth_bounds(i) < 0.0_r8) exit
       n_out = n_out + 1
    end do
    n_out = max(0, n_out - 1)     ! n boundaries -> n-1 layers

    if (n_out > 0) then
       if (vcoarsen_depth_bounds(1) < 0.0_r8) then
          call endrun(msg=subname//' ERROR: '//unset_msg//' '//errMsg(__FILE__, __LINE__))
       end if
       do i = 1, n_out
          if (vcoarsen_depth_bounds(i+1) <= vcoarsen_depth_bounds(i)) then
             call endrun(msg=subname//' ERROR: '//unset_msg//' '//errMsg(__FILE__, __LINE__))
          end if
       end do
    end if

    ! --- count the field lists -------------------------------------------
    n_soil_flds = 0
    do i = 1, vcoarsen_max_flds
       if (len_trim(vcoarsen_soil_flds(i)) == 0) exit
       n_soil_flds = n_soil_flds + 1
    end do
    n_band_flds = 0
    do i = 1, vcoarsen_max_flds
       if (len_trim(vcoarsen_band_flds(i)) == 0) exit
       n_band_flds = n_band_flds + 1
    end do

    if (n_soil_flds > 0 .and. n_out == 0) then
       call endrun(msg=subname//' ERROR: vcoarsen_soil_flds is set but '// &
            'vcoarsen_depth_bounds gives no layers '//errMsg(__FILE__, __LINE__))
    end if

    is_active = (n_soil_flds > 0 .or. n_band_flds > 0)
    if (.not. is_active) return

    if (masterproc) then
       write(iulog,*) subname,': collapsing ',n_soil_flds,' soil field(s) into ', &
            n_out,' layer(s), bounds (m): ',vcoarsen_depth_bounds(1:n_out+1)
       write(iulog,*) subname,': splitting ',n_band_flds,' radiation-band field(s)'
       if (n_out > 0) then
          if (vcoarsen_depth_bounds(n_out+1) > zisoi(nlevgrnd)) then
             write(iulog,*) subname,': WARNING deepest bound ', &
                  vcoarsen_depth_bounds(n_out+1),' m is below the bottom of the ', &
                  'ELM soil column (',zisoi(nlevgrnd),' m); that layer will only ', &
                  'average what exists above it'
          end if
       end if
    end if

    ! --- soil: allocate, replicate the interface profile, register -------
    if (n_soil_flds > 0) then
       nslots = n_soil_flds * n_out
       allocate(soil_out(bounds%begc:bounds%endc, nslots))
       soil_out(:,:) = spval

       allocate(zi_cols(bounds%begc:bounds%endc, nlevgrnd+1))
       do k = 1, nlevgrnd+1
          do c = bounds%begc, bounds%endc
             zi_cols(c, k) = zisoi(k-1)
          end do
       end do

       do i = 1, n_soil_flds
          call validate_soil_field(vcoarsen_soil_flds(i), units)
          do k = 1, n_out
             slot = (i-1)*n_out + k
             call make_layer_name(vcoarsen_soil_flds(i), k-1, fname)
             write(lname,'(a,a,f0.3,a,f0.3,a)') trim(vcoarsen_soil_flds(i)), &
                  ', thickness-weighted mean over ', vcoarsen_depth_bounds(k), &
                  ' to ', vcoarsen_depth_bounds(k+1), ' m depth'
             p1d => soil_out(:, slot)
             call hist_addfld1d(fname=trim(fname), units=trim(units), &
                  avgflag='A', long_name=trim(lname), &
                  ptr_col=p1d, l2g_scale_type='veg', default='inactive')
          end do
       end do
    end if

    ! --- radiation bands: allocate and register --------------------------
    if (n_band_flds > 0) then
       nslots = n_band_flds * numrad
       allocate(band_out(bounds%begp:bounds%endp, nslots))
       band_out(:,:) = spval

       do i = 1, n_band_flds
          call validate_band_field(vcoarsen_band_flds(i), units)
          do k = 1, numrad
             slot = (i-1)*numrad + k
             call make_band_name(vcoarsen_band_flds(i), k, fname)
             write(lname,'(a,a,a,a)') trim(vcoarsen_band_flds(i)), ', ', &
                  trim(band_label(k)), ' band'
             p1d => band_out(:, slot)
             call hist_addfld1d(fname=trim(fname), units=trim(units), &
                  avgflag='A', long_name=trim(lname), &
                  ptr_patch=p1d, c2l_scale_type='urbanf', default='inactive')
          end do
       end do
    end if

  end subroutine elm_vcoarsen_init

  !-----------------------------------------------------------------------
  subroutine elm_vcoarsen_update(bounds, surfalb_vars)
    !
    ! !DESCRIPTION:
    ! Recompute every collapsed field. Call once per timestep immediately
    ! before hist_update_hbuf, so the history buffers accumulate the values
    ! for this step.
    !
    ! surfalb_vars is passed in rather than taken from elm_instMod: that
    ! module uses controlMod, which in turn uses this one for the namelist,
    ! and the cycle would not compile.
    !
    ! !USES:
    use shr_vcoarsen_mod, only : shr_vcoarsen_avg_cols, shr_vcoarsen_select_index
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in) :: bounds
    type(surfalb_type), intent(in) :: surfalb_vars
    !
    ! !LOCAL VARIABLES:
    integer :: i, k, slot, ncol, npatch
    real(r8), pointer :: src2d(:,:)
    real(r8), allocatable :: tmp(:,:)
    integer,  allocatable :: nlev_max(:)
    character(len=*), parameter :: subname = 'elm_vcoarsen_update'
    !-----------------------------------------------------------------------

    if (.not. is_active) return

    ! --- soil profile ----------------------------------------------------
    if (n_soil_flds > 0) then
       ncol = bounds%endc - bounds%begc + 1
       allocate(tmp(ncol, n_out))
       do i = 1, n_soil_flds
          call get_soil_source(vcoarsen_soil_flds(i), src2d)
          call shr_vcoarsen_avg_cols(src2d(bounds%begc:bounds%endc, 1:nlevgrnd), &
               zi_cols(bounds%begc:bounds%endc, 1:nlevgrnd+1), ncol, nlevgrnd, &
               vcoarsen_depth_bounds(1:n_out+1), n_out, spval, tmp)
          do k = 1, n_out
             slot = (i-1)*n_out + k
             soil_out(bounds%begc:bounds%endc, slot) = tmp(1:ncol, k)
          end do
       end do
       deallocate(tmp)
    end if

    ! --- radiation bands -------------------------------------------------
    if (n_band_flds > 0) then
       npatch = bounds%endp - bounds%begp + 1
       allocate(nlev_max(npatch))
       nlev_max(:) = numrad
       do i = 1, n_band_flds
          call get_band_source(vcoarsen_band_flds(i), surfalb_vars, src2d)
          do k = 1, numrad
             slot = (i-1)*numrad + k
             call shr_vcoarsen_select_index( &
                  src2d(bounds%begp:bounds%endp, 1:numrad), npatch, numrad, &
                  k, nlev_max, spval, &
                  band_out(bounds%begp:bounds%endp, slot))
          end do
       end do
       deallocate(nlev_max)
    end if

  end subroutine elm_vcoarsen_update

  !-----------------------------------------------------------------------
  ! Private helpers
  !-----------------------------------------------------------------------

  subroutine get_soil_source(fname, src2d)
    !
    ! Bind to the live (begc:endc, 1:nlevgrnd) array behind a named soil field.
    ! ELM has no name->array registry, so this is an explicit table; keep it
    ! in step with validate_soil_field.
    !
    ! !USES:
    use ColumnDataType, only : col_es, col_ws
    !
    ! !ARGUMENTS:
    character(len=*) , intent(in)  :: fname
    real(r8), pointer, intent(out) :: src2d(:,:)
    !-----------------------------------------------------------------------

    select case (trim(fname))
    case ('TSOI')
       src2d => col_es%t_soisno
    case ('H2OSOI')
       src2d => col_ws%h2osoi_vol
    case default
       call endrun(msg='elm_vcoarsen: get_soil_source: unsupported field "'// &
            trim(fname)//'" '//errMsg(__FILE__, __LINE__))
    end select

  end subroutine get_soil_source

  !-----------------------------------------------------------------------
  subroutine validate_soil_field(fname, units)
    !
    ! Accept only fields whose thickness-weighted mean is meaningful, and
    ! return the unit string (the mean is unit-preserving).
    !
    character(len=*) , intent(in)  :: fname
    character(len=*) , intent(out) :: units
    !-----------------------------------------------------------------------

    select case (trim(fname))
    case ('TSOI')
       units = 'K'
    case ('H2OSOI')
       units = 'mm3/mm3'
    case default
       call endrun(msg='elm_vcoarsen: vcoarsen_soil_flds entry "'//trim(fname)// &
            '" is not supported. Add it to get_soil_source and '// &
            'validate_soil_field, and check first that a thickness-weighted '// &
            'MEAN is the right reduction for it (a per-layer mass needs a '// &
            'sum instead). '//errMsg(__FILE__, __LINE__))
    end select

  end subroutine validate_soil_field

  !-----------------------------------------------------------------------
  subroutine get_band_source(fname, surfalb_vars, src2d)
    !
    ! !ARGUMENTS:
    character(len=*)  , intent(in)  :: fname
    type(surfalb_type), intent(in)  :: surfalb_vars
    real(r8), pointer , intent(out) :: src2d(:,:)
    !-----------------------------------------------------------------------

    select case (trim(fname))
    case ('ALBD')
       src2d => surfalb_vars%albd_patch
    case ('ALBI')
       src2d => surfalb_vars%albi_patch
    case default
       call endrun(msg='elm_vcoarsen: get_band_source: unsupported field "'// &
            trim(fname)//'" '//errMsg(__FILE__, __LINE__))
    end select

  end subroutine get_band_source

  !-----------------------------------------------------------------------
  subroutine validate_band_field(fname, units)
    !
    ! !ARGUMENTS:
    character(len=*), intent(in)  :: fname
    character(len=*), intent(out) :: units
    !-----------------------------------------------------------------------

    select case (trim(fname))
    case ('ALBD', 'ALBI')
       units = 'proportion'
    case default
       call endrun(msg='elm_vcoarsen: vcoarsen_band_flds entry "'//trim(fname)// &
            '" is not supported. Add it to get_band_source and '// &
            'validate_band_field. '//errMsg(__FILE__, __LINE__))
    end select

  end subroutine validate_band_field

  !-----------------------------------------------------------------------
  subroutine make_layer_name(base, layer_idx, out_name)
    ! e.g. base='TSOI', layer_idx=0 => 'TSOI_0'  (0-based, as in EAM vcoarsen)
    character(len=*), intent(in)  :: base
    integer         , intent(in)  :: layer_idx
    character(len=*), intent(out) :: out_name
    character(len=8) :: idx_str
    !-----------------------------------------------------------------------

    write(idx_str,'(i0)') layer_idx
    out_name = trim(base)//'_'//trim(idx_str)

  end subroutine make_layer_name

  !-----------------------------------------------------------------------
  subroutine make_band_name(base, band_idx, out_name)
    ! e.g. base='ALBD', band_idx=1 => 'ALBD_vis'
    character(len=*), intent(in)  :: base
    integer         , intent(in)  :: band_idx
    character(len=*), intent(out) :: out_name
    !-----------------------------------------------------------------------

    out_name = trim(base)//'_'//trim(band_label(band_idx))

  end subroutine make_band_name

  !-----------------------------------------------------------------------
  function band_label(band_idx) result(lbl)
    ! ELM's numrad=2 solar bands, in order: 1 = visible, 2 = near-infrared.
    integer, intent(in) :: band_idx
    character(len=8)    :: lbl
    !-----------------------------------------------------------------------

    select case (band_idx)
    case (1)     ; lbl = 'vis'
    case (2)     ; lbl = 'nir'
    case default ; write(lbl,'(a,i0)') 'band', band_idx
    end select

  end function band_label

end module elmVcoarsenMod
