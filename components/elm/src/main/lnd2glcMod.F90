module lnd2glcMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle arrays used for exchanging data from land model to glc
  ! For now glc datais send and received on the lnd grid and decomposition.
  !
  ! The fields sent from the lnd component to the glc component via
  !  the coupler are labeled 's2x', or sno to coupler.
  ! The fields received by the lnd component from the glc component
  !  via the coupler are labeled 'x2s', or coupler to sno.
  ! 'Sno' is a misnomer in that the exchanged data are related to
  !  the ice beneath the snow, not the snow itself.  But by CESM convention,
  ! 'ice' refers to sea ice, not land ice.
  !
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use decompMod       , only : get_proc_bounds, bounds_type
  use domainMod       , only : ldomain
  use elm_varpar      , only : maxpatch_glcmec
  use elm_varctl      , only : iulog
  use elm_varcon      , only : spval, tfrz, namec
  use column_varcon   , only : col_itype_to_icemec_class
  use landunit_varcon , only : istice_mec, istsoil
  use abortutils      , only : endrun
  use WaterFluxType   , only : waterflux_type
  use TemperatureType , only : temperature_type
  use LandunitType    , only : lun_pp                
  use ColumnType      , only : col_pp
  use ColumnDataType  , only : col_es, col_wf  
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save

  ! land -> glc variables structure
  type, public :: lnd2glc_type
     real(r8), pointer :: tsrf_grc(:,:) => null()
     real(r8), pointer :: topo_grc(:,:) => null()
     real(r8), pointer :: qice_grc(:,:) => null()

   contains

     procedure, public  :: Init
     procedure, public  :: update_lnd2glc
     procedure, private :: InitAllocate
     procedure, private :: InitHistory

  end type lnd2glc_type

  ! !PUBLIC MEMBER FUNCTIONS:
  
  ! The following is public simply to support unit testing, and should not generally be
  ! called from outside this module.
  !
  ! Note that it is not a type-bound procedure, because it doesn't actually involve the
  ! lnd2glc_type. This suggests that perhaps it belongs in some other module.
  public :: bareland_normalization ! compute normalization factor for fluxes from the bare land portion of the grid cell
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(lnd2glc_type) :: this
    type(bounds_type), intent(in) :: bounds  

    call this%InitAllocate(bounds)
    call this%InitHistory(bounds)
    
  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize land variables required by glc
    !
    ! !USES:
    use elm_varcon , only : spval
    use histFileMod, only : hist_addfld1d
    !
    ! !ARGUMENTS:
    class(lnd2glc_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begg,endg 
    !------------------------------------------------------------------------

    begg = bounds%begg; endg = bounds%endg

    allocate(this%tsrf_grc(begg:endg,0:maxpatch_glcmec)) ; this%tsrf_grc(:,:)=0.0_r8
    allocate(this%topo_grc(begg:endg,0:maxpatch_glcmec)) ; this%topo_grc(:,:)=0.0_r8
    allocate(this%qice_grc(begg:endg,0:maxpatch_glcmec)) ; this%qice_grc(:,:)=0.0_r8

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !USES:
    use histFileMod, only : hist_addfld1d,hist_addfld2d 
    !
    ! !ARGUMENTS:
    class(lnd2glc_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    real(r8), pointer :: data2dptr(:,:)
    integer  :: begg, endg
    !---------------------------------------------------------------------

    begg = bounds%begg; endg = bounds%endg

    if (maxpatch_glcmec > 0) then
       this%qice_grc(begg:endg,0:maxpatch_glcmec) = spval
       ! For this and the following fields, set up a pointer to the field simply for the
       ! sake of changing the indexing, so that levels start with an index of 1, as is
       ! assumed by histFileMod - so levels go 1:(nec+1) rather than 0:nec
       data2dptr => this%qice_grc(:,0:maxpatch_glcmec)
       call hist_addfld2d (fname='QICE_FORC', units='mm/s', type2d='elevclas', &
            avgflag='A', long_name='qice forcing sent to GLC', &
            ptr_lnd=data2dptr, default='inactive')

       this%tsrf_grc(begg:endg,0:maxpatch_glcmec) = spval
       data2dptr => this%tsrf_grc(:,0:maxpatch_glcmec)
       call hist_addfld2d (fname='TSRF_FORC', units='K', type2d='elevclas', &
            avgflag='A', long_name='surface temperature sent to GLC', &
            ptr_lnd=data2dptr, default='inactive')

       this%topo_grc(begg:endg,0:maxpatch_glcmec) = spval
       data2dptr => this%topo_grc(:,0:maxpatch_glcmec)
       call hist_addfld2d (fname='TOPO_FORC', units='m', type2d='elevclas', &
            avgflag='A', long_name='topograephic height sent to GLC', &
            ptr_lnd=data2dptr, default='inactive')
    end if

  end subroutine InitHistory


  !------------------------------------------------------------------------------
  subroutine update_lnd2glc(this, bounds, num_do_smb_c, filter_do_smb_c, &
       temperature_vars, waterflux_vars, init)
    !
    ! !DESCRIPTION:
    ! Assign values to lnd2glc+
    !
    ! !ARGUMENTS:
    class(lnd2glc_type)    , intent(inout) :: this
    type(bounds_type)      , intent(in)    :: bounds  
    integer                , intent(in)    :: num_do_smb_c       ! number of columns in filter_do_smb_c
    integer                , intent(in)    :: filter_do_smb_c(:) ! column filter: columns where smb calculations are performed
    type(temperature_type) , intent(in)    :: temperature_vars
    type(waterflux_type)   , intent(in)    :: waterflux_vars
    logical                , intent(in)    :: init               ! if true=>only set a subset of fields
    !
    ! !LOCAL VARIABLES:
    integer  :: c, l, g, n, fc                   ! indices
    logical, allocatable :: fields_assigned(:,:) ! tracks whether fields have already been assigned for each index [begg:endg, 0:maxpatch_glcmec]
    real(r8) :: flux_normalization               ! factor by which fluxes should be normalized

    character(len=*), parameter :: subname = 'update_lnd2glc'
    !------------------------------------------------------------------------------

    ! Initialize to reasonable defaults

    this%qice_grc(bounds%begg : bounds%endg, :) = 0._r8
    this%tsrf_grc(bounds%begg : bounds%endg, :) = tfrz
    this%topo_grc(bounds%begg : bounds%endg, :) = 0._r8     
  
    ! Fill the lnd->glc data on the clm grid

    allocate(fields_assigned(bounds%begg:bounds%endg, 0:maxpatch_glcmec))
    fields_assigned(:,:) = .false.

    do fc = 1, num_do_smb_c
      c = filter_do_smb_c(fc)
      l = col_pp%landunit(c)
      g = col_pp%gridcell(c) 

      ! Set vertical index and a flux normalization, based on whether the column in question is glacier or vegetated.  
      if (lun_pp%itype(l) == istice_mec) then
         n = col_itype_to_icemec_class(col_pp%itype(c))
         flux_normalization = 1.0_r8
      else if (lun_pp%itype(l) == istsoil) then
         n = 0  !0-level index (bareland information)
         flux_normalization = bareland_normalization(c)
      else
         ! Other landunit types do not pass information in the lnd2glc fields.
         ! Note: for this to be acceptable, we need virtual vegetated columns in any grid
         ! cell that is made up solely of glacier plus some other special landunit (e.g.,
         ! glacier + lake) -- otherwise CISM wouldn't have any information for the non-
         ! glaciated portion of the grid cell.
         cycle
      end if

      ! Make sure we haven't already assigned the coupling fields for this point
      ! (this could happen, for example, if there were multiple columns in the
      ! istsoil landunit, which we aren't prepared to handle)
      if (fields_assigned(g,n)) then
         write(iulog,*) subname//' ERROR: attempt to assign coupling fields twice for the same index.'
         write(iulog,*) 'One possible cause is having multiple columns in the istsoil landunit,'
         write(iulog,*) 'which this routine cannot handle.'
         write(iulog,*) 'g, n = ', g, n
         call endrun(decomp_index=c, clmlevel=namec, msg=errMsg(__FILE__, __LINE__))
      end if

      ! Send surface temperature, topography, and SMB flux (qice) to coupler.
      ! t_soisno and glc_topo are valid even in initialization, so tsrf and topo
      ! are set here regardless of the value of init. But qflx_glcice is not valid
      ! until the run loop; thus, in initialization, we will use the default value
      ! for qice, as set above.
      fields_assigned(g,n) = .true.
      this%tsrf_grc(g,n) = col_es%t_soisno(c,1)
      this%topo_grc(g,n) = col_pp%glc_topo(c)
      if (.not. init) then
         this%qice_grc(g,n) = col_wf%qflx_glcice(c) * flux_normalization

         ! Check for bad values of qice
         if ( abs(this%qice_grc(g,n)) > 1.0_r8) then
            write(iulog,*) 'WARNING: qice out of bounds: g, n, qice =', g, n, this%qice_grc(g,n)
         end if
      end if

    end do

    deallocate(fields_assigned)
                
  end subroutine update_lnd2glc

  !-----------------------------------------------------------------------
  real(r8) function bareland_normalization(c)
    !
    ! !DESCRIPTION:
    ! Compute normalization factor for fluxes from the bare land portion of the grid
    ! cell. Fluxes should be multiplied by this factor before being sent to CISM.
    !
    ! The point of this is: CISM effectively has two land cover types: glaciated and
    ! bare. CLM, on the other hand, subdivides the bare land portion of the grid cell into
    ! multiple landunits. However, we currently don't do any sort of averaging of
    ! quantities computed in the different "bare land" landunits - instead, we simply send
    ! the values computed in the natural vegetated landunit - these fluxes (like SMB) are
    ! 0 in the other landunits. To achieve conservation, we need to normalize these
    ! natural veg. fluxes by the fraction of the "bare land" area accounted for by the
    ! natural veg. landunit.
    !
    ! For example, consider a grid cell that is:
    !   60% glacier_mec
    !   30% natural veg
    !   10% lake
    !
    ! According to CISM, this grid cell is 60% icesheet, 40% "bare land". Now suppose CLM
    ! has an SMB flux of 1m in the natural veg landunit. If we simply sent 1m of ice to
    ! CISM, conservation would be broken, since it would also apply 1m of ice to the 10%
    ! of the grid cell that CLM says is lake. So, instead, we must multiply the 1m of ice
    ! by (0.3/0.4), thus "spreading out" the SMB from the natural veg. landunit, so that
    ! 0.75m of ice is grown throughout the bare land portion of CISM.
    !
    ! Note: If the non-glaciated area of the grid cell is 0, then we arbitrarily return a
    ! normalization factor of 1.0, in order to avoid divide-by-zero errors.
    !
    ! Note: We currently aren't careful about how we would handle things if there are
    ! multiple columns within the vegetated landunit. If that possibility were introduced,
    ! this code - as well as the code in update_clm_s2x - may need to be reworked somewhat.
    !
    ! !USES:
    use subgridWeightsMod , only : get_landunit_weight
    !
    ! !ARGUMENTS:
    integer, intent(in) :: c  ! column index
    !
    ! !LOCAL VARIABLES:
    integer  :: t             ! topounit index
    real(r8) :: area_glacier  ! fractional area of the glacier_mec landunit in this grid cell
    real(r8) :: area_this_col ! fractional area of column c in the grid cell

    real(r8), parameter :: tol = 1.e-13_r8  ! tolerance for checking subgrid weight equality
    character(len=*), parameter :: subname = 'bareland_normalization'
    !-----------------------------------------------------------------------

    t = col_pp%topounit(c)
    
    area_glacier = get_landunit_weight(t, istice_mec)

    if (abs(area_glacier - 1.0_r8) < tol) then
       ! If the whole grid cell is glacier, then the normalization factor is arbitrary;
       ! set it to 1 so we don't do any normalization in this case
       bareland_normalization = 1.0_r8
    else
       area_this_col = col_pp%wttopounit(c)
       bareland_normalization = area_this_col / (1.0_r8 - area_glacier)
    end if

  end function bareland_normalization

end module lnd2glcMod

