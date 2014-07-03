module clm_glclnd

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle arrays used for exchanging data between glc and land model.
  ! Based on clm_atmlnd (but without mapping routines because glc data
  !  is send and received on the lnd decomposition, at least for now).
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
  use decompMod      , only : bounds_type
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use clm_varpar     , only : maxpatch_glcmec
  use clm_varctl     , only : iulog, glc_smb
  use abortutils     , only : endrun
  !
  ! !REVISION HISTORY:
  ! Created by William Lipscomb, Dec. 2007, based on clm_atmlnd.F90.
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save

  ! glc -> land variables structure
  type, public :: glc2lnd_type
     real(r8), pointer :: frac(:,:) => null()
     real(r8), pointer :: topo(:,:) => null()
     real(r8), pointer :: hflx(:,:) => null()
     real(r8), pointer :: icemask(:) => null()
  end type glc2lnd_type

  ! land -> glc variables structure
  type, public :: lnd2glc_type
     real(r8), pointer :: tsrf(:,:) => null()
     real(r8), pointer :: topo(:,:) => null()
     real(r8), pointer :: qice(:,:) => null()
  end type lnd2glc_type

  type (lnd2glc_type), public, target :: clm_s2x  ! s2x fields on clm grid
  type (glc2lnd_type), public, target :: clm_x2s  ! x2s fields on clm grid

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: init_glc2lnd_type
  public :: init_lnd2glc_type
  public :: update_clm_x2s
  public :: update_clm_s2x
  
  ! The following is public simply to support unit testing, and should not generally be
  ! called from outside this module
  public :: bareland_normalization ! compute normalization factor for fluxes from the bare land portion of the grid cell

  ! !PRIVATE MEMBER FUNCTIONS:
  private :: update_clm_x2s_icemask ! update icemask based on input from GLC
  private :: update_clm_x2s_fracs   ! update subgrid fractions based on input from GLC
  private :: update_clm_x2s_topo    ! update column-level topographic heights based on input from GLC

  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine init_glc2lnd_type(bounds, x2s)
    !
    ! !DESCRIPTION:
    ! Initialize glc variables required by the land
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    type (glc2lnd_type), intent(inout):: x2s
    !------------------------------------------------------------------------

    allocate(x2s%frac(bounds%begg:bounds%endg,0:maxpatch_glcmec))
    x2s%frac(bounds%begg:bounds%endg,0:maxpatch_glcmec)=0.0_r8
    allocate(x2s%topo(bounds%begg:bounds%endg,0:maxpatch_glcmec))
    x2s%topo(bounds%begg:bounds%endg,0:maxpatch_glcmec)=0.0_r8
    allocate(x2s%hflx(bounds%begg:bounds%endg,0:maxpatch_glcmec))
    x2s%hflx(bounds%begg:bounds%endg,0:maxpatch_glcmec)=0.0_r8
    allocate(x2s%icemask(bounds%begg:bounds%endg))
    x2s%icemask(bounds%begg:bounds%endg)=0.0_r8

  end subroutine init_glc2lnd_type

  !------------------------------------------------------------------------
  subroutine init_lnd2glc_type(bounds, s2x)
    !
    ! !DESCRIPTION:
    ! Initialize land variables required by glc
    
    use clm_varcon, only : tfrz
    ! !ARGUMENTS:

    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    type (lnd2glc_type), intent(inout):: s2x
    !------------------------------------------------------------------------

    allocate(s2x%tsrf(bounds%begg:bounds%endg,0:maxpatch_glcmec))
    s2x%tsrf(bounds%begg:bounds%endg,0:maxpatch_glcmec)=tfrz
    allocate(s2x%topo(bounds%begg:bounds%endg,0:maxpatch_glcmec))
    s2x%topo(bounds%begg:bounds%endg,0:maxpatch_glcmec)=0.0_r8
    allocate(s2x%qice(bounds%begg:bounds%endg,0:maxpatch_glcmec))
    s2x%qice(bounds%begg:bounds%endg,0:maxpatch_glcmec)=0.0_r8

  end subroutine init_lnd2glc_type

  !-----------------------------------------------------------------------
  subroutine update_clm_x2s(bounds)
    !
    ! !DESCRIPTION:
    ! Update values to derived-type CLM variables based on input from GLC (via the coupler)
    !
    ! icemask and topo are always updated (although note that this routine should only be
    ! called when create_glacier_mec_landunit is true, or some similar condition; this
    ! should be controlled in a conditional around the call to this routine); fracs are
    ! updated if glc_do_dynglacier is true
    !
    ! !USES:
    use clm_varctl , only : glc_do_dynglacier
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in) :: bounds ! bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'update_clm_x2s'
    !-----------------------------------------------------------------------

    call update_clm_x2s_icemask(bounds)

    if (glc_do_dynglacier) then
       call update_clm_x2s_fracs(bounds)
    end if

    call update_clm_x2s_topo(bounds)

  end subroutine update_clm_x2s

  !-----------------------------------------------------------------------
  subroutine update_clm_x2s_icemask(bounds)
    !
    ! !DESCRIPTION:
    ! Set CLM's internal icemask to be the same as icemask received from CISM via coupler
    !
    ! !USES:
    use domainMod , only : ldomain
    use clmtype   , only : grc, nameg
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in) :: bounds ! bounds
    !
    ! !LOCAL VARIABLES:
    integer :: g  ! grid cell index
    
    character(len=*), parameter :: subname = 'update_clm_x2s_icemask'
    !-----------------------------------------------------------------------
    
    do g = bounds%begg, bounds%endg
       grc%icemask(g)=clm_x2s%icemask(g)

       ! Ensure that icemask is a subset of glcmask. This is needed because we allocated
       ! memory based on glcmask, so it is a problem if the ice sheet tries to expand
       ! beyond the area defined by glcmask.
       if (grc%icemask(g) > 0._r8 .and. ldomain%glcmask(g) == 0._r8) then
          write(iulog,*) subname//' ERROR: icemask must be a subset of glcmask.'
          write(iulog,*) 'You can fix this problem by adding more grid cells'
          write(iulog,*) 'to the mask defined by the fglcmask file.'
          write(iulog,*) '(Change grid cells to 1 everywhere that CISM can operate.)'
          call endrun(decomp_index=g, clmlevel=nameg, msg=errMsg(__FILE__, __LINE__))
       end if
    end do

  end subroutine update_clm_x2s_icemask


  !-----------------------------------------------------------------------
  subroutine update_clm_x2s_fracs(bounds)
    !
    ! !DESCRIPTION:
    ! Update subgrid fractions based on input from GLC (via the coupler)
    !
    ! The weights updated here are some col%wtlunit and lun%wtgcell values
    !
    ! !USES:
    use clmtype           , only : grc, lun, col
    use clm_varcon        , only : istice_mec, ispval, col_itype_to_icemec_class
    use subgridWeightsMod , only : set_landunit_weight
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  ! bounds
    !
    ! !LOCAL VARIABLES:
    integer :: g,c                              ! indices
    real(r8):: area_ice_mec                     ! area of the ice_mec landunit
    integer :: l_ice_mec                        ! index of the ice_mec landunit
    integer :: icemec_class                     ! current icemec class (1..maxpatch_glcmec)
    logical :: frac_assigned(1:maxpatch_glcmec) ! whether clm_x2s%frac has been assigned for each elevation class
    logical :: error                            ! if an error was found
    
    character(len=*), parameter :: subname = 'update_clm_x2s_fracs'
    !-----------------------------------------------------------------------
    
    do g = bounds%begg, bounds%endg
       ! Values from GLC are only valid within the icemask, so we only update CLM's areas there
       if (grc%icemask(g) > 0._r8) then

          ! Set total icemec landunit area
          area_ice_mec = sum(clm_x2s%frac(g, 1:maxpatch_glcmec))
          call set_landunit_weight(g, istice_mec, area_ice_mec)

          ! If new landunit area is greater than 0, then update column areas
          ! (If new landunit area is 0, col%wtlunit is arbitrary, so we might as well keep the existing values)
          if (area_ice_mec > 0) then
             ! Determine index of the glc_mec landunit
             l_ice_mec = grc%landunit_indices(istice_mec, g)
             if (l_ice_mec == ispval) then
                write(iulog,*) subname//' ERROR: no ice_mec landunit found within the icemask, for g = ', g
                call endrun()
             end if
          
             frac_assigned(1:maxpatch_glcmec) = .false.
             do c = lun%coli(l_ice_mec), lun%colf(l_ice_mec)
                icemec_class = col_itype_to_icemec_class(col%itype(c))
                col%wtlunit(c) = clm_x2s%frac(g, icemec_class) / lun%wtgcell(l_ice_mec)
                frac_assigned(icemec_class) = .true.
             end do

             ! Confirm that all elevation classes that have non-zero area according to
             ! clm_x2s%frac have been assigned to a column in CLM's data structures
             error = .false.
             do icemec_class = 1, maxpatch_glcmec
                if (clm_x2s%frac(g, icemec_class) > 0._r8 .and. &
                     .not. frac_assigned(icemec_class)) then
                   error = .true.
                end if
             end do
             if (error) then
                write(iulog,*) subname//' ERROR: at least one glc_mec column has non-zero area from the coupler,'
                write(iulog,*) 'but there was no slot in memory for this column; g = ', g
                write(iulog,*) 'clm_x2s%frac(g, 1:maxpatch_glcmec) = ', &
                     clm_x2s%frac(g, 1:maxpatch_glcmec)
                write(iulog,*) 'frac_assigned(1:maxpatch_glcmec) = ', &
                     frac_assigned(1:maxpatch_glcmec)
                call endrun()
             end if  ! error
          end if  ! area_ice_mec > 0
       end if  ! grc%icemask(g) > 0
    end do  ! g

  end subroutine update_clm_x2s_fracs

  !-----------------------------------------------------------------------
  subroutine update_clm_x2s_topo(bounds)
    !
    ! !DESCRIPTION:
    ! Update column-level topographic heights based on input from GLC (via the coupler)
    !
    ! !USES:
    use clmtype    , only : grc, lun, col, cps
    use clm_varcon , only : istice_mec, col_itype_to_icemec_class
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in) :: bounds ! bounds
    !
    ! !LOCAL VARIABLES:
    integer :: c, l, g      ! indices
    integer :: icemec_class ! current icemec class (1..maxpatch_glcmec)

    character(len=*), parameter :: subname = 'update_clm_x2s_topo'
    !-----------------------------------------------------------------------
    
    ! It is tempting to use the do_smb_c filter here, since we only need glc_topo inside
    ! this filter. But the problem with using the filter is that this routine is called
    ! before the filters are updated to reflect the updated weights. As long as
    ! glacier_mec, natural veg and any other landunit within the smb filter are always
    ! active, regardless of their weights, this isn't a problem. But we don't want to
    ! build in assumptions that those rules will be in place regarding active flags.
    ! Other ways around this problem would be:
    ! (1) Use the inactive_and_active filter  - but we're trying to avoid use of that
    !     filter if possible, because it can be confusing
    ! (2) Call this topo update routine later in the driver loop, after filters have been
    !     updated  - but that leads to greater complexity in the driver loop.
    ! So it seems simplest just to take the minor performance hit of setting glc_topo
    ! over all columns, even those outside the do_smb_c filter.
    
    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       g = col%gridcell(c)

       ! Values from GLC are only valid within the icemask, so we only update CLM's topo values there
       if (grc%icemask(g) > 0._r8) then
          if (lun%itype(l) == istice_mec) then
             icemec_class = col_itype_to_icemec_class(col%itype(c))
          else
             ! If not on a glaciated column, assign topography to the bare-land value determined by GLC.
             icemec_class = 0
          end if

          cps%glc_topo(c) = clm_x2s%topo(g, icemec_class)
       end if
    end do

  end subroutine update_clm_x2s_topo


  !------------------------------------------------------------------------------
  subroutine update_clm_s2x(bounds, num_do_smb_c, filter_do_smb_c, init)
    !
    ! !DESCRIPTION:
    ! Update values sent to GLC (via the coupler)
    !
    ! !USES:
    use clmtype   , only : col, lun, ces, cps, cwf, namec
    use clm_varcon, only : istice_mec, istsoil
    use clm_varcon, only : spval, tfrz, col_itype_to_icemec_class
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type) , intent(in) :: bounds             ! bounds   
    integer           , intent(in) :: num_do_smb_c       ! number of columns in filter_do_smb_c
    integer           , intent(in) :: filter_do_smb_c(:) ! column filter: columns where smb calculations are performed
    logical           , intent(in) :: init               ! if true=>only set a subset of fields
    !
    ! !LOCAL VARIABLES:
    integer  :: c, l, g, n, fc                   ! indices
    logical, allocatable :: fields_assigned(:,:) ! tracks whether fields have already been assigned for each index [begg:endg, 0:maxpatch_glcmec]
    real(r8) :: flux_normalization               ! factor by which fluxes should be normalized

    character(len=*), parameter :: subname = 'update_clm_s2x'
    !------------------------------------------------------------------------------
        
    ! Initialize to reasonable defaults
    clm_s2x%qice(bounds%begg : bounds%endg, :) = 0._r8
    clm_s2x%tsrf(bounds%begg : bounds%endg, :) = tfrz
    clm_s2x%topo(bounds%begg : bounds%endg, :) = 0._r8     
  
    allocate(fields_assigned(bounds%begg:bounds%endg, 0:maxpatch_glcmec))
    fields_assigned(:,:) = .false.

    do fc = 1, num_do_smb_c
      c = filter_do_smb_c(fc)
      l = col%landunit(c)
      g = col%gridcell(c) 

      ! Set vertical index and a flux normalization, based on whether the column in question is glacier or vegetated.  
      if (lun%itype(l) == istice_mec) then
         n = col_itype_to_icemec_class(col%itype(c))
         flux_normalization = 1.0_r8
      else if (lun%itype(l) == istsoil) then
         n = 0  !0-level index (bareland information)
         flux_normalization = bareland_normalization(c)
      else
         ! Other landunit types do not pass information in the clm_s2x fields.
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
      clm_s2x%tsrf(g,n) = ces%t_soisno(c,1)
      clm_s2x%topo(g,n) = cps%glc_topo(c)
      if (.not. init) then
         clm_s2x%qice(g,n) = cwf%qflx_glcice(c) * flux_normalization
         ! Check for bad values of qice
         if ( abs(clm_s2x%qice(g,n)) > 1.0_r8) then
            write(iulog,*) 'WARNING: qice out of bounds: g, n, qice =', g, n, clm_s2x%qice(g,n)
         end if
      end if

    end do

    deallocate(fields_assigned)
                
  end subroutine update_clm_s2x

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
    use clmtype           , only : col, grc
    use clm_varcon        , only : istice_mec
    use subgridWeightsMod , only : get_landunit_weight
    !
    ! !ARGUMENTS:
    integer, intent(in) :: c  ! column index
    !
    ! !LOCAL VARIABLES:
    integer  :: g             ! grid cell index
    real(r8) :: area_glacier  ! fractional area of the glacier_mec landunit in this grid cell
    real(r8) :: area_this_col ! fractional area of column c in the grid cell

    real(r8), parameter :: tol = 1.e-13_r8  ! tolerance for checking subgrid weight equality
    character(len=*), parameter :: subname = 'bareland_normalization'
    !-----------------------------------------------------------------------

    g = col%gridcell(c)

    area_glacier = get_landunit_weight(g, istice_mec)

    if (abs(area_glacier - 1.0_r8) < tol) then
       ! If the whole grid cell is glacier, then the normalization factor is arbitrary;
       ! set it to 1 so we don't do any normalization in this case
       bareland_normalization = 1.0_r8
    else
       area_this_col = col%wtgcell(c)
       bareland_normalization = area_this_col / (1.0_r8 - area_glacier)
    end if

  end function bareland_normalization


end module clm_glclnd

