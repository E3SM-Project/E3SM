module seq_domain_mct

  use shr_kind_mod, only: R8=>shr_kind_r8, IN=>shr_kind_in
  use shr_kind_mod, only: CL=>shr_kind_cl
  use shr_sys_mod,  only: shr_sys_flush, shr_sys_abort
  use shr_mpi_mod,  only: shr_mpi_min, shr_mpi_max

  use mct_mod
  use seq_comm_mct
  use seq_infodata_mod
  use seq_map_mod     , only: seq_map_map
  use seq_map_type_mod, only: seq_map

  use component_type_mod

  implicit none
  private ! except
#include <mpif.h>
  save

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: seq_domain_check_moab
  public :: seq_domain_compare
  public :: seq_domain_areafactinit

  !--------------------------------------------------------------------------
  ! Public variables
  !--------------------------------------------------------------------------

  real(R8), parameter :: eps_tiny   = 1.0e-16_R8 ! roundoff eps
  real(R8), parameter :: eps_big    = 1.0e+02_R8 ! big eps
  real(R8), parameter :: eps_frac_samegrid = 1.0e-9_R8 ! epsilon for fractions for samegrid

  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------

  private :: seq_domain_check_grid
  private :: seq_domain_check_fracmask_moab

  !================================================================================
contains
  !================================================================================

  !===============================================================================
  !================================================================================
  ! MOAB-based domain checking subroutines
  !================================================================================

  !================================================================================
  ! subroutine seq_domain_check_moab
  !
  ! Purpose: Comprehensive domain validation using MOAB data structures
  !
  ! This routine validates grid consistency across all climate model components
  ! using MOAB tag data. It checks
  ! 1. Mask/fraction compatibility within each component
  ! NOT YET FOR
  ! 2. Coordinate consistency between components sharing grids
  ! 3. Fraction sum validation (land+ice+ocean = 1.0) on atmosphere grid
  !
  ! The validation ensures that:
  ! - Grid coordinates match within tolerance for same-grid components
  ! - Component fractions are physically consistent
  ! - Mask and fraction arrays are compatible
  ! - Domain mappings are correctly initialized

  subroutine seq_domain_check_moab( infodata, &
       samegrid_al, samegrid_ao, samegrid_ro, samegrid_lg)

    !-----------------------------------------------------------
    ! MOAB version of seq_domain_check
    ! Checks domain consistency using MOAB tag data instead of MCT data structures
    !
    use seq_comm_mct, only: mbaxid, mblxid, mboxid, mbixid, mbrxid
    use shr_moab_mod, only: mbGetnCells, mbGetCellTagVals
    use iso_c_binding, only: C_NULL_CHAR
    !
    ! Arguments
    !
    type (seq_infodata_type) , intent(inout) :: infodata
    logical                  , intent(in)    :: samegrid_al ! atm lnd grid same
    logical                  , intent(in)    :: samegrid_ao ! atm ocn grid same
    logical                  , intent(in)    :: samegrid_ro ! rof ocn grid same
    logical                  , intent(in)    :: samegrid_lg ! lnd glc grid same
    !
    ! Local variables
    !
    integer(IN) :: n            ! indicies
    integer(IN) :: mpicom_cplid
    !
    logical      :: atm_present              ! atm present flag
    logical      :: lnd_present              ! lnd present flag
    logical      :: ocn_present              ! ocn present flag
    logical      :: ice_present              ! ice present flag
    logical      :: glc_present              ! glc present flag
    logical      :: rof_present              ! rof present flag
    logical      :: ocnrof_prognostic        ! ocn rof prognostic flag
    integer(IN)  :: rcode                    ! error status
    integer(IN)  :: atmsize                  ! local  size of atm  grid
    integer(IN)  :: lndsize                  ! local  size of land grid
    integer(IN)  :: ocnsize                  ! local  size of ocn  grid
    integer(IN)  :: icesize                  ! local  size of ice  grid
    real(R8)     :: diff,dmaxo,dmaxi         ! difference tracker
    logical      :: iamroot                  ! local masterproc
    real(R8)     :: eps_frac                 ! epsilon for fractions
    real(R8)     :: eps_axmask               ! epsilon for masks, atm/lnd
    real(R8)     :: eps_axgrid               ! epsilon for grid coords, atm/lnd
    real(R8)     :: eps_axarea               ! epsilon for areas, atm/lnd
    real(R8)     :: eps_oimask               ! epsilon for masks, ocn/ice
    real(R8)     :: eps_oigrid               ! epsilon for grid coords, ocn/ice
    real(R8)     :: eps_oiarea               ! epsilon for areas, ocn/ice
    real(R8)     :: my_eps_frac              ! local eps_frac value
    !
    real(R8),allocatable :: fracl(:)         ! land fraction
    real(R8),allocatable :: fraco(:)         ! ocn  fraction
    real(R8),allocatable :: fraci(:)         ! ice  fraction
    real(R8),allocatable :: maskl(:)         ! land mask
    real(R8),allocatable :: maski(:)         ! ice  mask
    real(R8),allocatable :: masko(:)         ! ocn  mask
    !
    character(*),parameter :: F00 = "('(seq_domain_check_moab) ',4a)"
    character(*),parameter :: F01 = "('(seq_domain_check_moab) ',a,i6,a)"
    character(*),parameter :: F02 = "('(seq_domain_check_moab) ',a,g23.15)"
    character(*),parameter :: F0R = "('(seq_domain_check_moab) ',2A,2g23.15,A )"
    character(*),parameter :: subName = '(seq_domain_check_moab) '
    !-----------------------------------------------------------

    !--------------------------------------------------------------------------
    ! Section 1: Get coupler communication info and configuration parameters
    !
    ! Retrieve MPI communicator and root process flag for the coupler.
    ! Also extract component presence flags and tolerance parameters from
    ! infodata. The epsilon values control how strict the grid comparisons
    ! are for different attributes (mask, coordinates, areas).
    !--------------------------------------------------------------------------

    call seq_comm_setptrs(CPLID,iamroot=iamroot, mpicom=mpicom_cplid)

    call seq_infodata_GetData( infodata,      &
         lnd_present=lnd_present,             &
         ocn_present=ocn_present,             &
         ice_present=ice_present,             &
         glc_present=glc_present,             &
         atm_present=atm_present,             &
         rof_present=rof_present,             &
         ocnrof_prognostic=ocnrof_prognostic, &
         eps_frac=eps_frac,                   &
         eps_amask=eps_axmask,                &
         eps_agrid=eps_axgrid,                &
         eps_aarea=eps_axarea,                &
         eps_omask=eps_oimask,                &
         eps_ogrid=eps_oigrid,                &
         eps_oarea=eps_oiarea )

    ! Get sizes from MOAB apps
    if (atm_present .and. mbaxid >= 0) then
       atmsize = mbGetnCells(mbaxid)
    else
       atmsize = 0
    endif

    if (lnd_present .and. mblxid >= 0) then
       lndsize = mbGetnCells(mblxid)
    else
       lndsize = 0
    endif

    if (ocn_present .and. mboxid >= 0) then
       ocnsize = mbGetnCells(mboxid)
    else
       ocnsize = 0
    endif

    if (ice_present .and. mbixid >= 0) then
       icesize = mbGetnCells(mbixid)
    else
       icesize = 0
    endif

    !--------------------------------------------------------------------------
    ! Section 2: Validate mask/fraction consistency
    !
    ! Check mask/fraction consistency using MOAB tags (frac>0 requires mask>0)
    !--------------------------------------------------------------------------

    ! Check land domain
    if (atm_present .and. lnd_present .and. mbaxid >= 0 .and. mblxid >= 0) then
       if (iamroot) write(logunit,F00) ' --- checking land mask and frac (MOAB) ---'
       call seq_domain_check_fracmask_moab(mblxid)
    endif

    ! Check ocean domain
    if (atm_present .and. ocn_present .and. mbaxid >= 0 .and. mboxid >= 0) then
       if (iamroot) write(logunit,F00) ' --- checking ocean mask and frac (MOAB) ---'
       call seq_domain_check_fracmask_moab(mboxid)
    endif

    ! Check ice domain
    if (atm_present .and. ice_present .and. mbaxid >= 0 .and. mbixid >= 0) then
       if (iamroot) write(logunit,F00) ' --- checking ice maskand frac (MOAB) ---'
       call seq_domain_check_fracmask_moab(mbixid)
    endif

    !--------------------------------------------------------------------------
    !
    ! MOAB TODO: when samegrid_xy is true, check that the global sizes of the grids are the same.
    !
    ! MOAB TODO: check_grid
    !
    !  The MCT version of this first declared local Avs and mapped DOMAIN vars to those
    !  local Avs, specifcally: mask, lat, lon, area.  It mapped most grids to atm.
    !
    !  IF SAMEGRID for a  pair is true, the seq_domain_check_grid call would then check if the native value and
    !  copied or rearranged values from the other grid were the same within a tolerance.  This was done in parallel and
    !  the total number of diffs was gathered to node 0.  If more then 0, an error resulted.
    !
    !  THIS WON'T WORK IN MOAB because the way seq_map_map is implemented, the mesh values in the target app
    !  are REPLACED by the mapped versions (for both rearrange and actual map).  We don't want that.
    !  Also the mapping to atm happens even if samegrid is not true which would replace the real values with mapped values.
    !  To implement this correctly, would need to copy the domain vars to new tags which can be mapped
    !  without overwriting the real domain data and then compared.
    !
    !  The domain variables mapped are
    ! if (atm_present .and. lnd_present) l2a mappping done
    ! if (atm_present .and. ocn_present) o2a mapping done
    ! if (atm_present .and. ice_present) i2a mapping done
    ! if (atm_present .and. iac_present) z2a mapping done
    ! if (ice_present .and. ocn_present) i2o mapping done
    !
    ! checking and additional mapping:
    ! if( samegrid_lg)
    !   l2g mapping done, mask, lat, lon, area checked
    ! if (ocn_present .and. ice_present)  (assumes they are sme grid)
    !   using ice mapped to ocn, check mask, lat, lon, area
    ! if (atm_present .and. lnd_present .and. samegrid_al
    !   using lnd mapped to atm, check lat, lon, area
    ! if (atm_present .and. iac_present .and. samegrid_az)
    !   using iac mapped to atm, check lat, lon, area
    ! if (atm_present .and. ice_present .and. samegrid_ao) 
    !   using ice mapped to atm, check lat, lon, area
    ! if (atm_present .and. ocn_present .and. samegrid_ao)
    !   using ocn mapped to atm, check lat, lon, area
    !
    ! MOAB TODO: check fractions
    !
    !  After that checking, the values mapped to atm are used to check fraction consistency
    !  with fracl = frac from land-mapped-to-atm
    !       fraci = mask from seaice-mapped-to-atm
    !       fraco = mask from ocean-mapped-to-atm
    !  looping over local points
    !    if (lnd_present .and. ice_present) then
    !       (1. - fracl - fraci) is checked to be within my_eps_frac of 0.
    !        1 - fraci must be > eps_frac AND  fracl < eps_tiny
    !    if (lnd_present .and. ocn_present)
    !       (1. - fracl(n) - fraco(n) is checked to be within my_eps_frac of 0.
    !        1. - fraco(n)) > eps_frac AND  fracl(n) < eps_tiny
    !    error out if any of the above are violated.  Report the max difference otherwise
    !      (note the MCT version of max diff reporting does not do a gather)
    !  

    call shr_sys_flush(logunit)

  end subroutine seq_domain_check_moab

  !===============================================================================

  subroutine seq_domain_compare(dom1, dom2, mpicom, eps)

    !-----------------------------------------------------------

    ! Arguments

    type(mct_gGrid)  , intent(in) :: dom1
    type(mct_gGrid)  , intent(in) :: dom2
    integer(IN)      , intent(in) :: mpicom
    real(R8),optional, intent(in) :: eps    ! error condition for compare

    ! Local variables
    real(R8) :: leps
    character(*),parameter :: F00 = "('(seq_domain_compare) ',4a)"
    character(*),parameter :: F01 = "('(seq_domain_compare) ',a,i12,a)"
    character(*),parameter :: F02 = "('(seq_domain_compare) ',2a,g23.15)"
    character(*),parameter :: F0R = "('(seq_domain_compare) ',2A,2g23.15,A )"
    character(*),parameter :: subName = '(seq_domain_compare) '

    leps = eps_tiny
    if (present(eps)) then
       leps = eps
    endif

    call seq_domain_check_grid(dom1%data, dom2%data, 'mask', eps=leps, mpicom=mpicom)
    call seq_domain_check_grid(dom1%data, dom2%data, 'lat' , eps=leps, mpicom=mpicom)
    call seq_domain_check_grid(dom1%data, dom2%data, 'lon' , eps=leps, mpicom=mpicom)
    call seq_domain_check_grid(dom1%data, dom2%data, 'area', eps=leps, mpicom=mpicom)

  end subroutine seq_domain_compare

  !===============================================================================

  subroutine seq_domain_check_grid(dom1, dom2, attr, eps, mpicom, mask)

    !-----------------------------------------------------------

    ! Arguments

    type(mct_aVect) , intent(in) :: dom1
    type(mct_aVect) , intent(in) :: dom2
    character(len=*), intent(in) :: attr   ! grid attribute to compare
    real(R8)        , intent(in) :: eps    ! error condition for compare
    integer(IN)     , intent(in) :: mpicom
    real(R8)        , intent(in), optional :: mask(:)

    ! Local variables

    integer(in)       :: n,ndiff            ! indices
    integer(in)       :: npts1,npts2,npts   ! counters
    integer(in)       :: rcode              ! error code
    real(R8)          :: diff,max_diff      ! temporaries
    real(R8)          :: tot_diff           ! maximum diff across all pes
    integer(IN)       :: ier                ! error code
    real(R8), pointer :: data1(:)           ! temporaries
    real(R8), pointer :: data2(:)           ! temporaries
    real(R8), pointer :: lmask(:)           ! temporaries
    logical           :: iamroot            ! local masterproc

    character(*),parameter :: F00 = "('(seq_domain_check_grid) ',4a)"
    character(*),parameter :: F01 = "('(seq_domain_check_grid) ',a,i12,a)"
    character(*),parameter :: F02 = "('(seq_domain_check_grid) ',2a,g23.15)"
    character(*),parameter :: F0R = "('(seq_domain_check_grid) ',2A,2g23.15,A )"
    character(*),parameter :: subName = '(seq_domain_check_grid) '
    !-----------------------------------------------------------

    call seq_comm_setptrs(CPLID,iamroot=iamroot)

    npts1 = mct_aVect_lsize(dom1)
    npts2 = mct_aVect_lsize(dom2)
    npts  = npts1

    if (npts1 == npts2) then
       if (iamroot) write(logunit,F01) " the domain size is = ", npts
    else
       write(logunit,*) trim(subname)," domain size #1 = ", npts1
       write(logunit,*) trim(subname)," domain size #2 = ", npts2
       write(logunit,*) trim(subname)," ERROR: domain size mis-match"
       call shr_sys_abort(subName//" ERROR: domain size mis-match")
    end if

    allocate(data1(npts),stat=rcode)
    if(rcode /= 0) call shr_sys_abort(subname//' allocate data1')
    allocate(data2(npts),stat=rcode)
    if(rcode /= 0) call shr_sys_abort(subname//' allocate data2')
    allocate(lmask(npts),stat=rcode)
    if(rcode /= 0) call shr_sys_abort(subname//' allocate lmask')

    call mct_aVect_exportRAttr(dom1, trim(attr), data1, npts)
    call mct_aVect_exportRAttr(dom2, trim(attr), data2, npts)
    lmask = 1.0_R8
    if (present(mask)) then
       if (size(mask) /= npts) then
          call shr_sys_abort(subName//" ERROR: mask size mis-match")
       endif
       lmask = mask
    endif

    ! --- adjust lons to address wraparound issues, we're assuming degree here! ---

    if (trim(attr) == "lon") then
       do n = 1,npts
          if (data2(n) > data1(n)) then
             do while ( (data1(n)+360.0_R8) < (data2(n)+180.0_R8) ) ! longitude is periodic
                data1(n) = data1(n) + 360.0_R8
             end do
          else
             do while ( (data2(n)+360.0_R8) < (data1(n)+180.0_R8) ) ! longitude is periodic
                data2(n) = data2(n) + 360.0_R8
             end do
          endif
       enddo
    endif

    ! Only check consistency where mask is greater than zero, if mask is present

    max_diff = 0.0_R8
    ndiff = 0
    do n=1,npts
       if (lmask(n) > eps_tiny) then
          diff = abs(data1(n)-data2(n))
          max_diff = max(max_diff,diff)
          if (diff > eps) then
             write(logunit,150) n,data1(n),data2(n),diff,eps
             ndiff = ndiff + 1
          endif
       end if
    end do
150 format('seq_domain_check_grid - n:',I3,' d1:',F12.6,' d2:',F12.6,' diff:',F18.14,' eps:',F18.14)

    call mpi_reduce(max_diff,tot_diff,1,MPI_REAL8,MPI_MAX,0,mpicom,ier)
    if (iamroot) then
       write(logunit,F02) " maximum           difference for ",trim(attr),tot_diff
       write(logunit,F02) " maximum allowable difference for ",trim(attr),eps
       call shr_sys_flush(logunit)
    endif
    call mpi_barrier(mpicom,ier)

    if (ndiff > 0) then
       write(logunit,*) trim(subname)," ERROR: incompatible domain grid coordinates"
       call shr_sys_flush(logunit)
       call shr_sys_abort(subName//" incompatible domain grid coordinates")
    endif

    deallocate(data1,stat=rcode)
    if(rcode /= 0) call shr_sys_abort(subname//' deallocate data1')
    deallocate(data2,stat=rcode)
    if(rcode /= 0) call shr_sys_abort(subname//' deallocate data2')
    deallocate(lmask,stat=rcode)
    if(rcode /= 0) call shr_sys_abort(subname//' deallocate lmask')

  end subroutine seq_domain_check_grid

  !===============================================================================

  subroutine seq_domain_areafactinit(domain, mdl2drv, drv2mdl, &
       samegrid, mpicom, iamroot, comment)
    !-----------------------------------------------------------
    !
    ! Arguments
    !
    type(mct_gGrid)  , pointer             :: domain     ! component domain on component pes
    real(R8)         , pointer             :: mdl2drv(:) ! comp->cpl factor on component pes
    real(R8)         , pointer             :: drv2mdl(:) ! cpl->comp factor on component pes
    logical          , intent(in)          :: samegrid   ! true => two grids are same
    integer          , intent(in)          :: mpicom     ! mpi communicator on component pes
    logical          , intent(in)          :: iamroot
    character(len=*) , optional,intent(in) :: comment
    !
    ! Local variables
    !
    integer                :: j1,j2,m1,n,rcode
    integer                :: gridsize
    real(R8)               :: rmin1,rmax1,rmin,rmax
    real(R8)               :: rmask,rarea,raream
    character(cl)          :: lcomment
    character(len=*),parameter :: subName = '(seq_domain_areafactinit) '
    character(len=*),parameter :: F0R = "(2A,2g25.17,A )"
    !
    !-----------------------------------------------------------

    lcomment = ''
    if (present(comment)) lcomment = comment

    ! get sizes

    gridsize = mct_gGrid_lsize(domain)
    allocate(drv2mdl(gridsize),mdl2drv(gridsize),stat=rcode)
    if(rcode /= 0) call shr_sys_abort(subname//' allocate area correction factors')

    j1 = mct_gGrid_indexRA(domain,"area"    ,dieWith=subName)
    j2 = mct_gGrid_indexRA(domain,"aream"   ,dieWith=subName)
    m1 = mct_gGrid_indexRA(domain,"mask"    ,dieWith=subName)

    mdl2drv(:)=1.0_R8
    drv2mdl(:)=1.0_R8

    if (samegrid) then
       ! default 1.0
    else
       do n=1,gridsize
          rmask  = domain%data%rAttr(m1,n)
          rarea  = domain%data%rAttr(j1,n)
          raream = domain%data%rAttr(j2,n)
          if ( abs(rmask) >= 1.0e-06) then
             if (rarea * raream /= 0.0_R8) then
                mdl2drv(n) = rarea/raream
                drv2mdl(n) = 1.0_R8/mdl2drv(n)
                !if (mdl2drv(n) > 10.0 .or. mdl2drv(n) < 0.1) then
                !   write(logunit,*) trim(subname),' WARNING area,aream= ', &
                !      domain%data%rAttr(j1,n),domain%data%rAttr(j2,n),' in ',n,gridsize
                !endif
             else
                write(logunit,*) trim(subname),' ERROR area,aream= ', &
                     rarea,raream,' in ',n,gridsize
                call shr_sys_flush(logunit)
                call shr_sys_abort()
             endif
          endif
       enddo
    end if

    rmin1 = minval(mdl2drv)
    rmax1 = maxval(mdl2drv)
    call shr_mpi_min(rmin1,rmin,mpicom)
    call shr_mpi_max(rmax1,rmax,mpicom)
    if (iamroot) write(logunit,F0R) trim(subname),' : min/max mdl2drv ',rmin,rmax,trim(lcomment)

    rmin1 = minval(drv2mdl)
    rmax1 = maxval(drv2mdl)
    call shr_mpi_min(rmin1,rmin,mpicom)
    call shr_mpi_max(rmax1,rmax,mpicom)
    if (iamroot) write(logunit,F0R) trim(subname),' : min/max drv2mdl ',rmin,rmax,trim(lcomment)
    if (iamroot) call shr_sys_flush(logunit)

  end subroutine seq_domain_areafactinit

  !===============================================================================

  subroutine seq_domain_check_fracmask_moab(mbid)

    !-----------------------------------------------------------
    ! MOAB version of seq_domain_check_fracmask
    ! Checks that fraction and mask are consistent in MOAB mesh
    !
    use shr_moab_mod, only: mbGetnCells, mbGetCellTagVals
    use iso_c_binding, only: C_NULL_CHAR
    !
    ! Arguments
    !
    integer(IN), intent(in) :: mbid   ! MOAB app id

    ! Local variables
    integer(in) :: n,npts,ndiff
    integer(in) :: rcode
    real(R8), allocatable :: dmask(:)           ! temporaries
    real(R8), allocatable :: dfrac(:)           ! temporaries

    character(*),parameter :: F00 = "('(seq_domain_check_fracmask_moab) ',4a)"
    character(*),parameter :: F01 = "('(seq_domain_check_fracmask_moab) ',a,i12,a)"
    character(*),parameter :: F02 = "('(seq_domain_check_fracmask_moab) ',2a,g23.15)"
    character(*),parameter :: F0R = "('(seq_domain_check_fracmask_moab) ',2A,2g23.15,A )"
    character(*),parameter :: subName = '(seq_domain_check_fracmask_moab) '
    !-----------------------------------------------------------

    if (mbid < 0) return

    npts = mbGetnCells(mbid)
    if (npts <= 0) return

    allocate(dmask(npts),stat=rcode)
    if(rcode /= 0) call shr_sys_abort(subname//' allocate dmask')
    allocate(dfrac(npts),stat=rcode)
    if(rcode /= 0) call shr_sys_abort(subname//' allocate dfrac')

    call mbGetCellTagVals(mbid, 'mask'//C_NULL_CHAR, dmask, npts)
    call mbGetCellTagVals(mbid, 'frac'//C_NULL_CHAR, dfrac, npts)

    ndiff = 0
    do n = 1,npts
       if (abs(dfrac(n)) > eps_tiny .and. abs(dmask(n)) < eps_tiny) then
          ndiff = ndiff + 1
       endif
    enddo

    if (ndiff > 0) then
       write(logunit,*) trim(subname)," ERROR: incompatible domain mask and frac values"
       call shr_sys_flush(logunit)
       call shr_sys_abort(subName//" incompatible domain mask and frac values")
    endif

    deallocate(dmask,stat=rcode)
    if(rcode /= 0) call shr_sys_abort(subname//' deallocate dmask')
    deallocate(dfrac,stat=rcode)
    if(rcode /= 0) call shr_sys_abort(subname//' deallocate dfrac')

  end subroutine seq_domain_check_fracmask_moab

  !===============================================================================

end module seq_domain_mct
