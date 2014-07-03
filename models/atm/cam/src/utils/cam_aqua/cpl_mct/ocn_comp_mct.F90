module ocn_comp_mct

  use mct_mod
  use esmf
  use seq_flds_mod
  use seq_cdata_mod
  use seq_infodata_mod
  use seq_timemgr_mod
  use perf_mod

  use shr_kind_mod,      only: r8 => shr_kind_r8, shr_kind_cl
  use physconst,         only: pi
  use shr_sys_mod,       only: shr_sys_abort, shr_sys_flush

  use phys_grid,         only: get_ncols_p,get_rlat_all_p,get_rlon_all_p, &
                               get_area_all_p,ngcols, get_gcol_p, get_nlcols_p
  use ppgrid,            only: pcols,begchunk,endchunk
  use ocn_types,         only: ocn_out_t
  use ocn_comp
  use spmd_utils,        only : iam
  implicit none
  private
  save

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public ocn_init_mct         ! Initialization method
  public ocn_run_mct          ! Run method
  public ocn_final_mct        ! Finalization method

!--------------------------------------------------------------------------
! Private interfaces
!--------------------------------------------------------------------------

  private :: ocn_SetgsMap_mct
  private :: ocn_export_mct

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

  type(ocn_out_t), pointer :: ocn_out(:) 

  integer :: index_o2x_So_t


!===============================================================
contains
!===============================================================

  subroutine ocn_init_mct( EClock, cdata_o, x2o_o, o2x_o, NLFilename )

    !----------------------------------------------------------
    !
    ! Arguments
    !
    type(ESMF_Clock),             intent(in)    :: EClock
    type(seq_cdata),              intent(inout) :: cdata_o
    type(mct_aVect),              intent(inout) :: x2o_o, o2x_o
    character(len=*), optional,   intent(in)    :: NLFilename 
    !
    ! Local variables
    !	
    integer                                     :: OCNID	
    integer                                     :: mpicom_ocn       	
    type(mct_gsMap),              pointer       :: gsMap_ocn
    type(mct_gGrid),              pointer       :: dom_o
    type(seq_infodata_type),      pointer       :: infodata   
    integer ::  lsize

    !----------------------------------------------------------
    !	
    ! Set cdata pointers
    !
    call seq_cdata_setptrs(cdata_o,ID=OCNID, mpicom=mpicom_ocn, &
         gsMap=gsMap_ocn, dom=dom_o, infodata=infodata)
    ! 
    ! Initialize ocn model
    !
    call ocn_init( ocn_out )
    !
    ! Initialize MCT gsMap, domain and attribute vectors
    !
    call ocn_SetgsMap_mct( mpicom_ocn, OCNID, gsMap_ocn ) 	
    lsize = mct_gsMap_lsize(gsMap_ocn, mpicom_ocn)
    !
    ! Initialize mct domain
    !
    call ocn_domain_mct( lsize, gsMap_ocn, dom_o )
    !
    ! Inialize mct attribute vectors
    !
    call mct_aVect_init(x2o_o, rList=seq_flds_x2o_fields, lsize=lsize)
    call mct_aVect_zero(x2o_o)

    call mct_aVect_init(o2x_o, rList=seq_flds_o2x_fields, lsize=lsize)
    index_o2x_So_t      = mct_avect_indexra(o2x_o,'So_t')
    call mct_aVect_zero(o2x_o)
    !
    ! Create initial ocn export state
    !
    call ocn_export_mct( ocn_out, o2x_o )
    call seq_infodata_PutData( infodata, ocn_prognostic=.false., ocnrof_prognostic=.false.)


  end subroutine ocn_init_mct
  
!==========================================================================

  subroutine ocn_run_mct ( EClock, cdata_o, x2o_o, o2x_o)

    !----------------------------------------------------------
    type(ESMF_Clock)    , intent(in)    :: EClock
    type(seq_cdata)     , intent(in)    :: cdata_o
    type(mct_aVect)     , intent(inout) :: x2o_o
    type(mct_aVect)     , intent(inout) :: o2x_o
    !----------------------------------------------------------

    ! Run ocean model

    call t_startf ('ocn_run')
    call ocn_run( ocn_out )
    call t_stopf ('ocn_run')

    ! Extract export state

    call t_startf ('ocn_export')
    call ocn_export_mct (ocn_out, o2x_o )
    call t_stopf ('ocn_export')
    
  end subroutine ocn_run_mct

!==========================================================================

  subroutine ocn_final_mct( EClock, cdata_o, x2o_o, o2x_o)

    type(ESMF_Clock)            , intent(in)    :: EClock
    type(seq_cdata)             , intent(inout) :: cdata_o
    type(mct_aVect)             , intent(inout) :: x2o_o
    type(mct_aVect)             , intent(inout) :: o2x_o

     ! *********************    
     ! does nothing 
     ! *********************    
  end subroutine ocn_final_mct
  
!==========================================================================

  subroutine ocn_export_mct( ocn_out, o2x_o ) 

    !-------------------------------------------------------------------
    implicit none
    type(ocn_out_t), intent(in)    :: ocn_out(begchunk:endchunk) 
    type(mct_aVect), intent(inout) :: o2x_o

    integer :: i,c,ig        ! indices
    integer :: ncols         ! number of columns
    !-----------------------------------------------------------------------

    ig=1
    do c=begchunk, endchunk
       ncols = get_ncols_p(c)
       do i=1,ncols
          o2x_o%rAttr(index_o2x_So_t,ig) = ocn_out(c)%ts(i)      
          ig=ig+1
       end do
    end do

  end subroutine ocn_export_mct

!==========================================================================

  subroutine ocn_SetgsMap_mct( mpicom_ocn, OCNID, gsMap_ocn )
   
    !-------------------------------------------------------------------
    !
    ! Arguments
    !
    integer        , intent(in)  :: mpicom_ocn
    integer        , intent(in)  :: OCNID
    type(mct_gsMap), intent(out) :: gsMap_ocn
    !
    ! Local Variables
    !
    integer, allocatable :: gindex(:)
    integer :: i, startpoint, j, sizebuf, n, c, ncols, gcol, nlcols
    !-------------------------------------------------------------------

    ! Determine global seg map

    sizebuf=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       do i = 1,ncols
          sizebuf = sizebuf+1
       end do
    end do

    allocate(gindex(sizebuf))

    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       do i = 1,ncols
          n=n+1
          gindex(n) = get_gcol_p(c,i)
       end do
    end do
    nlcols = get_nlcols_p()
    call mct_gsMap_init( gsMap_ocn, gindex, mpicom_ocn, OCNID, nlcols, ngcols )

    deallocate(gindex)

  end subroutine ocn_SetgsMap_mct
  
!===============================================================================

  subroutine ocn_domain_mct( lsize, gsMap_o, dom_o )

    !-------------------------------------------------------------------
    !
    ! Arguments
    !
    integer        , intent(in)    :: lsize
    type(mct_gsMap), intent(inout) :: gsMap_o
    type(mct_ggrid), intent(inout) :: dom_o     
    !
    ! Local Variables
    !
    integer  :: n,j,i,c,ncols         ! indices	
    real(r8) :: lats(pcols)           ! array of global latitude indices
    real(r8) :: lons(pcols)           ! array of global longitude indices
    real(r8) :: area(pcols)           ! area in radians squared for each grid point
    real(r8), pointer :: data(:)      ! temporary
    integer , pointer :: idata(:)     ! temporary
    real(r8), parameter :: radtodeg = 180.0_r8/pi
    !-------------------------------------------------------------------
    !
    ! Initialize domain type
    !
    call mct_gGrid_init( GGrid=dom_o, CoordChars=trim(seq_flds_dom_coord), &
      OtherChars=trim(seq_flds_dom_other), lsize=lsize )
    !
    ! Allocate memory
    !
    allocate(data(lsize))
    !
    ! Determine global gridpoint number attribute, GlobGridNum, which is set automatically by MCT
    !
    call mct_gsMap_orderedPoints(gsMap_o, iam, idata)
    call mct_gGrid_importIAttr(dom_o,'GlobGridNum',idata,lsize)
    call mct_gGrid_importIAttr(dom_o,'GlobGridNum',idata,lsize)
    !
    ! Determine domain (numbering scheme is: West to East and South to North to South pole)
    ! Initialize attribute vector with special value
    !
    data(:) = -9999.0_R8 
    call mct_gGrid_importRAttr(dom_o,"lat"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_o,"lon"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_o,"area" ,data,lsize) 
    call mct_gGrid_importRAttr(dom_o,"aream",data,lsize) 
    data(:) = 0.0_R8     
    call mct_gGrid_importRAttr(dom_o,"mask" ,data,lsize) 
    call mct_gGrid_importRAttr(dom_o,"frac" ,data,lsize) 
    !
    ! Fill in correct values for domain components
    !
    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       call get_rlat_all_p(c, ncols, lats)
       do i=1,ncols
          n = n+1
          data(n) = lats(i)*radtodeg
       end do
    end do
    call mct_gGrid_importRAttr(dom_o,"lat",data,lsize) 

    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       call get_rlon_all_p(c, ncols, lons)
       do i=1,ncols
          n = n+1
          data(n) = lons(i)*radtodeg
       end do
    end do
    call mct_gGrid_importRAttr(dom_o,"lon",data,lsize) 

    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       call get_area_all_p(c, ncols, area)
       do i=1,ncols
          n = n+1
          data(n) = area(i) 
       end do
    end do
    call mct_gGrid_importRAttr(dom_o,"area",data,lsize) 

    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       do i=1,ncols
          n = n+1
          data(n) = 1._r8 ! mask for aqua_planet
       end do
    end do
    call mct_gGrid_importRAttr(dom_o,"mask",data,lsize) 

    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       do i=1,ncols
          n = n+1
          data(n) = 1._r8 ! frac for aqua_planet
       end do
    end do
    call mct_gGrid_importRAttr(dom_o,"frac",data,lsize) 

    deallocate(data)
    deallocate(idata)

  end subroutine ocn_domain_mct

end module ocn_comp_mct
