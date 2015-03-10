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

  public :: seq_domain_check
  public :: seq_domain_compare
  public :: seq_domain_areafactinit

!--------------------------------------------------------------------------
! Public variables
!--------------------------------------------------------------------------

  real(R8), parameter :: eps_tiny   = 1.0e-16_R8 ! roundoff eps
  real(R8), parameter :: eps_big    = 1.0e+02_R8 ! big eps
  real(R8), parameter :: eps_frac_samegrid = 1.0e-14_R8 ! epsilon for fractions for samegrid

!--------------------------------------------------------------------------
! Private interfaces
!--------------------------------------------------------------------------

  private :: seq_domain_check_grid

!================================================================================
contains
!================================================================================

!================================================================================

  subroutine seq_domain_check( infodata, &
       atm, ice, lnd, ocn, rof, glc, &
       samegrid_al, samegrid_ao, samegrid_ro)

    !-----------------------------------------------------------
    ! Uses
    !
    use prep_atm_mod, only: prep_atm_get_mapper_Fi2a
    use prep_atm_mod, only: prep_atm_get_mapper_Fl2a
    use prep_atm_mod, only: prep_atm_get_mapper_Fo2a
    use prep_lnd_mod, only: prep_lnd_get_mapper_Fa2l
    use prep_ocn_mod, only: prep_ocn_get_mapper_SFi2o
    use prep_glc_mod, only: prep_glc_get_mapper_SFl2g
    !
    ! Arguments
    !
    type (seq_infodata_type) , intent(inout) :: infodata
    type(component_type)     , intent(in)    :: atm
    type(component_type)     , intent(in)    :: ice
    type(component_type)     , intent(in)    :: lnd
    type(component_type)     , intent(in)    :: ocn
    type(component_type)     , intent(in)    :: rof
    type(component_type)     , intent(in)    :: glc
    logical                  , intent(in)    :: samegrid_al ! atm lnd grid same
    logical                  , intent(in)    :: samegrid_ao ! atm ocn grid same
    logical                  , intent(in)    :: samegrid_ro ! rof ocn grid same
    !
    ! Local variables
    !
    type(seq_map)   , pointer :: mapper_i2a ! inout needed for lower methods
    type(seq_map)   , pointer :: mapper_i2o ! inout needed for lower methods
    type(seq_map)   , pointer :: mapper_o2a !
    type(seq_map)   , pointer :: mapper_l2g !
    type(seq_map)   , pointer :: mapper_a2l !
    type(seq_map)   , pointer :: mapper_l2a !
                                            !
    type(mct_gGrid) , pointer :: atmdom_a   ! atm domain
    type(mct_gGrid) , pointer :: icedom_i   ! ice domain
    type(mct_gGrid) , pointer :: lnddom_l   ! lnd domain
    type(mct_gGrid) , pointer :: ocndom_o   ! ocn domain
    type(mct_gGrid) , pointer :: glcdom_g   ! glc domain
                                            !
    type(mct_gsMap) , pointer :: gsMap_a    ! atm global seg map 
    type(mct_gsMap) , pointer :: gsMap_i    ! ice global seg map 
    type(mct_gsMap) , pointer :: gsMap_l    ! lnd global seg map 
    type(mct_gsMap) , pointer :: gsMap_o    ! ocn global seg map 
    type(mct_gsMap) , pointer :: gsMap_r    ! ocn global seg map 
    type(mct_gsMap) , pointer :: gsMap_g    ! glc global seg map 
    !
    type(mct_gGrid) :: lnddom_a              ! lnd domain info on atm decomp
    type(mct_gGrid) :: lnddom_g              ! lnd domain info on glc decomp
    type(mct_gGrid) :: icedom_a              ! ice domain info on atm decomp (all grids same)
    type(mct_gGrid) :: ocndom_a              ! ocn domain info on atm decomp (all grids same)
    type(mct_gGrid) :: icedom_o              ! ocn domain info on ocn decomp (atm/ocn grid different)
    !
    real(R8), pointer :: fracl(:)            ! land fraction on atm decomp 
    real(R8), pointer :: fraco(:)            ! ocn  fraction on atm decomp 
    real(R8), pointer :: fraci(:)            ! ice  fraction on atm decomp 
    real(R8), pointer :: maskl(:)            ! land mask on atm decomp (all grids same)
    real(R8), pointer :: maski(:)            ! ice  mask on atm decomp (all grids same)
    real(R8), pointer :: masko(:)            ! ocn  mask on atm decomp (all grids same)
    !
    integer(IN) :: n, kl, ko, ki             ! indicies
    integer(IN) :: k1,k2,k3                  ! indicies
    !
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
    integer(IN)  :: glcsize                  ! local  size of glc  grid
    integer(IN)  :: gatmsize                 ! global size of atm  grid
    integer(IN)  :: glndsize                 ! global size of land grid
    integer(IN)  :: gocnsize                 ! global size of ocn  grid
    integer(IN)  :: grofsize                 ! global size of ocn  grid
    integer(IN)  :: gicesize                 ! global size of ice  grid
    integer(IN)  :: gglcsize                 ! global size of glc  grid
    integer(IN)  :: npts                     ! local size temporary
    integer(IN)  :: ier                      ! error code
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
    real(R8)     :: rmin1,rmax1,rmin,rmax    ! local min max computation
    !
    real(R8),allocatable :: mask (:)         ! temporary real vector, domain mask
    !
    character(*),parameter :: F00 = "('(seq_domain_check) ',4a)"
    character(*),parameter :: F01 = "('(seq_domain_check) ',a,i6,a)"
    character(*),parameter :: F02 = "('(seq_domain_check) ',a,g23.15)"
    character(*),parameter :: F0R = "('(seq_domain_check) ',2A,2g23.15,A )"
    character(*),parameter :: subName = '(seq_domain_check) '
    !-----------------------------------------------------------

    mapper_i2a => prep_atm_get_mapper_Fi2a()
    mapper_i2o => prep_ocn_get_mapper_SFi2o()
    mapper_o2a => prep_atm_get_mapper_Fo2a()
    mapper_l2g => prep_glc_get_mapper_SFl2g()
    mapper_a2l => prep_lnd_get_mapper_Fa2l()
    mapper_l2a => prep_atm_get_mapper_Fl2a()

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

    ! Get info

    gsmap_a  => component_get_gsmap_cx(atm) ! gsmap_ax
    atmdom_a => component_get_dom_cx(atm)   ! dom_ax
    atmsize  = mct_avect_lsize(atmdom_a%data)
    gatmsize = mct_gsMap_gsize(gsMap_a)

    if (atm_present .and. lnd_present) then
       gsmap_l  => component_get_gsmap_cx(lnd) ! gsmap_lx
       lnddom_l => component_get_dom_cx(lnd)   ! dom_lx
       lndsize  = mct_avect_lsize(lnddom_l%data)
       glndsize = mct_gsMap_gsize(gsMap_l) 

       if (samegrid_al .and. gatmsize /= glndsize) then
          write(logunit,*) subname,' error: global atmsize = ',&
               gatmsize,' global lndsize= ',glndsize
          call shr_sys_flush(logunit)
          call shr_sys_abort(subname//' atm and lnd grid must have the same global size')
       end if
       if (iamroot) write(logunit,F00) ' --- checking land maskfrac ---'
       call seq_domain_check_fracmask(lnddom_l%data)
       call mct_gGrid_init(oGGrid=lnddom_a, iGGrid=lnddom_l, lsize=atmsize)
       call mct_aVect_zero(lnddom_a%data)
       call seq_map_map(mapper_l2a, lnddom_l%data, lnddom_a%data, norm=.false.)
       allocate(maskl(atmsize),stat=rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' allocate maskl')
       allocate(fracl(atmsize),stat=rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' allocate fracl')
       call mct_aVect_exportRAttr(lnddom_a%data, 'mask', maskl, atmsize)
       call mct_aVect_exportRAttr(lnddom_a%data, 'frac', fracl, atmsize)
    endif

    if (atm_present .and. ocn_present) then
       gsmap_o  => component_get_gsmap_cx(ocn) ! gsmap_ox
       ocndom_o => component_get_dom_cx(ocn)   ! dom_ox
       ocnsize  = mct_avect_lsize(ocndom_o%data)
       gocnsize = mct_gsMap_gsize(gsMap_o)

       if (samegrid_ao .and. gatmsize /= gocnsize) then
          write(logunit,*) subname,' error: global atmsize = ',gatmsize,' global ocnsize= ',gocnsize
          call shr_sys_flush(logunit)
          call shr_sys_abort(subname//' atm and ocn grid must have the same global size')
       end if
       if (iamroot) write(logunit,F00) ' --- checking ocean maskfrac ---'
       call seq_domain_check_fracmask(ocndom_o%data)
       call mct_gGrid_init(oGGrid=ocndom_a, iGGrid=ocndom_o, lsize=atmsize)
       call mct_aVect_zero(ocndom_a%data)
       call seq_map_map(mapper_o2a, ocndom_o%data, ocndom_a%data, norm=.false.)
       allocate(masko(atmsize),stat=rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' allocate masko')
       allocate(fraco(atmsize),stat=rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' allocate fraco')
       call mct_aVect_exportRAttr(ocndom_a%data, 'mask', masko, atmsize)
       if (samegrid_ao) then
          call mct_aVect_exportRattr(ocndom_a%data, 'frac', fraco, atmsize)
       else
          call mct_aVect_exportRattr(ocndom_a%data, 'mask', fraco, atmsize)
       endif
    endif
   
    if (atm_present .and. ice_present) then
       gsmap_i  => component_get_gsmap_cx(ice) ! gsmap_ix
       icedom_i => component_get_dom_cx(ice)   ! dom_ix
       icesize  = mct_avect_lsize(icedom_i%data)
       gicesize = mct_gsMap_gsize(gsMap_i) 

       if (samegrid_ao .and. gatmsize /= gicesize) then
          write(logunit,*) subname,' error: global atmsize = ',&
               gatmsize,' global icesize= ',gicesize
          call shr_sys_flush(logunit)
          call shr_sys_abort(subname//' atm and ice grid must have the same global size')
       end if
       if (iamroot) write(logunit,F00) ' --- checking ice maskfrac ---'
       call seq_domain_check_fracmask(icedom_i%data)
       call mct_gGrid_init(oGGrid=icedom_a, iGGrid=icedom_i, lsize=atmsize)
       call mct_aVect_zero(icedom_a%data)
       call seq_map_map(mapper_i2a, icedom_i%data, icedom_a%data, norm=.false.)
       allocate(maski(atmsize),stat=rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' allocate maski')
       allocate(fraci(atmsize),stat=rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' allocate fraci')
       call mct_aVect_exportRAttr(icedom_a%data, 'mask', maski, atmsize)
       if (samegrid_ao) then
          call mct_aVect_exportRattr(icedom_a%data, 'frac', fraci, atmsize)
       else
          call mct_aVect_exportRattr(icedom_a%data, 'mask', fraci, atmsize)
       endif
    endif

    if (lnd_present .and. glc_present) then
       gsmap_l  => component_get_gsmap_cx(lnd) ! gsmap_lx
       lnddom_l => component_get_dom_cx(lnd)   ! dom_lx
       lndsize  = mct_avect_lsize(lnddom_l%data)
       glndsize = mct_gsMap_gsize(gsMap_l) 

       gsmap_g  => component_get_gsmap_cx(glc) ! gsmap_gx
       glcdom_g => component_get_dom_cx(glc)   ! dom_gx
       glcsize  = mct_avect_lsize(glcdom_g%data)
       gglcsize = mct_gsMap_gsize(gsMap_g) 

       if (gglcsize /= glndsize) then
          write(logunit,*) subname,' error: global glcsize = ',gglcsize,' global lndsize= ',glndsize
          call shr_sys_flush(logunit)
          call shr_sys_abort(subname//' glc and lnd grid must have the same global size')
       end if
       if (iamroot) write(logunit,F00) ' --- checking glc maskfrac ---'
       call seq_domain_check_fracmask(glcdom_g%data)
       if (iamroot) write(logunit,F00) ' --- checking lnd maskfrac ---'
       call seq_domain_check_fracmask(lnddom_l%data)
       call mct_gGrid_init(oGGrid=lnddom_g, iGGrid=lnddom_l, lsize=glcsize)
       call mct_aVect_zero(lnddom_g%data)
       call seq_map_map(mapper_l2g, lnddom_l%data, lnddom_g%data, norm=.false.)
       if (iamroot) write(logunit,F00) ' --- checking glc/lnd domains ---'
       npts = glcsize
       allocate(mask(npts),stat=rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' allocate mask')
       call mct_aVect_getRAttr(lnddom_g%data,"mask",mask,rcode)
       where (mask < eps_axmask) mask = 0.0_R8
       call seq_domain_check_grid(glcdom_g%data, lnddom_g%data, 'mask', eps=eps_axmask, mpicom=mpicom_cplid, mask=mask)
       call seq_domain_check_grid(glcdom_g%data, lnddom_g%data, 'lat' , eps=eps_axgrid, mpicom=mpicom_cplid, mask=mask)
       call seq_domain_check_grid(glcdom_g%data, lnddom_g%data, 'lon' , eps=eps_axgrid, mpicom=mpicom_cplid, mask=mask)
       call seq_domain_check_grid(glcdom_g%data, lnddom_g%data, 'area', eps=eps_axarea, mpicom=mpicom_cplid, mask=mask)
       deallocate(mask,stat=rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' deallocate mask')
    endif

    if (ice_present .and. ocn_present) then
       gsmap_i  => component_get_gsmap_cx(ice) ! gsmap_ix
       icedom_i => component_get_dom_cx(ice)   ! dom_ix
       icesize  = mct_avect_lsize(icedom_i%data)
       gicesize = mct_gsMap_gsize(gsMap_i) 

       gsmap_o  => component_get_gsmap_cx(ocn) ! gsmap_ox
       ocndom_o => component_get_dom_cx(ocn)   ! dom_ox
       ocnsize  = mct_avect_lsize(ocndom_o%data)
       gocnsize = mct_gsMap_gsize(gsMap_o)

       if (gocnsize /= gicesize) then
          write(logunit,*) subname,' error: global ocnsize = ',gocnsize,' global icesize= ',gicesize
          call shr_sys_flush(logunit)
          call shr_sys_abort(subname//' ocean and ice grid must have the same global size')
       endif
       call mct_gGrid_init(oGGrid=icedom_o, iGGrid=icedom_i, lsize=ocnsize)
       call mct_aVect_zero(icedom_o%data)
       call seq_map_map(mapper_i2o, icedom_i%data, icedom_o%data, norm=.false.)
    end if

    if (rof_present .and. ocnrof_prognostic .and. samegrid_ro) then
       gsmap_r  => component_get_gsmap_cx(glc) ! gsmap_gx
       grofsize = mct_gsMap_gsize(gsMap_r)

       if (gocnsize /= grofsize) then
          write(logunit,*) subname,' error: global ocnsize = ',gocnsize,' global rofsize= ',grofsize
          call shr_sys_flush(logunit)
          call shr_sys_abort(subname//' ocean and rof grid must have the same global size')
       endif
    end if

    !------------------------------------------------------------------------------
    ! Check ice/ocean grid consistency
    !------------------------------------------------------------------------------

     if (ocn_present .and. ice_present) then
!    if (samegrid_oi) then       ! doesn't yet exist

       npts = ocnsize
       allocate(mask(npts),stat=rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' allocate mask')

       if (iamroot) write(logunit,F00) ' --- checking ocn/ice domains ---'
       call seq_domain_check_grid(ocndom_o%data, icedom_o%data,'mask', eps=eps_oigrid, mpicom=mpicom_cplid)
       call mct_aVect_getRAttr(ocndom_o%data,"mask",mask,rcode)
       where (mask < eps_oimask) mask = 0.0_R8

       call seq_domain_check_grid(ocndom_o%data, icedom_o%data,'lat' , eps=eps_oigrid, mpicom=mpicom_cplid, mask=mask)
       call seq_domain_check_grid(ocndom_o%data, icedom_o%data,'lon' , eps=eps_oigrid, mpicom=mpicom_cplid, mask=mask)
       call seq_domain_check_grid(ocndom_o%data, icedom_o%data,'area', eps=eps_oiarea, mpicom=mpicom_cplid, mask=mask)

       deallocate(mask,stat=rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' deallocate mask')

!    endif
     endif

    !------------------------------------------------------------------------------
    ! Check atm/lnd grid consistency
    !------------------------------------------------------------------------------

    if (atm_present .and. lnd_present .and. samegrid_al) then
       if (iamroot) write(logunit,F00) ' --- checking atm/land domains ---'
       call seq_domain_check_grid(atmdom_a%data, lnddom_a%data, 'lat' , eps=eps_axgrid, mpicom=mpicom_cplid, mask=maskl)
       call seq_domain_check_grid(atmdom_a%data, lnddom_a%data, 'lon' , eps=eps_axgrid, mpicom=mpicom_cplid, mask=maskl)
       call seq_domain_check_grid(atmdom_a%data, lnddom_a%data, 'area', eps=eps_axarea, mpicom=mpicom_cplid, mask=maskl)
    endif

    !------------------------------------------------------------------------------
    ! Check atm/ocn and atm/ice grid consistency (if samegrid)
    !------------------------------------------------------------------------------

    if (atm_present .and. ice_present .and. samegrid_ao) then
       if (iamroot) write(logunit,F00) ' --- checking atm/ice domains ---'
       call seq_domain_check_grid(atmdom_a%data, icedom_a%data, 'lat' , eps=eps_axgrid, mpicom=mpicom_cplid, mask=maski)
       call seq_domain_check_grid(atmdom_a%data, icedom_a%data, 'lon' , eps=eps_axgrid, mpicom=mpicom_cplid, mask=maski)
       call seq_domain_check_grid(atmdom_a%data, icedom_a%data, 'area', eps=eps_axarea, mpicom=mpicom_cplid, mask=maski)
    endif

    if (atm_present .and. ocn_present .and. samegrid_ao) then
       if (iamroot) write(logunit,F00) ' --- checking atm/ocn domains ---'
       call seq_domain_check_grid(atmdom_a%data, ocndom_a%data, 'lat' , eps=eps_axgrid, mpicom=mpicom_cplid, mask=masko)
       call seq_domain_check_grid(atmdom_a%data, ocndom_a%data, 'lon' , eps=eps_axgrid, mpicom=mpicom_cplid, mask=masko)
       call seq_domain_check_grid(atmdom_a%data, ocndom_a%data, 'area', eps=eps_axarea, mpicom=mpicom_cplid, mask=masko)
    endif

    !------------------------------------------------------------------------------
    ! Check consistency of land fraction with ocean mask on grid
    !------------------------------------------------------------------------------

    my_eps_frac = eps_frac
    if (samegrid_ao) my_eps_frac = eps_frac_samegrid
    if (.not. samegrid_al) my_eps_frac = eps_big

    if (iamroot) write(logunit,F00) ' --- checking fractions in domains ---'
    dmaxi = 0.0_R8
    dmaxo = 0.0_R8
    do n = 1,atmsize
       if (atm_present .and. lnd_present .and. ice_present) then
          diff = abs(1._R8 - fracl(n) - fraci(n))
          dmaxi = max(diff,dmaxi)
          if (diff > my_eps_frac) then
             write(logunit,*)'inconsistency between land fraction and sea ice fraction'
             write(logunit,*)'n= ',n,' fracl= ',fracl(n),' fraci= ',fraci(n),' sum= ',fracl(n)+fraci(n)
             call shr_sys_flush(logunit)
             call shr_sys_abort(subname//' inconsistency between land fraction and sea ice fraction')
          end if
          if ((1._R8-fraci(n)) > eps_frac .and. fracl(n) < eps_tiny) then
             write(logunit,*)'inconsistency between land mask and sea ice mask'
             write(logunit,*)'n= ',n,' fracl= ',fracl(n),' fraci= ',fraci(n)
             call shr_sys_flush(logunit)
             call shr_sys_abort(subname//'  inconsistency between land mask and sea ice mask')
          end if
       endif
       if (atm_present .and. lnd_present .and. ocn_present) then
          diff = abs(1._R8 - fracl(n) - fraco(n))
          dmaxo = max(diff,dmaxo)
          if (diff > my_eps_frac) then
             write(logunit,*)'inconsistency between land fraction and ocn land fraction'
             write(logunit,*)'n= ',n,' fracl= ',fracl(n),' fraco= ',fraco(n),' sum= ',fracl(n)+fraco(n)
             call shr_sys_flush(logunit)
             call shr_sys_abort(subname//'  inconsistency between land fraction and ocn land fraction')
          end if
          if ((1._R8-fraco(n)) > eps_frac .and. fracl(n) < eps_tiny) then
             write(logunit,*)'inconsistency between land mask and ocn land mask'
             write(logunit,*)'n= ',n,' fracl= ',fracl(n),' fraco= ',fraco(n)
             call shr_sys_flush(logunit)
             call shr_sys_abort(subname//'  inconsistency between land mask and ocn land mask')
          end if
       endif
    end do 
    if (iamroot) then
       write(logunit,F02) ' maximum           difference for ofrac sum ',dmaxo
       write(logunit,F02) ' maximum           difference for ifrac sum ',dmaxi
       write(logunit,F02) ' maximum allowable difference for  frac sum ',my_eps_frac
       write(logunit,F02) ' maximum allowable tolerance for valid frac ',eps_frac
       call shr_sys_flush(logunit)
    endif

    !------------------------------------------------------------------------------
    ! Clean up allocated memory
    !------------------------------------------------------------------------------

    if (atm_present .and. lnd_present) then
       deallocate(fracl,stat=rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' deallocate fracl')
       deallocate(maskl,stat=rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' deallocate maskl')
       call mct_gGrid_clean(lnddom_a, rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' clean lnddom_a')
    endif

    if (atm_present .and. ocn_present) then
       deallocate(fraco,stat=rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' deallocate fraco')
       deallocate(masko,stat=rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' deallocate masko')
       call mct_gGrid_clean(ocndom_a, rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' clean ocndom_a')
    endif

    if (atm_present .and. ice_present) then
       deallocate(fraci,stat=rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' deallocate fraci')
       deallocate(maski,stat=rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' deallocate maski')
       call mct_gGrid_clean(icedom_a, rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' clean icedom_o')
    endif

    if (ocn_present .and. ice_present) then
       call mct_gGrid_clean(icedom_o, rcode)
       if(rcode /= 0) call shr_sys_abort(subname//' clean icedom_o')
    endif

    call shr_sys_flush(logunit)

  end subroutine seq_domain_check

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
  
  subroutine seq_domain_check_fracmask(dom1)
   
    !-----------------------------------------------------------

    ! Arguments

    type(mct_aVect) , intent(in) :: dom1

    ! Local variables
    integer(in) :: n,npts,ndiff
    integer(in) :: rcode
    real(R8), pointer :: dmask(:)           ! temporaries
    real(R8), pointer :: dfrac(:)           ! temporaries

    character(*),parameter :: F00 = "('(seq_domain_check_fracmask) ',4a)"
    character(*),parameter :: F01 = "('(seq_domain_check_fracmask) ',a,i12,a)"
    character(*),parameter :: F02 = "('(seq_domain_check_fracmask) ',2a,g23.15)"
    character(*),parameter :: F0R = "('(seq_domain_check_fracmask) ',2A,2g23.15,A )"
    character(*),parameter :: subName = '(seq_domain_check_fracmask) '
    !-----------------------------------------------------------

    npts = mct_aVect_lsize(dom1)

    allocate(dmask(npts),stat=rcode)
    if(rcode /= 0) call shr_sys_abort(subname//' allocate dmask')
    allocate(dfrac(npts),stat=rcode)
    if(rcode /= 0) call shr_sys_abort(subname//' allocate dfrac')

    call mct_aVect_exportRAttr(dom1, 'mask', dmask, npts)
    call mct_aVect_exportRAttr(dom1, 'frac', dfrac, npts)

    ndiff = 0
    do n = 1,npts
       if (abs(dfrac(n)) > eps_tiny .and. abs(dmask(n)) < eps_tiny) then
!debug            write(logunit,*)'n= ',n,' dfrac= ',dfrac(n),' dmask= ',dmask(n)
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

 end subroutine seq_domain_check_fracmask

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
      !debug            write(logunit,*)'n= ',n,' data1= ',data1(n),' data2= ',data2(n),' diff= ',diff, ' eps= ',eps
	     ndiff = ndiff + 1
	  endif
       end if
    end do

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
    integer                :: gridsize,m2dsize,d2msize
    real(R8)               :: rmin1,rmax1,rmin,rmax
    real(R8)               :: rmask,rarea,raream
    character(cl)          :: lcomment
    character(len=*),parameter :: subName = '(seq_domain_areafactinit) '
    character(len=*),parameter :: F0R = "(2A,2g23.15,A )"
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

end module seq_domain_mct



