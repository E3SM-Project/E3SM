module seq_flux_mct
  
  use shr_kind_mod,      only: r8 => shr_kind_r8, in=>shr_kind_in
  use shr_sys_mod,       only: shr_sys_abort
  use shr_flux_mod,      only: shr_flux_atmocn
  use shr_orb_mod,       only: shr_orb_params, shr_orb_cosz, shr_orb_decl
  use shr_mct_mod

  use mct_mod
  use seq_flds_mod
  use seq_comm_mct
  use seq_cdata_mod
  use seq_infodata_mod

  implicit none
  private 	
  save

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public seq_flux_init_mct
  public seq_flux_initexch_mct
  public seq_flux_ocnalb_mct
  public seq_flux_atmocn_mct
  public seq_flux_atmocnexch_mct

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

  real(r8), pointer     :: lats(:)    ! latitudes  (degrees)
  real(r8), pointer     :: lons(:)    ! longitudes (degrees)
  integer(in),allocatable ::  mask(:) ! ocn domain mask: 0 <=> inactive cell
  integer(in),allocatable :: emask(:) ! ocn mask on exchange grid decomp

  real(r8), allocatable ::  uocn (:)  ! ocn velocity, zonal
  real(r8), allocatable ::  vocn (:)  ! ocn velocity, meridional
  real(r8), allocatable ::  tocn (:)  ! ocean temperature
  real(r8), allocatable ::  zbot (:)  ! atm level height
  real(r8), allocatable ::  ubot (:)  ! atm velocity, zonal     
  real(r8), allocatable ::  vbot (:)  ! atm velocity, meridional
  real(r8), allocatable ::  thbot(:)  ! atm potential T
  real(r8), allocatable ::  shum (:)  ! atm specific humidity
  real(r8), allocatable ::  dens (:)  ! atm density
  real(r8), allocatable ::  tbot (:)  ! atm bottom surface T
  real(r8), allocatable ::  sen  (:)  ! heat flux: sensible 
  real(r8), allocatable ::  lat  (:)  ! heat flux: latent   
  real(r8), allocatable ::  lwup (:)  ! lwup over ocean
  real(r8), allocatable ::  evap (:)  ! water flux: evaporation
  real(r8), allocatable ::  taux (:)  ! wind stress, zonal
  real(r8), allocatable ::  tauy (:)  ! wind stress, meridional
  real(r8), allocatable ::  tref (:)  ! diagnostic:  2m ref T
  real(r8), allocatable ::  qref (:)  ! diagnostic:  2m ref Q
  real(r8), allocatable :: duu10n(:)  ! diagnostic: 10m wind speed squared

  real(r8), allocatable ::  ustar(:)  ! saved ustar
  real(r8), allocatable ::  re   (:)  ! saved re
  real(r8), allocatable ::  ssq  (:)  ! saved sq

  ! Conversion from degrees to radians

  real(r8),parameter :: const_pi      = SHR_CONST_PI       ! pi
  real(r8),parameter :: const_deg2rad = const_pi/180.0_R8  ! deg to rads

  ! Coupler field indices

  integer :: index_a2x_Sa_z    
  integer :: index_a2x_Sa_u    
  integer :: index_a2x_Sa_v    
  integer :: index_a2x_Sa_tbot 
  integer :: index_a2x_Sa_ptem 
  integer :: index_a2x_Sa_shum 
  integer :: index_a2x_Sa_dens 
  integer :: index_o2x_So_t      
  integer :: index_o2x_So_u
  integer :: index_o2x_So_v
  integer :: index_xao_So_tref    
  integer :: index_xao_So_qref    
  integer :: index_xao_So_avsdr   
  integer :: index_xao_So_avsdf   
  integer :: index_xao_So_anidr   
  integer :: index_xao_So_anidf   
  integer :: index_xao_Faox_taux  
  integer :: index_xao_Faox_tauy   
  integer :: index_xao_Faox_lat   
  integer :: index_xao_Faox_sen   
  integer :: index_xao_Faox_evap  
  integer :: index_xao_Faox_lwup
  integer :: index_xao_So_ustar
  integer :: index_xao_So_re   
  integer :: index_xao_So_ssq  
  integer :: index_xao_So_duu10n 
  integer :: index_xao_So_u10

  character(len=16) :: fluxsetting = 'unknown'
  character(len=*),parameter  :: fluxsetting_atmocn = 'atmocn'
  character(len=*),parameter  :: fluxsetting_exchange = 'exchange'

  !--- for exchange grid ---
  type(mct_rearr) :: Re_a2e, Re_e2a, Re_o2e, Re_e2o  ! atm/ocn/exch rearrangers
  type(mct_sMat ) :: sMata2o, sMato2a                ! decomp sMat 
  type(mct_gsMap) :: gsmap_ae, gsmap_oe              ! gsmaps for atm/ocn on exch grid
  integer(in)     :: nloc_a2o,nloc_o2a,nloc_o,nloc_a,nloc_ae,nloc_oe 

!===============================================================================
contains
!===============================================================================

  subroutine seq_flux_init_mct( cdata, fractions)

    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(seq_cdata), intent(in)  :: cdata
    type(mct_aVect), intent(in)  :: fractions
    !
    ! Local variables
    !
    integer(in)                     :: nloc
    type(mct_gGrid), pointer        :: dom
    type(mct_gsMap), pointer        :: gsMap
    integer                         :: mpicom
    integer                         :: kx,kr   ! fractions indices
    integer                         :: ier
    real(r8), pointer     ::  rmask(:)  ! ocn domain mask
    character(*),parameter :: subName =   '(seq_flux_init_mct) '
    !-----------------------------------------------------------------------

    ! Set cdata pointers
    call seq_cdata_setptrs(cdata, gsMap=gsMap, dom=dom, mpicom=mpicom)
    nloc = mct_gsMap_lsize(gsMap, mpicom)

    ! Input fields atm
    allocate( zbot(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate zbot',ier)
    allocate( ubot(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate ubot',ier)
    allocate( vbot(nloc))
    if(ier/=0) call mct_die(subName,'allocate vbot',ier)
    allocate(thbot(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate thbot',ier)
    allocate(shum(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate shum',ier)
    allocate(dens(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate dens',ier)
    allocate(tbot(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate tbot',ier)
    allocate(ustar(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate ustar',ier)
    allocate(re(nloc), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate re',ier)
    allocate(ssq(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate ssq',ier)
    allocate( uocn(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate uocn',ier)
    allocate( vocn(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate vocn',ier)
    allocate( tocn(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate tocn',ier)

    ! Output fields 
    allocate(sen (nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate sen',ier)
    allocate(lat (nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate lat',ier)
    allocate(evap(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate evap',ier)
    allocate(lwup(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate lwup',ier)
    allocate(taux(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate taux',ier)
    allocate(tauy(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate tauy',ier)
    allocate(tref(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate tref',ier)
    allocate(qref(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate qref',ier)
    allocate(duu10n(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate duu10n',ier)

    ! Grid fields
    allocate( lats(nloc),stat=ier )
    if(ier/=0) call mct_die(subName,'allocate lats',ier)
    allocate( lons(nloc),stat=ier )
    if(ier/=0) call mct_die(subName,'allocate lons',ier)
    allocate( emask(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate emask',ier)
    allocate(mask(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate mask',ier)
    
    ! Get lat, lon, mask, which is time-invariant
    allocate(rmask(nloc),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate rmask',ier)
    call mct_gGrid_exportRAttr(dom, 'lat' , lats , nloc) 
    call mct_gGrid_exportRAttr(dom, 'lon' , lons , nloc) 
    call mct_gGrid_exportRAttr(dom, 'mask', rmask, nloc)
!tcx, want to mask properly, but applying this changes answers to roundoff for some reason
!    kx = mct_aVect_indexRA(fractions,"ofrac")
!    mask = 0
!    where (fractions%rAttr(kx,:) > 0.0_r8) mask = nint(rmask)
    mask = nint(rmask)
    deallocate(rmask)
    emask = mask

    fluxsetting = trim(fluxsetting_atmocn)

  end subroutine seq_flux_init_mct

!===============================================================================

  subroutine seq_flux_initexch_mct(cdata_a,cdata_o)

    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(seq_cdata), intent(in)  :: cdata_a
    type(seq_cdata), intent(in)  :: cdata_o
    !
    ! Local variables
    !

    integer(in) :: kw,ka,ko,iw,ia,io,n
    type(mct_gGrid), pointer        :: dom_a
    type(mct_gGrid), pointer        :: dom_o
    type(mct_gsMap), pointer        :: gsmap_a
    type(mct_gsMap), pointer        :: gsmap_o
    integer(in)                     :: mpicom_a, mpicom_o
    integer(in)                     :: ID_a, ID_o
    character(len=128)              :: strat
    integer                         :: ier
    integer                         :: mytask
    integer(in)                     :: kmsk  ! field indices
    character(len=128)              :: ConfigFileName  ! config file to read
    character(len=128)              :: MapLabel        ! map name
    character(len=128)              :: MapTypeLabel    ! map type
    character(len=256)              :: fileName
    character(len=1)                :: maptype
    character(len=3)                :: Smaptype
    type(mct_aVect)                 :: avdom_oe
    type(mct_list)                  :: sort_keys
    character(*),parameter :: subName =   '(seq_flux_initexch_mct) '
    !-----------------------------------------------------------------------

    !--- Set cdata pointers
    call seq_cdata_setptrs(cdata_a, dom=dom_a, gsmap=gsmap_a, mpicom=mpicom_a, ID=ID_a)
    call seq_cdata_setptrs(cdata_o, dom=dom_o, gsmap=gsmap_o, mpicom=mpicom_o, ID=ID_o)

    call shr_mpi_commrank(mpicom_o,mytask)

    !--- Get mapping file info
    do n = 1,2
       ConfigFileName = "seq_maps.rc"
       if (n == 1) then
          MapLabel = "atm2ocn_fmapname:"
          MapTypeLabel = "atm2ocn_fmaptype:"
       elseif (n == 2) then
          MapLabel = "ocn2atm_fmapname:"
          MapTypeLabel = "ocn2atm_fmaptype:"
       else
          call shr_sys_abort(trim(subname)//' do error1')
       endif

       call I90_allLoadF(ConfigFileName,0,mpicom_o,ier)
       if(ier /= 0) then
          write(logunit,*) trim(subname)," Cant find config file = ",ConfigFileName
          call mct_die("mct_sMapP_initnc","File Not Found")
       endif
       call i90_label(trim(MapLabel),ier)
       if(ier /= 0) then
          write(logunit,*) trim(subname)," Cant find label = ",MapLabel
          call mct_die("mct_sMapP_initnc","Label Not Found")
       endif
       call i90_gtoken(fileName,ier)
       if(ier /= 0) then
          write(logunit,*) trim(subname)," Error reading token = ",fileName
          call mct_die("mct_sMapP_initnc","Error on read")
       endif
       write(logunit,*) trim(subname)," map weight filename = ",trim(fileName)
       call i90_label(trim(MapTypeLabel),ier)
       call i90_gtoken(maptype,ier)

       call i90_Release(ier)

       !--- hardwire decomposition to gsmap_o
       if (n == 1) then
          Smaptype = "src"
          call shr_mct_sMatReaddnc(sMata2o, gsmap_a, gsmap_o, Smaptype, &
             filename=fileName, mytask=mytask, mpicom=mpicom_o)
       elseif (n == 2) then
          Smaptype = "dst"
          call shr_mct_sMatReaddnc(sMato2a, gsmap_o, gsmap_a, Smaptype, &
             filename=fileName, mytask=mytask, mpicom=mpicom_o)
       else
          call shr_sys_abort(trim(subname)//' do error2')
       endif

    enddo

    !--- the two mapping files must have their local indices in identical order
    !--- sort the global indices as a starting point

    call mct_list_init(sort_keys,'grow:gcol')
    call mct_sMat_SortPermute(sMata2o,sort_keys)
    call mct_list_clean(sort_keys)
    call mct_list_init(sort_keys,'gcol:grow')
    call mct_sMat_SortPermute(sMato2a,sort_keys)
    call mct_list_clean(sort_keys)

    !--- now check that they are sorted properly

    nloc_a2o= mct_sMat_lsize(sMata2o)
    nloc_o2a= mct_sMat_lsize(sMato2a)

    if (nloc_a2o /= nloc_o2a) then
       write(logunit,*) trim(subname),' ERROR: sMat sizes',nloc_a2o,nloc_o2a
       call shr_sys_abort(trim(subname)//' ERROR in sMat sizes')
    endif
    ko = mct_sMat_indexIA(sMata2o,'grow')    ! local row (dst) index
    ka = mct_sMat_indexIA(sMato2a,'gcol')    ! local column (src) index
    do n = 1,nloc_a2o
       io = sMata2o%data%iAttr(ko,n)
       ia = sMato2a%data%iAttr(ka,n)
       if (io /= ia) then
          write(logunit,*) trim(subname),' ERROR: sMat indices1 ',io,ia
          call shr_sys_abort(trim(subname)//' ERROR in sMat indices1')
       endif
    enddo
    ko = mct_sMat_indexIA(sMata2o,'gcol')    ! local column (src) index
    ka = mct_sMat_indexIA(sMato2a,'grow')    ! local row (dst) index
    do n = 1,nloc_a2o
       io = sMata2o%data%iAttr(ko,n)
       ia = sMato2a%data%iAttr(ka,n)
       if (io /= ia) then
          write(logunit,*) trim(subname),' ERROR: sMat indices2 ',io,ia
          call shr_sys_abort(trim(subname)//' ERROR in sMat indices2')
       endif
    enddo

    !--- instantiate/create/compute various datatypes

    call mct_sMat_2XgsMap(sMata2o, gsmap_ae, 0, mpicom_a, ID_a)
    call mct_sMat_2YgsMap(sMata2o, gsmap_oe, 0, mpicom_a, ID_a)
    call mct_rearr_init(gsmap_a, gsmap_ae, mpicom_a, Re_a2e)
    call mct_rearr_init(gsmap_ae,gsmap_a,  mpicom_a, Re_e2a)
    call mct_rearr_init(gsmap_o, gsmap_oe, mpicom_a, Re_o2e)
    call mct_rearr_init(gsmap_oe,gsmap_o,  mpicom_a, Re_e2o)
    call mct_sMat_g2lMat(sMata2o,gsmap_ae,'column',mpicom_a)
    call mct_sMat_g2lMat(sMata2o,gsmap_oe,'row',   mpicom_a)
    call mct_sMat_g2lMat(sMato2a,gsmap_ae,'row',   mpicom_a)
    call mct_sMat_g2lMat(sMato2a,gsmap_oe,'column',mpicom_a)

    nloc_a  = mct_gsmap_lsize(gsmap_a,mpicom_a)
    nloc_o  = mct_gsmap_lsize(gsmap_o,mpicom_a)
    nloc_ae = mct_gsmap_lsize(gsmap_ae,mpicom_a)
    nloc_oe = mct_gsmap_lsize(gsmap_oe,mpicom_a)

    call mct_gsmap_clean(gsmap_ae)
    call mct_gsmap_clean(gsmap_oe)

    ! Input fields atm
    allocate( emask(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate emask',ier)
    allocate( zbot(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate zbot',ier)
    allocate( ubot(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate ubot',ier)
    allocate( vbot(nloc_a2o))
    if(ier/=0) call mct_die(subName,'allocate vbot',ier)
    allocate(thbot(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate thbot',ier)
    allocate(shum(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate shum',ier)
    allocate(dens(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate dens',ier)
    allocate(tbot(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate tbot',ier)
    allocate(ustar(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate ustar',ier)
    allocate(re(nloc_a2o), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate re',ier)
    allocate(ssq(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate ssq',ier)
    allocate( uocn(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate uocn',ier)
    allocate( vocn(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate vocn',ier)
    allocate( tocn(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate tocn',ier)

    ! Output fields 
    allocate(sen (nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate sen',ier)
    allocate(lat (nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate lat',ier)
    allocate(evap(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate evap',ier)
    allocate(lwup(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate lwup',ier)
    allocate(taux(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate taux',ier)
    allocate(tauy(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate tauy',ier)
    allocate(tref(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate tref',ier)
    allocate(qref(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate qref',ier)
    allocate(duu10n(nloc_a2o),stat=ier)
    if(ier/=0) call mct_die(subName,'allocate duu10n',ier)

    ! set emask

    call mct_avect_init(avdom_oe,dom_o%data,lsize=nloc_oe)
    call mct_rearr_rearrange(dom_o%data, avdom_oe, Re_o2e, VECTOR=mct_usevector, ALLTOALL=mct_usealltoall)
    ko = mct_sMat_indexIA(sMata2o,'lrow')    ! local dst index
    kmsk = mct_aVect_indexRA(avdom_oe,"mask",dieWith=subName)
    do n = 1,nloc_a2o
       io = sMata2o%data%iAttr(ko,n)
       emask(n) = nint(avdom_oe%rAttr(kmsk,io))
       if (emask(n) == 0) then
          write(logunit,*) trim(subname),' ERROR: weights use masked ocean value'
          call shr_sys_abort(trim(subname)//' ERROR: weights use masked ocean value')
       endif
    enddo

    call mct_aVect_clean(avdom_oe)

    fluxsetting = trim(fluxsetting_exchange)

  end subroutine seq_flux_initexch_mct

!===============================================================================

  subroutine seq_flux_ocnalb_mct( cdata_o, xao_o, fractions_o)

    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(seq_cdata), intent(in)    :: cdata_o
    type(mct_aVect), intent(inout) :: xao_o
    type(mct_aVect), intent(inout) :: fractions_o
    !
    ! Local variables
    !
    type(seq_infodata_type),pointer :: infodata
    type(mct_gGrid), pointer        :: dom_o
    integer(in) :: n,i          ! indices
    real(r8)    :: rlat         ! gridcell latitude in radians
    real(r8)    :: rlon         ! gridcell longitude in radians
    real(r8)    :: cosz         ! Cosine of solar zenith angle
    real(r8)    :: eccen        ! Earth orbit eccentricity
    real(r8)    :: mvelpp       ! Earth orbit
    real(r8)    :: lambm0       ! Earth orbit
    real(r8)    :: obliqr       ! Earth orbit
    real(r8)    :: delta        ! Solar declination angle  in radians
    real(r8)    :: eccf         ! Earth orbit eccentricity factor
    real(r8)    :: calday       ! calendar day including fraction, at 0e
    real(r8)    :: nextsw_cday  ! calendar day of next atm shortwave
    real(r8)    :: anidr        ! albedo: near infrared, direct
    real(r8)    :: avsdr        ! albedo: visible      , direct
    real(r8)    :: anidf        ! albedo: near infrared, diffuse
    real(r8)    :: avsdf        ! albedo: visible      , diffuse
    integer(in) :: ID           ! comm ID
    integer(in) :: ier          ! error code
    logical     :: flux_albav   ! flux avg option
    integer(in) :: kx,kr        ! fractions indices
    integer(in) :: klat,klon,kmsk  ! field indices
    logical     :: update_alb   ! was albedo updated
    logical,save:: first_call = .true. 
    !
    real(R8),parameter :: albdif = 0.06_R8 ! 60 deg reference albedo, diffuse
    real(R8),parameter :: albdir = 0.07_R8 ! 60 deg reference albedo, direct 
    character(*),parameter :: subName =   '(seq_flux_ocnalb_mct) '
    !
    !-----------------------------------------------------------------------

    ! Determine indices

    update_alb = .false.

    call seq_cdata_setptrs(cdata_o, infodata=infodata, ID=ID)
    call seq_infodata_GetData(infodata,flux_albav=flux_albav)

    if (first_call) then
       index_xao_So_anidr  = mct_aVect_indexRA(xao_o,'So_anidr')
       index_xao_So_anidf  = mct_aVect_indexRA(xao_o,'So_anidf')
       index_xao_So_avsdr  = mct_aVect_indexRA(xao_o,'So_avsdr')
       index_xao_So_avsdf  = mct_aVect_indexRA(xao_o,'So_avsdf')

       call seq_cdata_setptrs(cdata_o, dom=dom_o)
       nloc_o  = mct_ggrid_lsize(dom_o)
       klat = mct_gGrid_indexRA(dom_o,"lat" ,dieWith=subName)
       klon = mct_gGrid_indexRA(dom_o,"lon" ,dieWith=subName)
       allocate( lats(nloc_o),stat=ier )
       if(ier/=0) call mct_die(subName,'allocate lats',ier)
       allocate( lons(nloc_o),stat=ier )
       if(ier/=0) call mct_die(subName,'allocate lons',ier)
       do n = 1,nloc_o
          lats(n) = dom_o%data%rAttr(klat,n)
          lons(n) = dom_o%data%rAttr(klon,n)
       enddo
       first_call = .false.
    endif

    if (flux_albav) then
       do n=1,nloc_o   
          anidr = albdir
          avsdr = albdir
          anidf = albdif
          avsdf = albdif

          ! Albedo is now function of latitude (will be new implementation)
          !rlat = const_deg2rad * lats(n)
          !anidr = 0.069_R8 - 0.011_R8 * cos(2._R8 * rlat)
          !avsdr = anidr
          !anidf = anidr
          !avsdf = anidr

          xao_o%rAttr(index_xao_So_avsdr,n) = avsdr
          xao_o%rAttr(index_xao_So_anidr,n) = anidr
          xao_o%rAttr(index_xao_So_avsdf,n) = avsdf
          xao_o%rAttr(index_xao_So_anidf,n) = anidf

       end do
       update_alb = .true.

    else
       ! Solar declination 
       ! Will only do albedo calculation if nextsw_cday is not -1.
       
       call seq_infodata_GetData(infodata,nextsw_cday=nextsw_cday,orb_eccen=eccen, &
          orb_mvelpp=mvelpp, orb_lambm0=lambm0, orb_obliqr=obliqr)
       if (nextsw_cday >= -0.5_r8) then
          calday = nextsw_cday
          call shr_orb_decl(calday, eccen, mvelpp,lambm0, obliqr, delta, eccf)
          ! Compute albedos 
          do n=1,nloc_o
             rlat = const_deg2rad * lats(n)
             rlon = const_deg2rad * lons(n)
             cosz = shr_orb_cosz( calday, rlat, rlon, delta )
             if (cosz  >  0.0_R8) then !--- sun hit --
                anidr = (.026_R8/(cosz**1.7_R8 + 0.065_R8)) +   &
                        (.150_R8*(cosz         - 0.100_R8 ) *   &
                                 (cosz         - 0.500_R8 ) *   &
                                 (cosz         - 1.000_R8 )  )
                avsdr = anidr
                anidf = albdif
                avsdf = albdif
             else !--- dark side of earth ---
                anidr = 1.0_R8
                avsdr = 1.0_R8
                anidf = 1.0_R8
                avsdf = 1.0_R8
             end if

             xao_o%rAttr(index_xao_So_avsdr,n) = avsdr
             xao_o%rAttr(index_xao_So_anidr,n) = anidr
             xao_o%rAttr(index_xao_So_avsdf,n) = avsdf
             xao_o%rAttr(index_xao_So_anidf,n) = anidf

          end do   ! nloc_o
          update_alb = .true.
       endif    ! nextsw_cday
    end if   ! flux_albav

    !--- update current ifrad/ofrad values if albedo was updated

    if (update_alb) then
       kx = mct_aVect_indexRA(fractions_o,"ifrac")
       kr = mct_aVect_indexRA(fractions_o,"ifrad")
       fractions_o%rAttr(kr,:) = fractions_o%rAttr(kx,:)
       kx = mct_aVect_indexRA(fractions_o,"ofrac")
       kr = mct_aVect_indexRA(fractions_o,"ofrad")
       fractions_o%rAttr(kr,:) = fractions_o%rAttr(kx,:)
    endif
       
    end subroutine seq_flux_ocnalb_mct

!===============================================================================

  subroutine seq_flux_atmocnexch_mct( cdata_a, cdata_o, a2x, o2x, xao_a, xao_o, fractions_a, fractions_o)

    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(seq_cdata), intent(in)    :: cdata_a
    type(seq_cdata), intent(in)    :: cdata_o
    type(mct_aVect), intent(in)    :: a2x  
    type(mct_aVect), intent(in)    :: o2x
    type(mct_aVect), intent(inout) :: xao_a
    type(mct_aVect), intent(inout) :: xao_o
    type(mct_aVect), intent(in)    :: fractions_a
    type(mct_aVect), intent(in)    :: fractions_o
    !
    ! Local variables
    !
    type(seq_infodata_type),pointer :: infodata
    type(mct_aVect) :: a2x_e
    type(mct_aVect) :: o2x_e
    type(mct_aVect) :: xaop_ae
    type(mct_aVect) :: xaop_oe
    type(mct_aVect) :: xaop_a
    type(mct_aVect) :: xaop_o
    type(mct_gsmap),pointer :: gsmap_a
    type(mct_gsmap),pointer :: gsmap_o
    type(mct_aVect) :: fractions_oe

    integer(in) :: kw,ka,ko,iw,ia,io,kf
    integer(in) :: n,i          ! indices
    logical     :: dead_comps   ! .true.  => dead components are used
    integer(in) :: index_tref  
    integer(in) :: index_qref  
    integer(in) :: index_duu10n
    integer(in) :: index_ustar 
    integer(in) :: index_ssq   
    integer(in) :: index_re    
    integer(in) :: index_u10   
    integer(in) :: index_taux  
    integer(in) :: index_tauy  
    integer(in) :: index_lat   
    integer(in) :: index_sen   
    integer(in) :: index_evap  
    integer(in) :: index_lwup  
    integer(in) :: index_sumwt
    integer(in) :: atm_nx,atm_ny,ocn_nx,ocn_ny
    real(r8)    :: wt
    character(len=256) :: fldlist  ! subset of xao fields
    !
    character(*),parameter :: subName =   '(seq_flux_atmocnexch_mct) '
    !
    !-----------------------------------------------------------------------

    if (trim(fluxsetting) /= trim(fluxsetting_exchange)) then
       call shr_sys_abort(trim(subname)//' ERROR with init')
    endif

    call seq_cdata_setptrs(cdata_o, infodata=infodata)
    call seq_cdata_setptrs(cdata_a, gsmap=gsmap_a)
    call seq_cdata_setptrs(cdata_o, gsmap=gsmap_o)
    !-----------------------------------------------------

    ! Update ocean surface fluxes 
    ! Must fabricate "reasonable" data (using dead components)

    call seq_infodata_GetData(infodata, dead_comps=dead_comps, &
        atm_nx=atm_nx, atm_ny=atm_ny, &
        ocn_nx=ocn_nx, ocn_ny=ocn_ny)

    if (dead_comps) then
       do n = 1,nloc_a2o
          tocn(n) = 290.0_R8 ! ocn temperature            ~ Kelvin
          uocn(n) =   0.0_R8 ! ocn velocity, zonal        ~ m/s
          vocn(n) =   0.0_R8 ! ocn velocity, meridional   ~ m/s
          zbot(n) =  55.0_R8 ! atm height of bottom layer ~ m
          ubot(n) =   0.0_R8 ! atm velocity, zonal        ~ m/s
          vbot(n) =   2.0_R8 ! atm velocity, meridional   ~ m/s
          thbot(n)= 301.0_R8 ! atm potential temperature  ~ Kelvin
          shum(n) = 1.e-2_R8 ! atm specific humidity      ~ kg/kg
          dens(n) =   1.0_R8 ! atm density                ~ kg/m^3
          tbot(n) = 300.0_R8 ! atm temperature            ~ Kelvin
       enddo
    else        

       !--- instantiate exchange grid aVects
       call mct_AVect_init(a2x_e, a2x, nloc_ae)
       call mct_AVect_zero(a2x_e)
       call mct_AVect_init(o2x_e, o2x, nloc_oe)
       call mct_AVect_zero(o2x_e)

       !--- rearrange a2x and o2x into exchange grid

       call mct_rearr_rearrange(a2x, a2x_e, Re_a2e, VECTOR=mct_usevector, ALLTOALL=mct_usealltoall)
       call mct_rearr_rearrange(o2x, o2x_e, Re_o2e, VECTOR=mct_usevector, ALLTOALL=mct_usealltoall)

       !--- extract fields from a2x and o2x (_e) into local arrays on exchange grid

       ko = mct_sMat_indexIA(sMata2o,'lrow')    ! local row index
       ka = mct_sMat_indexIA(sMata2o,'lcol')    ! local column index

       do n = 1,nloc_a2o
          io = sMata2o%data%iAttr(ko,n)
          ia = sMata2o%data%iAttr(ka,n)
          zbot(n) = a2x_e%rAttr(index_a2x_Sa_z   ,ia)
          ubot(n) = a2x_e%rAttr(index_a2x_Sa_u   ,ia)
          vbot(n) = a2x_e%rAttr(index_a2x_Sa_v   ,ia)
          thbot(n)= a2x_e%rAttr(index_a2x_Sa_ptem,ia)
          shum(n) = a2x_e%rAttr(index_a2x_Sa_shum,ia)
          dens(n) = a2x_e%rAttr(index_a2x_Sa_dens,ia)
          tbot(n) = a2x_e%rAttr(index_a2x_Sa_tbot,ia)
          tocn(n) = o2x_e%rAttr(index_o2x_So_t   ,io)   
          uocn(n) = o2x_e%rAttr(index_o2x_So_u   ,io)
          vocn(n) = o2x_e%rAttr(index_o2x_So_v   ,io)
       enddo
       call mct_aVect_clean(a2x_e)
       call mct_aVect_clean(o2x_e)
    end if

    call shr_flux_atmocn (nloc_a2o , zbot , ubot, vbot, thbot, &
                          shum , dens , tbot, uocn, vocn , &
                          tocn , emask, sen , lat , lwup , &
                          evap , taux , tauy, tref, qref , &
                          duu10n,ustar, re  , ssq , missval = 0.0_r8 )

    !--- create temporary aVects on exchange, atm, or ocn decomp as needed

    fldlist = trim(seq_flds_xao_states)//":"//trim(seq_flds_xao_fluxes)//":sumwt"
    call mct_aVect_init(xaop_ae,rList=trim(fldlist),lsize=nloc_ae)
    call mct_aVect_zero(xaop_ae)
    call mct_aVect_init(xaop_oe,rList=trim(fldlist),lsize=nloc_oe)
    call mct_aVect_zero(xaop_oe)
    call mct_aVect_init(xaop_a, rList=trim(fldlist),lsize=nloc_a)
    call mct_aVect_zero(xaop_a)
    call mct_aVect_init(xaop_o, rList=trim(fldlist),lsize=nloc_o)
    call mct_aVect_zero(xaop_o)

    index_tref   = mct_aVect_indexRA(xaop_ae,"So_tref")
    index_qref   = mct_aVect_indexRA(xaop_ae,"So_qref")
    index_duu10n = mct_aVect_indexRA(xaop_ae,"So_duu10n")
    index_ustar  = mct_aVect_indexRA(xaop_ae,"So_ustar")
    index_ssq    = mct_aVect_indexRA(xaop_ae,"So_ssq")
    index_re     = mct_aVect_indexRA(xaop_ae,"So_re")
    index_u10    = mct_aVect_indexRA(xaop_ae,"So_u10")
    index_taux   = mct_aVect_indexRA(xaop_ae,"Faox_taux")
    index_tauy   = mct_aVect_indexRA(xaop_ae,"Faox_tauy")
    index_lat    = mct_aVect_indexRA(xaop_ae,"Faox_lat")
    index_sen    = mct_aVect_indexRA(xaop_ae,"Faox_sen")
    index_evap   = mct_aVect_indexRA(xaop_ae,"Faox_evap")
    index_lwup   = mct_aVect_indexRA(xaop_ae,"Faox_lwup")
    index_sumwt  = mct_aVect_indexRA(xaop_ae,"sumwt")

    !--- aggregate ocean values locally based on exchange grid decomp

    ko = mct_sMat_indexIA(sMata2o,'lrow')    ! local row index
    ka = mct_sMat_indexIA(sMata2o,'lcol')    ! local column index
    kw = mct_sMat_indexRA(sMata2o,'weight')  ! weight index

    do n = 1,nloc_a2o
       io = sMata2o%data%iAttr(ko,n)
       ia = sMata2o%data%iAttr(ka,n)
       wt = sMata2o%data%rAttr(kw,n)
       xaop_oe%rAttr(index_sen   ,io) = xaop_oe%rAttr(index_sen   ,io) + sen(n) * wt
       xaop_oe%rAttr(index_lat   ,io) = xaop_oe%rAttr(index_lat   ,io) + lat(n) * wt
       xaop_oe%rAttr(index_taux  ,io) = xaop_oe%rAttr(index_taux  ,io) + taux(n)* wt
       xaop_oe%rAttr(index_tauy  ,io) = xaop_oe%rAttr(index_tauy  ,io) + tauy(n)* wt
       xaop_oe%rAttr(index_evap  ,io) = xaop_oe%rAttr(index_evap  ,io) + evap(n)* wt
       xaop_oe%rAttr(index_tref  ,io) = xaop_oe%rAttr(index_tref  ,io) + tref(n)* wt
       xaop_oe%rAttr(index_qref  ,io) = xaop_oe%rAttr(index_qref  ,io) + qref(n)* wt
       xaop_oe%rAttr(index_ustar ,io) = xaop_oe%rAttr(index_ustar ,io) + ustar(n)*wt   ! friction velocity
       xaop_oe%rAttr(index_re    ,io) = xaop_oe%rAttr(index_re    ,io) + re(n)  * wt   ! reynolds number
       xaop_oe%rAttr(index_ssq   ,io) = xaop_oe%rAttr(index_ssq   ,io) + ssq(n) * wt   ! s.hum. saturation at Ts
       xaop_oe%rAttr(index_lwup  ,io) = xaop_oe%rAttr(index_lwup  ,io) + lwup(n)* wt   
       xaop_oe%rAttr(index_duu10n,io) = xaop_oe%rAttr(index_duu10n,io) + duu10n(n)*wt  
       xaop_oe%rAttr(index_u10   ,io) = xaop_oe%rAttr(index_u10   ,io) + sqrt(duu10n(n))*wt
       xaop_oe%rAttr(index_sumwt ,io) = xaop_oe%rAttr(index_sumwt ,io) + wt
    enddo

    !--- aggregate atm values locally based on exchange grid decomp

    ko = mct_sMat_indexIA(sMato2a,'lcol')    ! local column index
    ka = mct_sMat_indexIA(sMato2a,'lrow')    ! local row index
    kw = mct_sMat_indexRA(sMato2a,'weight')  ! weight index
    kf = mct_aVect_indexRA(fractions_o,"ofrac")

    !--- to apply fraction corrections, the indexing must be correct so rearrange
    call mct_avect_init(fractions_oe,fractions_o,lsize=nloc_oe)
    call mct_rearr_rearrange(fractions_o, fractions_oe, Re_o2e, VECTOR=mct_usevector, ALLTOALL=mct_usealltoall)
    do n = 1,nloc_o2a
       io = sMato2a%data%iAttr(ko,n)
       ia = sMato2a%data%iAttr(ka,n)
!tcx   wt = sMato2a%data%rAttr(kw,n)
       wt = sMato2a%data%rAttr(kw,n) * fractions_oe%rAttr(kf,io)
       xaop_ae%rAttr(index_sen   ,ia) = xaop_ae%rAttr(index_sen   ,ia) + sen(n) * wt
       xaop_ae%rAttr(index_lat   ,ia) = xaop_ae%rAttr(index_lat   ,ia) + lat(n) * wt
       xaop_ae%rAttr(index_taux  ,ia) = xaop_ae%rAttr(index_taux  ,ia) + taux(n)* wt
       xaop_ae%rAttr(index_tauy  ,ia) = xaop_ae%rAttr(index_tauy  ,ia) + tauy(n)* wt
       xaop_ae%rAttr(index_evap  ,ia) = xaop_ae%rAttr(index_evap  ,ia) + evap(n)* wt
       xaop_ae%rAttr(index_tref  ,ia) = xaop_ae%rAttr(index_tref  ,ia) + tref(n)* wt
       xaop_ae%rAttr(index_qref  ,ia) = xaop_ae%rAttr(index_qref  ,ia) + qref(n)* wt
       xaop_ae%rAttr(index_ustar ,ia) = xaop_ae%rAttr(index_ustar ,ia) + ustar(n)*wt   ! friction velocity
       xaop_ae%rAttr(index_re    ,ia) = xaop_ae%rAttr(index_re    ,ia) + re(n)  * wt   ! reynolds number
       xaop_ae%rAttr(index_ssq   ,ia) = xaop_ae%rAttr(index_ssq   ,ia) + ssq(n) * wt   ! s.hum. saturation at Ts
       xaop_ae%rAttr(index_lwup  ,ia) = xaop_ae%rAttr(index_lwup  ,ia) + lwup(n)* wt   
       xaop_ae%rAttr(index_duu10n,ia) = xaop_ae%rAttr(index_duu10n,ia) + duu10n(n)*wt  
       xaop_ae%rAttr(index_u10   ,ia) = xaop_ae%rAttr(index_u10   ,ia) + sqrt(duu10n(n))*wt
       xaop_ae%rAttr(index_sumwt ,ia) = xaop_ae%rAttr(index_sumwt ,ia) + wt
    enddo

    call mct_aVect_clean(fractions_oe)

    !--- rearrange and sum from exchange grid to gsmap_a and gsmap_o decomps

    call mct_rearr_rearrange(xaop_ae, xaop_a, Re_e2a, sum=.true., VECTOR=mct_usevector, ALLTOALL=mct_usealltoall)
    call mct_rearr_rearrange(xaop_oe, xaop_o, Re_e2o, sum=.true., VECTOR=mct_usevector, ALLTOALL=mct_usealltoall)

    !--- normalize by sum of wts associated with mapping

    do n = 1,nloc_a
       wt = xaop_a%rAttr(index_sumwt,n)
       if (wt /= 0.0_r8) then 
          wt = 1.0_r8/wt
       else
          wt = 1.0_r8
       endif
       xaop_a%rAttr(:,n) = xaop_a%rAttr(:,n) * wt
    enddo

    do n = 1,nloc_o
       wt = xaop_o%rAttr(index_sumwt,n)
       if (wt /= 0.0_r8) then 
          wt = 1.0_r8/wt
       else
          wt = 1.0_r8
       endif
       xaop_o%rAttr(:,n) = xaop_o%rAttr(:,n) * wt
    enddo

    !--- copy subset of fields to xao_a and xao_o and clean up

    call mct_avect_clean(xaop_ae)
    call mct_avect_clean(xaop_oe)

    call mct_avect_copy(xaop_a, xao_a)
    call mct_avect_copy(xaop_o, xao_o)

    call mct_avect_clean(xaop_a)
    call mct_avect_clean(xaop_o)

  end subroutine seq_flux_atmocnexch_mct

!===============================================================================

  subroutine seq_flux_atmocn_mct( cdata, a2x, o2x, xao)

    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(seq_cdata), intent(in)    :: cdata
    type(mct_aVect), intent(in)    :: a2x  
    type(mct_aVect), intent(in)    :: o2x
    type(mct_aVect), intent(inout) :: xao
    !
    ! Local variables
    !
    type(seq_infodata_type),pointer :: infodata
    integer(in) :: n,i          ! indices
    real(r8)    :: rlat         ! gridcell latitude in radians
    real(r8)    :: rlon         ! gridcell longitude in radians
    real(r8)    :: cosz         ! Cosine of solar zenith angle
    real(r8)    :: eccen        ! Earth orbit eccentricity
    real(r8)    :: mvelpp       ! Earth orbit
    real(r8)    :: lambm0       ! Earth orbit
    real(r8)    :: obliqr       ! Earth orbit
    real(r8)    :: delta        ! Solar declination angle  in radians
    real(r8)    :: eccf         ! Earth orbit eccentricity factor
    real(r8)    :: calday       ! calendar day including fraction, at 0e
    real(r8)    :: nextsw_cday  ! calendar day of next atm shortwave
    real(r8)    :: anidr        ! albedo: near infrared, direct
    real(r8)    :: avsdr        ! albedo: visible      , direct
    real(r8)    :: anidf        ! albedo: near infrared, diffuse
    real(r8)    :: avsdf        ! albedo: visible      , diffuse
    integer(in) :: nloc         ! number of gridcells
    integer(in) :: ID           ! comm ID
    logical     :: flux_albav   ! flux avg option
    logical     :: dead_comps   ! .true.  => dead components are used
    logical     :: first_call = .true.
    !
    real(R8),parameter :: albdif = 0.06_R8 ! 60 deg reference albedo, diffuse
    real(R8),parameter :: albdir = 0.07_R8 ! 60 deg reference albedo, direct 
    character(*),parameter :: subName =   '(seq_flux_atmocn_mct) '
    !
    !-----------------------------------------------------------------------

    if (first_call) then
       index_xao_So_tref   = mct_aVect_indexRA(xao,'So_tref')
       index_xao_So_qref   = mct_aVect_indexRA(xao,'So_qref')
       index_xao_So_ustar  = mct_aVect_indexRA(xao,'So_ustar')  
       index_xao_So_re     = mct_aVect_indexRA(xao,'So_re')  
       index_xao_So_ssq    = mct_aVect_indexRA(xao,'So_ssq')
       index_xao_So_u10    = mct_aVect_indexRA(xao,'So_u10')
       index_xao_So_duu10n = mct_aVect_indexRA(xao,'So_duu10n')
       index_xao_Faox_taux = mct_aVect_indexRA(xao,'Faox_taux')
       index_xao_Faox_tauy = mct_aVect_indexRA(xao,'Faox_tauy')  
       index_xao_Faox_lat  = mct_aVect_indexRA(xao,'Faox_lat')   
       index_xao_Faox_sen  = mct_aVect_indexRA(xao,'Faox_sen')   
       index_xao_Faox_evap = mct_aVect_indexRA(xao,'Faox_evap')   
       index_xao_Faox_lwup = mct_aVect_indexRA(xao,'Faox_lwup')  
       
       index_a2x_Sa_z      = mct_aVect_indexRA(a2x,'Sa_z')
       index_a2x_Sa_u      = mct_aVect_indexRA(a2x,'Sa_u')
       index_a2x_Sa_v      = mct_aVect_indexRA(a2x,'Sa_v')
       index_a2x_Sa_tbot   = mct_aVect_indexRA(a2x,'Sa_tbot')
       index_a2x_Sa_ptem   = mct_aVect_indexRA(a2x,'Sa_ptem')
       index_a2x_Sa_shum   = mct_aVect_indexRA(a2x,'Sa_shum')
       index_a2x_Sa_dens   = mct_aVect_indexRA(a2x,'Sa_dens')
       
       index_o2x_So_t      = mct_aVect_indexRA(o2x,'So_t')
       index_o2x_So_u      = mct_aVect_indexRA(o2x,'So_u')
       index_o2x_So_v      = mct_aVect_indexRA(o2x,'So_v')
       first_call = .false.
    end if
       
    if (trim(fluxsetting) /= trim(fluxsetting_atmocn)) then
       call shr_sys_abort(trim(subname)//' ERROR with init')
    endif

    call seq_cdata_setptrs(cdata, infodata=infodata, ID=ID)
    call seq_infodata_GetData(infodata,flux_albav=flux_albav)

    nloc = mct_aVect_lsize(xao)

    ! Update ocean surface fluxes 
    ! Must fabricate "reasonable" data (using dead components)

    call seq_infodata_GetData(infodata, dead_comps=dead_comps)

    emask = mask
    if (dead_comps) then
       do n = 1,nloc
          mask(n) =   1      ! ocn domain mask            ~ 0 <=> inactive cell
          tocn(n) = 290.0_R8 ! ocn temperature            ~ Kelvin
          uocn(n) =   0.0_R8 ! ocn velocity, zonal        ~ m/s
          vocn(n) =   0.0_R8 ! ocn velocity, meridional   ~ m/s
          zbot(n) =  55.0_R8 ! atm height of bottom layer ~ m
          ubot(n) =   0.0_R8 ! atm velocity, zonal        ~ m/s
          vbot(n) =   2.0_R8 ! atm velocity, meridional   ~ m/s
          thbot(n)= 301.0_R8 ! atm potential temperature  ~ Kelvin
          shum(n) = 1.e-2_R8 ! atm specific humidity      ~ kg/kg
          dens(n) =   1.0_R8 ! atm density                ~ kg/m^3
          tbot(n) = 300.0_R8 ! atm temperature            ~ Kelvin
       enddo
    else	
       do n = 1,nloc
          if (mask(n) /= 0) then	
             zbot(n) = a2x%rAttr(index_a2x_Sa_z   ,n)
             ubot(n) = a2x%rAttr(index_a2x_Sa_u   ,n)
             vbot(n) = a2x%rAttr(index_a2x_Sa_v   ,n)
             thbot(n)= a2x%rAttr(index_a2x_Sa_ptem,n)
             shum(n) = a2x%rAttr(index_a2x_Sa_shum,n)
             dens(n) = a2x%rAttr(index_a2x_Sa_dens,n)
             tbot(n) = a2x%rAttr(index_a2x_Sa_tbot,n)
             tocn(n) = o2x%rAttr(index_o2x_So_t   ,n)   
             uocn(n) = o2x%rAttr(index_o2x_So_u   ,n)
             vocn(n) = o2x%rAttr(index_o2x_So_v   ,n)
             !--- mask missing atm or ocn data
             if (dens(n) < 1.0e-12 .or. tocn(n) < 1.0) then
                emask(n) = 0
                !write(logunit,*) 'aoflux tcx1',n,dens(n),tocn(n)
             endif
          end if
       enddo
    end if

    call shr_flux_atmocn (nloc , zbot , ubot, vbot, thbot, &
                          shum , dens , tbot, uocn, vocn , &
                          tocn , emask, sen , lat , lwup , &
                          evap , taux , tauy, tref, qref , &
! missval should not be needed if flux calc consistent with mrgx2a fraction
!                         duu10n,ustar, re  , ssq, missval = 0.0_r8 )
                          duu10n,ustar, re  , ssq)

    do n = 1,nloc
       if (mask(n) /= 0) then	
          xao%rAttr(index_xao_Faox_sen ,n) = sen(n)
          xao%rAttr(index_xao_Faox_lat ,n) = lat(n)
          xao%rAttr(index_xao_Faox_taux,n) = taux(n)
          xao%rAttr(index_xao_Faox_tauy,n) = tauy(n)
          xao%rAttr(index_xao_Faox_evap,n) = evap(n)
          xao%rAttr(index_xao_So_tref  ,n) = tref(n)
	  xao%rAttr(index_xao_So_qref  ,n) = qref(n)
          xao%rAttr(index_xao_So_ustar ,n) = ustar(n)  ! friction velocity
          xao%rAttr(index_xao_So_re    ,n) = re(n)     ! reynolds number
          xao%rAttr(index_xao_So_ssq   ,n) = ssq(n)    ! s.hum. saturation at Ts
          xao%rAttr(index_xao_Faox_lwup,n) = lwup(n)   
          xao%rAttr(index_xao_So_duu10n,n) = duu10n(n)  
          xao%rAttr(index_xao_So_u10   ,n) = sqrt(duu10n(n))  
       end if
    enddo

  end subroutine seq_flux_atmocn_mct

!===============================================================================

end module seq_flux_mct
