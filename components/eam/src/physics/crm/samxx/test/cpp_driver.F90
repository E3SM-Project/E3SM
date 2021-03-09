
program driver
  use crmdims
  use params, only: crm_iknd, crm_lknd
  use params_kind, only: crm_rknd
  use cpp_interface_mod, only: crm
  use crm_input_module
  use crm_output_module
  use crm_state_module
  use crm_rad_module
  use dmdf
#if HAVE_MPI
  use mpi
#endif
  use iso_c_binding, only: c_bool, c_double
  use gator_mod, only: gator_init, gator_finalize
  implicit none
  integer :: ncrms
  type(crm_input_type)         :: crm_input
  type(crm_output_type)        :: crm_output
  type(crm_state_type)         :: crm_state
  type(crm_rad_type)           :: crm_rad
  integer          , parameter :: plev   = PLEV
  character(len=64), parameter :: fname_in = 'input.nc'
  real(crm_rknd), pointer, contiguous  :: lat0  (:)
  real(crm_rknd), pointer, contiguous  :: long0 (:)
  real(crm_rknd), pointer, contiguous  :: dt_gl (:)
  integer                      :: icrm, ierr, nranks, rank, myTasks_beg, myTasks_end, irank
  logical                      :: masterTask
  real(crm_rknd), allocatable :: read_crm_input_zmid       (:,:)
  real(crm_rknd), allocatable :: read_crm_input_zint       (:,:)
  real(crm_rknd), allocatable :: read_crm_input_tl         (:,:)
  real(crm_rknd), allocatable :: read_crm_input_ql         (:,:)
  real(crm_rknd), allocatable :: read_crm_input_qccl       (:,:)
  real(crm_rknd), allocatable :: read_crm_input_qiil       (:,:)
  real(crm_rknd), allocatable :: read_crm_input_pmid       (:,:)
  real(crm_rknd), allocatable :: read_crm_input_pint       (:,:)
  real(crm_rknd), allocatable :: read_crm_input_pdel       (:,:)
  real(crm_rknd), allocatable :: read_crm_input_ul         (:,:)
  real(crm_rknd), allocatable :: read_crm_input_vl         (:,:)
  real(crm_rknd), allocatable :: read_crm_state_u_wind     (:,:,:,:)
  real(crm_rknd), allocatable :: read_crm_state_v_wind     (:,:,:,:)
  real(crm_rknd), allocatable :: read_crm_state_w_wind     (:,:,:,:)
  real(crm_rknd), allocatable :: read_crm_state_temperature(:,:,:,:)
  real(crm_rknd), allocatable :: read_crm_state_qt         (:,:,:,:)
  real(crm_rknd), allocatable :: read_crm_state_qp         (:,:,:,:)
  real(crm_rknd), allocatable :: read_crm_state_qn         (:,:,:,:)
  real(crm_rknd), allocatable :: read_crm_rad_qrad         (:,:,:,:)
  real(crm_rknd), allocatable :: read_crm_rad_temperature  (:,:,:,:)
  real(crm_rknd), allocatable :: read_crm_rad_qv           (:,:,:,:)
  real(crm_rknd), allocatable :: read_crm_rad_qc           (:,:,:,:)
  real(crm_rknd), allocatable :: read_crm_rad_qi           (:,:,:,:)
  real(crm_rknd), allocatable :: read_crm_rad_cld          (:,:,:,:)
  real(crm_rknd), allocatable :: crm_clear_rh(:,:)
  integer       , allocatable :: gcolp(:)
  character(len=64) :: fprefix = 'cpp_output'
  integer(8) :: t1, t2, tr

#if HAVE_MPI
  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world,nranks,ierr)
  call mpi_comm_rank(mpi_comm_world,rank,ierr)
#else
  nranks = 1
  rank = 0
#endif
  call distr_indices(NCRMS,nranks,rank,myTasks_beg,myTasks_end)
  ncrms = myTasks_end - myTasks_beg + 1
  masterTask = rank == 0

  call gator_init()

  if (masterTask) then
    write(*,*) "File   : ", trim(fname_in)
    write(*,*) "Samples: ", ncrms
    write(*,*) "crm_nx : ", crm_nx
    write(*,*) "crm_ny : ", crm_ny
    write(*,*) "crm_dx : ", crm_dx
    write(*,*) "crm_dt : ", crm_dt
    write(*,*) "plev   : ", plev 
  endif

  ! Allocate model data
  call crm_input%initialize (           ncrms,plev)
  call crm_output_initialize(crm_output,ncrms,plev)
  ! These are normally allocated by pbuf, so we have to do it explicitly
  allocate( crm_state%u_wind     (ncrms,crm_nx,crm_ny,crm_nz) )
  allocate( crm_state%v_wind     (ncrms,crm_nx,crm_ny,crm_nz) )
  allocate( crm_state%w_wind     (ncrms,crm_nx,crm_ny,crm_nz) )
  allocate( crm_state%temperature(ncrms,crm_nx,crm_ny,crm_nz) )
  allocate( crm_state%qt         (ncrms,crm_nx,crm_ny,crm_nz) )
  allocate( crm_state%qp         (ncrms,crm_nx,crm_ny,crm_nz) )
  allocate( crm_state%qn         (ncrms,crm_nx,crm_ny,crm_nz) )
  allocate( crm_rad%qrad         (ncrms,crm_nx_rad,crm_ny_rad,crm_nz) )
  allocate( crm_rad%temperature  (ncrms,crm_nx_rad,crm_ny_rad,crm_nz) )
  allocate( crm_rad%qv           (ncrms,crm_nx_rad,crm_ny_rad,crm_nz) )
  allocate( crm_rad%qc           (ncrms,crm_nx_rad,crm_ny_rad,crm_nz) )
  allocate( crm_rad%qi           (ncrms,crm_nx_rad,crm_ny_rad,crm_nz) )
  allocate( crm_rad%cld          (ncrms,crm_nx_rad,crm_ny_rad,crm_nz) )
  allocate( lat0                 (ncrms) )
  allocate( long0                (ncrms) )
  allocate( dt_gl                (ncrms) )
  allocate( gcolp                (ncrms) )
  allocate( crm_clear_rh         (ncrms,crm_nz) )
  
  ! Allocate transposed arrays because this is the storage format in netcdf
  allocate( read_crm_input_zmid       (plev  ,ncrms))
  allocate( read_crm_input_zint       (plev+1,ncrms))
  allocate( read_crm_input_tl         (plev  ,ncrms))
  allocate( read_crm_input_ql         (plev  ,ncrms))
  allocate( read_crm_input_qccl       (plev  ,ncrms))
  allocate( read_crm_input_qiil       (plev  ,ncrms))
  allocate( read_crm_input_pmid       (plev  ,ncrms))
  allocate( read_crm_input_pint       (plev+1,ncrms))
  allocate( read_crm_input_pdel       (plev  ,ncrms))
  allocate( read_crm_input_ul         (plev  ,ncrms))
  allocate( read_crm_input_vl         (plev  ,ncrms))
  allocate( read_crm_state_u_wind     (crm_nx    ,crm_ny    ,crm_nz,ncrms) )
  allocate( read_crm_state_v_wind     (crm_nx    ,crm_ny    ,crm_nz,ncrms) )
  allocate( read_crm_state_w_wind     (crm_nx    ,crm_ny    ,crm_nz,ncrms) )
  allocate( read_crm_state_temperature(crm_nx    ,crm_ny    ,crm_nz,ncrms) )
  allocate( read_crm_state_qt         (crm_nx    ,crm_ny    ,crm_nz,ncrms) )
  allocate( read_crm_state_qp         (crm_nx    ,crm_ny    ,crm_nz,ncrms) )
  allocate( read_crm_state_qn         (crm_nx    ,crm_ny    ,crm_nz,ncrms) )
  allocate( read_crm_rad_qrad         (crm_nx_rad,crm_ny_rad,crm_nz,ncrms) )
  allocate( read_crm_rad_temperature  (crm_nx_rad,crm_ny_rad,crm_nz,ncrms) )
  allocate( read_crm_rad_qv           (crm_nx_rad,crm_ny_rad,crm_nz,ncrms) )
  allocate( read_crm_rad_qc           (crm_nx_rad,crm_ny_rad,crm_nz,ncrms) )
  allocate( read_crm_rad_qi           (crm_nx_rad,crm_ny_rad,crm_nz,ncrms) )
  allocate( read_crm_rad_cld          (crm_nx_rad,crm_ny_rad,crm_nz,ncrms) )

  ! Read in the samples to drive the code
  if (masterTask) then
    write(*,*) 'Reading the data'
  endif
  call dmdf_read( dt_gl                      , fname_in , trim("dt_gl            ") , myTasks_beg , myTasks_end , .true.  , .false. )
  call dmdf_read( lat0                       , fname_in , trim("latitude0        ") , myTasks_beg , myTasks_end , .false. , .false. )
  call dmdf_read( long0                      , fname_in , trim("longitude0       ") , myTasks_beg , myTasks_end , .false. , .false. )
  call dmdf_read( read_crm_input_zmid        , fname_in , trim("in_zmid          ") , myTasks_beg , myTasks_end , .false. , .false. )
  call dmdf_read( read_crm_input_zint        , fname_in , trim("in_zint          ") , myTasks_beg , myTasks_end , .false. , .false. )
  call dmdf_read( read_crm_input_tl          , fname_in , trim("in_tl            ") , myTasks_beg , myTasks_end , .false. , .false. )
  call dmdf_read( read_crm_input_ql          , fname_in , trim("in_ql            ") , myTasks_beg , myTasks_end , .false. , .false. )
  call dmdf_read( read_crm_input_qccl        , fname_in , trim("in_qccl          ") , myTasks_beg , myTasks_end , .false. , .false. )
  call dmdf_read( read_crm_input_qiil        , fname_in , trim("in_qiil          ") , myTasks_beg , myTasks_end , .false. , .false. )
  call dmdf_read( read_crm_input_pmid        , fname_in , trim("in_pmid          ") , myTasks_beg , myTasks_end , .false. , .false. )
  call dmdf_read( read_crm_input_pint        , fname_in , trim("in_pint          ") , myTasks_beg , myTasks_end , .false. , .false. )
  call dmdf_read( read_crm_input_pdel        , fname_in , trim("in_pdel          ") , myTasks_beg , myTasks_end , .false. , .false. )
  call dmdf_read( read_crm_input_ul          , fname_in , trim("in_ul            ") , myTasks_beg , myTasks_end , .false. , .false. )
  call dmdf_read( read_crm_input_vl          , fname_in , trim("in_vl            ") , myTasks_beg , myTasks_end , .false. , .false. )
  call dmdf_read( read_crm_state_u_wind      , fname_in , trim("state_u_wind     ") , myTasks_beg , myTasks_end , .false. , .false. )
  call dmdf_read( read_crm_state_v_wind      , fname_in , trim("state_v_wind     ") , myTasks_beg , myTasks_end , .false. , .false. )
  call dmdf_read( read_crm_state_w_wind      , fname_in , trim("state_w_wind     ") , myTasks_beg , myTasks_end , .false. , .false. )
  call dmdf_read( read_crm_state_temperature , fname_in , trim("state_temperature") , myTasks_beg , myTasks_end , .false. , .false. )
  call dmdf_read( read_crm_state_qt          , fname_in , trim("state_qt         ") , myTasks_beg , myTasks_end , .false. , .false. )
  call dmdf_read( read_crm_state_qp          , fname_in , trim("state_qp         ") , myTasks_beg , myTasks_end , .false. , .false. )
  call dmdf_read( read_crm_state_qn          , fname_in , trim("state_qn         ") , myTasks_beg , myTasks_end , .false. , .false. )
  call dmdf_read( read_crm_rad_qrad          , fname_in , trim("rad_qrad         ") , myTasks_beg , myTasks_end , .false. , .false. )
  call dmdf_read( read_crm_rad_temperature   , fname_in , trim("rad_temperature  ") , myTasks_beg , myTasks_end , .false. , .false. )
  call dmdf_read( read_crm_rad_qv            , fname_in , trim("rad_qv           ") , myTasks_beg , myTasks_end , .false. , .false. )
  call dmdf_read( read_crm_rad_qc            , fname_in , trim("rad_qc           ") , myTasks_beg , myTasks_end , .false. , .false. )
  call dmdf_read( read_crm_rad_qi            , fname_in , trim("rad_qi           ") , myTasks_beg , myTasks_end , .false. , .false. )
  call dmdf_read( read_crm_rad_cld           , fname_in , trim("rad_cld          ") , myTasks_beg , myTasks_end , .false. , .false. )
  call dmdf_read( crm_input%ps               , fname_in , trim("in_ps            ") , myTasks_beg , myTasks_end , .false. , .false. )
  call dmdf_read( crm_input%phis             , fname_in , trim("in_phis          ") , myTasks_beg , myTasks_end , .false. , .false. )
  call dmdf_read( crm_input%ocnfrac          , fname_in , trim("in_ocnfrac       ") , myTasks_beg , myTasks_end , .false. , .false. )
  call dmdf_read( crm_input%tau00            , fname_in , trim("in_tau00         ") , myTasks_beg , myTasks_end , .false. , .false. )
  call dmdf_read( crm_input%wndls            , fname_in , trim("in_wndls         ") , myTasks_beg , myTasks_end , .false. , .false. )
  call dmdf_read( crm_input%bflxls           , fname_in , trim("in_bflxls        ") , myTasks_beg , myTasks_end , .false. , .false. )
  call dmdf_read( crm_input%fluxu00          , fname_in , trim("in_fluxu00       ") , myTasks_beg , myTasks_end , .false. , .false. )
  call dmdf_read( crm_input%fluxv00          , fname_in , trim("in_fluxv00       ") , myTasks_beg , myTasks_end , .false. , .false. )
  call dmdf_read( crm_input%fluxt00          , fname_in , trim("in_fluxt00       ") , myTasks_beg , myTasks_end , .false. , .false. )
  call dmdf_read( crm_input%fluxq00          , fname_in , trim("in_fluxq00       ") , myTasks_beg , myTasks_end , .false. , .false. )
  call dmdf_read( crm_output%subcycle_factor   , fname_in , trim("out_subcycle_factor") , myTasks_beg , myTasks_end , .false. , .true.  )

  do icrm = 1 , ncrms
    crm_input%zmid       (icrm,:)     = read_crm_input_zmid       (:    ,icrm)                       
    crm_input%zint       (icrm,:)     = read_crm_input_zint       (:    ,icrm)                       
    crm_input%tl         (icrm,:)     = read_crm_input_tl         (:    ,icrm)                       
    crm_input%ql         (icrm,:)     = read_crm_input_ql         (:    ,icrm)                       
    crm_input%qccl       (icrm,:)     = read_crm_input_qccl       (:    ,icrm)                       
    crm_input%qiil       (icrm,:)     = read_crm_input_qiil       (:    ,icrm)                       
    crm_input%pmid       (icrm,:)     = read_crm_input_pmid       (:    ,icrm)                       
    crm_input%pint       (icrm,:)     = read_crm_input_pint       (:    ,icrm)                       
    crm_input%pdel       (icrm,:)     = read_crm_input_pdel       (:    ,icrm)                       
    crm_input%ul         (icrm,:)     = read_crm_input_ul         (:    ,icrm)                       
    crm_input%vl         (icrm,:)     = read_crm_input_vl         (:    ,icrm)                       
    crm_state%u_wind     (icrm,:,:,:) = read_crm_state_u_wind     (:,:,:,icrm) 
    crm_state%v_wind     (icrm,:,:,:) = read_crm_state_v_wind     (:,:,:,icrm) 
    crm_state%w_wind     (icrm,:,:,:) = read_crm_state_w_wind     (:,:,:,icrm) 
    crm_state%temperature(icrm,:,:,:) = read_crm_state_temperature(:,:,:,icrm) 
    crm_state%qt         (icrm,:,:,:) = read_crm_state_qt         (:,:,:,icrm) 
    crm_state%qp         (icrm,:,:,:) = read_crm_state_qp         (:,:,:,icrm) 
    crm_state%qn         (icrm,:,:,:) = read_crm_state_qn         (:,:,:,icrm) 
    crm_rad%qrad         (icrm,:,:,:) = read_crm_rad_qrad         (:,:,:,icrm) 
    crm_rad%temperature  (icrm,:,:,:) = read_crm_rad_temperature  (:,:,:,icrm) 
    crm_rad%qv           (icrm,:,:,:) = read_crm_rad_qv           (:,:,:,icrm) 
    crm_rad%qc           (icrm,:,:,:) = read_crm_rad_qc           (:,:,:,icrm) 
    crm_rad%qi           (icrm,:,:,:) = read_crm_rad_qi           (:,:,:,icrm) 
    crm_rad%cld          (icrm,:,:,:) = read_crm_rad_cld          (:,:,:,icrm) 
  enddo

  if (masterTask) then
    write(*,*) 'Running the CRM'
  endif

#if HAVE_MPI
  call mpi_barrier(mpi_comm_world,ierr)
#endif
  if (masterTask) then
    call system_clock(t1)
  endif

  ! NOTE - the crm_output%tkew variable is a diagnostic quantity that was 
  ! recently added for the 2020 INCITE simulations, so if you get a build error
  ! here you might need to remove this argument

  ! Run the code
  call crm(ncrms, ncrms, dt_gl(1), plev, crm_input%bflxls, crm_input%wndls, crm_input%zmid, crm_input%zint, &
           crm_input%pmid, crm_input%pint, crm_input%pdel, crm_input%ul, crm_input%vl, &
           crm_input%tl, crm_input%qccl, crm_input%qiil, crm_input%ql, crm_input%tau00, &
           crm_state%u_wind, crm_state%v_wind, crm_state%w_wind, crm_state%temperature, &
           crm_state%qt, crm_state%qp, crm_state%qn, crm_rad%qrad, crm_rad%temperature, &
           crm_rad%qv, crm_rad%qc, crm_rad%qi, crm_rad%cld, crm_output%subcycle_factor, &
           crm_output%prectend, crm_output%precstend, crm_output%cld, crm_output%cldtop, &
           crm_output%gicewp, crm_output%gliqwp, crm_output%mctot, crm_output%mcup, crm_output%mcdn, &
           crm_output%mcuup, crm_output%mcudn, crm_output%qc_mean, crm_output%qi_mean, crm_output%qs_mean, &
           crm_output%qg_mean, crm_output%qr_mean, crm_output%mu_crm, crm_output%md_crm, crm_output%eu_crm, &
           crm_output%du_crm, crm_output%ed_crm, crm_output%flux_qt, crm_output%flux_u, crm_output%flux_v, &
           ! crm_output%fluxsgs_qt, crm_output%tkez, crm_output%tkesgsz, crm_output%tkz, crm_output%flux_qp, & 
           crm_output%fluxsgs_qt, crm_output%tkez, crm_output%tkew, crm_output%tkesgsz, crm_output%tkz, crm_output%flux_qp, &
           crm_output%precflux, crm_output%qt_trans, crm_output%qp_trans, crm_output%qp_fall, crm_output%qp_evp, &
           crm_output%qp_src, crm_output%qt_ls, crm_output%t_ls, crm_output%jt_crm, crm_output%mx_crm, crm_output%cltot, &
           crm_output%clhgh, crm_output%clmed, crm_output%cllow, &
           crm_output%sltend, crm_output%qltend, crm_output%qcltend, crm_output%qiltend, &
           crm_output%tk, crm_output%tkh, crm_output%qcl, crm_output%qci, crm_output%qpl, crm_output%qpi, &
           crm_output%z0m, crm_output%taux, crm_output%tauy, crm_output%precc, crm_output%precl, crm_output%precsc, &
           crm_output%precsl, crm_output%prec_crm, crm_clear_rh, &
           lat0, long0, gcolp, 2, logical(.true.,c_bool) , 2._c_double , logical(.true.,c_bool) )


#if HAVE_MPI
  call mpi_barrier(mpi_comm_world,ierr)
#endif
  if (masterTask) then
    call system_clock(t2,tr)
    write(*,*) "Elapsed Time: " , real(t2-t1,8) / real(tr,8)
  endif

  if (masterTask) then
    write(*,*) 'Writing output data'
  endif
  ! dmdf_write(dat,rank,fprefix,vname       ,first,last) !For scalar values
  ! dmdf_write(dat,rank,fprefix,vname,dnames,first,last) !For array values
  do irank = 0 , nranks-1
    if (irank == rank) then
      do icrm = 1 , ncrms
        call dmdf_write( crm_state%u_wind         (icrm,:,:,:) , 1 , fprefix , trim('state_u_wind        ') , (/'crm_nx','crm_ny','crm_nz'/)             , .true.  , .false. )
        call dmdf_write( crm_state%v_wind         (icrm,:,:,:) , 1 , fprefix , trim('state_v_wind        ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
        call dmdf_write( crm_state%w_wind         (icrm,:,:,:) , 1 , fprefix , trim('state_w_wind        ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
        call dmdf_write( crm_state%temperature    (icrm,:,:,:) , 1 , fprefix , trim('state_temperature   ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
        call dmdf_write( crm_state%qt             (icrm,:,:,:) , 1 , fprefix , trim('state_qt            ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
        call dmdf_write( crm_state%qp             (icrm,:,:,:) , 1 , fprefix , trim('state_qp            ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
        call dmdf_write( crm_state%qn             (icrm,:,:,:) , 1 , fprefix , trim('state_qn            ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
        call dmdf_write( crm_output%qcl           (icrm,:,:,:) , 1 , fprefix , trim('output_qcl          ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
        call dmdf_write( crm_output%qci           (icrm,:,:,:) , 1 , fprefix , trim('output_qci          ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
        call dmdf_write( crm_output%qpl           (icrm,:,:,:) , 1 , fprefix , trim('output_qpl          ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
        call dmdf_write( crm_output%qpi           (icrm,:,:,:) , 1 , fprefix , trim('output_qpi          ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
        call dmdf_write( crm_output%tk            (icrm,:,:,:) , 1 , fprefix , trim('output_tk           ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
        call dmdf_write( crm_output%tkh           (icrm,:,:,:) , 1 , fprefix , trim('output_tkh          ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
        call dmdf_write( crm_output%prec_crm      (icrm,:,:)   , 1 , fprefix , trim('output_prec_crm     ') , (/'crm_nx','crm_ny'/)                      , .false. , .false. )
        call dmdf_write( crm_output%wvar          (icrm,:,:,:) , 1 , fprefix , trim('output_wvar         ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
        call dmdf_write( crm_output%aut           (icrm,:,:,:) , 1 , fprefix , trim('output_aut          ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
        call dmdf_write( crm_output%acc           (icrm,:,:,:) , 1 , fprefix , trim('output_acc          ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
        call dmdf_write( crm_output%evpc          (icrm,:,:,:) , 1 , fprefix , trim('output_evpc         ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
        call dmdf_write( crm_output%evpr          (icrm,:,:,:) , 1 , fprefix , trim('output_evpr         ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
        call dmdf_write( crm_output%mlt           (icrm,:,:,:) , 1 , fprefix , trim('output_mlt          ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
        call dmdf_write( crm_output%sub           (icrm,:,:,:) , 1 , fprefix , trim('output_sub          ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
        call dmdf_write( crm_output%dep           (icrm,:,:,:) , 1 , fprefix , trim('output_dep          ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
        call dmdf_write( crm_output%con           (icrm,:,:,:) , 1 , fprefix , trim('output_con          ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
        call dmdf_write( crm_output%cltot         (icrm)       , 1 , fprefix , trim('output_cltot        ')                                              , .false. , .false. )
        call dmdf_write( crm_output%cllow         (icrm)       , 1 , fprefix , trim('output_cllow        ')                                              , .false. , .false. )
        call dmdf_write( crm_output%clmed         (icrm)       , 1 , fprefix , trim('output_clmed        ')                                              , .false. , .false. )
        call dmdf_write( crm_output%clhgh         (icrm)       , 1 , fprefix , trim('output_clhgh        ')                                              , .false. , .false. )
        call dmdf_write( crm_output%precc         (icrm)       , 1 , fprefix , trim('output_precc        ')                                              , .false. , .false. )
        call dmdf_write( crm_output%precl         (icrm)       , 1 , fprefix , trim('output_precl        ')                                              , .false. , .false. )
        call dmdf_write( crm_output%precsc        (icrm)       , 1 , fprefix , trim('output_precsc       ')                                              , .false. , .false. )
        call dmdf_write( crm_output%precsl        (icrm)       , 1 , fprefix , trim('output_precsl       ')                                              , .false. , .false. )
        call dmdf_write( crm_output%cldtop        (icrm,:)     , 1 , fprefix , trim('output_cldtop       ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%qc_mean       (icrm,:)     , 1 , fprefix , trim('output_qc_mean      ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%qi_mean       (icrm,:)     , 1 , fprefix , trim('output_qi_mean      ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%qs_mean       (icrm,:)     , 1 , fprefix , trim('output_qs_mean      ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%qg_mean       (icrm,:)     , 1 , fprefix , trim('output_qg_mean      ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%qr_mean       (icrm,:)     , 1 , fprefix , trim('output_qr_mean      ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%sltend        (icrm,:)     , 1 , fprefix , trim('output_sltend       ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%qltend        (icrm,:)     , 1 , fprefix , trim('output_qltend       ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%qcltend       (icrm,:)     , 1 , fprefix , trim('output_qcltend      ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%qiltend       (icrm,:)     , 1 , fprefix , trim('output_qiltend      ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%cld           (icrm,:)     , 1 , fprefix , trim('output_cld          ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%gicewp        (icrm,:)     , 1 , fprefix , trim('output_gicewp       ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%gliqwp        (icrm,:)     , 1 , fprefix , trim('output_gliqwp       ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%mctot         (icrm,:)     , 1 , fprefix , trim('output_mctot        ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%mcup          (icrm,:)     , 1 , fprefix , trim('output_mcup         ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%mcdn          (icrm,:)     , 1 , fprefix , trim('output_mcdn         ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%mcuup         (icrm,:)     , 1 , fprefix , trim('output_mcuup        ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%mcudn         (icrm,:)     , 1 , fprefix , trim('output_mcudn        ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%mu_crm        (icrm,:)     , 1 , fprefix , trim('output_mu_crm       ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%md_crm        (icrm,:)     , 1 , fprefix , trim('output_md_crm       ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%du_crm        (icrm,:)     , 1 , fprefix , trim('output_du_crm       ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%eu_crm        (icrm,:)     , 1 , fprefix , trim('output_eu_crm       ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%ed_crm        (icrm,:)     , 1 , fprefix , trim('output_ed_crm       ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%jt_crm        (icrm)       , 1 , fprefix , trim('output_jt_crm       ')                                              , .false. , .false. )
        call dmdf_write( crm_output%mx_crm        (icrm)       , 1 , fprefix , trim('output_mx_crm       ')                                              , .false. , .false. )
        call dmdf_write( crm_output%flux_qt       (icrm,:)     , 1 , fprefix , trim('output_flux_qt      ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%fluxsgs_qt    (icrm,:)     , 1 , fprefix , trim('output_fluxsgs_qt   ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%tkez          (icrm,:)     , 1 , fprefix , trim('output_tkez         ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%tkesgsz       (icrm,:)     , 1 , fprefix , trim('output_tkesgsz      ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%tkz           (icrm,:)     , 1 , fprefix , trim('output_tkz          ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%flux_u        (icrm,:)     , 1 , fprefix , trim('output_flux_u       ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%flux_v        (icrm,:)     , 1 , fprefix , trim('output_flux_v       ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%flux_qp       (icrm,:)     , 1 , fprefix , trim('output_flux_qp      ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%precflux      (icrm,:)     , 1 , fprefix , trim('output_precflux     ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%qt_ls         (icrm,:)     , 1 , fprefix , trim('output_qt_ls        ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%qt_trans      (icrm,:)     , 1 , fprefix , trim('output_qt_trans     ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%qp_trans      (icrm,:)     , 1 , fprefix , trim('output_qp_trans     ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%qp_fall       (icrm,:)     , 1 , fprefix , trim('output_qp_fall      ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%qp_src        (icrm,:)     , 1 , fprefix , trim('output_qp_src       ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%qp_evp        (icrm,:)     , 1 , fprefix , trim('output_qp_evp       ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%t_ls          (icrm,:)     , 1 , fprefix , trim('output_t_ls         ') , (/'nlev'/)                                 , .false. , .false. )
        call dmdf_write( crm_output%prectend      (icrm)       , 1 , fprefix , trim('output_prectend     ')                                              , .false. , .false. )
        call dmdf_write( crm_output%precstend     (icrm)       , 1 , fprefix , trim('output_precstend    ')                                              , .false. , .false. )
        call dmdf_write( crm_output%taux          (icrm)       , 1 , fprefix , trim('output_taux         ')                                              , .false. , .false. )
        call dmdf_write( crm_output%tauy          (icrm)       , 1 , fprefix , trim('output_tauy         ')                                              , .false. , .false. )
        call dmdf_write( crm_output%z0m           (icrm)       , 1 , fprefix , trim('output_z0m          ')                                              , .false. , .false. )
        call dmdf_write( crm_output%subcycle_factor (icrm)       , 1 , fprefix , trim('output_subcycle_factor')                                              , .false. , .false. )
        call dmdf_write( crm_rad%qrad             (icrm,:,:,:) , 1 , fprefix , trim('rad_qrad            ') , (/'crm_nx_rad','crm_ny_rad','crm_nz    '/) , .false. , .false. )
        call dmdf_write( crm_rad%temperature      (icrm,:,:,:) , 1 , fprefix , trim('rad_temperature     ') , (/'crm_nx_rad','crm_ny_rad','crm_nz    '/) , .false. , .false. )
        call dmdf_write( crm_rad%qv               (icrm,:,:,:) , 1 , fprefix , trim('rad_qv              ') , (/'crm_nx_rad','crm_ny_rad','crm_nz    '/) , .false. , .false. )
        call dmdf_write( crm_rad%qc               (icrm,:,:,:) , 1 , fprefix , trim('rad_qc              ') , (/'crm_nx_rad','crm_ny_rad','crm_nz    '/) , .false. , .false. )
        call dmdf_write( crm_rad%qi               (icrm,:,:,:) , 1 , fprefix , trim('rad_qi              ') , (/'crm_nx_rad','crm_ny_rad','crm_nz    '/) , .false. , .false. )
        call dmdf_write( crm_rad%cld              (icrm,:,:,:) , 1 , fprefix , trim('rad_cld             ') , (/'crm_nx_rad','crm_ny_rad','crm_nz    '/) , .false. , .false. )
        call dmdf_write( crm_clear_rh             (icrm,:)     , 1 , fprefix , trim('crm_clear_rh        ') , (/'crm_nz'/)                               , .false. , .true.  )
      enddo
    endif
#if HAVE_MPI
    call mpi_barrier(mpi_comm_world,ierr)
#endif
  enddo

  call gator_finalize()
#if HAVE_MPI
  call mpi_finalize(ierr)
#endif


contains


  subroutine distr_indices(nTasks,nThreads,myThreadID,myTasks_beg,myTasks_end)
    implicit none
    integer, intent(in   ) :: nTasks
    integer, intent(in   ) :: nThreads
    integer, intent(in   ) :: myThreadID
    integer, intent(  out) :: myTasks_beg
    integer, intent(  out) :: myTasks_end
    real :: nper
    nper = real(nTasks)/nThreads
    myTasks_beg = nint( nper* myThreadID    )+1
    myTasks_end = nint( nper*(myThreadID+1) )
  end subroutine distr_indices


end program driver


