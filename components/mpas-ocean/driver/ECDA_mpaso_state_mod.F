!BOP
!   !MODULE: ecda_mpaso_state_mod
!   !INTERFACE:
module ECDA_mpaso_state_mod
!   !DESCRIPTION:
!   ! Define the state variable type and initialize the state variabel for ODA-EnKF
!   !
!   !REVISION HISTORY:
!    May  2022 Y. Liu <liu6@tamu.edu> initial release
!   !USES:
  use shr_kind_mod, only : r8=> shr_kind_r8
  use ecda_shr_type_mod
  use seq_comm_mct,  only : logunit


  implicit none
  public :: ocn_DAstate_init, ocn_DAstate_register,set_state_type_value, &
            average_state_data_calculate,average_state_data_initial, &
            IAU_state_data_initial,relax_inflation, &
            ocn_DAstate_allocate,ocn_DAstate_deallocate,&
            mpaso_DAstate_allocate,mpaso_DAstate_register,mpas_relax_inflation
 
  integer, parameter, public :: n_OCN_state =  4    ! Currently using T,S
  integer, parameter, public :: OCN_TEMP_index = 1     ! index of ocean temeprature for CELL
  integer, parameter, public :: OCN_SALT_index = 2     ! index of ocean salinity for ncell
  integer, parameter, public :: OCN_VEL_index = 3      ! index of ocean velocity
  integer, parameter, public :: OCN_SSH_index = 4      ! index of ocean SSH


  real(r8),public            ::   Mx_da_increm(4) = (/5.0,5.0,2.0,1.0/)              ! set the maxium analysis increment 
  integer, parameter, public :: OCN_state_dim(n_OCN_state) = (/2,2,2,1/)   ! dimesion of ocean temeprature for CELL
  
  integer,public             ::  steps_IAU, IAU_update

  integer,public             ::  my_inst, ensemble_size
  integer, public            ::  allCellsSolve, DACellsSolve,DACellsSolve_halo
  integer, public            ::  allEdgesSolve, DAEdgesSolve,DAEdgesSolve_halo
  integer, public            ::  DA_nvert                  ! total vertical layers

  !    variable for Data assimilation
  integer,public :: mpicom_ocn_local           ! mpaso model communicator
  integer,public :: mpicom_da_all,mpicom_da  !comunicator
  integer,public :: iam_da_all, iam_da       ! comunicater ID
  logical,public :: iamin_da_all,iamin_da
  logical,public :: iamroot_da_all,iamroot_da ! root


  integer,public             :: da_its,da_ite,da_jts,da_jte            !DA subdomain grid range without halo
  integer,public             :: da_ibs,da_ibe,da_jbs,da_jbe            !DA subdomain grid range with halo
  integer,public             :: da_Cell_nts,da_Cell_nte, da_Cell_nbs,da_Cell_nbe           ! for mpaso, [da_Cell_nbs,da_Cell_nbe]=[da_Cell_nts da_Cell_nte]+halo
  integer,public             :: da_Edge_nts,da_Edge_nte, da_Edge_nbs,da_Edge_nbe           ! 

    integer,public             :: ibs, ibe, jbs, jbe, &
                                ids, ide, jds, jde, &
                                its, ite, jts, jte, &
                                iue, iuw, jvn, jvs

  
  real(r8),public            ::  RR=0.0                                   ! the default maxium localization cutoff 2*RR=0 
  logical,public             ::  average_obs = .FALSE.
  character(len=10),public   ::  YMDH ='2010010100'
  
  integer :: ocnDALogUnit     ! unit number for ocn DA log
!  
  real(r8),allocatable,public          :: h_da(:,:), lon_da(:,:), lat_da(:,:),z_r_da(:,:,:)
  real(r8),allocatable,public    :: DA_lonCell(:), DA_latCell(:),DA_lonEdge(:), DA_latEdge(:), DA_depth(:)
  integer, allocatable,public    :: DA_maxLevelCell(:), DA_indexToCellID(:), DA_indexToEdgeID(:) 
  type(ECDA_state_type),public         ::  OCN_DAstate_data(n_OCN_state) 
contains

!==============================================================================
!BOP
!   !ROUTINE: ocn_DAstate_init
!   pre-initialize the OCN state for DA
!   !INTERFACE:

    subroutine ocn_DAstate_init()
!EOP

!BOC
        implicit none
        ! --- Local variables ---
        character(*), parameter :: subName =   '(ocn_DAstate_init) '
        integer :: ivar
        
        OCN_DAstate_data(OCN_TEMP_index)%name = 'OCN_temp'
        OCN_DAstate_data(OCN_SALT_index)%name = 'OCN_salt'
        OCN_DAstate_data(OCN_VEL_index)%name = 'OCN_vel'
        OCN_DAstate_data(OCN_SSH_index)%name = 'OCN_zeta'
        do ivar =1, n_OCN_state
           OCN_DAstate_data(ivar)%state_dim = OCN_state_dim(ivar)
           OCN_DAstate_data(ivar)%required = .false.
           OCN_DAstate_data(ivar)%update = .false.
           OCN_DAstate_data(ivar)%averaged = .false.
           OCN_DAstate_data(ivar)%r2p_fac = 0.0         ! no R2P by default, @YL, not useful 
        end do       

        !-- set the range of states
        ! temperature 
        OCN_DAstate_data(OCN_TEMP_index)%mn_value = -1.97  
        OCN_DAstate_data(OCN_TEMP_index)%mx_value = 40.0 
        ! salinity
        OCN_DAstate_data(OCN_SALT_index)%mn_value = 0.001
        OCN_DAstate_data(OCN_SALT_index)%mx_value = 60.0
        ! ssh
        OCN_DAstate_data(OCN_SALT_index)%mn_value = -2.0
        OCN_DAstate_data(OCN_SALT_index)%mx_value = 2.0
        ! vel
        OCN_DAstate_data(OCN_VEl_index)%mn_value = -2.0
        OCN_DAstate_data(OCN_VEL_index)%mx_value = 2.0

   end subroutine ocn_DAstate_init


!==============================================================================
!BOP
!   !ROUTINE: mpaso_DAstate_register
!   register the mpas_o state for DA when it is required
!   active IAU_increm wen the state is updated for model simulation
!   
!   !INTERFACE:

  subroutine mpaso_DAstate_register()
!EOP

!BOC
        implicit none
        ! --- Local variables ---
        character(*), parameter :: subName =   '(ocn_DAstate_register) '
        integer :: ivar,ierr

        do ivar = 1, n_OCN_state
           if (.not.  OCN_DAstate_data(ivar)%required) CYCLE
           select case (OCN_DAstate_data(ivar)%state_dim)
              case(1)
                     if(OCN_DAstate_data(ivar)%averaged) allocate(OCN_DAstate_data(ivar)%ave_value1d(allCellsSolve),stat=ierr)
                     if(OCN_DAstate_data(ivar)%update) allocate(OCN_DAstate_data(ivar)%IAU_increm1D(allCellsSolve),stat=ierr)
              case(2) 
                     if(OCN_DAstate_data(ivar)%averaged) allocate(OCN_DAstate_data(ivar)%ave_value2d(allCellsSolve,DA_nvert),stat=ierr)
                     if(OCN_DAstate_data(ivar)%update) allocate(OCN_DAstate_data(ivar)%IAU_increm2D(allCellsSolve,DA_nvert),stat=ierr)
              CASE DEFAULT
                  write(logunit,*) trim(subname),'state dimension:',OCN_DAstate_data(ivar)%state_dim,'to be develop'
            end select

        end do

   end subroutine mpaso_DAstate_register

!==============================================================================
!BOP
!   !ROUTINE: ocn_DAstate_register
!   register the OCN state for DA when it is required
!   active IAU_increm wen the state is updated for model simulation
!   
!   !INTERFACE:

  subroutine ocn_DAstate_register()
!EOP

!BOC
        implicit none
        ! --- Local variables ---
        character(*), parameter :: subName =   '(ocn_DAstate_register) '
        integer :: ivar,ierr

        do ivar = 1, n_OCN_state
           if (.not.  OCN_DAstate_data(ivar)%required) CYCLE

           ! allocate 2D state variable ensemble 
           if (OCN_DAstate_data(ivar)%state_dim == 3) then
              if(OCN_DAstate_data(ivar)%averaged) allocate(OCN_DAstate_data(ivar)%ave_value(its:ite,jts:jte,DA_nvert),stat=ierr)
              if(OCN_DAstate_data(ivar)%update)allocate (OCN_DAstate_data(ivar)%IAU_increm(ibs:ibe,jbs:jbe,DA_nvert), stat=ierr)
           endif

           ! allocate 2D state variable ensemble 
           if (OCN_DAstate_data(ivar)%state_dim == 2) then
              if(OCN_DAstate_data(ivar)%averaged) allocate(OCN_DAstate_data(ivar)%ave_value2d(its:ite,jts:jte),stat=ierr)
              if(OCN_DAstate_data(ivar)%update)allocate (OCN_DAstate_data(ivar)%IAU_increm2d(ibs:ibe,jbs:jbe), stat=ierr)
           endif

        end do

 end subroutine ocn_DAstate_register

!==============================================================================
!BOP
!   !ROUTINE: ocn_DAstate_allocation
!    allocation the OCN prior,posterior for four DA
!   active IAU_increm when the state is updated for model simulation
!   
!   !INTERFACE:

  subroutine mpaso_DAstate_allocate()
!EOP

!BOC
        implicit none
        ! --- Local variables ---
        character(*), parameter :: subName =   '(mpaso_DAstate_deallocate) '
        integer :: da_nbs,da_nbe
        integer :: ivar,ierr


        do ivar = 1, n_OCN_state
           if(.not. OCN_DAstate_data(ivar)%required) CYCLE
           !@YL, need add the edge/cell settings for state type, now only for cell
           da_nbs=da_Cell_nbs
           da_nbe=da_Cell_nbe
           select case (OCN_DAstate_data(ivar)%state_dim)
                 case(1)
                      allocate(OCN_DAstate_data(ivar)%prior1d(da_nbs:da_nbe,ensemble_size),&
                          OCN_DAstate_data(ivar)%posterior1d(da_nbs:da_nbe,ensemble_size),stat=ierr)
                      if(OCN_DAstate_data(ivar)%averaged) &
                          allocate(OCN_DAstate_data(ivar)%prior1d_ave(da_nbs:da_nbe,ensemble_size),stat=ierr)

                 case(2)
                      allocate(OCN_DAstate_data(ivar)%prior2d(da_nbs:da_nbe,da_nvert,ensemble_size),&
                          OCN_DAstate_data(ivar)%posterior2d(da_nbs:da_nbe,da_nvert,ensemble_size),stat=ierr)
                      if(OCN_DAstate_data(ivar)%averaged) &
                          allocate(OCN_DAstate_data(ivar)%prior2d_ave(da_nbs:da_nbe,da_nvert,ensemble_size),stat=ierr)

                 CASE DEFAULT
                      write(logunit,*)trim(subname),'More dimension to develop'
                      return
           end select     
        enddo
        !write(logunit,*) trim(subname),' Finished'
     end subroutine mpaso_DAstate_allocate

!==============================================================================
!==============================================================================
!BOP
!   !ROUTINE: ocn_DAstate_allocation
!    allocation the OCN prior,posterior for four DA
!   active IAU_increm when the state is updated for model simulation
!   
!   !INTERFACE:

  subroutine ocn_DAstate_allocate()
!EOP

!BOC
        implicit none
        ! --- Local variables ---
        character(*), parameter :: subName =   '(ocn_DAstate_deallocate) '
        integer :: ivar,ierr

        do ivar = 1, n_OCN_state
           if(.not. OCN_DAstate_data(ivar)%required) CYCLE
           if(OCN_DAstate_data(ivar)%state_dim == 2) then
              allocate (OCN_DAstate_data(ivar)%prior2d(da_ibs:da_ibe,da_jbs:da_jbe,ensemble_size), &
                    OCN_DAstate_data(ivar)%posterior2d(da_ibs:da_ibe,da_jbs:da_jbe,ensemble_size),&
                    OCN_DAstate_data(ivar)%KKlev(1), stat=ierr)
              if(average_obs) &
               allocate(OCN_DAstate_data(ivar)%prior2d_ave(da_ibs:da_ibe,da_jbs:da_jbe,ensemble_size),stat=ierr)                  
           endif
           
           if(OCN_DAstate_data(ivar)%state_dim == 3) then
             allocate (OCN_DAstate_data(ivar)%prior(da_ibs:da_ibe,da_jbs:da_jbe,DA_nvert,ensemble_size), &
                        OCN_DAstate_data(ivar)%posterior(da_ibs:da_ibe,da_jbs:da_jbe,DA_nvert,ensemble_size),&
                        OCN_DAstate_data(ivar)%KKlev(DA_nvert), stat=ierr)
              if(average_obs) &
               allocate(OCN_DAstate_data(ivar)%prior_ave(da_ibs:da_ibe,da_jbs:da_jbe,DA_nvert,ensemble_size),stat=ierr)                  
           endif

        enddo
        !write(logunit,*) trim(subname),' Finished'
     end subroutine ocn_DAstate_allocate

!==============================================================================
!BOP
!   !ROUTINE: ocn_DAstate_deallocation
!    deallocation the OCN prior,posterior for after DA
!   active IAU_increm wen the state is updated for model simulation
!   
!   !INTERFACE:

  subroutine ocn_DAstate_deallocate()
!EOP

!BOC
        implicit none
        ! --- Local variables ---
        character(*), parameter :: subName =   '(ocn_DAstate_deallocate) '
        integer :: ivar,ierr

        do ivar = 1, n_OCN_state
          select case (OCN_DAstate_data(ivar)%state_dim)
                 case(1)
                         deallocate(OCN_DAstate_data(ivar)%prior1d, OCN_DAstate_data(ivar)%posterior1d, &
                                 OCN_DAstate_data(ivar)%prior1d_ave,stat=ierr)
                 case(2)
                         deallocate(OCN_DAstate_data(ivar)%prior2d,OCN_DAstate_data(ivar)%posterior2d,&
                          OCN_DAstate_data(ivar)%prior2d_ave,OCN_DAstate_data(ivar)%KKlev, stat=ierr)
                 case(3)
                         deallocate(OCN_DAstate_data(ivar)%prior,OCN_DAstate_data(ivar)%posterior,&
                          OCN_DAstate_data(ivar)%prior_ave,OCN_DAstate_data(ivar)%KKlev, stat=ierr)
          CASE DEFAULT
                  return
          end select     
        enddo
        !write(logunit,*) trim(subname),' Finished'
 end subroutine ocn_DAstate_deallocate


!==============================================================================

!BOP
! ! Rountine set_state_type_value(obs_data,...)
! ! set the obs configures

    subroutine set_state_type_value(ID,type_name, state_dim,required, update,averaged,r2p_fac)

! EOP 
      integer, intent(in)                                    ::      ID
      character(*),optional,intent(in)                       ::      type_name
      integer, optional,intent(in)                           ::      state_dim
      logical, optional,intent(in)                           ::      required, &
                                                                     update, &
                                                                     averaged
      real(r8),optional,intent(in)                           ::      r2p_fac

      if(present(type_name)) OCN_DAstate_data(ID)%name    = type_name
      if(present(state_dim)) OCN_DAstate_data(ID)%state_dim = state_dim
      if(present(required)) OCN_DAstate_data(ID)%required = required
      if(present(update)) OCN_DAstate_data(ID)%update = update
      if(present(averaged)) OCN_DAstate_data(ID)%averaged = averaged
      if(present(r2p_fac)) OCN_DAstate_data(ID)%r2p_fac = r2p_fac

    end subroutine set_state_type_value



!==============================================================================
!BOP
! ! Rountine average_state_data_initial()
! ! set the obs configures

    subroutine average_state_data_initial()

! EOP 
      integer :: ivar
      do ivar=1,n_OCN_state
         if( OCN_DAstate_data(ivar)%required .and. OCN_DAstate_data(ivar)%averaged) then
            if(OCN_DAstate_data(ivar)%state_dim == 3) OCN_DAstate_data(ivar)%ave_value =0.0
            if(OCN_DAstate_data(ivar)%state_dim == 2) OCN_DAstate_data(ivar)%ave_value2d =0.0
            if(OCN_DAstate_data(ivar)%state_dim == 1) OCN_DAstate_data(ivar)%ave_value1d =0.0
            OCN_DAstate_data(ivar)%ave_num =  0
         endif
      enddo

    end subroutine average_state_data_initial

!==============================================================================
!BOP
! ! Rountine average_state_data_calculate()
! ! set the obs configures

    subroutine average_state_data_calculate()

! EOP 
      integer :: ivar
      do ivar=1,n_OCN_state
         if( OCN_DAstate_data(ivar)%required .and. OCN_DAstate_data(ivar)%averaged) then
            if(OCN_DAstate_data(ivar)%state_dim == 3)&
               OCN_DAstate_data(ivar)%ave_value = OCN_DAstate_data(ivar)%ave_value/ OCN_DAstate_data(ivar)%ave_num
            if(OCN_DAstate_data(ivar)%state_dim == 2) &
               OCN_DAstate_data(ivar)%ave_value2d =OCN_DAstate_data(ivar)%ave_value2d/ OCN_DAstate_data(ivar)%ave_num
            if(OCN_DAstate_data(ivar)%state_dim == 1) &
               OCN_DAstate_data(ivar)%ave_value1d =OCN_DAstate_data(ivar)%ave_value1d/ OCN_DAstate_data(ivar)%ave_num
            OCN_DAstate_data(ivar)%ave_num =  0
         endif
      enddo

    end subroutine average_state_data_calculate
 
!==============================================================================
!BOP
! ! Rountine IAU_state_data_initial()
! ! set the obs configures

    subroutine IAU_state_data_initial()

! EOP 
      integer :: ivar
      do ivar=1,n_OCN_state
         if(.not. OCN_DAstate_data(ivar)%update) return
         if(OCN_DAstate_data(ivar)%state_dim == 3) OCN_DAstate_data(ivar)%IAU_increm =0.0
         if(OCN_DAstate_data(ivar)%state_dim == 2) OCN_DAstate_data(ivar)%IAU_increm2D =0.0
         if(OCN_DAstate_data(ivar)%state_dim == 1) OCN_DAstate_data(ivar)%IAU_increm1D =0.0
      enddo

    end subroutine IAU_state_data_initial

 
!==============================================================================
!BOP
!  ! routine relax2prior_inflation()
!   relaxiation inflation to prior spread
    subroutine mpas_relax_inflation()
! EOP
    use sEnKF_mod,     only :inflation_RTPS
    
!BOC
        implicit none
        ! --- Local variables ---
        character(*), parameter :: subName =   '(mpas_relax_inflation) '
        integer                 :: i,j,k,ivar
        integer                 :: da_nbs,da_nbe
        real(r8)                :: RTP_rate
        real(r8)                :: prior(ensemble_size),analysis(ensemble_size)

        !write(logunit,*) trim(subname),'inflation of RTPS'
        flush(logunit)
        do ivar  = 1, n_OCN_state
           if(.not. OCN_DAstate_data(ivar)%update) cycle
           RTP_rate=OCN_DAstate_data(ivar)%r2p_fac
           !write(logunit,*) trim(subname),trim(OCN_DAstate_data(ivar)%name),'RTPS rate',RTP_rate
           if(RTP_rate <= tiny(RTP_rate)) cycle 
           da_nbs=da_Cell_nbs
           da_nbe=da_cell_nbe
           do i=da_nbs,da_nbe
               if (OCN_DAstate_data(ivar)%state_dim==1) then
                     prior=OCN_DAstate_data(ivar)%prior1d(i,:)
                     analysis=OCN_DAstate_data(ivar)%posterior1d(i,:)
                     call inflation_RTPS(analysis,prior,RTP_rate,ensemble_size)
                     OCN_DAstate_data(ivar)%posterior1d(i,:)=analysis
                 endif
                 if (OCN_DAstate_data(ivar)%state_dim==2) then
                    do k=1,DA_nvert
                       prior=OCN_DAstate_data(ivar)%prior2d(i,k,:)
                       analysis=OCN_DAstate_data(ivar)%posterior2d(i,k,:)
                       call inflation_RTPS(analysis,prior,RTP_rate,ensemble_size)
                       OCN_DAstate_data(ivar)%posterior2d(i,k,:)=analysis
                    end do
                 endif
           end do ! i
        end do ! ivar
        !write(logunit,*) trim(subname),' finished'
        flush(logunit)
    end subroutine mpas_relax_inflation
!==============================================================================
!BOP
!  ! routine relax2prior_inflation()
!   relaxiation inflation to prior spread
    subroutine relax_inflation()
! EOP
    use sEnKF_mod,     only :inflation_RTPS
    
!BOC
        implicit none
        ! --- Local variables ---
        character(*), parameter :: subName =   '(relax2prior_inflation) '
        integer                 :: i,j,k,ivar
        real(r8)                :: RTP_rate
        real(r8)                :: prior(ensemble_size),analysis(ensemble_size)

        !write(logunit,*) trim(subname),'inflation of RTPS'
        flush(logunit)
        do ivar  = 1, n_OCN_state
           if(.not. OCN_DAstate_data(ivar)%update) cycle
           RTP_rate=OCN_DAstate_data(ivar)%r2p_fac
           !write(logunit,*) trim(subname),trim(OCN_DAstate_data(ivar)%name),'RTPS rate',RTP_rate
           if(RTP_rate <= tiny(RTP_rate)) cycle 
           do i=da_ibs,da_ibe
              do j=da_jbs,da_jbe
                 if (OCN_DAstate_data(ivar)%state_dim==2) then
                     prior=OCN_DAstate_data(ivar)%prior2d(i,j,:)
                     analysis=OCN_DAstate_data(ivar)%posterior2d(i,j,:)
                     call inflation_RTPS(analysis,prior,RTP_rate,ensemble_size)
                     OCN_DAstate_data(ivar)%posterior2d(i,j,:)=analysis
                 endif
                 if (OCN_DAstate_data(ivar)%state_dim==3) then
                    do k=1,DA_nvert
                       prior=OCN_DAstate_data(ivar)%prior(i,j,k,:)
                       analysis=OCN_DAstate_data(ivar)%posterior(i,j,k,:)
                       call inflation_RTPS(analysis,prior,RTP_rate,ensemble_size)
                       OCN_DAstate_data(ivar)%posterior(i,j,k,:)=analysis
                    end do
                 endif
              end do ! 
           end do ! i
        end do ! ivar
        !write(logunit,*) trim(subname),'inflation of RTPS finished'
        flush(logunit)
    end subroutine relax_inflation


end module ECDA_mpaso_state_mod


        




                 
