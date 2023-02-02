!BOP
!   !MODULE: ecda_mpaso_obs_mod
!   !INTERFACE:
module ECDA_mpaso_obs_mod
!   !DESCRIPTION:
!   ! Define the obs variable type and corresponding obsope
!   !
!   !REVISION HISTORY:
!    Jul  2019 Y. Liu <liu6@tamu.edu> initial release
!   !USES:

  use shr_kind_mod, only : r8=> shr_kind_r8
  use ecda_shr_type_mod
  use ECDA_mpaso_state_mod 
  use seq_comm_mct,  only : logunit
  use shr_file_mod,      only : shr_file_getUnit,     &
                                  shr_file_freeUnit,    &
                                  shr_file_setIO,       &
                                  shr_file_getLogUnit,  &
                                  shr_file_setLogUnit,  &
                                  shr_file_getLogLevel, &
                                  shr_file_setLogLevel

  implicit none
  public :: ocn_obs_init, set_obs_type_value, ocn_obs_type_register, &
            observation_deallocate, phys2ij_mercator
  
  integer, parameter, public :: max_OCN_obs_type =  10    ! Currently using T,S,U,V,zeta 
  
  integer, parameter, public :: O_rea_TEMP_ID = 1     ! index of reanalysis ocn temperature
  integer, parameter, public :: O_rea_SALT_ID = 2     ! index of reanalysis ocn salinity
  integer, parameter, public :: O_rea_VEL_ID = 3     ! index of reanalysis ocn velocity
  integer, parameter, public :: O_rea_SSH_ID = 4     ! index of reanalysis ocn surface height


  integer, parameter, public :: O_GRID_MODEL     = 1        ! instrument type for reanalysis with the same grids as ocean model 
  integer, parameter, public :: O_GRID_NONMODEL =  2         ! grid observation but not the model grid
  integer, parameter, public :: O_profile =  3         ! grid observation but not the model grid

  integer, public :: NUM_OCN_obs_type
!  type(ECDA_obs_type),public         ::  OCN_obs_data(max_OCN_obs_type) 
  type(ECDA_obs_type),allocatable,public         ::  OCN_obs_data(:) 

contains

!==============================================================================
!BOP
!   !ROUTINE: ocn_obs_init
!   initialize the OCN obs for DA
!   define the observations
!   !INTERFACE:

    subroutine ocn_obs_init()
!EOP

!BOC
        implicit none
        ! --- Local variables ---
        character(*), parameter :: subName =   '(ocn_obs_init) '
        integer :: iob,i,j,k,ierr
        integer :: inst_type
!        logical :: update_states(n_OCN_state),required_states(n_OCN_state)
!        allocate (OCN_obs_data(max_OCN_obs_type)        
        do iob=1,NUM_OCN_obs_type
           OCN_obs_data(iob)%assimilate=.False.    
           allocate (OCN_obs_data(iob)%required_states(n_OCN_state), stat=ierr)
           allocate (OCN_obs_data(iob)%update_states(n_OCN_state), stat=ierr)
           OCN_obs_data(iob)%required_states = .FALSE.
           OCN_obs_data(iob)%update_states   = .FALSE.
           OCN_obs_data(iob)%grid_type    = O_GRID_MODEL       ! set model grid as default
        enddo

     end subroutine ocn_obs_init
    
!==============================================================================
!BOP
! ! Rountine set_obs_type_value(obs_data,...)
! ! set the obs configures
 
    subroutine set_obs_type_value(ID,obs_name,filename,dm_name,inst_type,assimilate, required_states, update_states)
      
! EOP 
       integer, intent(in)                                    ::      ID
       character(*),optional,intent(in)                       ::      obs_name
       character(*),optional,intent(in)                       ::      filename
       character(*),optional,intent(in)                       ::      dm_name
       integer, optional,intent(in)                           ::      inst_type
       logical, optional,intent(in)                           ::      assimilate, &
                                                                      required_states(n_OCN_state), &
                                                                      update_states(n_OCN_state)
                                                       
      if(present(obs_name)) OCN_obs_data(ID)%obs_name    = obs_name
      if(present(filename)) OCN_obs_data(ID)%filename = filename
      if(present(dm_name)) OCN_obs_data(ID)%dm_name = dm_name
      if(present(inst_type)) OCN_obs_data(ID)%inst_type = inst_type
      if(present(assimilate)) OCN_obs_data(ID)%assimilate = assimilate
      if(present(required_states)) OCN_obs_data(ID)%required_states = required_states
      if(present(update_states)) OCN_obs_data(ID)%update_states = update_states

    end subroutine set_obs_type_value
    
!==============================================================================

!BOP
!   !ROUTINE: ocn_obs_type_register
!   register the ocean observation type  for DA
!
!   !INTERFACE:

  subroutine ocn_obs_type_register()
!EOP
!BOC
        implicit none
        ! --- Local variables ---
        character(*), parameter :: subName =   '(ocn_obs_type_register) '
        character(len = filename_len) :: nml_inpfile = 'ocean_da.namelist'
        integer :: i,j,k,iob,iunit,rc,ierr
        logical :: exists
        logical :: obs_avalaible(max_OCN_obs_type)=.False.
        character(len = type_name_len) ::  assim_obs_types(max_OCN_obs_type) = 'null'
        namelist /obs_kind_nml/assim_obs_types
        
        NUM_OCN_obs_type=0
        ! -read namelist for observation type        
        inquire(file=trim(nml_inpfile), exist=exists)
        if (.not. exists) return
        iunit = shr_file_getUnit()
        open(unit = iunit, file = trim(nml_inpfile), action = 'read')
        read(unit = iunit, nml = obs_kind_nml, iostat = ierr)
        close(unit = iunit)
        call shr_file_freeUnit( iunit )
        !- loop to get total number of assim. obs type
        do i = 1, max_OCN_obs_type      
           if(assim_obs_types(i) == 'null' .or. len_trim(assim_obs_types(i)) == 0 ) exit
           NUM_OCN_obs_type = i
        end do
        if(NUM_OCN_obs_type == 0) return

        !- initialize the observation variabel
        allocate(OCN_obs_data(NUM_OCN_obs_type),stat=ierr)
        call ocn_obs_init()

        !- specific obs type  to assim, read the coresponding namelist file for
        !assim. setting
        do i= 1, NUM_OCN_obs_type  
              OCN_obs_data(i)%obs_name=assim_obs_types(i)   
              OCN_obs_data(i)%assimilate= .true.
              call read_obs_assim_namelist(i,exists)   ! if ierr is not true, error, need revisit @ YL
               RR=max(RR,OCN_obs_data(i)%RR_xy)         ! find maxium localization cutoff for  grid setting
                ! active related state variable for DA
              do k=1,n_OCN_state
                 if( OCN_obs_data(i)%required_states(k).or. OCN_obs_data(i)%update_states(k)) &
                     OCN_DAstate_data(k)%required = .true.
              end do
        end do
        ! Falling off the end is an error
         write(logunit,*) trim(subname), 'finish','NUM_OCN_obs_type',NUM_OCN_obs_type

     end subroutine ocn_obs_type_register


!==============================================================================
!BOP
    !ROUTINe: read_obs_assim_namelist
    !real namelist for assimilation setting of a observation
!   !INTERFACE:
!  subroutine read_obs_assim_namelist(ID,exists)
  subroutine read_obs_assim_namelist(ID,exists)
! @YL need refine for define all values
!EOP
    integer, intent(in)   :: ID
    logical, intent(out) :: exists
!EOC
    character(*), parameter :: subName =   '(read_obs_assim_namelist) '
    character(len=type_name_len)         ::   obs_name              ! state variable  name,ex : ocn_temp,ocn_salt...
    integer                              ::   obs_type            ! observation ID  name,ex : e.g. O_rea_TEMP_ID...
    character(len=filename_len)          ::   filename='null'          ! the obs file name or part
    character(len=var_name_len)          ::   varname='varname'           ! the obs variable name in ncfile 
    character(len=var_name_len)          ::   uncert_name='null'    ! the obs variable name in ncfile 
    character(len=var_name_len)          ::   lonname='longitude'           ! the lon coordinator name in ncfile 
    character(len=var_name_len)          ::   latname='latitude'            ! the lat coordinator name in ncfile 
    character(len=var_name_len)          ::   vertical_name= 'depth'         ! vertical
    character(len=var_name_len)          ::   dm_name= 'records'
    character(len=10)                    ::   init_obs_YMDH='2010010112'     ! the inital time to observations
    integer                              ::   uncert_type = 1       ! variance 0 ,std 1 default, not valid when err_var0 >0 applies
    integer                              ::   grid_type = 1      !  netcdf file, 1,2 are grid file, 1: model grid, 2 other resolution, need read coordinator information, 3: station data, (lon,lat,lev, value,error) 
    integer                              ::   obs_dm = 2                                ! data dimension in file, only used for grid file
    integer                              ::   start_layer = 1                           ! data start  from (vertical layer) for 3D,only used for grid file
    integer                              ::   nlayer                                 ! total vertical layer for 3D,only used for grid file
    integer                              ::   ilayer                                 ! pick vertical layer ,only used for grid file
    integer                              ::   inst_type =0                             ! instrument types are defined by platform class (e.g. MOORING, DROP) and instrument type (XBT, CTD, ...), not active now @YL
    logical                              ::   required_states(n_OCN_state)=.false.   ! the required state variables for obsope 
    logical                              ::   update_states(n_OCN_state)=.false.     ! the updated  state variables for current variable 
    real(r8)                             ::   skip_sfc(n_OCN_state)=0                ! skip surface for ssh lik observation 
    real(r8)                             ::   RR_xy=10.0, RR_z=200.0                 !
    real(r8)                             ::   Mx_PO_dist=10.0                        ! maxium distance allowed between observation and forecast
    real(r8)                             ::   Mx_misfit=25                           ! maxium rate allowed for inflate  observation uncertainty
    real(r8)                             ::   err_var0 = -1.0                          ! specify constant observation error variance (>0). When <=0, input from observation file, default.
    real(r8)                             ::   adj_scale= 1.0               !some variable need adjust scale

    integer :: ierr
    integer :: iunit
    real(r8):: obs_err_var0
    namelist /obs_assim_nml/obs_name,obs_type,filename,inst_type,varname,lonname,latname,vertical_name,&
                           grid_type, obs_dm,uncert_name, start_layer, nlayer, ilayer,&
                           required_states,update_states,skip_sfc,RR_xy,RR_z,Mx_PO_dist,Mx_misfit,err_var0, &
                           init_obs_YMDH,dm_name,uncert_type,adj_scale

  !  write(logunit,*) trim(subname),' obs namelist: ',trim(filename),',',trim(lonname),',',trim(varname),',',trim(vertical_name),',',trim(latname) &
  !                    ,',', trim(uncert_name)
    inquire(file=trim(OCN_obs_data(ID)%obs_name)//'.nml', exist=exists)
    if (.not. exists) then 
        write(logunit,*) trim(subname),'error: can not find mamelist file: ',trim(OCN_obs_data(ID)%obs_name)//'.nml'
       return
    endif
    iunit = shr_file_getUnit()
    open(unit = iunit, file=trim(OCN_obs_data(ID)%obs_name)//'.nml', action = 'read')
    read(unit = iunit, nml = obs_assim_nml, iostat = ierr)
    close(unit = iunit)
    call shr_file_freeUnit(iunit)
    if(trim(obs_name) .ne. trim(OCN_obs_data(ID)%obs_name))then
        write(logunit,*) trim(subname),'obs_name do name match ',trim(obs_name),trim(OCN_obs_data(ID)%obs_name)
        return    
    endif
   ! write(logunit,*) trim(subname),' obs namelist: ',trim(filename),',',trim(lonname),',',trim(varname),',',trim(vertical_name),',',trim(latname) &
   !                   ,',', trim(uncert_name)
        

    allocate (OCN_obs_data(ID)%required_states(n_OCN_state), stat=ierr)
    allocate (OCN_obs_data(ID)%update_states(n_OCN_state), stat=ierr)
    allocate (OCN_obs_data(ID)%skip_sfc(n_OCN_state), stat=ierr)

    OCN_obs_data(ID)%obs_type  = obs_type
    OCN_obs_data(ID)%filename  = filename
    OCN_obs_data(ID)%inst_type  = inst_type
    OCN_obs_data(ID)%dm_name   = dm_name 
    OCN_obs_data(ID)%varname   = varname
    OCN_obs_data(ID)%uncert_name   = uncert_name
    OCN_obs_data(ID)%init_obs_YMDH   = init_obs_YMDH
    
    OCN_obs_data(ID)%required_states = required_states
    OCN_obs_data(ID)%update_states = update_states
    OCN_obs_data(ID)%grid_type = grid_type
    OCN_obs_data(ID)%obs_dm = obs_dm
    OCN_obs_data(ID)%lonname = lonname
    OCN_obs_data(ID)%latname = latname
    OCN_obs_data(ID)%vertical_name = vertical_name
    OCN_obs_data(ID)%dm_name = dm_name
    OCN_obs_data(ID)%skip_sfc = skip_sfc
    OCN_obs_data(ID)%RR_xy=RR_xy
    OCN_obs_data(ID)%RR_z=RR_z
    OCN_obs_data(ID)%err_var0=err_var0 
    OCN_obs_data(ID)%Mx_PO_dist=Mx_PO_dist
    OCN_obs_data(ID)%Mx_misfit=Mx_misfit
    OCN_obs_data(ID)%uncert_type=uncert_type
    OCN_obs_data(ID)%adj_scale=adj_scale
    if(grid_type .ne. 3) then
        read(init_obs_YMDH(1:4),'(i)') OCN_obs_data(ID)%init_obs_YR
        read(init_obs_YMDH(5:6),'(i)') OCN_obs_data(ID)%init_obs_MO
        read(init_obs_YMDH(7:8),'(i)') OCN_obs_data(ID)%init_obs_DY
        read(init_obs_YMDH(9:10),'(i)') OCN_obs_data(ID)%init_obs_HR
    endif
   if(iamroot_da_all)then
      write(logunit,*) trim(subname),' request/update states: ',required_states,update_states
      write(logunit,*) trim(subname),' obs namelist: ',trim( OCN_obs_data(ID)%filename),',',trim( OCN_obs_data(ID)%lonname),',',trim( OCN_obs_data(ID)%varname),',',trim( OCN_obs_data(ID)%vertical_name)
      write(logunit,*) trim(subname),'finsh obs namelist: ',trim(OCN_obs_data(ID)%obs_name)//'.nml'
    end if 
  end subroutine read_obs_assim_namelist
!==============================================================================
!BOP
!   !ROUTINE: mpaso_obs_operator
!   register the ocean observation type  for DA
!
!   !INTERFACE:
  subroutine mpaso_obs_operator(ID,iob,fobs,obs,obs_var, icell_edge)

!EOP
    integer, intent(in)   :: ID,iob
    real(r8),intent(out)  :: obs,obs_var
    real(r8),intent(out)  :: fobs(ensemble_size)
    integer, intent(out)  ::  icell_edge     ! use for localization @YL

!BOC
    character(*), parameter :: subName =   '(mpaso_obs_operator) '
    integer :: da_nbs,da_nbe
    real(r8),allocatable :: diff_lon(:),diff_lat(:),distance(:)
    integer  :: state_ID
    integer  :: ierr, ilocation(1)
    


    ! only for cell elements now @YL
    da_nbs=da_Cell_nbs
    da_nbe=da_Cell_nbe

    allocate(diff_lon(da_nbs:da_nbe),diff_lat(da_nbs:da_nbe),distance(da_nbs:da_nbe), stat=ierr)


    obs      =   OCN_obs_data(ID)%values(iob)
    obs_var  =   OCN_obs_data(ID)%err_var(iob)    

    ! assume all the grid is profile, first to find the nearest grid, no interptation @YL
    ! when use EAKF sequential, obs operateor works on updated model state ( here posterior), the prior do not change during analysis step
    ! and use to calculate analysis increment and inflation

!= select corresponding state ID
    select case(OCN_obs_data(ID)%obs_type)
      case(O_rea_TEMP_ID)
         state_ID=OCN_TEMP_index
      case(O_rea_SALT_ID)
         state_ID=OCN_SALT_index
      case(O_rea_VEL_ID)
         state_ID=OCN_VEL_index
      case(O_rea_SSH_ID)
         state_ID=OCN_SSH_index
      CASE DEFAULT
         write(logunit,*) trim(subname),'observation type:', OCN_obs_data(ID)%obs_type,'not support now'
    end select
    
!=find the closest cell ( only for cell now)   @YL 
    diff_lon =  abs(OCN_obs_data(ID)%lon(iob)- DA_lonCell)
    diff_lat =  abs(OCN_obs_data(ID)%lat(iob)- DA_latCell)
    distance = diff_lon + diff_lat     ! quick fix @YL
    ilocation=minloc(distance)+ da_nbs-1   ! system think  minloc return with a array. Define an arry  
    icell_edge= ilocation(1)
    ! write(logunit,*) trim(subname),'min distance',minval(distance), distance(icell_edge)
    ! write(logunit,*) trim(subname), &
    ! OCN_obs_data(ID)%lon(iob),OCN_obs_data(ID)%lat(iob),DA_lonCell(icell_edge),DA_latCell(icell_edge)

    ! assume the surface observation for now 
    if(OCN_DAstate_data(state_ID)%state_dim==2)  fobs= OCN_DAstate_data(state_ID)%prior2d(icell_edge,1,:)
    if(OCN_DAstate_data(state_ID)%state_dim==1)  fobs= OCN_DAstate_data(state_ID)%prior1d(icell_edge,:)
!= spatail anf vertical intep need @YL
    
    deallocate(diff_lon,diff_lat,distance)

 
  end subroutine mpaso_obs_operator


!==============================================================================
!BOP
!   !ROUTINE: OCN_obs_operator
!   register the ocean observation type  for DA
!
!   !INTERFACE:

  subroutine OCN_obs_operator(ID,iob,fobs,obs,obs_var)
!EOP
    integer, intent(in)   :: ID,iob
    real(r8),intent(out)  :: obs,obs_var
    real(r8),intent(out)  :: fobs(ensemble_size)

!BOC
    character(*), parameter :: subName =   '(OCN_obs_operator) '
   
    obs      =   OCN_obs_data(ID)%values(iob)
    obs_var  =   OCN_obs_data(ID)%err_var(iob)    
    select case(OCN_obs_data(ID)%grid_type)
      case(1)
         call obs_operator_O_reanalysis(ID,iob,fobs)
      case(2)
         call obs_operator_O_profile(ID,iob,fobs)
      case(3)
         call obs_operator_O_profile(ID,iob,fobs)   
      CASE DEFAULT
         write(logunit,*) trim(subname),'grid type:', OCN_obs_data(ID)%grid_type,'not support now'   
    end select
!   write(logunit,*) trim(subname),'obs:',obs,obs_var,' forec obs: ',fobs   

  end subroutine OCN_obs_operator
!==============================================================================
!BOP
!   !ROUTINE: obs_operator_O_profile
!   obs operator for obs with on model grid. no interpotation need 
!
!   !INTERFACE:

  subroutine obs_operator_O_profile(ID,iob,fobs)
!EOP
    integer, intent(in)   :: ID,iob
    real(r8),intent(out)  :: fobs(ensemble_size)

!BOC
    character(*), parameter :: subName =   '(obs_operator_O_profile) '
    integer  :: state_ID
    integer  :: i , j, k
    integer  :: i0 , j0, k0 , i1,j1,k1
    real(r8) :: ri,rj,rk,ai,aj,ak,depth_diff
    integer  :: kk(1)
! here, analysis grid is consistent with model grids, so values of rx,ry,rz are integer. otherwise, intepotation need. 
! when use EAKF sequential, obs operateor works on updated model state ( here posterior), the prior do not change during analysis step
! and use to calculate analysis increment and inflation

!= select corresponding state ID

    select case(OCN_obs_data(ID)%obs_type)
      case(O_rea_TEMP_ID)
         state_ID=OCN_TEMP_index
      case(O_rea_SALT_ID)
         state_ID=OCN_SALT_index
!      case(O_rea_UVEL_ID)
!         state_ID=OCN_UVEL_index
!      case(O_rea_VVEL_ID)
!         state_ID=OCN_VVEL_index
!      case(O_rea_SSH_ID)
!         state_ID=OCN_SSH_index
      CASE DEFAULT
         write(logunit,*) trim(subname),'observation type:', OCN_obs_data(ID)%obs_type,'not support now'   
    end select

!=   for spatial interpotation, implement for 2 D variable ( zeta here @YL) 
   if(OCN_DAstate_data(state_ID)%state_dim==2) then 
      i=floor(OCN_obs_data(ID)%rx(iob))       !i,j should within the range of current grids
      j=floor(OCN_obs_data(ID)%ry(iob))
      if ( i< da_ibs .or. i>=da_ibe .or. j< da_jbs .or. j>= da_jbe) then
        write(logunit,*) trim(subname), " terror out of range",da_ibs,da_ibe,i, da_jbs,da_jbe,j
        fobs=0.0
        return
      endif 
    !k0=floor(OCN_obs_data(ID)%rz(iob))
      ri= OCN_obs_data(ID)%rx(iob)-i
      rj= OCN_obs_data(ID)%ry(iob)-j
      ai=1-ri
      aj=1-rj
      fobs= OCN_DAstate_data(state_ID)%posterior2d(i,j,:)*ai*aj &
           +OCN_DAstate_data(state_ID)%posterior2d(i+1,j,:)*ri*aj &
           +OCN_DAstate_data(state_ID)%posterior2d(i,j+1,:)*ai*rj &
           +OCN_DAstate_data(state_ID)%posterior2d(i+1,j+1,:)*ri*rj 
       
   endif
    
!=   no spatial interpotation for 3d variable, only vertical
  
   if(OCN_DAstate_data(state_ID)%state_dim==3) then 
      i=nint(OCN_obs_data(ID)%rx(iob))
      j=nint(OCN_obs_data(ID)%ry(iob))
!= vertical interpotation
      kk=minloc(abs(z_r_da(i,j,:)-OCN_obs_data(ID)%depth(iob)))
      k=kk(1)
      depth_diff=z_r_da(i,j,k)-OCN_obs_data(ID)%depth(iob) 
      if(abs(depth_diff) <tiny(z_r_da)) then 
         fobs= OCN_DAstate_data(state_ID)%posterior(i,j,k,:)
         return
      endif
      if(depth_diff > 0.0) then
         k0=k-1
      else
         k0=k+1
      endif
      rk=abs(depth_diff/(z_r_da(i,j,k)-z_r_da(i,j,k0))) 
      ak=1-rk
      fobs= OCN_DAstate_data(state_ID)%posterior(i,j,k,:)*ak &
           +OCN_DAstate_data(state_ID)%posterior(i,j,k0,:)*rk 
   end if
!   write(logunit,*) trim(subname),'obs:',i,j,OCN_obs_data(ID)%values(iob),OCN_obs_data(ID)%err_var(iob),' forec obs: ',fobs   
  end SUBROUTINE obs_operator_O_profile
!==============================================================================
!BOP
!==============================================================================
!BOP
!   !ROUTINE: obs_operator_O_reanalysis
!   obs operator for obs with on model grid. no interpotation need 
!
!   !INTERFACE:

  subroutine obs_operator_O_reanalysis(ID,iob,fobs)
!EOP
    integer, intent(in)   :: ID,iob
    real(r8),intent(out)  :: fobs(ensemble_size)

!BOC
    character(*), parameter :: subName =   '(obs_operator_O_reanalysis) '
    integer  :: state_ID
    integer  :: i , j, k
    integer  :: i0 , j0, k0 , i1,j1,k1
    real(r8) :: ai,aj,ak
! here, analysis grid is consistent with model grids, so values of rx,ry,rz are integer. otherwise, intepotation need. 
! when use EAKF sequential, obs operateor works on updated model state ( here posterior), the prior do not change during analysis step
! and use to calculate analysis increment and inflation

!   for interpotation @YL  
    i0=floor(OCN_obs_data(ID)%rx(iob))
    j0=floor(OCN_obs_data(ID)%ry(iob))
    k0=floor(OCN_obs_data(ID)%rz(iob))
    ai= OCN_obs_data(ID)%rx(iob)-i0
    aj= OCN_obs_data(ID)%ry(iob)-j0
    ak= OCN_obs_data(ID)%rz(iob)-k0
    
    i=nint(OCN_obs_data(ID)%rx(iob))
    j=nint(OCN_obs_data(ID)%ry(iob))
    
    select case(OCN_obs_data(ID)%obs_type)
      case(O_rea_TEMP_ID)
         state_ID=OCN_TEMP_index
      case(O_rea_SALT_ID)
         state_ID=OCN_SALT_index
!      case(O_rea_UVEL_ID)
!         state_ID=OCN_UVEL_index
!      case(O_rea_VVEL_ID)
!         state_ID=OCN_VVEL_index
!      case(O_rea_SSH_ID)
!         state_ID=OCN_SSH_index
      CASE DEFAULT
         write(logunit,*) trim(subname),'observation type:', OCN_obs_data(ID)%obs_type,'not support now'   
    end select


    
   
    if(OCN_DAstate_data(state_ID)%state_dim==2) then 
      ! assume the 2D state variabel is on surface layer for now @YL 
        fobs= OCN_DAstate_data(state_ID)%posterior2d(i,j,:) 
    endif
    if(OCN_DAstate_data(state_ID)%state_dim==3) then
        k=nint(OCN_obs_data(ID)%rz(iob))
        fobs= OCN_DAstate_data(state_ID)%posterior(i,j,k,:)
    end if
  end SUBROUTINE obs_operator_O_reanalysis
!==============================================================================
!BOP
!   !ROUTINE:observation_deallocate 
!   ! after analysis step ,deallocate  variabels in OCN_obs_data
!   !   
!   
  SUBROUTINE observation_deallocate()
!EOP
!BOC
        implicit none
        ! --- Local variables ---
        character(*), parameter :: subName =   '(observation_deallocate) '
        integer                 :: ierr, iobs
        do iobs=1,NUM_OCN_obs_type
         !  if(.not. OCN_obs_data(iobs)%assimilate) cycle
           OCN_obs_data(iobs)%nobs=0     ! initial observation numbers
           deallocate(OCN_obs_data(iobs)%rx,OCN_obs_data(iobs)%ry,OCN_obs_data(iobs)%rz,stat=ierr)
           deallocate(OCN_obs_data(iobs)%lon,OCN_obs_data(iobs)%lat,OCN_obs_data(iobs)%depth,stat=ierr)
           deallocate(OCN_obs_data(iobs)%values,OCN_obs_data(iobs)%err_var,OCN_obs_data(iobs)%qc_flag,stat=ierr)
        enddo
        !write(logunit,*) trim(subname),'Finished'

  end SUBROUTINE observation_deallocate
!==============================================================================
!BOP
!   !ROUTINE: phys2ij_mercator
!   observation lat and lon on model grid. 
!   for  Mercator projection, dlon is even but the dlat is not
!   not consider the periodic boundary condition @YL
!   but  latitude may not even
!
!   !INTERFACE:
      SUBROUTINE phys2ij_mercator(rlon,rlat,ri,rj)
!EOP
      IMPLICIT NONE
      real(r8),INTENT(IN) :: rlon,rlat
      REAL(r8),INTENT(OUT) :: ri,rj
      integer :: i , j,jj(1),ii(1) 
      real(r8)  ::  dlon,dlat,dri,drj   

!-    lontigude, even on degree  
      dlon= lon_da(da_ibs+1,da_jbs)-lon_da(da_ibs,da_jbs)
      ri= (rlon-lon_da(da_ibs,da_jbs))/dlon+da_ibs
      i=nint(ri)

!-   latitude , no even
      jj= minloc(abs(rlat-lat_da(i,:)))
      j= jj(1)+ da_jbs-1  ! find nearest J, da_jbs is start point
      if(rlat<lat_da(i,j)) j=j-1     ! make sure rlat locate @ [j, j+1]
      drj = rlat-lat_da(i,j)
      if(abs(drj) <tiny(drj))then
        rj=j
      else
        dlat=lat_da(i,j+1)-lat_da(i,j)
        rj=j+drj/dlat
      endif     !
       
       RETURN
      END SUBROUTINE phys2ij_mercator

!==============================================================================
!BOP
!   !ROUTINE: com_distll_1
!   copy from Letf  (Takemasa MIYOSHI)
!   use the Vincenty's formulae
!   DISTANCE BETWEEN TWO POINTS (LONa,LATa)-(LONb,LATb)
      SUBROUTINE com_distll_1(alon,alat,blon,blat,dist, chg_unit)

      IMPLICIT NONE
      REAL(r8),INTENT(IN) :: alon
      REAL(r8),INTENT(IN) :: alat
      REAL(r8),INTENT(IN) :: blon
      REAL(r8),INTENT(IN) :: blat
      REAL(r8),INTENT(OUT) :: dist
      REAL(r8),optional,INTENT(in) :: chg_unit 

      REAL(r8),PARAMETER :: r180=1.0d0/180.0d0
      REAL(r8),PARAMETER :: re=6371.3d3        ! earth radius
      REAL(r8) :: lon1,lon2,lat1,lat2
      REAL(r8) :: cosd
      
      ! in common situation, need change from degree to radius, chg_unit= pi /180
      if( present(chg_unit)) then 
         lon1 = alon * chg_unit
         lon2 = blon * chg_unit
         lat1 = alat * chg_unit
         lat2 = blat * chg_unit
      else
         lon1 = alon 
         lon2 = blon
         lat1 = alat
         lat2 = blat 
      end if 

      cosd = SIN(lat1)*SIN(lat2) + COS(lat1)*COS(lat2)*COS(lon2-lon1)
      cosd = MIN( 1.d0,cosd)
      cosd = MAX(-1.d0,cosd)
      dist = ACOS( cosd ) * re

      RETURN
      END SUBROUTINE com_distll_1


!-----------------------------------------------------------------------

end module ECDA_mpaso_obs_mod


             


        




                 
