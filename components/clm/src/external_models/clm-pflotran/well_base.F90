module Well_Base_class
#ifdef WELL_CLASS

  use PFLOTRAN_Constants_module
  use WellSpec_Base_class

  use Connection_module
  use Hydrostatic_Common_module


  implicit none

  private

#include "petsc/finclude/petscsys.h"

  type, public :: well_base_type 
    !PetscInt :: id                            ! well id is not needed - use coupler id
    character(len=MAXWORDLENGTH) :: name       ! same as associated coupler name
    PetscMPIInt :: comm                        ! well communicator
    PetscMPIInt :: group                       ! well group 
    PetscMPIInt :: cntr_rank                   ! rank where the controlling connection is located
    PetscInt, pointer  :: w_rank_conn(:)       ! number of well conns in each well rank
    PetscInt, pointer  :: disp_rank_conn(:)    ! conns stride for each well rank
    PetscReal :: z_pw_ref                      ! well elevation where the reference pressure is defined
    PetscInt  :: iwconn_ref                    ! index of reference connection in w_conn_z(:)
    PetscInt :: well_num_conns                 ! number of connection for the entire well - not only portion local to the rank
    PetscBool :: cntrl_conn                    ! true if the local well segment contains the control connection  
    PetscInt  :: cntrl_lcell_id                ! if(cntrl_conn) = local (non-ghosted) cell id of control connection, otherwise -999  
    PetscReal, pointer :: conn_factors(:)      ! well connection factors
    PetscInt, pointer  :: conn_drill_dir(:)    ! connection drilling directions
    PetscInt, pointer  :: conn_status(:)       ! connection status (can be open, closed or auto)
    PetscInt, pointer  :: w_conn_order(:)      ! connection order for ascending elevation
    PetscInt, pointer  :: conn_l2w(:)          ! map the list of local well conns to well conns 
    PetscReal, pointer :: w_conn_z(:)          ! all well connection elevations ordered by ascending z
    class(one_dim_grid_type), pointer :: fine_grid !vertical finer grid for hydrostatic pressure computation - order for ascending z
    class(well_spec_base_type), pointer :: spec  !well_spec pointer
    type(connection_set_type), pointer :: connection_set !pointer to well connection_set
  contains  ! add here type-bound procedure 
    procedure, public :: PrintMsg => PrintBase
    procedure, public :: ConnInit => WellBaseConnInit
    procedure, public :: ExplUpdate => BaseExplUpdate
    procedure, public :: PrintOutputHeader => PrintOutputHeaderWellBase
    procedure, public :: ExplRes => WellBaseExplRes
    procedure, public :: ExplJDerivative => WellBaseExplJDerivative
    procedure, public :: InitRun => WellBaseInitRun
    procedure, public :: InitTimeStep => WellBaseInitTimeStep
    procedure, public :: Output => BaseOutput  
    procedure, public :: Setup
    procedure, public :: WellFactorUpdate
    procedure, public :: OneDimGridVarsSetup => WellBase1DGridVarsSetup
    procedure, public :: DataOutput => BaseDataOutput
    procedure  :: InitWellZRefCntrlConn
    procedure  :: WellConnSort
  end type  well_base_type

  public :: WellBaseInit, WellBaseConnInit, BaseWellStrip, &
            WellBaseInitRun, WellBaseInitTimeStep

contains

! ************************************************************************** !

subroutine PrintBase(this)

  implicit none

  class(well_base_type) :: this

  write(*,*) "Well Base Printing message"

end subroutine PrintBase

! ************************************************************************** !

subroutine PrintOutputHeaderWellBase(this,output_option,file_unit)

  use Output_Aux_module

  implicit none

  class(well_base_type) :: this
  type(output_option_type), intent(in) :: output_option
  PetscInt, intent(in) :: file_unit


  write(*,*) "PrintOutputHeaderWellBase must be extended"
  stop
 
end subroutine PrintOutputHeaderWellBase

! ************************************************************************** !

subroutine WellBaseInit(this,well_spec,option)
  ! 
  ! Initializes variables/objects in base well class
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 05/18/16
  ! 

  use Option_module

  implicit none

  class(well_base_type) :: this
  class(well_spec_base_type), pointer :: well_spec
  type(option_type) :: option

  this%name='';
  this%comm=0;                         
  this%group=0;
  this%cntr_rank=0;                       
  this%z_pw_ref = 0.0;
  this%iwconn_ref = -111
  this%well_num_conns = 0;
  this%cntrl_conn = PETSC_FALSE 
  this%cntrl_lcell_id = -999

  nullify(this%w_rank_conn);   
  nullify(this%disp_rank_conn);
  nullify(this%conn_factors);  
  nullify(this%conn_drill_dir); 
  nullify(this%conn_status); 
  nullify(this%w_conn_order); 
  nullify(this%conn_l2w); 
  nullify(this%w_conn_z);
  nullify(this%fine_grid);

  nullify(this%spec);
  this%spec => well_spec

  nullify(this%connection_set);

end subroutine WellBaseInit

! *************************************************************************** !

subroutine Setup(this,connection_set,grid,option)
  ! 
  ! Initilize well connection arrays
  ! Determin well controll connection and rank
  ! Reordering well connection for ascending elevation 
  ! 
  ! Author: Paolo Orsini - OpenGoSim
  ! Date: 6/23/2015

  use Connection_module
  use Grid_module
  use Option_module

  implicit none

  class(well_base_type) :: this
  type(connection_set_type), pointer :: connection_set
  type(grid_type), pointer :: grid 
  type(option_type) :: option

  PetscReal :: min_z, max_z
  PetscReal :: datum(3)

  call this%ConnInit(connection_set%num_connections,option)

  call this%InitWellZRefCntrlConn(grid,connection_set,option)

  call this%WellConnSort(grid,connection_set,option)

  if (connection_set%num_connections > 0) then
    !create well one-dimensional fine grid for hydrostatic computation
    min_z = this%w_conn_z(1) - 1.0d0                  !1m buffer
    max_z = this%w_conn_z(this%well_num_conns) +1.0d0 !1m buffer
    !max_z = max(grid%z_max_global,datum(Z_DIRECTION))+1.d0 ! add 1m buffer
    !min_z = min(grid%z_min_global,datum(Z_DIRECTION))-1.d0
    datum = 0.0d0
    datum(Z_DIRECTION) = this%z_pw_ref  

    this%fine_grid => CreateOneDimGrid(min_z,max_z,datum)

    call this%OneDimGridVarsSetup(option)
  end if

end subroutine Setup

! ************************************************************************** !

subroutine WellBase1DGridVarsSetup(this,option)

  use Option_module

  implicit none

  class(well_base_type) :: this
  type(option_type) :: option 

  print *, "WellBase1DGridVarsSetup must be extended"
  stop  

end subroutine WellBase1DGridVarsSetup

! ************************************************************************** !

subroutine WellBaseConnInit(this,num_connections,option)
  ! 
  ! Allocate and initilize well_base connections arrays
  ! 
  ! Author: Paolo Orsini - OpenGoSim
  ! Date: 6/23/2015

  use Option_module

  implicit none
  
  class(well_base_type) :: this
  PetscInt, intent(in) :: num_connections 
  type(option_type) :: option  
  
  if(num_connections > 0) then

    allocate(this%conn_factors(num_connections)); 
    allocate(this%conn_drill_dir(num_connections)); 
    allocate(this%conn_status(num_connections));

    select case(this%spec%well_fact_itype)
      case(WELL_FACTOR_CONST)
        this%conn_factors = this%spec%const_well_fact
      case(WELL_FACTOR_PEACEMAN) 
        ! do nothing: the well connection factors will be initialize
        ! witht from other well variables
    end select

    if(this%spec%input_const_drill_dir) then
      this%conn_drill_dir = this%spec%const_drill_dir
    else
      option%io_buffer = "WellConnInit: " // &
               "wells with mix drilling directions not implemented"
      call printErrMsg(option)
    end if 

    select case(this%spec%status)
      case(WELL_STATUS_OPEN)
        this%conn_status = CONN_STATUS_OPEN
      case(WELL_STATUS_CLOSE)
        this%conn_status = CONN_STATUS_CLOSE
    end select

  end if

end subroutine WellBaseConnInit

! *************************************************************************** !
subroutine InitWellZRefCntrlConn(this,grid,connection_set,option)
  ! 
  ! - If not given in the input it computes z_pw_ref where pw_ref is defined 
  !   as the elevation of the shallower well connection
  !   Note: in case of a perfect horizontal well segment (same elevations)
  !         the first local cell id is taken for max_z. The same situation
  !         between more well ranks, returns the rank with smaller id 
  ! - If z_pw_ref is given in the input it identifies the well connection 
  !   where pw_ref is defined. Required for implicit wells (pw_ref sol. var)
  ! 
  ! Author: Paolo Orsini - OpenGoSim
  ! Date: 6/13/2015
  ! 

  use Connection_module
  use Grid_module
  use Option_module

  implicit none

  class(well_base_type) :: this
  type(connection_set_type), pointer :: connection_set
  type(grid_type), pointer :: grid 
  type(option_type) :: option

  PetscReal :: max_z !Zmax
  PetscReal :: max_z_snd(2), max_z_rcv(2)
  PetscInt :: output_rank
  PetscInt :: iconn, local_id, ghosted_id
  PetscInt :: lcell_id
  PetscInt :: ierr
  PetscMPIInt :: w_myrank, w_cntr_myrank 

  print *, "InitZref ", "myrank =", option%myrank
  print *, "InitZref ", "num_connections =", connection_set%num_connections
  ! do nothing for empty wells
  if( connection_set%num_connections == 0 ) return

  call MPI_Comm_rank(this%comm, w_myrank, ierr )

  ! all segments belonging to a well have the same input_z_pw_ref
  if( .not.this%spec%input_z_pw_ref ) then
    max_z = -1.0d20 ! initialize to a large negative number
    do iconn=1,connection_set%num_connections
      local_id = connection_set%id_dn(iconn);
      ghosted_id = grid%nL2G(local_id);
      if(grid%z(ghosted_id) > max_z) then
        max_z = grid%z(ghosted_id)
        lcell_id = local_id  
      end if
    end do
    max_z_snd(1) = max_z 
    !max_z_snd(2) = option%myrank
    max_z_snd(2) = w_myrank
    call MPI_ALLREDUCE(max_z_snd, max_z_rcv, 1, MPI_2DOUBLE_PRECISION, &
                       MPI_MAXLOC,this%comm, ierr)
    this%z_pw_ref = max_z_rcv(1)
    !output_rank = max_z_rcv(2) ! casting from double to interger (see MPI)
    w_cntr_myrank = int(max_z_rcv(2)) ! available to all ranks
    this%cntr_rank = w_cntr_myrank

#ifdef WELL_DEBUG
      print *, "InitWellZRef well_rank", w_myrank, &
               "option-myrank=",option%myrank
      !print *, "InitWellZRefCntrlConn - option-myrank", option%myrank    
#endif

    if(w_myrank == w_cntr_myrank ) then  
      !this%cntrl_lcell_id = local_id
      this%cntrl_lcell_id = lcell_id
      this%cntrl_conn = PETSC_TRUE

#ifdef WELL_DEBUG
      !print *, "contrl cell id =",this%cntrl_lcell_id,"well myrank =", w_myrank
      !can't call Bcast here - not invoked in ranks where this%cntrl_conn==false
      !call MPI_Bcast ( this%cntrl_lcell_id,1,MPI_INTEGER,w_myrank, &
      !                 this%comm, ierr )
#endif

    end if

    call MPI_Bcast ( this%cntrl_lcell_id,1,MPI_INTEGER,w_cntr_myrank, &
                     this%comm, ierr )
#ifdef WELL_DEBUG
      print *, "after Bcast cntrl_cell_id =",this%cntrl_lcell_id, &
               "well myrank =", w_myrank
      print *, "rcv max_z =",max_z_rcv(1),"rcv rank_max =",w_cntr_myrank, &
               "myrank =", option%myrank
      print *, "contrl cell id =",this%cntrl_lcell_id,"myrank =", option%myrank    
#endif
   

  else
   ! CAN BE IMPLEMENTED later, this is another option to assign Zref
   ! - loop over the well connections
   ! - compute minimum distance between the well grid block and Zref, 
   !   record lcell_id for minmimum distance
   ! - MPI_Allreduce (...MPI_MINLOC..) with dist_min
   ! - identify well segment with minimum distance, and corresponding rank 
   !   determine which segement contain the control connection  
   option%io_buffer = "InitWellZRefCntrlConn: " // &
             "Well Z_ref from inut file not yet implemented"
   call printErrMsg(option)
   
  end if
  
end subroutine InitWellZRefCntrlConn

! *************************************************************************** !
subroutine WellConnSort(this,grid,connection_set,option)
  ! 
  ! - Sort well connections for growing elevations, to facilitate the 
  !   conputation of the well hydorstatic pressure corrections 
  ! - create an array of well connection elevations for each patch
  ! - allocate arrays related to well connection orders 
  !
  ! Author: Paolo Orsini - OpenGoSim
  ! Date: 30/6/2015
  ! 

  use Connection_module
  use Grid_module
  use Option_module

  implicit none

  class(well_base_type) :: this
  type(connection_set_type), pointer :: connection_set
  type(grid_type), pointer :: grid 
  type(option_type) :: option

  PetscInt :: iconn, irank
  PetscInt :: ghosted_id, local_id
  PetscMPIInt :: w_myrank
  PetscInt :: well_num_conns
  PetscInt :: bub_step
  PetscInt :: ierr 
  PetscInt :: iconn_ref 
  PetscReal, pointer :: conns_z_snd(:)
  PetscReal :: tmp_real
  PetscInt :: tmp_int
  PetscInt, pointer :: w_rank_conn(:) ! each rank has its own num_connections 
  PetscInt, pointer :: disp_rank_conn(:) ! conns displacement  
  PetscInt :: num_w_ranks

  if( connection_set%num_connections == 0 ) return


  call MPI_Comm_size( this%comm, num_w_ranks, ierr)

  allocate(this%w_rank_conn(num_w_ranks))
  this%w_rank_conn = -1
  allocate(this%disp_rank_conn(num_w_ranks)) 
  this%disp_rank_conn = -999


  call MPI_Allgather(connection_set%num_connections, 1, MPI_INTEGER, &
                    this%w_rank_conn, 1, MPI_INTEGER, this%comm, ierr);

  this%disp_rank_conn(1) = 0
  do irank=0,num_w_ranks-1
    if(irank > 0) then
      this%disp_rank_conn(irank+1) = this%disp_rank_conn(irank) + &
                                     this%w_rank_conn(irank)
    end if
  end do

  call MPI_Comm_rank(this%comm, w_myrank, ierr )

#ifdef WELL_DEBUG
  print *,"Well_myrank", w_myrank
#endif
  
  allocate(this%conn_l2w(connection_set%num_connections))
  do iconn=1,connection_set%num_connections
    ! note that indices in this%conn_l2w are 1 based, 
    ! while those in disp_rank_conn are zero based
    this%conn_l2w(iconn) = this%disp_rank_conn(w_myrank+1) + iconn
  end do

#ifdef WELL_DEBUG
  !print *,"w_rank_array", w_rank_conn(1:num_w_ranks)
  print *,"w_rank_array", this%w_rank_conn(1:num_w_ranks) 
  !print *,"disp_array", disp_rank_conn(1:num_w_ranks)
  print *,"disp_array", this%disp_rank_conn(1:num_w_ranks)
  print *,"conn_l2w_array", this%conn_l2w(1:connection_set%num_connections)
#endif

  call MPI_ALLREDUCE(connection_set%num_connections, well_num_conns, 1, &
                     MPI_INTEGER, MPI_SUM,this%comm, ierr)  


#ifdef WELL_DEBUG
  print *,"myrank", option%myrank," num of well connections =", well_num_conns
#endif
  
  ! each rank need to know the number of connection on the entire well 
  this%well_num_conns = well_num_conns

  allocate(this%w_conn_order(well_num_conns)) 
  do iconn=1,well_num_conns
    this%w_conn_order(iconn) = iconn
  end do

  allocate(this%w_conn_z(well_num_conns)) 
  this%w_conn_z = 0.0d0

  allocate(conns_z_snd(connection_set%num_connections))
  conns_z_snd = -1.d30


  do iconn=1,connection_set%num_connections
    local_id = connection_set%id_dn(iconn);
    ghosted_id = grid%nL2G(local_id);
    conns_z_snd(iconn) = grid%z(ghosted_id)
    if(this%cntrl_conn) then
      if(this%cntrl_lcell_id == local_id) then
        iconn_ref = iconn 
#ifdef WELL_DEBUG
        print *,"w_myrank_zref", w_myrank
        print *,"iconn_ref_z", grid%z(ghosted_id)
#endif          
      end if
    end if
  end do

  !return for debugging - from here

#ifdef WELL_DEBUG
  print *,"Before Well Conn MPI_Allgatherv" 
#endif
  !MPI_Allgatherv because each well segment can have a different number of conns
  call MPI_Allgatherv(conns_z_snd, connection_set%num_connections, &
                   MPI_DOUBLE_PRECISION, this%w_conn_z, this%w_rank_conn, &
                   this%disp_rank_conn,MPI_DOUBLE_PRECISION, this%comm, ierr)
#ifdef WELL_DEBUG
  print *,"after Well Conn MPI_Allgatherv" 
#endif

#ifdef WELL_DEBUG
  print *,"w_conn_z_array", this%w_conn_z(1:well_num_conns)
#endif

  !return for debugging
  !return


  ! perform well conns sorting for ascending Z - simple bubble sort
  ! the num_conns of a single well is usually less than 100 
  ! reorder this%w_conn_z and create order array

  do bub_step=1,well_num_conns
    do  iconn = 1, well_num_conns - bub_step
      if(this%w_conn_z(iconn) > this%w_conn_z(iconn + 1)) then
        tmp_real = this%w_conn_z(iconn)
        tmp_int = this%w_conn_order(iconn)
        this%w_conn_z(iconn) = this%w_conn_z(iconn + 1)
        this%w_conn_z(iconn + 1) = tmp_real
        this%w_conn_order(iconn) = this%w_conn_order(iconn + 1)
        this%w_conn_order(iconn + 1) = tmp_int
      end if
    end do
  end do

#ifdef WELL_DEBUG
  print *,"sorted w_conn_z_array", this%w_conn_z(1:well_num_conns)
  print *,"order_array", this%w_conn_order(1:well_num_conns)
#endif

  ! this is the rank containing the controlling connection
  if(this%cntrl_conn) then
    this%iwconn_ref = this%w_conn_order( this%conn_l2w(iconn_ref) )
#ifdef WELL_DEBUG
    print *,"iwconn_ref", this%iwconn_ref
#endif
  end if
   
  call MPI_Bcast ( this%iwconn_ref,1,MPI_INTEGER,this%cntr_rank, &
                      this%comm, ierr )


  deallocate(conns_z_snd); nullify(conns_z_snd)

end subroutine WellConnSort

! *************************************************************************** !

subroutine WellFactorUpdate(this,grid,connection_set,material_auxvars,option)
  ! 
  ! Computes the well factor for each well connection
  !
  ! Author: Paolo Orsini - OpenGoSim
  ! Date: 05/25/2016
  ! 

  use Connection_module
  use Grid_module
  use Option_module
  use Material_Aux_class, only : material_auxvar_type, perm_xx_index, &
                                 perm_yy_index, perm_zz_index  

  implicit none

  class(well_base_type) :: this
  type(connection_set_type), pointer :: connection_set
  type(grid_type), pointer :: grid 
  ! material_auxvars(:) should be declared as a class, because material_auxvar_type
  ! is a class (it has bound-procedure, and it is designed to be extendable)
  ! however, there is a bug reported in for gfortran <= 4.9 for which 
  ! class arrays transfers do not behave corretly (array values differ between the 
  ! calling program and the called functions). See link below   
  ! https://wiki.ucar.edu/display/ccsm/Fortran+Compiler+Bug+List#FortranCompilerBugList-Fortran2003
  ! The work around suggested is: 
  ! "if possible, declare the offending variable using 'type' rather than 'class'"
  ! In here it seems to work, however id the intel compiler won't accept this, there
  ! a second option: passing a scalar object "material_auxvar" 
  ! instead of "material_auxvars(:)". This will requiring taking the connection 
  ! loop outside.
  !class(material_auxvar_type), intent(in) :: material_auxvars(:)
  type(material_auxvar_type), intent(in) :: material_auxvars(:)
  !class(material_auxvar_type) :: material_auxvars(:)
  !class(material_auxvar_type), pointer :: material_auxvars(:)
  type(option_type) :: option

  PetscInt :: iconn, local_id, ghosted_id

  PetscReal :: dx,dy,dz
  PetscReal :: dx1,dx2,dh,k1,k2,r0  

  if (this%spec%well_fact_itype == WELL_FACTOR_CONST) then
    this%conn_factors = this%spec%const_well_fact
  else  
    do iconn=1,connection_set%num_connections
      local_id = connection_set%id_dn(iconn);
      ghosted_id = grid%nL2G(local_id);
      if (associated(this%spec%dxyz_const)) then
         dx = this%spec%dxyz_const(1)
         dy = this%spec%dxyz_const(2)
         dz = this%spec%dxyz_const(3)         
      else 
         dx = grid%structured_grid%dx(ghosted_id)
         dy = grid%structured_grid%dy(ghosted_id)
         dz = grid%structured_grid%dz(ghosted_id)
      end if  
      select case(this%conn_drill_dir(iconn))
        case(X_DIRECTION) 
          dx1 = dy
          dx2 = dz
          dh = dx
          k1 = material_auxvars(ghosted_id)%permeability(perm_yy_index)
          k2 = material_auxvars(ghosted_id)%permeability(perm_zz_index)
        case(Y_DIRECTION)
          dx1 = dx
          dx2 = dz
          dh = dy
          k1 = material_auxvars(ghosted_id)%permeability(perm_xx_index)        
          k2 = material_auxvars(ghosted_id)%permeability(perm_zz_index)
        case(Z_DIRECTION)
          dx1 = dx
          dx2 = dy
          dh = dz
          k1 = material_auxvars(ghosted_id)%permeability(perm_xx_index)        
          k2 = material_auxvars(ghosted_id)%permeability(perm_yy_index)
      end select

#ifdef WELL_DEBUG
      print *,"permx idx = ", perm_xx_index
      print *,"permx idy = ", perm_yy_index
      print *,"permx idz = ", perm_zz_index
      print *,  material_auxvars(ghosted_id)%permeability(perm_xx_index)
      print *,  material_auxvars(ghosted_id)%permeability(perm_yy_index)
      print *,  material_auxvars(ghosted_id)%permeability(perm_zz_index)
#endif     

      r0 = (dx1**2.d0 * (k2/k1)**0.5d0 + dx2**2.d0 * (k1/k2)**0.5d0)**0.5d0 * &
           0.28d0 / ((k2/k1)**0.25d0 + (k1/k2)**0.25d0)
  
      this%conn_factors(iconn) = 2.0d0 * PI * dh * dsqrt(k1*k2) * &
                         this%spec%theta_frac / &
                         ( dlog(r0/this%spec%radius) + this%spec%skin_factor )

#ifdef WELL_DEBUG
      print *,"conn_id = ", iconn," conn fact = ", this%conn_factors(iconn)
#endif
    end do !end loop on well connections

  end if !end if well factor type


end subroutine WellFactorUpdate

! *************************************************************************** !

!subroutine BaseExplUpdate(this,grid,ss_fluxes,option)
subroutine BaseExplUpdate(this,grid,option)
  ! 
  ! - Update FlowEnergy well vars
  ! - Perform a limit on well checks 
  ! - Update well control variable in case of switch when a limit is reached
  !
  ! Author: Paolo Orsini (OGS)
  ! Date: 6/03/2016
  ! 

  use Grid_module
  use Option_module

  implicit none

  class(well_base_type) :: this
  type(grid_type), pointer :: grid
  type(option_type) :: option

  print *, "Well => BaseExplUpdate must be extended"
  stop  

end subroutine BaseExplUpdate

! ************************************************************************** !

subroutine WellBaseExplJDerivative(this,iconn,ghosted_id,isothermal, &
                                   energy_equation_index,option,Jac)
  ! 
  ! Computes the well derivatives terms for the jacobian
  ! 
  ! Author: Paolo Orsini
  ! Date: 6/06/16
  ! 

  use Option_module

  implicit none

  class(well_base_type) :: this
  PetscInt :: iconn
  PetscInt :: ghosted_id 
  PetscBool :: isothermal
  PetscInt :: energy_equation_index
  type(option_type) :: option
  !PetscReal :: Jac(option%nflowdof,option%nflowdof)
  PetscReal :: Jac(:,:)

  print *, "WellFlowEnergyExplJDerivative must be extended"
  stop  
  
end subroutine WellBaseExplJDerivative

! ************************************************************************** !

subroutine WellBaseInitRun(this,grid,material_auxvars,output_option,option)
  ! 
  ! Initialise well for a run
  ! 
  ! Author: Paolo Orsini
  ! Date: 4/08/16
  ! 

  use Grid_module
  use Material_Aux_class, only : material_auxvar_type
  use Output_Aux_module
  use Option_module

  implicit none

  class(well_base_type) :: this
  type(grid_type), pointer :: grid
  type(material_auxvar_type), intent(in) :: material_auxvars(:)
  type(output_option_type), intent(in) :: output_option
  type(option_type) :: option

  PetscMPIInt :: cur_w_myrank
  PetscInt :: ierr
  character(len=MAXSTRINGLENGTH) :: wfile_name

                                              !can get rid of connection_set   
  call this%WellFactorUpdate(grid,this%connection_set, &
                             material_auxvars,option)
     
  ! create well outputfile - should be moved into a well class
  ! For now open files to print the well variables by default 
  ! TODO: add to well_spec user options to control well printing
  call MPI_Comm_rank(this%comm, cur_w_myrank, ierr )  
  if (this%cntr_rank == cur_w_myrank ) then
    wfile_name = trim(option%global_prefix) // "_" // &
                      trim(this%name) // ".tec" 
    open(unit=IUNIT_TEMP,file=wfile_name)
    call this%PrintOutputHeader(output_option,IUNIT_TEMP)
    close(unit=IUNIT_TEMP)
  end if


end subroutine WellBaseInitRun

! ************************************************************************** !

subroutine WellBaseInitTimeStep(this,grid,material_auxvars,option)
  ! 
  ! Initialise well time step
  ! 
  ! Author: Paolo Orsini
  ! Date: 4/08/16
  ! 

  use Grid_module
  use Material_Aux_class, only : material_auxvar_type
  use Option_module

  implicit none
  
  class(well_base_type) :: this
  type(grid_type), pointer :: grid
  type(material_auxvar_type), intent(in) :: material_auxvars(:)
  type(option_type) :: option

  if (option%update_flow_perm) then
    call this%WellFactorUpdate(grid,this%connection_set, &
                               material_auxvars,option)
  end if

end subroutine WellBaseInitTimeStep

! ************************************************************************** !

subroutine WellBaseExplRes(this,iconn,ss_flow_vol_flux,isothermal, &
                                ghosted_id, dof,option,res)
  ! 
  ! Compute residual term for a TOilIms Water injector
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 06/06/16
  ! 
  use Option_module
  
  implicit none

  class(well_base_type) :: this
  PetscInt :: iconn
  PetscBool :: isothermal
  PetscInt :: ghosted_id, dof
  type(option_type) :: option
  !PetscReal :: Res(1:option%nflowdof)
  !PetscReal :: ss_flow_vol_flux(1:option%nphase)
  PetscReal :: Res(:)
  PetscReal :: ss_flow_vol_flux(:)


  print *, "WellBaseExplRes must be extended"
  stop  


end subroutine  WellBaseExplRes

! ************************************************************************** !

subroutine BaseOutput(this,output_file_unit,output_option,option)
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 05/18/16
  ! 
  use Option_module
  use Output_Aux_module

  implicit none

  class(well_base_type) :: this
  PetscInt, intent(in) :: output_file_unit
  type(output_option_type), intent(in) :: output_option
  type(option_type) :: option

  print *, "Well BaseOutput must be extended"
  stop  

end subroutine BaseOutput

! ************************************************************************** !

subroutine BaseDataOutput(this,grid,src_name,option)
  !
  ! Write well pressure and perforated grid lock profile
  ! Overwrites previous file - currently for debugging
  ! TO DO - should add control at which time step to print the profiles 
  !
  ! Author: Paolo Orsini (OpenGoSim)  
  ! Date : 6/20/2016

  use Grid_module
  use Option_module
  
  implicit none

  class(well_base_type) :: this
  type(grid_type), pointer :: grid
  character(len=MAXWORDLENGTH) :: src_name
  type(option_type) :: option

  print *, "Well BaseDataOutput must be extended"
  stop  

end subroutine BaseDataOutput

!*****************************************************************************!

subroutine BaseWellStrip(well)
  !
  ! Strip well_flow and all its parent members
  !
  ! Author: Paolo Orsini (OpenGoSim)  
  ! Date : 6/14/2016
  !

  use Utility_module, only : DeallocateArray 

  implicit none

  class(well_base_type) :: well

  call DeallocateArray(well%w_rank_conn)
  call DeallocateArray(well%disp_rank_conn)
  call DeallocateArray(well%conn_factors)
  call DeallocateArray(well%conn_drill_dir)
  call DeallocateArray(well%conn_status)
  call DeallocateArray(well%w_conn_order)
  call DeallocateArray(well%conn_l2w)
  call DeallocateArray(well%w_conn_z)

  call DestroyOneDimGrid(well%fine_grid)
  nullify(well%fine_grid)

  !these are pointer only 
  nullify(well%spec)
  nullify(well%connection_set)  
 
end subroutine BaseWellStrip

!*****************************************************************************!

#endif   
end module Well_Base_class
!end of WELL_CLASS


