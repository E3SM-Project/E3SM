module e4d_vars

  implicit none
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscmat.h90"
#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petscviewer.h90"  
#include "petsc/finclude/petscksp.h"
#include "petsc/finclude/petscksp.h90"


  logical :: tank_flag = .false.

  integer, dimension(:), allocatable :: e4d_ranks,pf_e4d_ranks
  integer :: mpi_comm_grp,mpi_e4d_grp,mpi_pfe4d_grp
  integer :: my_wrank,my_pfe4d_rank,n_pfe4drank             !!my mpi rank
  integer :: tn_rank                                        !!number of processes 

  character*40 :: mshfile                                  !!file containing the mesh options
  character*40 :: efile                                    !!survey configuration file
  character*40 :: mapfile
  character*40 :: sigfile                                  !!baseline conductivity file
  character*40 :: list_file                                !!conductivity list file
  character*40 :: csrv_file                                !!survey file
  character*40 :: ccond_file                               !!conductivity file
  character*40 :: log_file                                 !!e4d_log file
  character*40 :: fmn_file                                 !!archies parameters file name

  integer :: my_rank                                       !!my mpi rank
  integer :: ierr                                          !!generall error
  integer :: n_rank                                        !!number of processes 
  integer :: E4D_COMM                                      !!E4D_COMMUNICATOR
  integer :: PFE4D_MASTER_COMM                             !!PF to E4D MASTER COMMUNICATOR
  integer :: ne,tne                                        !!num my electrodes, total electrodes
  integer :: nsig                                          !!number of elements, sigma values
  integer :: nm                                            !!number of measurements
  integer :: nnodes                                        !!number of nodes
  integer :: nelem                                         !!number of elements (same as nsig)
  integer :: nfaces                                        !!number of faces
  integer :: nvals                                         !!number of non-zeros in coupling mat
  integer :: my_ne                                         !!number of electrodes I'm assigned
  integer :: nmy_drows                                     !!number of data in my assembly vector
  integer :: nmap
  integer :: ntime                                         !!number of e4d simulation times
  integer :: mode
  integer, dimension(1) :: i_zpot                          !!ghost node index for tank sims 

  integer :: pfnx                                          !!number of pf cells in x dim
  integer :: pfny                                          !!number of pf cells in y dim
  integer :: pfnz                                          !!number of pf cells in z dim

  real :: gw_sig                                           !!groundwater electrical conductivity
  real :: sw_sig                                           !!surface water electrical condctivity
  real :: FF                                               !!formation factor
  real :: Cbeg,Cend,etm                                    !!timing variables
  real*8 :: e4d_time, pf_time
  
  integer, dimension(:,:), allocatable :: map_inds
  integer, dimension(:,:), allocatable :: s_conf           !!abmn survey configuration
  integer, dimension(:,:), allocatable :: eind             !!electrode assignments
  integer, dimension(:,:), allocatable :: jind             !!element map assignments
  integer, dimension(:), allocatable :: nbounds,zones      !!node boundaries and element zones
  integer, dimension(:,:), allocatable :: elements         !!elements connections
  integer, dimension(:,:), allocatable :: faces            !!face connections
  integer, dimension(:), allocatable :: e_nods             !!indices of electrode nodes
  integer, dimension(:), allocatable :: rows,cols
  integer, dimension(:), allocatable :: trows,tcols
  integer, dimension(:), allocatable :: A_map              !!coupling matrix mapping vector
  integer, dimension(:), allocatable :: S_map              !!Sigma mapping vector
  integer, dimension(:), allocatable :: my_drows           !!rows of my data assemble vector

  real, dimension(:,:), allocatable :: e_pos
  real, dimension(:,:), allocatable :: nodes               !!node positions
  real, dimension(:,:), allocatable :: poles               !!pole solutions
  real, dimension(:,:), allocatable :: FMN                 !!static Archie's parameters
  real, dimension(:), allocatable :: pf_porosity           !!pflotran porosity
  real, dimension(:), allocatable :: pf_tracer             !!pflotran tracer solution
  real, dimension(:), allocatable :: pf_saturation         !!pflotran saturation solution
  real, dimension(:), allocatable :: pf_saturation_0       !!pflotran saturation solution at time 0
  real, dimension(:), allocatable :: pf_temperature        !!pflotran temperature
  real, dimension(:), allocatable :: sigma                 !!element conductivities
  real, dimension(:), allocatable :: dpred                 !!simulated data vector
  real, dimension(:), allocatable :: dobs                  !!observed data
  real, dimension(:), allocatable :: sd                    !!observed data standard deviations
  real, dimension(:), allocatable :: my_dvals              !!values in my data assembly vector
  real, dimension(:), allocatable :: map
  real, dimension(:), allocatable :: base_sigma            !!baseline element conductivity
  real, dimension(:), allocatable :: ffac                  !!formation factor

  real, dimension(:), allocatable :: pfxcb                 !!pf cell boundaries in x
  real, dimension(:), allocatable :: pfycb                 !!pf cell boundaries in y
  real, dimension(:), allocatable :: pfzcb                 !!pf cell boundaries in z
 
  

  PetscInt, dimension(:), allocatable :: d_nnz             !!petsc prealloc vec (diag blocks)
  PetscReal, dimension(:), allocatable :: delA
  Mat :: A,Ai
  PetscErrorCode :: perr
  MatType :: tp
  PetscInt :: prn(1),pcn(1)
  PetscInt :: d_nz,o_nz
  Vec :: psol
  Vec :: X
  Vec :: B
  KSP :: KS,KSi
  PC :: P
  logical :: nzero_flag=.true.
  PetscScalar, pointer :: vloc(:)
  Vec :: pflotran_tracer_vec_mpi
  Vec :: pflotran_tracer_vec_seq
  Vec :: pflotran_saturation_vec_mpi
  Vec :: pflotran_saturation_vec_seq
  Vec :: pflotran_temperature_vec_mpi
  Vec :: pflotran_temperature_vec_seq
  VecScatter :: pflotran_scatter
  PetscInt :: pflotran_vec_size
  character(len=32) :: pflotran_group_prefix

contains
!_________________________________________________________________
subroutine elog(com,i1,i2)
  implicit none
  integer :: com,i1,i2
  logical :: exst
  integer :: d1,d2,d3
  real :: v1,v2,v3

  select case (com)
     
  case(0)
     inquire(file='e4d.inp',exist=exst)
     if (.not. exst) then
        open(13,file=trim(log_file),status='old',action='write', &
             position='append')
        write(13,*) 'Cannot find the input file e4d.inp'
        write(*,*) 'Aborting E4D'
        i2=-1
        close(13)
        return
     else
        open(13,file=trim(log_file),status='old',action='write',position='append')
        write(13,*) 'INITIALIZING E4D: FOUND e4d.inp'
        close(13)
        i2=0
     end if

     case(1)
        open(13,file=trim(log_file),status='old',action='write', &
             position='append')
        if (i1 .ne. 0) then
           write(13,*) "There was a problem reading the mesh file name in e4d.inp"
           write(13,*) "Aborting E4D"
           i2=-1
        else
           i2=0
           write(13,*) "The specified mesh file is: ",trim(mshfile)
           inquire(file=trim(mshfile),exist=exst)
           if (.not. exst) then
              write(13,*) "Cannot find the mesh file: ",trim(mshfile)
              write(*,*) "Aborting E4D"
              i2=-1
           end if
        end if
        close(13)
        return

     case(2)
        open(13,file=trim(log_file),status='old',action='write', &
             position='append')
        if (i1 .ne. 0) then
           write(13,*) "There was a problem reading the survey file name in e4d.inp"
           write(13,*) "Aborting E4D"
           i2=-1
        else
           i2=0
           write(13,*) "The specified survey file is: ",trim(efile)
           inquire(file=trim(efile),exist=exst)
           if (.not. exst) then
              write(13,*) "Cannot find the survey file: ",trim(efile)
              write(*,*) "Aborting E4D"
              i2=-1
           end if
        end if
        close(13)
        return

     case(3)
        open(13,file=trim(log_file),status='old',action='write', &
             position='append')
        if (i1 .ne. 0) then
           write(13,*) "There was a problem reading the conductivity list file name in e4d.inp"
           write(13,*) "Aborting E4D"
           i2=-1
        else
           i2=0
           write(13,*) "The specified conductivity list file is: ",trim(sigfile)
           inquire(file=trim(sigfile),exist=exst)
           if (.not. exst) then
              write(13,*) "Cannot find the survey file: ",trim(sigfile)
              write(*,*) "Aborting E4D"
              i2=-1
           end if
        end if
        close(13)
        return

    case(4)
        open(13,file=trim(log_file),status='old',action='write', &
             position='append')
        if (i1 .ne. 0) then
           write(13,*) "There was a problem reading the FMN file name in e4d.inp"
           write(13,*) "Aborting E4D"
           i2=-1
        else
           i2=0
           write(13,*) "The specified FMN (Archies parameter) file is: ",trim(fmn_file)
           inquire(file=trim(fmn_file),exist=exst)
           if (.not. exst) then
              write(13,*) "Cannot find the FMN file: ",trim(fmn_file)
              write(*,*) "Aborting E4D"
              i2=-1
           end if
        end if
        close(13)
        return

     case(5)
        open(13,file=trim(log_file),status='old',action='write', &
             position='append')
        if (i1 .ne. 0) then
           write(13,*) "There was a problem reading the number of electrodes in: ",trim(efile)
           write(*,*) "Aborting E4D"
        else
           write(13,*) "Number of electrodes: ",ne
        end if
        close(13)
        return

     case(6)
        open(13,file=trim(log_file),status='old',action='write',position='append')
        write(13,*) "There was a problem reading the parameters for electrode: ",i1
        write(*,*) "Aborting E4D"
        close(13)
        return
        
     case(7)
        open(13,file=trim(log_file),status='old',action='write',position='append')
        write(13,*) "E4D internal mesh translation file: ",mshfile(1:i1)//".trn"
        inquire(file=mshfile(1:i1)//".trn",exist=exst)
        if (.not. exst) then
           write(13,*) "Cannot find the mesh translation file: ",mshfile(1:i1)//".trn"
           write(*,*) "Aborting E4D"
           i2=-1
        else
           i2=0
        end if
        close(13)
        return

     case(8)
        open(13,file=trim(log_file),status='old',action='write',position='append')
        write(13,*) "There was a problem reading the internal mesh translation values"
        write(13,*) "Aborting E4D"
        close(13)
        return

     case(9)
        open(13,file=trim(log_file),status='old',action='write', &
             position='append')
        if (i1 .ne. 0) then
           write(13,*) "There was a problem reading the number of measurements"
           write(13,*) "in the survey file: ",trim(efile)
           write(13,*) "Aborting E4D"
           i2=-1
        else
           i2=0
           write(13,*) "The number of measurements per survey is: ",nm
        end if
        close(13)
        return
 
     case(10)
        open(13,file=trim(log_file),status='old',action='write',position='append')
        write(13,*) 'There wase a problem reading measurement number :',i1
        write(13,*) 'Aborting E4D'
        close(13)
        return

     case(11)
        open(13,file=trim(log_file),status='old',action='write', &
             position='append')
        if (i1 .ne. 0) then
           write(13,*) 'There was a problem reading the first line of the '
           write(13,*) 'conductivity file :',trim(sigfile)
           write(13,*) 'The first line of the conductivity file must contain the '
           write(13,*) 'following parameters: '
           write(13,*) 'Number_of_values Formation_Factor Cond_surface_water Cond_groundwater'
           write(13,*) 'Aborting E4D'
 
        else
           write(13,*) 'Number of conductivity values: ',i2    
        end if
        close(13)
        return

     case(12)
        open(13,file=trim(log_file),status='old',action='write',position='append')
        write(13,*) 'There was a problem reading the conductivity for element: ',i1
        write(13,*) 'Aborting E4D'
        close(13)
        return

     case(13)
        open(13,file=trim(log_file),status='old',action='write', &
             position='append')
        if (i1 .ne. 0) then
           write(13,*) 'There was a problem reading the number of mapping values'
           write(13,*) 'in the mapping file: ',trim(mapfile)
           write(13,*) 'Aborting E4D'
           i2=-1
        else
           write(13,*) "Number of mapping values: ",nmap
           i2=0
        end if
        close(13)
        return

     case(14)
        open(13,file=trim(log_file),status='old',action='write',position='append')
        write(13,*) "There was a problem reading mapping value: ",i1
        write(13,*) "Aborting E4D"
        close(13)
        return

     case(15)
        open(13,file=trim(log_file),status='old',action='write',position='append')
        inquire(file=mshfile(1:i1)//'.node',exist=exst)
        if (.not. exst) then
           write(13,*) 'Could not find the mesh node file: ',mshfile(1:i1)//'.node'
           write(13,*) 'Aborting E4D'
           i2=-1
        else
           write(13,*) "E4D mesh node file: ",mshfile(1:i1)//'.node'
           i2=0
        end if
        close(13)
        return

     case(16)
        open(13,file=trim(log_file),status='old',action='write',position='append')
        write(13,*) 'There was a problem reading the first line of the node file.'
        write(13,*) 'Aborting E4D'
        close(13)
        return

     case(17)
        open(13,file=trim(log_file),status='old',action='write',position='append')
        write(13,*) 'There was a problem reading node number: ',i1
        write(13,*) 'Aborting E4D'
        close(13)
        return
        
    case(18)
        open(13,file=trim(log_file),status='old',action='write',position='append')
        inquire(file=mshfile(1:i1)//'.ele',exist=exst)
        if (.not. exst) then
           write(13,*) 'Could not find the mesh node file: ',mshfile(1:i1)//'.ele'
           write(13,*) 'Aborting E4D'
           i2=-1
        else
           write(13,*) "E4D mesh node file: ",mshfile(1:i1)//'.ele'
           i2=0
        end if
        close(13)
        return
        
     case(19)
        open(13,file=trim(log_file),status='old',action='write',position='append')
        write(13,*) 'There was a problem reading the first line of the element file.'
        write(13,*) 'Aborting E4D'
        close(13)
        return

     case(20)
        open(13,file=trim(log_file),status='old',action='write',position='append')
        write(13,*) 'There was a problem reading element number: ',i1
        write(13,*) 'Aborting E4D'
        close(13)
        return

     case(21)
        open(13,file=trim(log_file),status='old',action='write',position='append')
        write(13,*) 
        write(13,*) 'Done setting up forward run ......'
        write(13,*) '    Number of nodes: ',nnodes
        write(13,*) '    Number of elements: ',nelem
        !write(13,*) '    Minimum initial conductivity: ',minval(base_sigma)
        !write(13,*)  '    Maximum intitial conductivity ',maxval(base_sigma)
        close(13)

     case(22)
        open(13,file=trim(log_file),status='old',action='write',position='append')
        write(13,*)
        write(13,*) 'mcomm = ',i1
        write(13,*) 'waiting for pflotran solution'
        close(13)

     case(23)
        open(13,file=trim(log_file),status='old',action='write', &
             position='append')
        if (i1 .ne. 0) then
           write(13,*) "There was a problem reading the list file name in e4d.inp"
           write(13,*) "Aborting E4D"
           i2=-1
        else
           i2=0
           write(13,*) "The specified conductivity list file is: ",trim(list_file)
           inquire(file=trim(list_file),exist=exst)
           if (.not. exst) then
              write(13,*) "Cannot find the list file: ",trim(list_file)
              write(*,*) "Aborting E4D"
              i2=-1
           end if
        end if
        close(13)
        return

     case(24)
        open(13,file=trim(log_file),status='old',action='write', &
             position='append')
        if (i1 .ne. 0) then
           write(13,*) 'There was a problem reading the first line of the list file.'
           write(13,*) 'The first line of the list file must contain the '
           write(13,*) 'Number_of_E4D_times. '
           write(13,*) 'Aborting E4D'
           i2=-1
        else
           write(13,*) "Number of E4D Times: ",ntime
           i2=0
        end if
        close(13)
        return
        
     case(25)
        open(13,file=trim(log_file),status='old',action='write', &
             position='append')
        if (i1 .ne. 0) then
           write(13,*) "There was a problem reading list file line: ",i2
           write(13,*) "Each list file line must contain the following: "
           write(13,*) "E4D_time Survey_file_name Conductivity_file_name"
           write(13,*) "Aborting E4D"
        else
           write(13,"(I5,F12.0,A42,A42)") i2,e4d_time,csrv_file,ccond_file  
        end if
        close(13)
        return
        
     case(26)
        open(13,file=trim(log_file),status='old',action='write',position='append')
        inquire(file=trim(csrv_file),exist=exst)
        if (.not. exst) then
           write(13,*) 'Cannot open the file: ',csrv_file
           write(13,*) 'which is listed on line: ',i2
           write(13,*) 'of the list file: ',list_file
           write(13,*) 'Aborting E4D'
           i2=-1
           close(13)
           return
        end if
        inquire(file=trim(ccond_file),exist=exst)
        if (.not. exst) then
           write(13,*) 'Cannot open the file: ',ccond_file
           write(13,*) 'which is listed on line: ',i2
           write(13,*) 'of the list file: ',list_file
           write(13,*) 'Aborting E4D'
           i2=-1
           close(13)
           return
        end if
        i2=0
        close(13)
        return
        
     case(27)
        open(13,file=trim(log_file),status='old',action='write', &
             position='append')
        if (i1 .ne. 0) then
           write(13,*) 'There was a problem reading the number of electrodes in: ',trim(csrv_file)
           write(13,*) 'Aborting E4D'
           i1=-1
        elseif (i2 .ne. ne) then
           write(13,*) 'The number of electrodes in file: ',trim(csrv_file)
           write(13,*) 'is: ',i2
           write(13,*) 'The number of electrodes in the baseline survey file is: ',ne
           write(13,*) 'Each survey geometry must be equivalent.'
           write(13,*) 'Aborting E4D'
           i1=-1
        end if
        close(13)
        return

     case(28)
        open(13,file=trim(log_file),status='old',action='write',position='append')
        write(13,*) 'There was a problem reading electrode number: ',i1
        write(13,*) 'in file: ',trim(csrv_file)
        write(13,*) 'Make sure the electrodes are specified exactly as in the baseline survey file.'
        write(13,*) 'Aborting E4D'
        close(13)
        return

     case(29)
        open(13,file=trim(log_file),status='old',action='write', &
             position='append')
        if (i1 .ne. 0) then
           write(13,*) 'There was a problem reading the number of electrodes in: ',trim(csrv_file)
           write(13,*) 'Aborting E4D'
           i1=-1
        elseif (i2 .ne. nm) then
           write(13,*) 'The number of measurements in file: ',trim(csrv_file)
           write(13,*) 'is: ',i2
           write(13,*) 'The number of measurements in the baseline survey file is: ',nm
           write(13,*) 'Each survey geometry must be equivalent.'
           write(13,*) 'Aborting E4D'
           i1=-1
        end if
        close(13)
        return

     case(30)
        open(13,file=trim(log_file),status='old',action='write',position='append')
        write(13,*) 'There was a problem reading measurement number: ',i1
        write(13,*) 'in file: ',trim(csrv_file)
        write(13,*) 'Make sure the a,b,m,n is specified exactly as in the baseline survey file.'
        write(13,*) 'Aborting E4D'
        close(13)
        return

     case(31)
        open(13,file=trim(log_file),status='old',action='write',position='append')
        write(13,*) 'There was a problem reading the number of conductivity values in: ',trim(ccond_file)
        write(13,*) 'Aborting E4D'
        close(13)
        return

        
     case(32)
        open(13,file=trim(log_file),status='old',action='write',position='append')
        write(13,*) 'The number of conductivity values specified in: ',trim(ccond_file)
        write(13,*) 'is: ',i1
        write(13,*) 'The number of conductivity values specifed in the elment file is'
        write(13,*) 'is: ',i2
        write(13,*) 'Aborting E4D'
        close(13)
        return

     case(33)
        open(13,file=trim(log_file),status='old',action='write',position='append')
        write(13,*) 'There was a problem reading conductivity value: ',i1
        write(13,*) 'in file: ',trim(ccond_file)
        write(13,*) 'Aborting E4D'
        close(13)
        return

     case(34)
        open(13,file=trim(log_file),status='old',action='write',position='append')
        write(13,*) 'The number of elements in the element file is: ',i2
        write(13,*) 'The number of values in the conductivity files is: ',i1
        write(13,*) 'Aborting E4D'
        close(13)
        return

     case(35)
        open(13,file=trim(log_file),status='old',action='write',position='append')
        write(13,*) "Executing E4D Simulation for Time: ",e4d_time
        close(13)
        
     case(36)
        open(13,file=trim(log_file),status='old',action='write',position='append')
        write(13,*) "Received Time: ",pf_time, " from PFLOTRAN"
        close(13)
     

     case(37)
        inquire(file='pf_mesh.txt',exist=exst)
        if(.not.exst) then
           open(13,file=trim(log_file),status='old',action='write',position='append')
           write(*,*) "E4D couldn't find the pflotran mesh description file pf_mesh.txt"
           write(*,*) "E4D is aborting ..."
           write(13,*) "E4D couldn't find the pflotran mesh description file pf_mesh.txt"
           write(13,*) "Aborting ..."
           i1=-1
           close(13)
           return
        else
           i1=0
           return
        end if
           

     case(38)
        open(13,file=trim(log_file),status='old',action='write', &
             position='append')
        if (i1 .ne. 0) then
           write(13,*) 'There was a problem reading the first line of the '
           write(13,*) 'FMN file :',trim(fmn_file)
           write(13,*) 'The first line of the FMN file must contain the '
           write(13,*) 'number of FMN values, which must be the number of mesh elements.'
           write(13,*) 'Aborting E4D'
 
        else
           write(13,*) 'Number of FMN values: ',i2
        end if
        !if(i2 .ne. nelem) then
        !   write(13,*) 'The number of FMN values must equal the number of elements'
        !   write(13,*) 'Aborting ...'
        !   i1=-1
        !end if
        close(13)
        return

     case(39)
        open(13,file=trim(log_file),status='old',action='write', &
             position='append')
        write(13,*) 'There was a problem reading row number: ',i2
        write(13,*) 'in the FMN file: ',trim(fmn_file)
        write(13,*) 'Aborting E4D'
        close(13)
        return

     case(40)
        open(13,file=trim(log_file),status='old',action='write', &
             position='append')
        if (i1 .ne. 0) then
           write(13,*) "There was a problem reading the mode in e4d.inp"
           write(13,*) "Aborting E4D"
           i2=-1
        else
           if(mode==33) then
              write(13,*) "Running in tank simulation mode"
           end if
           i2 = 0
        end if
        close(13)
        return

     end select


        
end subroutine elog
!_________________________________________________________________	
 
end module e4d_vars
