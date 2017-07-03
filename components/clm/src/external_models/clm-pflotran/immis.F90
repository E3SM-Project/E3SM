module Immis_module
  
  use Immis_Aux_module
  use Global_Aux_module

  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "petsc/finclude/petscsys.h"
  
!#include "include/petscf90.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
  ! It is VERY IMPORTANT to make sure that the above .h90 file gets included.
  ! Otherwise some very strange things will happen and PETSc will give no
  ! indication of what the problem is.
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscmat.h90"
#include "petsc/finclude/petscdm.h"
#include "petsc/finclude/petscdm.h90"
!#ifdef USE_PETSC216
!#include "petsc/finclude/petscsles.h"
!#endif
#include "petsc/finclude/petscsnes.h"
#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petscsysdef.h"
#include "petsc/finclude/petscis.h"
#include "petsc/finclude/petscis.h90"
#include "petsc/finclude/petsclog.h"
#include "petsc/finclude/petscerror.h"

! Cutoff parameters
  PetscReal, parameter :: formeps   = 1.D-4
  PetscReal, parameter :: eps = 1.D-8
  PetscReal, parameter :: dfac = 1D-8
  PetscReal, parameter :: floweps   = 1.D-24
!  PetscReal, parameter :: satcuteps = 1.D-5
  PetscReal, parameter :: zerocut =0.D0  !1D-8
  

  PetscInt, parameter :: jh2o=1, jco2=2

  public ImmisResidual,ImmisJacobian, &
         ImmisUpdateFixedAccumulation,ImmisTimeCut, &
         ImmisSetup,ImmisUpdateReason, &
         ImmisMaxChange,ImmisUpdateSolution, &
         ImmisGetTecplotHeader,ImmisInitializeTimestep, &
         ImmisUpdateAuxVars, &
         ImmisComputeMassBalance, &
         ImmisDestroy

contains

! ************************************************************************** !

subroutine ImmisTimeCut(realization)
  ! 
  ! Resets arrays for time step cut
  ! 
  ! Author: Chuan Lu
  ! Date: 9/13/08
  ! 
 
  use Realization_Subsurface_class
  use Option_module
  use Field_module
 
  implicit none
  
  type(realization_subsurface_type) :: realization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  
  PetscReal, pointer :: xx_p(:),yy_p(:)
  PetscErrorCode :: ierr
  PetscInt :: local_id

  option => realization%option
  field => realization%field

end subroutine ImmisTimeCut

! ************************************************************************** !

subroutine ImmisSetup(realization)
  ! 
  ! Author: Chuan Lu
  ! Date: 9/13/08
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Output_Aux_module
  use co2_span_wagner_module
  use co2_sw_module
  use co2_span_wagner_spline_module
   
  type(realization_subsurface_type) :: realization
  
  type(patch_type), pointer :: cur_patch
  type(output_variable_list_type), pointer :: list

  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call ImmisSetupPatch(realization)
    cur_patch => cur_patch%next
  enddo

  list => realization%output_option%output_snap_variable_list
  call ImmisSetPlotVariables(list)
  list => realization%output_option%output_obs_variable_list
  call ImmisSetPlotVariables(list)

end subroutine ImmisSetup

! ************************************************************************** !

subroutine ImmisSetupPatch(realization)
  ! 
  ! Creates arrays for auxiliary variables
  ! 
  ! Author: Chuan Lu
  ! Date: 10/1/08
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Coupler_module
  use Connection_module
  use Grid_module
 
  implicit none
  
  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type),pointer :: patch
  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink

  PetscInt :: ghosted_id, iconn, sum_connection, ipara
  type(Immis_auxvar_type), pointer :: auxvars(:)
  type(Immis_auxvar_type), pointer :: auxvars_bc(:)  
  type(Immis_auxvar_type), pointer :: auxvars_ss(:)  
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid

  patch%aux%Immis => ImmisAuxCreate()
  
!  option%io_buffer = 'Before Immis can be run, the thc_parameter object ' // &
!                     'must be initialized with the proper variables ' // &
!                     'ImmisAuxCreate() is called anyhwere.'
!  call printErrMsg(option)
! print *,' ims setup get Aux', option%nphase, size(patch%saturation_function_array)     
! immis_parameters create *********************************************
! Sir
  allocate(patch%aux%Immis%Immis_parameter%sir(option%nphase, &
                                  size(patch%saturation_function_array)))
                                
  do ipara = 1, size(patch%saturation_function_array)
    patch%aux%Immis%immis_parameter%sir(:,patch% &
        saturation_function_array(ipara)%ptr%id) = &
      patch%saturation_function_array(ipara)%ptr%Sr(:)
  enddo

! dencpr  
  allocate(patch%aux%Immis%Immis_parameter%dencpr(size(patch%material_property_array)))
  do ipara = 1, size(patch%material_property_array)
    patch%aux%Immis%Immis_parameter%dencpr(iabs(patch% &
        material_property_array(ipara)%ptr%internal_id)) = &
      patch%material_property_array(ipara)%ptr%rock_density*option%scale*&
      patch%material_property_array(ipara)%ptr%specific_heat
  enddo
! ckwet
  allocate(patch%aux%Immis%Immis_parameter%ckwet(size(patch%material_property_array)))
  do ipara = 1, size(patch%material_property_array)
    patch%aux%Immis%Immis_parameter%ckwet(iabs(patch% &
        material_property_array(ipara)%ptr%internal_id)) = &
      patch%material_property_array(ipara)%ptr%thermal_conductivity_wet*option%scale
  enddo
! immis_parameters create_end *****************************************

! allocate auxvar data structures for all grid cells  
  allocate(auxvars(grid%ngmax))
  ! print *,' ims setup get Aux alloc', grid%ngmax
  do ghosted_id = 1, grid%ngmax
    call ImmisAuxVarInit(auxvars(ghosted_id),option)
  enddo
  patch%aux%Immis%auxvars => auxvars
  patch%aux%Immis%num_aux = grid%ngmax
  ! print *,' ims setup get Aux init'

  allocate(patch%aux%Immis%delx(option%nflowdof, grid%ngmax))
  allocate(patch%aux%Immis%res_old_AR(grid%nlmax,option%nflowdof))
  allocate(patch%aux%Immis%res_old_FL(ConnectionGetNumberInList(patch%grid%&
           internal_connection_set_list),option%nflowdof))
  ! print *,' ims setup allocate app array'
   ! count the number of boundary connections and allocate
  ! auxvar data structures for them  
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    sum_connection = sum_connection + &
                     boundary_condition%connection_set%num_connections
    boundary_condition => boundary_condition%next
  enddo
  allocate(auxvars_bc(sum_connection))
  ! print *,' ims setup get AuxBc alloc', sum_connection
  do iconn = 1, sum_connection
    call ImmisAuxVarInit(auxvars_bc(iconn),option)
  enddo
  patch%aux%Immis%auxvars_bc => auxvars_bc
  patch%aux%Immis%num_aux_bc = sum_connection
  
 ! Allocate source /sink  
  source_sink => patch%source_sink_list%first
  sum_connection = 0    
  do 
    if (.not.associated(source_sink)) exit
    sum_connection = sum_connection + &
                     source_sink%connection_set%num_connections
    source_sink => source_sink%next
  enddo
  allocate(auxvars_ss(sum_connection))
  do iconn = 1, sum_connection
    call ImmisAuxVarInit(auxvars_ss(iconn),option)
  enddo
  patch%aux%Immis%auxvars_ss => auxvars_ss
  patch%aux%Immis%num_aux_ss = sum_connection
  
  option%flow%numerical_derivatives = PETSC_TRUE

end subroutine ImmisSetupPatch

! ************************************************************************** !

subroutine ImmisComputeMassBalance(realization,mass_balance)
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  ! 

  use Realization_Subsurface_class
  use Patch_module

  type(realization_subsurface_type) :: realization
! PetscReal :: mass_balance(realization%option%nflowspec,realization%option%nphase)
  PetscReal :: mass_balance(realization%option%nflowspec,1)
  
  type(patch_type), pointer :: cur_patch
  
  mass_balance = 0.d0
  
  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call ImmisComputeMassBalancePatch(realization,mass_balance)
    cur_patch => cur_patch%next
  enddo

end subroutine ImmisComputeMassBalance

! ************************************************************************** !

subroutine ImmisComputeMassBalancePatch(realization,mass_balance)
  ! 
  ! Initializes mass balance
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/19/08
  ! 
 
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Field_module
  use Grid_module
  use Material_Aux_class
 
  implicit none
  
  type(realization_subsurface_type) :: realization
! PetscReal :: mass_balance(realization%option%nflowspec,realization%option%nphase)
  PetscReal :: mass_balance(realization%option%nflowspec,1)

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(immis_auxvar_type), pointer :: immis_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)

  PetscErrorCode :: ierr
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: iphase
  PetscInt :: ispec_start, ispec_end, ispec

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field
  material_auxvars => realization%patch%aux%Material%auxvars

  immis_auxvars => patch%aux%immis%auxvars

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    ! mass = volume * saturation * density
    do iphase = 1, option%nphase
      ispec = iphase
      mass_balance(ispec,1) = mass_balance(ispec,1) + &
          immis_auxvars(ghosted_id)%auxvar_elem(0)%den(iphase)* &
          immis_auxvars(ghosted_id)%auxvar_elem(0)%sat(iphase)* &
          material_auxvars(ghosted_id)%porosity* &
          material_auxvars(ghosted_id)%volume
    enddo
  enddo

end subroutine ImmisComputeMassBalancePatch

! ************************************************************************** !

subroutine ImmisZeroMassBalDeltaPatch(realization)
  ! 
  ! Zeros mass balance delta array
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/19/08
  ! 
 
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Grid_module
 
  implicit none
  
  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)

  PetscInt :: iconn

  option => realization%option
  patch => realization%patch

  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss

#ifdef COMPUTE_INTERNAL_MASS_FLUX
  do iconn = 1, patch%aux%Immis%num_aux
    patch%aux%Global%auxvars(iconn)%mass_balance_delta = 0.d0
  enddo
#endif

  ! Intel 10.1 on Chinook reports a SEGV if this conditional is not
  ! placed around the internal do loop - geh
  if (patch%aux%Immis%num_aux_bc > 0) then
    do iconn = 1, patch%aux%Immis%num_aux_bc
      global_auxvars_bc(iconn)%mass_balance_delta = 0.d0
    enddo
  endif
  
  if (patch%aux%Immis%num_aux_ss > 0) then
    do iconn = 1, patch%aux%Immis%num_aux_ss
      global_auxvars_ss(iconn)%mass_balance_delta = 0.d0
    enddo
  endif

end subroutine ImmisZeroMassBalDeltaPatch

! ************************************************************************** !

  function  ImmisInitGuessCheck(realization)
  ! 
  ! Immisinitguesscheckpatch:
  ! 
  ! Author: Chuan Lu
  ! Date: 12/10/07
  ! 
 
  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  
  PetscInt ::  ImmisInitGuessCheck
  type(realization_subsurface_type) :: realization
  type(option_type), pointer :: option
  type(patch_type), pointer :: cur_patch
  PetscInt :: ipass, ipass0
  PetscErrorCode :: ierr    

  option => realization%option
  ipass = 1
  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    ipass= ImmisInitGuessCheckPatch(realization)
    if (ipass<=0)then
      exit 
    endif
    cur_patch => cur_patch%next
  enddo

  call MPI_Barrier(option%mycomm,ierr)
  if (option%mycommsize >1)then
      call MPI_Allreduce(ipass,ipass0,ONE_INTEGER_MPI,MPIU_INTEGER, &
                         MPI_SUM,option%mycomm,ierr)
      if (ipass0 < option%mycommsize) ipass=-1
   endif
   ImmisInitGuessCheck =ipass

end function ImmisInitGuessCheck

! ************************************************************************** !

subroutine ImmisUpdateReasonPatch(reason,realization)
  ! 
  ! Immisinitguesscheckpatch:
  ! 
  ! Author: Chuan Lu
  ! Date: 10/10/08
  ! 
   use Realization_Subsurface_class
   use Patch_module
   use Field_module
   use Option_module
   use Grid_module

  implicit none
 
  PetscInt, intent(out) :: reason
  type(realization_subsurface_type) :: realization  
  type(patch_type),pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(option_type), pointer :: option 
  PetscReal, pointer :: xx_p(:), yy_p(:) 
  PetscInt :: n,n0,re
  PetscInt :: re0
  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field  
  patch => realization%patch
  grid => patch%grid

  re=1
 
  if (re>0)then
     call VecGetArrayReadF90(field%flow_xx, xx_p, ierr);CHKERRQ(ierr)
     call VecGetArrayF90(field%flow_yy, yy_p, ierr);CHKERRQ(ierr)

     do n = 1,grid%nlmax
!**** clu-Ignore inactive cells with inactive materials **************
        if (associated(patch%imat)) then
           if (patch%imat(grid%nL2G(n)) <= 0) cycle
        endif
        n0=(n-1)* option%nflowdof
  
! ******** Too huge change in pressure ****************     
!geh: I don't believe that this code is being used.  Therefore, I will add an
!     error message and let someone sort the use of option%dpmxe later
        option%io_buffer = 'option%dpmxe and option%dtmpmxe needs to be ' // &
          'refactored in ImmisUpdateReasonPatch'
        call printErrMsg(option)
!geh        if (dabs(xx_p(n0 + 1)- yy_p(n0 + 1))> (10.0D0 * option%dpmxe))then
           re=0; print *,'huge change in p', xx_p(n0 + 1), yy_p(n0 + 1)
           exit
!geh        endif

! ******** Too huge change in temperature ****************
!geh        if (dabs(xx_p(n0 + 2)- yy_p(n0 + 2))> (10.0D0 * option%dtmpmxe))then
           re=0; print *,'huge change in T', xx_p(n0 + 2), yy_p(n0 + 2)
           exit
!geh        endif
 
! ******* Check 0<=sat/con<=1 **************************
           if (xx_p(n0 + 3) > 1.D0)then
              re=0; exit
           endif
           if (xx_p(n0 + 3) < 0.)then
              re=0; exit
           endif
     end do
  
    if (re<=0) print *,'Sat out of Region at: ',n,xx_p(n0+1:n0+3)
    call VecRestoreArrayReadF90(field%flow_xx, xx_p, ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(field%flow_yy, yy_p, ierr);CHKERRQ(ierr)

   endif
  
end subroutine ImmisUpdateReasonPatch

! ************************************************************************** !

subroutine ImmisUpdateReason(reason, realization)
  ! 
  ! ImmisUpdateAuxVars: Updates the auxiliary variables associated with
  ! the Richards problem
  ! 
  ! Author: Chuan Lu
  ! Date: 10/10/08
  ! 

  use Realization_Subsurface_class
  use Patch_module
  implicit none

  type(realization_subsurface_type) :: realization
  
  type(patch_type), pointer :: cur_patch
  PetscInt :: reason

  PetscInt :: re, re0
  PetscErrorCode :: ierr

  re = 1
  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call ImmisUpdateReasonPatch(re, realization)
    if (re<=0)then
      exit 
    endif
    cur_patch => cur_patch%next
  enddo

  call MPI_Barrier(realization%option%mycomm,ierr)
  
  if (realization%option%mycommsize >1)then
     call MPI_Allreduce(re,re0,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM, &
          realization%option%mycomm,ierr)
     if (re0<realization%option%mycommsize) re=0
  endif
  reason=re
  
  if (reason<=0 .and. realization%option%myrank ==0) print *,'Sat or Con out of Region', re
end subroutine ImmisUpdateReason

! ************************************************************************** !

  function ImmisInitGuessCheckPatch(realization)
  ! 
  ! Author: Chuan Lu
  ! Date: 10/10/08
  ! 
   
    use co2_span_wagner_module
     
    use Realization_Subsurface_class
    use Patch_module
    use Field_module
    use Grid_module
    use Option_module
    implicit none
    
    PetscInt :: ImmisInitGuessCheckPatch 
    type(realization_subsurface_type) :: realization
    type(grid_type), pointer :: grid
    type(patch_type), pointer :: patch
    type(option_type), pointer :: option
    type(field_type), pointer :: field
      
    PetscInt :: local_id, ghosted_id, ipass
    PetscReal, pointer :: xx_p(:)
    PetscErrorCode :: ierr


    patch => realization%patch
    grid => patch%grid
    option => realization%option
    field => realization%field
    
    call VecGetArrayReadF90(field%flow_xx,xx_p, ierr);CHKERRQ(ierr)
    
    ipass=1
    do local_id = 1, grid%nlmax
       ghosted_id = grid%nL2G(local_id)
       !geh - Ignore inactive cells with inactive materials
       if (associated(patch%imat)) then
          if (patch%imat(ghosted_id) <= 0) cycle
       endif
      
!   insure zero liquid sat not passed to ptran (no effect on pflow)
       if (xx_p((local_id-1)*option%nflowdof+3) < 0.D0)xx_p((local_id-1)*option%nflowdof+3) = zerocut
       if (xx_p((local_id-1)*option%nflowdof+3) > 1.D0)xx_p((local_id-1)*option%nflowdof+3) = 1.D0 - zerocut
    
!   check if p,T within range of table  
       if (xx_p((local_id-1)*option%nflowdof+1)< p0_tab*1D6 &
            .or. xx_p((local_id-1)*option%nflowdof+1)>(ntab_p*dp_tab + p0_tab)*1D6)then
          ipass=-1; exit  
       endif
       if (xx_p((local_id-1)*option%nflowdof+2)< t0_tab -273.15D0 &
            .or. xx_p((local_id-1)*option%nflowdof+2)>ntab_t*dt_tab + t0_tab-273.15D0)then
          ipass=-1; exit
       endif
    enddo

    call VecRestoreArrayReadF90(field%flow_xx,xx_p, ierr);CHKERRQ(ierr)
    ImmisInitGuessCheckPatch = ipass
  end function ImmisInitGuessCheckPatch

! ************************************************************************** !

subroutine ImmisUpdateAuxVars(realization)
  ! 
  ! Updates the auxiliary variables associated with
  ! the Immis problem
  ! 
  ! Author: Chuan Lu
  ! Date: 10/10/08
  ! 

  use Realization_Subsurface_class
  use Patch_module

  type(realization_subsurface_type) :: realization
  
  type(patch_type), pointer :: cur_patch
  
  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call ImmisUpdateAuxVarsPatch(realization)
    cur_patch => cur_patch%next
  enddo

end subroutine ImmisUpdateAuxVars

! ************************************************************************** !

subroutine ImmisUpdateAuxVarsPatch(realization)
  ! 
  ! Updates the auxiliary variables associated with
  ! the Immis problem
  ! 
  ! Author: Chuan Lu
  ! Date: 12/10/07
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Field_module
  use Option_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Material_module
  
  implicit none

  type(realization_subsurface_type) :: realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  type(Immis_auxvar_type), pointer :: auxvars(:)
  type(Immis_auxvar_type), pointer :: auxvars_bc(:)
  type(Immis_auxvar_type), pointer :: auxvars_ss(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)

  PetscInt :: ghosted_id, local_id, istart, iend, sum_connection, idof, iconn
  PetscReal, pointer :: xx_loc_p(:), icap_loc_p(:)
  PetscReal :: xxbc(realization%option%nflowdof)
  PetscErrorCode :: ierr
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field
  
  auxvars => patch%aux%Immis%auxvars
  auxvars_bc => patch%aux%Immis%auxvars_bc
  auxvars_ss => patch%aux%Immis%auxvars_ss

  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss

  
  call VecGetArrayF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%icap_loc,icap_loc_p,ierr);CHKERRQ(ierr)

  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    iend = ghosted_id*option%nflowdof
    istart = iend-option%nflowdof+1
    if (.not. associated(patch%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr))then
       print*, 'error!!! saturation function not allocated', ghosted_id,icap_loc_p(ghosted_id)
    endif
   
    call ImmisAuxVarCompute_NINC(xx_loc_p(istart:iend), &
                       auxvars(ghosted_id)%auxvar_elem(0), &
                       patch%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr, &
                       realization%fluid_properties,option)

 ! update global variables
    if (associated(global_auxvars))then
      global_auxvars(ghosted_id)%pres(:) = auxvars(ghosted_id)%auxvar_elem(0)%pres -&
               auxvars(ghosted_id)%auxvar_elem(0)%pc(:)
      global_auxvars(ghosted_id)%temp = auxvars(ghosted_id)%auxvar_elem(0)%temp
      global_auxvars(ghosted_id)%sat(:) = auxvars(ghosted_id)%auxvar_elem(0)%sat(:)
      global_auxvars(ghosted_id)%den(:) = auxvars(ghosted_id)%auxvar_elem(0)%den(:)
      global_auxvars(ghosted_id)%den_kg(:) = auxvars(ghosted_id)%auxvar_elem(0)%den(:) &
                                          * auxvars(ghosted_id)%auxvar_elem(0)%avgmw(:)
      
!     print *,'immis: update_global:den ',global_auxvars(ghosted_id)%den_kg(:), &
!                               auxvars(ghosted_id)%auxvar_elem(0)%avgmw(:)
    else
      print *,'Not associated global for IMS'
    endif
  enddo

  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
      do idof = 1, option%nflowdof
        select case(boundary_condition%flow_condition%itype(idof))
        case(DIRICHLET_BC)
!         xxbc(:) = boundary_condition%flow_aux_real_var(:,iconn)
          xxbc(idof) = boundary_condition%flow_aux_real_var(idof,iconn)
        case(HYDROSTATIC_BC)
!         xxbc(1) = boundary_condition%flow_aux_real_var(1,iconn)
!         xxbc(2:option%nflowdof) = &
!           xx_loc_p((ghosted_id-1)*option%nflowdof+2:ghosted_id*option%nflowdof)
          xxbc(MPH_PRESSURE_DOF) = boundary_condition%flow_aux_real_var(MPH_PRESSURE_DOF,iconn)
          if (idof >= MPH_PRESSURE_DOF) then
            xxbc(idof) = xx_loc_p((ghosted_id-1)*option%nflowdof+idof)
          endif
        case(NEUMANN_BC,ZERO_GRADIENT_BC)
!         xxbc(:) = xx_loc_p((ghosted_id-1)*option%nflowdof+1:ghosted_id*option%nflowdof)
          xxbc(idof) = xx_loc_p((ghosted_id-1)*option%nflowdof+idof)
        end select
      enddo
 
      call ImmisAuxVarCompute_NINC(xxbc,auxvars_bc(sum_connection)%auxvar_elem(0), &
              patch%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr, &
                         realization%fluid_properties, option)

      if (associated(global_auxvars_bc)) then
        global_auxvars_bc(sum_connection)%pres(:)= auxvars_bc(sum_connection)%auxvar_elem(0)%pres -&
                     auxvars(ghosted_id)%auxvar_elem(0)%pc(:)
        global_auxvars_bc(sum_connection)%temp=auxvars_bc(sum_connection)%auxvar_elem(0)%temp
        global_auxvars_bc(sum_connection)%sat(:)=auxvars_bc(sum_connection)%auxvar_elem(0)%sat(:)
        !    global_auxvars(ghosted_id)%sat_store = 
        global_auxvars_bc(sum_connection)%den(:)=auxvars_bc(sum_connection)%auxvar_elem(0)%den(:)
        global_auxvars_bc(sum_connection)%den_kg = auxvars_bc(sum_connection)%auxvar_elem(0)%den(:) &
              * auxvars_bc(sum_connection)%auxvar_elem(0)%avgmw(:)
  !     global_auxvars(ghosted_id)%den_kg_store
  !     global_auxvars(ghosted_id)%mass_balance
  !     global_auxvars(ghosted_id)%mass_balance_delta
      endif

    enddo
    boundary_condition => boundary_condition%next
  enddo


! source/sinks
  source_sink => patch%source_sink_list%first
  sum_connection = 0    
  do 
    if (.not.associated(source_sink)) exit
    cur_connection_set => source_sink%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle

      call ImmisAuxVarCopy(auxvars(ghosted_id)%auxvar_elem(0), &
                          auxvars_ss(sum_connection)%auxvar_elem(0),option)
      call GlobalAuxVarCopy(global_auxvars(ghosted_id), &
                          global_auxvars_ss(sum_connection),option)

    enddo
    source_sink => source_sink%next
  enddo


  call VecRestoreArrayF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%icap_loc,icap_loc_p,ierr);CHKERRQ(ierr)
  
  patch%aux%Immis%auxvars_up_to_date = PETSC_TRUE

end subroutine ImmisUpdateAuxVarsPatch

! ************************************************************************** !

subroutine ImmisInitializeTimestep(realization)
  ! 
  ! Update data in module prior to time step
  ! 
  ! Author: Chuan Lu
  ! Date: 10/12/08
  ! 

  use Realization_Subsurface_class
  
  implicit none
  
  type(realization_subsurface_type) :: realization

  call ImmisUpdateFixedAccumulation(realization)

end subroutine ImmisInitializeTimestep

! ************************************************************************** !

subroutine ImmisUpdateSolution(realization)
  ! 
  ! Updates data in module after a successful time step
  ! 
  ! Author: Chuan Lu
  ! Date: 10/13/08
  ! 

  use Realization_Subsurface_class
  use Field_module
  use Patch_module
  
  implicit none
  
  type(realization_subsurface_type) :: realization

  type(field_type), pointer :: field
  type(patch_type), pointer :: cur_patch
  PetscErrorCode :: ierr
  PetscViewer :: viewer
  
  field => realization%field
  
  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
!     call ImmisUpdateSolutionPatch(realization)
    cur_patch => cur_patch%next
  enddo

! make room for hysteric s-Pc-kr

end subroutine ImmisUpdateSolution

! ************************************************************************** !

subroutine ImmisUpdateSolutionPatch(realization)
  ! 
  ! Updates mass balance
  ! 
  ! Author: PCL
  ! Date: 11/18/11
  ! 

  use Realization_Subsurface_class
    
  implicit none
  
  type(realization_subsurface_type) :: realization

  if (realization%option%compute_mass_balance_new) then
    call ImmisUpdateMassBalancePatch(realization)
  endif

end subroutine ImmisUpdateSolutionPatch

! ************************************************************************** !

subroutine ImmisUpdateMassBalancePatch(realization)
  ! 
  ! Updates mass balance
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/19/08
  ! 
 
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Grid_module
 
  implicit none
  
  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)

  PetscInt :: iconn

  option => realization%option
  patch => realization%patch

  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss

  ! Intel 10.1 on Chinook reports a SEGV if this conditional is not
  ! placed around the internal do loop - geh
  if (patch%aux%immis%num_aux_bc > 0) then
    do iconn = 1, patch%aux%immis%num_aux_bc
      global_auxvars_bc(iconn)%mass_balance = &
        global_auxvars_bc(iconn)%mass_balance + &
        global_auxvars_bc(iconn)%mass_balance_delta*option%flow_dt
    enddo
  endif
  
  if (patch%aux%immis%num_aux_ss > 0) then
    do iconn = 1, patch%aux%immis%num_aux_ss
      global_auxvars_ss(iconn)%mass_balance = &
        global_auxvars_ss(iconn)%mass_balance + &
        global_auxvars_ss(iconn)%mass_balance_delta*option%flow_dt
    enddo
  endif

end subroutine ImmisUpdateMassBalancePatch

! ************************************************************************** !

subroutine ImmisUpdateFixedAccumulation(realization)
  ! 
  ! Updates the fixed portion of the
  ! accumulation term
  ! 
  ! Author: Chuan Lu
  ! Date: 10/12/08
  ! 

  use Realization_Subsurface_class
  use Patch_module

  type(realization_subsurface_type) :: realization
  
  type(patch_type), pointer :: cur_patch
  
  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call ImmisUpdateFixedAccumPatch(realization)
    cur_patch => cur_patch%next
  enddo

end subroutine ImmisUpdateFixedAccumulation

! ************************************************************************** !

subroutine ImmisUpdateFixedAccumPatch(realization)
  ! 
  ! Updates the fixed portion of the
  ! accumulation term
  ! 
  ! Author: Chuan Lu
  ! Date: 10/12/08
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module

  implicit none
  
  type(realization_subsurface_type) :: realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(Immis_parameter_type), pointer :: immis_parameter
  type(Immis_auxvar_type), pointer :: auxvars(:)

  PetscInt :: ghosted_id, local_id, istart, iend
  PetscReal, pointer :: xx_p(:), icap_loc_p(:)
  PetscReal, pointer :: porosity_loc_p(:), tortuosity_loc_p(:), volume_p(:), &
                        ithrm_loc_p(:), accum_p(:)
                          
  PetscErrorCode :: ierr
  
  call ImmisUpdateAuxVarsPatch(realization) 
  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid


  immis_parameter => patch%aux%Immis%immis_parameter
  auxvars => patch%aux%Immis%auxvars
    
  call VecGetArrayReadF90(field%flow_xx,xx_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%icap_loc,icap_loc_p,ierr);CHKERRQ(ierr)
!geh: refactor  call VecGetArrayF90(field%porosity_loc,porosity_loc_p,ierr)
!geh: refactor  call VecGetArrayF90(field%tortuosity_loc,tortuosity_loc_p,ierr)
!geh: refactor  call VecGetArrayF90(field%volume,volume_p,ierr)
  call VecGetArrayF90(field%ithrm_loc,ithrm_loc_p,ierr);CHKERRQ(ierr)

  call VecGetArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    iend = local_id*option%nflowdof
    istart = iend-option%nflowdof+1

    call ImmisAccumulation(auxvars(ghosted_id)%auxvar_elem(0), &
                              porosity_loc_p(ghosted_id), &
                              volume_p(local_id), &
                              immis_parameter%dencpr(int(ithrm_loc_p(ghosted_id))), &
                              option,ZERO_INTEGER, accum_p(istart:iend)) 
  enddo

  call VecRestoreArrayReadF90(field%flow_xx,xx_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%icap_loc,icap_loc_p,ierr);CHKERRQ(ierr)
!geh refactor  call VecRestoreArrayF90(field%porosity_loc,porosity_loc_p,ierr)
!geh refactor  call VecRestoreArrayF90(field%tortuosity_loc,tortuosity_loc_p,ierr)
!geh refactor  call VecRestoreArrayF90(field%volume,volume_p,ierr)
  call VecRestoreArrayF90(field%ithrm_loc,ithrm_loc_p,ierr);CHKERRQ(ierr)

  call VecRestoreArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)

#if 0
!  call ImmisNumericalJacobianTest(field%flow_xx,realization)
#endif

end subroutine ImmisUpdateFixedAccumPatch

! ************************************************************************** !

subroutine ImmisAccumulation(auxvar,por,vol,rock_dencpr,option,iireac,Res)
  ! 
  ! Computes the non-fixed portion of the accumulation
  ! term for the residual
  ! 
  ! Author: Chuan Lu
  ! Date: 10/12/08
  ! 

  use Option_module
  
  implicit none

  type(Immis_auxvar_elem_type) :: auxvar
  type(option_type) :: option
  PetscReal Res(1:option%nflowdof) 
  PetscReal vol,por,rock_dencpr
     
  PetscInt :: ispec, np, iireac
  PetscReal :: porXvol, mol(option%nflowspec), eng
  
! if (present(ireac)) iireac=ireac

  porXvol = por*vol
  mol=0.d0; eng=0.D0
  do np = 1, option%nphase
    mol(np) = mol(np) + auxvar%sat(np) * auxvar%den(np)
    eng = eng + auxvar%sat(np) * auxvar%den(np) * auxvar%u(np)
  enddo
  mol = mol * porXvol
 ! if (option%use_isothermal == PETSC_FALSE) &
  eng = eng * porXvol + (1.d0 - por) * vol * rock_dencpr * auxvar%temp 
 
! Reaction terms here
! Note if iireac >0, then it is the node global index
! if (option%run_coupled == PETSC_TRUE .and. iireac>0) then
!H2O
!   mol(1)= mol(1) - option%flow_dt * option%rtot(iireac,1)
!   mol(2)= mol(2) - option%flow_dt * option%rtot(iireac,2)
! endif
  
! if (option%use_isothermal) then
!   Res(1:option%nflowdof) = mol(:)
! else
    Res(1:option%nphase) = mol(:)
    Res(option%nflowdof) = eng
! endif
end subroutine ImmisAccumulation

! ************************************************************************** !

subroutine ImmisSourceSink(mmsrc,nsrcpara,psrc,tsrc,hsrc,auxvar,isrctype,Res, &
                           qsrc_phase,energy_flag,option)
  ! 
  ! Computes source/sink
  ! 
  ! Author: Chuan Lu
  ! Date: 10/12/08
  ! 

  use Option_module
  
  use EOS_Water_module
! use Gas_EOS_module  
  use co2eos_module
  use co2_span_wagner_spline_module, only: sw_prop
  use co2_sw_module, only: co2_sw_interp
  use co2_span_wagner_module 
  
  implicit none

  type(Immis_auxvar_elem_type) :: auxvar
  type(option_type) :: option
  PetscReal :: Res(1:option%nflowdof) 
  PetscReal, pointer :: mmsrc(:)
! PetscReal :: mmsrc(option%nflowspec), psrc(option%nphase),tsrc,hsrc 
  PetscReal :: psrc(option%nphase),tsrc,hsrc 
  PetscInt :: isrctype
  PetscInt :: nsrcpara
  PetscBool :: energy_flag
  PetscReal :: qsrc_phase(:) ! volumetric rate of injection/extraction for each phase
     
  PetscReal, allocatable :: msrc(:)
! PetscReal :: msrc(option%nflowspec),dw_kg, dw_mol,dddt,dddp
  PetscReal :: dw_kg, dw_mol,dddt,dddp
  PetscReal :: enth_src_h2o, enth_src_co2 
  PetscReal :: rho, fg, dfgdp, dfgdt, eng, dhdt, dhdp, visc, dvdt, dvdp, xphi
  PetscReal :: ukvr, v_darcy, dq, dphi
  PetscReal :: well_status, well_diameter
  PetscReal :: pressure_bh, well_factor, pressure_max, pressure_min
  PetscReal :: well_inj_water, well_inj_co2
  PetscInt :: np
  PetscInt :: iflag
  PetscErrorCode :: ierr
  
  Res = 0.D0
  allocate(msrc(nsrcpara))
  msrc = mmsrc(1:nsrcpara)

! if (present(ireac)) iireac=ireac
! if (energy_flag) then
!   Res(option%nflowdof) = Res(option%nflowdof) + hsrc * option%flow_dt   
! endif         

  qsrc_phase = 0.d0
 
  select case(isrctype)
    case(MASS_RATE_SS)
      msrc(1) =  msrc(1) / FMWH2O
      msrc(2) =  msrc(2) / FMWCO2
      if (msrc(1) /= 0.d0) then ! H2O injection
        call EOSWaterDensity(tsrc,auxvar%pres,dw_kg,dw_mol,ierr)
        call EOSWaterEnthalpy(tsrc,auxvar%pres,enth_src_h2o,ierr)
        ! J/kmol -> whatever units
        enth_src_h2o = enth_src_h2o * option%scale
!           units: dw_mol [mol/dm^3]; dw_kg [kg/m^3]
!           qqsrc = qsrc1/dw_mol ! [kmol/s (mol/dm^3 = kmol/m^3)]
!       Res(jh2o) = Res(jh2o) + msrc(1)*(1.d0-csrc)*option%flow_dt
!       Res(jco2) = Res(jco2) + msrc(1)*csrc*option%flow_dt
        Res(jh2o) = Res(jh2o) + msrc(1)*option%flow_dt
        Res(jco2) = Res(jco2) + msrc(1)*option%flow_dt
        if (energy_flag) Res(option%nflowdof) = Res(option%nflowdof) + &
          msrc(1)*enth_src_h2o*option%flow_dt
          
        ! store volumetric rate for ss_fluid_fluxes()
        qsrc_phase(1) = msrc(1)/dw_mol
      endif  
    
      if (msrc(2) > 0.d0) then ! CO2 injection
!       call printErrMsg(option,"concentration source not yet implemented in Immis")
        if (option%co2eos == EOS_SPAN_WAGNER) then
         !  span-wagner
          rho = auxvar%den(jco2)*FMWCO2  
          select case(option%itable)  
            case(0,1,2,4,5)
              if (option%itable >=4) then
              call co2_sw_interp(auxvar%pres*1.D-6,&
                tsrc,rho,dddt,dddp,fg,dfgdp,dfgdt, &
                eng,enth_src_co2,dhdt,dhdp,visc,dvdt,dvdp,option%itable)
              else
                iflag = 1
                call co2_span_wagner(auxvar%pres*1.D-6,&
                  tsrc+273.15D0,rho,dddt,dddp,fg,dfgdp,dfgdt, &
                  eng,enth_src_co2,dhdt,dhdp,visc,dvdt,dvdp,iflag,option%itable)
              endif 
            case(3) 
              call sw_prop(tsrc,auxvar%pres*1.D-6,rho, &
                     enth_src_co2, eng, fg)
          end select     

         !  units: rho [kg/m^3]; csrc1 [kmol/s]
          enth_src_co2 = enth_src_co2 * FMWCO2
          
          ! store volumetric rate for ss_fluid_fluxes()
          ! qsrc_phase [m^3/sec] = msrc [kmol/sec] / [kg/m^3] * [kg/kmol]  
          qsrc_phase(2) = msrc(2)/auxvar%den(jco2)
          
        else if (option%co2eos == EOS_MRK)then
! MRK eos [modified version from  Kerrick and Jacobs (1981) and Weir et al. (1996).]
          call CO2(tsrc,auxvar%pres, rho,fg, xphi,enth_src_co2)
          qsrc_phase(2) = msrc(2)/auxvar%den(jco2)
          enth_src_co2 = enth_src_co2*FMWCO2*option%scale
      else
         call printErrMsg(option,'pflow Immis ERROR: Need specify CO2 EOS')
      endif
              
      Res(jco2) = Res(jco2) + msrc(2)*option%flow_dt
      if (energy_flag) Res(option%nflowdof) = Res(option%nflowdof) + msrc(2) * &
           enth_src_co2 *option%flow_dt
       endif


    case(WELL_SS) ! production well
     !if node pessure is lower than the given extraction pressure, shut it down
    ! Flow term
!  well parameter explaination
!   1. well status. 1 injection; -1 production; 0 shut in
!                   2 rate controled injection (same as rate_ss, with max pressure control, not completed yet) 
!                  -2 rate controled production(not implemented for now) 
!
!   2. well factor [m^3],  the effective permeability [m^2/s]
!   3. bottomhole pressure:  [Pa]
!   4. max pressure: [Pa]
!   5. min pressure: [Pa]   
!   6. preferred mass flux of water [kg/s]
!   7. preferred mass flux of Co2 [kg/s]
!   8. well diameter, not used now
!   9. skin factor, not used now

      well_status = msrc(1)
      well_factor = msrc(2)
      pressure_bh = msrc(3)
      pressure_max = msrc(4)
      pressure_min = msrc(5)
      well_inj_water = msrc(6)
      well_inj_co2 = msrc(7)
    
!     if (pressure_min < 0D0) pressure_min = 0D0 !not limited by pressure lower bound   

    ! production well (well status = -1)
      if ( dabs(well_status + 1D0) < 1D-1) then 
        if (auxvar%pres > pressure_min) then
          Dq = well_factor 
          do np = 1, option%nphase
            dphi = auxvar%pres - auxvar%pc(np) - pressure_bh
            if (dphi>=0.D0) then ! outflow only
              ukvr = auxvar%kvr(np)
              if (ukvr<1e-20) ukvr=0D0
              v_darcy=0D0
              if (ukvr*Dq>floweps) then
                v_darcy = Dq * ukvr * dphi
                ! store volumetric rate for ss_fluid_fluxes()
                qsrc_phase(1) = -1.d0*v_darcy
                Res(1) = Res(1) - v_darcy* auxvar%den(np)* &
!                 auxvar%xmol((np-1)*option%nflowspec+1)*
                  option%flow_dt
                Res(2) = Res(2) - v_darcy* auxvar%den(np)* &
!                 auxvar%xmol((np-1)*option%nflowspec+2)*
                  option%flow_dt
                if (energy_flag) Res(3) = Res(3) - v_darcy * auxvar%den(np)* &
                  auxvar%h(np)*option%flow_dt
              ! print *,'produce: ',np,v_darcy
              endif
            endif
          enddo
        endif
      endif 
     !print *,'well-prod: ',  auxvar%pres,psrc(1), res
    ! injection well (well status = 2)
      if ( dabs(well_status - 2.D0) < 1.D-1) then 

        call EOSWaterDensity(tsrc,auxvar%pres,dw_kg,dw_mol,ierr)
        call EOSWaterEnthalpy(tsrc,auxvar%pres,enth_src_h2o,ierr)
        ! J/kmol -> whatever units
        enth_src_h2o = enth_src_h2o * option%scale

        Dq = msrc(2) ! well parameter, read in input file
                      ! Take the place of 2nd parameter 
        ! Flow term
        if ( auxvar%pres < pressure_max)then  
          do np = 1, option%nphase
            dphi = pressure_bh - auxvar%pres + auxvar%pc(np)
            if (dphi>=0.D0) then ! outflow only
              ukvr = auxvar%kvr(np)
              v_darcy=0.D0
              if (ukvr*Dq>floweps) then
                v_darcy = Dq * ukvr * dphi
                ! store volumetric rate for ss_fluid_fluxes()
                qsrc_phase(1) = v_darcy
                Res(1) = Res(1) + v_darcy* auxvar%den(np)* &
!                 auxvar%xmol((np-1)*option%nflowspec+1) * option%flow_dt
!                 (1.d0-csrc) * option%flow_dt
                  option%flow_dt
                Res(2) = Res(2) + v_darcy* auxvar%den(np)* &
!                 auxvar%xmol((np-1)*option%nflowspec+2) * option%flow_dt
!                 csrc * option%flow_dt
                  option%flow_dt
!               if (energy_flag) Res(3) = Res(3) + v_darcy*auxvar%den(np)* &
!                 auxvar%h(np)*option%flow_dt
                if (energy_flag) Res(3) = Res(3) + v_darcy*auxvar%den(np)* &
                  enth_src_h2o*option%flow_dt
                
!               print *,'inject: ',np,v_darcy
              endif
            endif
          enddo
        endif
      endif    
    case default
      print *,'Unrecognized Source/Sink condition: ', isrctype 
  end select      
  deallocate(msrc)
      
end subroutine ImmisSourceSink

! ************************************************************************** !

subroutine ImmisFlux(auxvar_up,por_up,tor_up,sir_up,dd_up,perm_up,Dk_up, &
                        auxvar_dn,por_dn,tor_dn,sir_dn,dd_dn,perm_dn,Dk_dn, &
                        area,dist_gravity,upweight, &
                        option,vv_darcy,Res)
  ! 
  ! Computes the internal flux terms for the residual
  ! 
  ! Author: Chuan Lu
  ! Date: 10/12/08
  ! 
  use Option_module                              
  
  implicit none
  
  type(Immis_auxvar_elem_type) :: auxvar_up, auxvar_dn
  type(option_type) :: option
  PetscReal :: sir_up(:), sir_dn(:)
  PetscReal :: por_up, por_dn
  PetscReal :: tor_up, tor_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: Dk_up, Dk_dn
  PetscReal :: vv_darcy(:),area
  PetscReal :: Res(1:option%nflowdof) 
  PetscReal :: dist_gravity  ! distance along gravity vector
     
  PetscInt :: ispec, np, ind
  PetscReal :: fluxm(option%nflowspec),fluxe,q, v_darcy
  PetscReal :: uh,uxmol(1:option%nflowspec),ukvr,difff,diffdp, DK,Dq
  PetscReal :: upweight,density_ave,cond,gravity,dphi
     
  Dq = (perm_up * perm_dn)/(dd_up*perm_dn + dd_dn*perm_up)
!  diffdp = (por_up *tor_up * por_dn*tor_dn) / (dd_dn*por_up*tor_up + dd_up*por_dn*tor_dn)*area
  
  fluxm = 0.D0
  fluxe = 0.D0
  vv_darcy =0.D0 
  
! Flow term
  do np = 1, option%nphase
!   if (auxvar_up%sat(np) > sir_up(np) .or. auxvar_dn%sat(np) > sir_dn(np)) then
    if ((auxvar_up%kvr(np) + auxvar_dn%kvr(np)) > eps) then
      upweight= dd_dn/(dd_up+dd_dn)
      if (auxvar_up%sat(np) <eps) then
        upweight=0.d0
      else if (auxvar_dn%sat(np) <eps) then
        upweight=1.d0
      endif
      density_ave = upweight*auxvar_up%den(np) + (1.D0-upweight)*auxvar_dn%den(np)
        
      gravity = (upweight*auxvar_up%den(np) * auxvar_up%avgmw(np) + &
             (1.D0-upweight)*auxvar_dn%den(np) * auxvar_dn%avgmw(np)) &
             * dist_gravity

      dphi = auxvar_up%pres - auxvar_dn%pres &
             - auxvar_up%pc(np) + auxvar_dn%pc(np) &
             + gravity

      v_darcy = 0.D0
      ukvr = 0.D0
      uh = 0.D0
      uxmol = 0.D0

      ! note uxmol only contains one phase xmol
      if (dphi >= 0.D0) then
        ukvr = auxvar_up%kvr(np)
        ! if (option%use_isothermal == PETSC_FALSE)&
        uh = auxvar_up%h(np)
      else
        ukvr = auxvar_dn%kvr(np)
        ! if (option%use_isothermal == PETSC_FALSE)&
        uh = auxvar_dn%h(np)
      endif
   

      if (ukvr > floweps) then
        v_darcy = Dq * ukvr * dphi
        vv_darcy(np) = v_darcy
        q = v_darcy * area
        fluxm(np)=fluxm(np) + q * density_ave
        ! if (option%use_isothermal == PETSC_FALSE)&
        fluxe = fluxe + q*density_ave*uh
      endif
    endif

#if 0 
! Diffusion term   
! Note : average rule may not be correct  
    if ((auxvar_up%sat(np) > eps) .and. (auxvar_dn%sat(np) > eps)) then
      difff = diffdp * 0.25D0*(auxvar_up%sat(np) + auxvar_dn%sat(np))* &
             (auxvar_up%den(np) + auxvar_dn%den(np))
      do ispec=1, option%nflowspec
        ind = ispec + (np-1)*option%nflowspec
        fluxm(ispec) = fluxm(ispec) + difff * .5D0 * &
                (auxvar_up%diff(ind) + auxvar_dn%diff(ind))* &
                (auxvar_up%xmol(ind) - auxvar_dn%xmol(ind))
      enddo
    endif
#endif
  enddo

! conduction term
  !if (option%use_isothermal == PETSC_FALSE) then     
  Dk = (Dk_up * Dk_dn) / (dd_dn*Dk_up + dd_up*Dk_dn)
  cond = Dk*area*(auxvar_up%temp-auxvar_dn%temp)
  fluxe = fluxe + cond
 ! end if

  !if (option%use_isothermal)then
  !   Res(1:option%nflowdof) = fluxm(:) * option%flow_dt
 ! else
  Res(1:option%nphase) = fluxm(:) * option%flow_dt
  Res(option%nflowdof) = fluxe * option%flow_dt
 ! end if
 ! note: Res is the flux contribution, for node 1 R = R + Res_FL
 !                                              2 R = R - Res_FL  

end subroutine ImmisFlux

! ************************************************************************** !

subroutine ImmisBCFlux(ibndtype,auxvars,auxvar_up,auxvar_dn, &
     por_dn,tor_dn,sir_dn,dd_up,perm_dn,Dk_dn, &
     area,dist_gravity,option,vv_darcy,Res)
  ! 
  ! Computes the  boundary flux terms for the residual
  ! 
  ! Author: Chuan Lu
  ! Date: 10/12/08
  ! 
  use Option_module
  
  implicit none
  
  PetscInt :: ibndtype(:)
  type(Immis_auxvar_elem_type) :: auxvar_up, auxvar_dn
  type(option_type) :: option
  PetscReal :: dd_up, sir_dn(:)
  PetscReal :: auxvars(:) ! from aux_real_var array
  PetscReal :: por_dn,perm_dn,Dk_dn,tor_dn
  PetscReal :: vv_darcy(:), area
  PetscReal :: Res(1:option%nflowdof) 
  
  PetscReal :: dist_gravity  ! distance along gravity vector
          
  PetscInt :: ispec, np
  PetscReal :: fluxm(option%nflowspec),fluxe,q,density_ave, v_darcy
  PetscReal :: uh,uxmol(1:option%nflowspec),ukvr,diff,diffdp,DK,Dq
  PetscReal :: upweight,cond,gravity,dphi

  fluxm = 0.d0
  fluxe = 0.d0

  do np = 1, option%nphase
    v_darcy = 0.d0
    density_ave = 0.d0
    q = 0.d0
    ukvr=0.D0
    select case(ibndtype(MPH_PRESSURE_DOF))
    case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC)
      Dq = perm_dn / dd_up

      ! only consider phase that exists, this also deals with IPHASE=1,2
!     if (auxvar_up%sat(np) > sir_dn(np) .or. auxvar_dn%sat(np) > sir_dn(np)) then
      if ((auxvar_up%kvr(np) + auxvar_dn%kvr(np)) > eps) then

        ! upweight according to saturation?
        upweight = 1.D0
        if (auxvar_up%sat(np) < eps) then 
          upweight = 0.d0
        else if (auxvar_dn%sat(np) < eps) then 
          upweight = 1.d0
        endif
        density_ave = upweight*auxvar_up%den(np) + (1.D0-upweight)*auxvar_dn%den(np)

        gravity = (upweight*auxvar_up%den(np) * auxvar_up%avgmw(np) + &
             (1.D0-upweight)*auxvar_dn%den(np) * auxvar_dn%avgmw(np)) &
             * dist_gravity

        ! calculate the pressure gradient
        dphi = auxvar_up%pres - auxvar_dn%pres &
             - auxvar_up%pc(np) + auxvar_dn%pc(np) &
             + gravity

        ! upwind rel perm by the pressure gradient
        if (dphi>=0.D0) then
          ukvr = auxvar_up%kvr(np)
        else
          ukvr = auxvar_dn%kvr(np)
        endif

        if (ukvr*Dq>floweps) then
          v_darcy = Dq * ukvr * dphi
        endif
      endif

    case(NEUMANN_BC) ! fixed by etc, 1/23/2012
      ! only consider phase that exists, this also deals with IPHASE=1,2
      if (auxvar_up%sat(np) > sir_dn(np) .or. auxvar_dn%sat(np) > sir_dn(np)) then

        ! upwind by imposed Neumann velocity
        if (dabs(auxvars(MPH_PRESSURE_DOF)) > floweps) then
          v_darcy = auxvars(MPH_PRESSURE_DOF)
          if (v_darcy > 0.d0) then
            density_ave = auxvar_up%den(np)
          else
            density_ave = auxvar_dn%den(np)
          endif
        endif
      end if
     end select

     q = v_darcy * area
     vv_darcy(np) = v_darcy
     uh=0.D0
     uxmol=0.D0

     if (v_darcy >= 0.D0) then
        !if (option%use_isothermal == PETSC_FALSE)&
         uh = auxvar_up%h(np)
        ! uxmol(:)=auxvar_up%xmol((np-1)*option%nflowspec+1 : np * option%nflowspec)
     else
         !if (option%use_isothermal == PETSC_FALSE)&
        uh = auxvar_dn%h(np)
         ! uxmol(:)=auxvar_dn%xmol((np-1)*option%nflowspec+1 : np * option%nflowspec)
     endif

        fluxm(np) = fluxm(np) + q*density_ave ! *uxmol(ispec)

      !if (option%use_isothermal == PETSC_FALSE) &
      fluxe = fluxe + q*density_ave*uh
 !print *,'FLBC', ibndtype(1),np, ukvr, v_darcy, uh, uxmol
   enddo

#if 0 
    ! Diffusion term   
  select case(ibndtype(3))
  case(DIRICHLET_BC) 
     !      if (auxvar_up%sat > eps .and. auxvar_dn%sat > eps) then
     !        diff = diffdp * 0.25D0*(auxvar_up%sat+auxvar_dn%sat)*(auxvar_up%den+auxvar_dn%den)
        do np = 1, option%nphase
          if (auxvar_up%sat(np)>eps .and. auxvar_dn%sat(np)>eps)then
              diff =diffdp * 0.25D0*(auxvar_up%sat(np)+auxvar_dn%sat(np))*&
                    (auxvar_up%den(np)+auxvar_up%den(np))
           do ispec = 1, option%nflowspec
              fluxm(ispec) = fluxm(ispec) + diff * auxvar_dn%diff((np-1)* option%nflowspec+ispec)* &
                   (auxvar_up%xmol((np-1)* option%nflowspec+ispec) &
                   -auxvar_dn%xmol((np-1)* option%nflowspec+ispec))
           enddo
          endif         
        enddo

  end select
#endif
  ! Conduction term
! if (option%use_isothermal == PETSC_FALSE) then
    select case(ibndtype(2))
    case(DIRICHLET_BC)
       Dk =  Dk_dn / dd_up
       cond = Dk*area*(auxvar_up%temp - auxvar_dn%temp) 
       fluxe = fluxe + cond
    case(NEUMANN_BC)
       fluxe = fluxe + auxvars(2)*area*option%scale
       ! from W to MW, Added by Satish Karra 10/19/11
    case(ZERO_GRADIENT_BC)
       ! No change in fluxe	
    end select
! end if

  Res(1:option%nphase)=fluxm(:)* option%flow_dt
  Res(option%nflowdof)=fluxe * option%flow_dt

end subroutine ImmisBCFlux

! ************************************************************************** !

subroutine ImmisResidual(snes,xx,r,realization,ierr)
  ! 
  ! Computes the residual equation
  ! 
  ! Author: Chuan Lu
  ! Date: 10/10/08
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Discretization_module
  use Field_module
  use Option_module
  use Grid_module 

  implicit none

  SNES :: snes
  Vec :: xx
  Vec :: r
  type(realization_subsurface_type) :: realization
  PetscErrorCode :: ierr
  
  type(discretization_type), pointer :: discretization
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(patch_type), pointer :: cur_patch
  PetscInt :: ichange  

  field => realization%field
  grid => realization%patch%grid
  option => realization%option
  discretization => realization%discretization
  
 
!  call DiscretizationGlobalToLocal(discretization,xx,field%flow_xx_loc,NFLOWDOF)
 ! check initial guess -----------------------------------------------
  ierr = ImmisInitGuessCheck(realization)
  if (ierr<0)then
    !ierr = PETSC_ERR_ARG_OUTOFRANGE
    if (option%myrank==0) print *,'table out of range: ',ierr
    call SNESSetFunctionDomainError() 
    return
  endif 
  ! end check ---------------------------------------------------------

  ! Communication -----------------------------------------
  ! These 3 must be called before ImmisUpdateAuxVars()
  call DiscretizationGlobalToLocal(discretization,xx,field%flow_xx_loc,NFLOWDOF)
  call DiscretizationLocalToLocal(discretization,field%icap_loc,field%icap_loc,ONEDOF)

!geh refactor  call DiscretizationLocalToLocal(discretization,field%perm_xx_loc,field%perm_xx_loc,ONEDOF)
!geh refactor  call DiscretizationLocalToLocal(discretization,field%perm_yy_loc,field%perm_yy_loc,ONEDOF)
!geh refactor  call DiscretizationLocalToLocal(discretization,field%perm_zz_loc,field%perm_zz_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%ithrm_loc,field%ithrm_loc,ONEDOF)
  
  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call ImmisResidualPatch(snes,xx,r,realization,ierr)
    cur_patch => cur_patch%next
  enddo

end subroutine ImmisResidual

! ************************************************************************** !

subroutine ImmisResidualPatch(snes,xx,r,realization,ierr)
  ! 
  ! Computes the residual equation at patch level
  ! 
  ! Author: Chuan Lu
  ! Date: 10/10/08
  ! 

  use Connection_module
  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module  
  use Field_module
  use Debug_module
  
  implicit none

  SNES, intent(in) :: snes
  Vec, intent(inout) :: xx
  Vec, intent(out) :: r
  type(realization_subsurface_type) :: realization

  PetscErrorCode :: ierr
  PetscInt :: i, iphase, jn
  PetscInt :: ip1, ip2
  PetscInt :: local_id, ghosted_id, local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn

  PetscReal, pointer ::accum_p(:)

  PetscReal, pointer :: r_p(:), porosity_loc_p(:), volume_p(:), &
               xx_loc_p(:), xx_p(:), yy_p(:),&
               tortuosity_loc_p(:),&
               perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)
                          
               
  PetscReal, pointer :: icap_loc_p(:), ithrm_loc_p(:)

  PetscInt :: icap_up, icap_dn, ithrm_up, ithrm_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: dd, f_up, f_dn, ff
  PetscReal :: perm_up, perm_dn
  PetscReal :: D_up, D_dn  ! "Diffusion" constants at upstream, downstream faces.
  PetscReal :: dw_kg, dw_mol,dddt,dddp
  PetscReal :: tsrc1, qsrc1, csrc1, enth_src_h2o, enth_src_co2 , hsrc1
  PetscReal :: rho, fg, dfgdp, dfgdt, eng, dhdt, dhdp, visc, dvdt, dvdp, xphi
  PetscReal :: upweight
  PetscReal :: Res(realization%option%nflowdof), v_darcy(realization%option%nphase)
  PetscReal :: xxbc(realization%option%nflowdof)
! PetscReal :: msrc(1:realization%option%nflowspec)
  PetscReal :: psrc(1:realization%option%nphase)
  PetscViewer :: viewer


  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(immis_type), pointer :: immis
  type(Immis_parameter_type), pointer :: immis_parameter
  
  type(Immis_auxvar_type), pointer :: auxvars(:), auxvars_bc(:), auxvars_ss(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscReal, pointer :: msrc(:)

  PetscBool :: enthalpy_flag
  PetscInt :: ng
  PetscInt :: iconn, idof, istart, iend
  PetscInt :: nsrcpara
  PetscInt :: sum_connection
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity
  PetscReal :: ss_flow_vol_flux(realization%option%nphase)
  
  character(len=MAXSTRINGLENGTH) :: string

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  immis => patch%aux%Immis
  immis_parameter => patch%aux%Immis%immis_parameter
  auxvars => patch%aux%Immis%auxvars
  auxvars_bc => patch%aux%Immis%auxvars_bc
  auxvars_ss => patch%aux%Immis%auxvars_ss
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss
  

 ! call ImmisUpdateAuxVarsPatchNinc(realization)
  ! override flags since they will soon be out of date  
 ! patch%ImmisAux%auxvars_up_to_date = PETSC_FALSE 
 
  if (option%compute_mass_balance_new) then
    call ImmisZeroMassBalDeltaPatch(realization)
  endif

! now assign access pointer to local variables
  call VecGetArrayF90(field%flow_xx_loc, xx_loc_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(r, r_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
 
! call VecGetArrayF90(field%flow_yy,yy_p,ierr)
!geh refactor  call VecGetArrayF90(field%porosity_loc, porosity_loc_p, ierr)
!geh refactor  call VecGetArrayF90(field%tortuosity_loc, tortuosity_loc_p, ierr)
!geh refactor  call VecGetArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
!geh refactor  call VecGetArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
!geh refactor  call VecGetArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
!geh refactor  call VecGetArrayF90(field%volume, volume_p, ierr)
  call VecGetArrayF90(field%ithrm_loc, ithrm_loc_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%icap_loc, icap_loc_p, ierr);CHKERRQ(ierr)
 
 
! Multiphase flash calculation is more expensive, so calculate once per iteration
#if 1
  ! Pertubations for aux terms --------------------------------
  do ng = 1, grid%ngmax
    if (grid%nG2L(ng) < 0) cycle
    if (associated(patch%imat)) then
      if (patch%imat(ng) <= 0) cycle
    endif
        
    istart = (ng-1)*option%nflowdof + 1; iend = istart - 1 + option%nflowdof
    call ImmisAuxVarCompute_Ninc(xx_loc_p(istart:iend),auxvars(ng)%auxvar_elem(0), &
      patch%saturation_function_array(int(icap_loc_p(ng)))%ptr, &
      realization%fluid_properties,option)

    if (option%flow%numerical_derivatives) then
      patch%aux%Immis%delx(1,ng) = xx_loc_p((ng-1)*option%nflowdof+1)*dfac * 1.D-3
        patch%aux%Immis%delx(2,ng) = xx_loc_p((ng-1)*option%nflowdof+2)*dfac
 
      if (xx_loc_p((ng-1)*option%nflowdof+3) <=0.9)then
        patch%aux%Immis%delx(3,ng) = dfac*xx_loc_p((ng-1)*option%nflowdof+3) 
      else
        patch%aux%Immis%delx(3,ng) = -dfac*xx_loc_p((ng-1)*option%nflowdof+3) 
      endif
           
      if (patch%aux%Immis%delx(3,ng) < 1D-12 .and. patch%aux%Immis%delx(3,ng) >= 0.D0) &
        patch%aux%Immis%delx(3,ng) = 1D-12
      if (patch%aux%Immis%delx(3,ng) > -1D-12 .and. patch%aux%Immis%delx(3,ng) < 0.D0) &
        patch%aux%Immis%delx(3,ng) = -1D-12
        
      if ((patch%aux%Immis%delx(3,ng)+xx_loc_p((ng-1)*option%nflowdof+3))>1.D0) then
        patch%aux%Immis%delx(3,ng) = (1.D0-xx_loc_p((ng-1)*option%nflowdof+3))*1D-6
      endif
      if (( patch%aux%Immis%delx(3,ng)+xx_loc_p((ng-1)*option%nflowdof+3))<0.D0)then
        patch%aux%Immis%delx(3,ng) = xx_loc_p((ng-1)*option%nflowdof+3)*1D-6
      endif
      call ImmisAuxVarCompute_Winc(xx_loc_p(istart:iend),patch%aux%Immis%delx(:,ng), &
          auxvars(ng)%auxvar_elem(1:option%nflowdof), &
          patch%saturation_function_array(int(icap_loc_p(ng)))%ptr, &
          realization%fluid_properties,option)
    endif
  enddo
#endif

   patch%aux%Immis%res_old_AR=0.D0
   patch%aux%Immis%res_old_FL=0.D0
   r_p = 0.d0

#if 1
  ! Accumulation terms ------------------------------------
  r_p = - accum_p

  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    iend = local_id*option%nflowdof
    istart = iend-option%nflowdof+1
    call ImmisAccumulation(auxvars(ghosted_id)%auxvar_elem(0),porosity_loc_p(ghosted_id), &
                              volume_p(local_id), &
                              immis_parameter%dencpr(int(ithrm_loc_p(ghosted_id))), &
                              option,ONE_INTEGER,Res) 
    r_p(istart:iend) = r_p(istart:iend) + Res(1:option%nflowdof)
    !print *,'REs, acm: ', res
    patch%aux%Immis%res_old_AR(local_id, :) = Res(1:option%nflowdof)
  enddo
#endif

#if 1
  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sink_list%first 
  sum_connection = 0
  do 
    if (.not.associated(source_sink)) exit
    !print *, 'RES s/s begin'
    ! check whether enthalpy dof is included
  !  if (source_sink%flow_condition%num_sub_conditions > 3) then
      enthalpy_flag = PETSC_TRUE
   ! else
   !   enthalpy_flag = PETSC_FALSE
   ! endif
      
    if (associated(source_sink%flow_condition%pressure)) then
      psrc(:) = source_sink%flow_condition%pressure%dataset%rarray(:)
    endif
!    qsrc1 = source_sink%flow_condition%pressure%dataset%rarray(1)
    tsrc1 = source_sink%flow_condition%temperature%dataset%rarray(1)
    csrc1 = source_sink%flow_condition%concentration%dataset%rarray(1)
    if (enthalpy_flag) hsrc1 = source_sink%flow_condition%enthalpy%dataset%rarray(1)
!     hsrc1=0D0
!     qsrc1 = qsrc1 / FMWH2O ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
!     csrc1 = csrc1 / FMWCO2
!     msrc(1) = qsrc1; msrc(2) =csrc1
!     msrc(:) = psrc(:)
!     msrc(1) = msrc(1) / FMWH2O
!     msrc(2) = msrc(2) / FMWCO2
      
!clu add
    select case(source_sink%flow_condition%itype(1))
      case(MASS_RATE_SS)
      msrc => source_sink%flow_condition%rate%dataset%rarray
      nsrcpara = 2
    case(WELL_SS)
      msrc => source_sink%flow_condition%well%dataset%rarray
      nsrcpara = 7 + option%nflowspec 
     
!     print *,'src/sink: ',nsrcpara,msrc
    case default
      print *, 'IMS mode does not support source/sink type: ', source_sink%flow_condition%itype(1)
      stop  
    end select

!clu end change

    cur_connection_set => source_sink%connection_set
    do iconn = 1, cur_connection_set%num_connections      
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle
      sum_connection = sum_connection + 1
      call ImmisSourceSink(msrc,nsrcpara,psrc,tsrc1,hsrc1,auxvars(ghosted_id)%auxvar_elem(0),&
            source_sink%flow_condition%itype(1),Res, &
            ss_flow_vol_flux, &
            enthalpy_flag, option)
      if (associated(patch%ss_flow_fluxes)) then
        patch%ss_flow_fluxes(:,sum_connection) = Res(:)/option%flow_dt
      endif
      if (associated(patch%ss_flow_vol_fluxes)) then
        patch%ss_flow_vol_fluxes(:,sum_connection) = ss_flow_vol_flux/option%flow_dt
      endif
      if (option%compute_mass_balance_new) then
        global_auxvars_ss(sum_connection)%mass_balance_delta(:,1) = &
          global_auxvars_ss(sum_connection)%mass_balance_delta(:,1) - &
            Res(:)/option%flow_dt
      endif
 
      r_p((local_id-1)*option%nflowdof + jh2o) = &
           r_p((local_id-1)*option%nflowdof + jh2o) - Res(jh2o)
      r_p((local_id-1)*option%nflowdof + jco2) = &
           r_p((local_id-1)*option%nflowdof + jco2) - Res(jco2)
      patch%aux%Immis%res_old_AR(local_id,jh2o) = &
           patch%aux%Immis%res_old_AR(local_id,jh2o) - Res(jh2o)
      patch%aux%Immis%res_old_AR(local_id,jco2) = &
           patch%aux%Immis%res_old_AR(local_id,jco2) - Res(jco2)
      if (enthalpy_flag) then
        r_p( local_id*option%nflowdof) = r_p(local_id*option%nflowdof) - Res(option%nflowdof)
           patch%aux%Immis%res_old_AR(local_id,option%nflowdof) = &
             patch%aux%Immis%res_old_AR(local_id,option%nflowdof) - Res(option%nflowdof)
      endif 
  !  else if (qsrc1 < 0.d0) then ! withdrawal
  !  endif
    enddo
    source_sink => source_sink%next
  enddo
#endif

#if 1
  ! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection_set
        
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif

      if (ghosted_id<=0) then
        print *, "Wrong boundary node index... STOP!!!"
        stop
      endif

      ithrm_dn = int(ithrm_loc_p(ghosted_id))
      D_dn = immis_parameter%ckwet(ithrm_dn)

      ! for now, just assume diagonal tensor
      perm_dn = perm_xx_loc_p(ghosted_id)*abs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id)*abs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id)*abs(cur_connection_set%dist(3,iconn))
      ! dist(0,iconn) = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3,iconn) = vector(3) - unit vector
      distance_gravity = cur_connection_set%dist(0,iconn) * &
                         dot_product(option%gravity, &
                                     cur_connection_set%dist(1:3,iconn))

      icap_dn = int(icap_loc_p(ghosted_id))  
! Then need fill up increments for BCs
       do idof = 1, option%nflowdof
         select case(boundary_condition%flow_condition%itype(idof))
           case(DIRICHLET_BC)
             xxbc(idof) = boundary_condition%flow_aux_real_var(idof,iconn)
           case(HYDROSTATIC_BC)
             xxbc(MPH_PRESSURE_DOF) = boundary_condition%flow_aux_real_var(MPH_PRESSURE_DOF,iconn)
             if (idof>=MPH_TEMPERATURE_DOF)then
               xxbc(idof) = xx_loc_p((ghosted_id-1)*option%nflowdof+idof)
             endif
           case(NEUMANN_BC, ZERO_GRADIENT_BC)
          ! solve for pb from Darcy's law given qb /= 0
             xxbc(idof) = xx_loc_p((ghosted_id-1)*option%nflowdof+idof)
         end select
      enddo

 
      call ImmisAuxVarCompute_Ninc(xxbc,auxvars_bc(sum_connection)%auxvar_elem(0),&
         patch%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr,&
         realization%fluid_properties, option)

      call ImmisBCFlux(boundary_condition%flow_condition%itype, &
         boundary_condition%flow_aux_real_var(:,iconn), &
         auxvars_bc(sum_connection)%auxvar_elem(0), &
         auxvars(ghosted_id)%auxvar_elem(0), &
         porosity_loc_p(ghosted_id), &
         tortuosity_loc_p(ghosted_id), &
         immis_parameter%sir(:,icap_dn), &
         cur_connection_set%dist(0,iconn),perm_dn,D_dn, &
         cur_connection_set%area(iconn), &
         distance_gravity,option, &
         v_darcy,Res)

      if (option%compute_mass_balance_new) then
        ! contribution to boundary
        do iphase = 1, option%nphase
          global_auxvars_bc(sum_connection)%mass_balance_delta(iphase,iphase) = &
            global_auxvars_bc(sum_connection)%mass_balance_delta(iphase,iphase) - &
            Res(iphase)/option%flow_dt
        enddo
      endif

      patch%boundary_velocities(:,sum_connection) = v_darcy(:)
      if (associated(patch%boundary_flow_fluxes)) then
        patch%boundary_flow_fluxes(:,sum_connection) = Res(:)
      endif        
      iend = local_id*option%nflowdof
      istart = iend-option%nflowdof+1
      r_p(istart:iend) = r_p(istart:iend) - Res(1:option%nflowdof)
      patch%aux%Immis%res_old_AR(local_id,1:option%nflowdof) = &
        patch%aux%Immis%res_old_AR(local_id,1:option%nflowdof) - Res(1:option%nflowdof)
    enddo
    boundary_condition => boundary_condition%next
  enddo
#endif

#if 1
  ! Interior Flux Terms -----------------------------------
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0  
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1

      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping   

      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id_up) <= 0 .or.  &
            patch%imat(ghosted_id_dn) <= 0) cycle
      endif

      fraction_upwind = cur_connection_set%dist(-1,iconn)
      distance = cur_connection_set%dist(0,iconn)
      ! distance = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3,iconn) = vector(3) - unit vector
      distance_gravity = distance * &
                         dot_product(option%gravity, &
                                     cur_connection_set%dist(1:3,iconn))
      dd_up = distance*fraction_upwind
      dd_dn = distance-dd_up ! should avoid truncation error
      ! upweight could be calculated as 1.d0-fraction_upwind
      ! however, this introduces ever so slight error causing pflow-overhaul not
      ! to match pflow-orig.  This can be changed to 1.d0-fraction_upwind
      upweight = dd_dn/(dd_up+dd_dn)
        
      ! for now, just assume diagonal tensor
      perm_up = perm_xx_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(3,iconn))

      perm_dn = perm_xx_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(3,iconn))

      ithrm_up = int(ithrm_loc_p(ghosted_id_up))
      ithrm_dn = int(ithrm_loc_p(ghosted_id_dn))
      icap_up = int(icap_loc_p(ghosted_id_up))
      icap_dn = int(icap_loc_p(ghosted_id_dn))
   
      D_up = immis_parameter%ckwet(ithrm_up)
      D_dn = immis_parameter%ckwet(ithrm_dn)

      call ImmisFlux(auxvars(ghosted_id_up)%auxvar_elem(0),porosity_loc_p(ghosted_id_up), &
                tortuosity_loc_p(ghosted_id_up),immis_parameter%sir(:,icap_up), &
                dd_up,perm_up,D_up, &
                auxvars(ghosted_id_dn)%auxvar_elem(0),porosity_loc_p(ghosted_id_dn), &
                tortuosity_loc_p(ghosted_id_dn),immis_parameter%sir(:,icap_dn), &
                dd_dn,perm_dn,D_dn, &
                cur_connection_set%area(iconn),distance_gravity, &
                upweight,option,v_darcy,Res)

      patch%internal_velocities(:,sum_connection) = v_darcy(:)
      if (associated(patch%internal_flow_fluxes)) then
        patch%internal_flow_fluxes(:,sum_connection) = Res(:)
      endif      
      patch%aux%Immis%res_old_FL(sum_connection,1:option%nflowdof)= Res(1:option%nflowdof)
 
      if (local_id_up>0) then
        iend = local_id_up*option%nflowdof
        istart = iend-option%nflowdof+1
        r_p(istart:iend) = r_p(istart:iend) + Res(1:option%nflowdof)
      endif
   
      if (local_id_dn>0) then
        iend = local_id_dn*option%nflowdof
        istart = iend-option%nflowdof+1
        r_p(istart:iend) = r_p(istart:iend) - Res(1:option%nflowdof)
      endif

    enddo
    cur_connection_set => cur_connection_set%next
  enddo    
#endif

! adjust residual to R/dt
  select case (option%idt_switch) 
  case(1) 
    r_p(:) = r_p(:)/option%flow_dt
  case(-1)
    if (option%flow_dt>1.D0) r_p(:) = r_p(:)/option%flow_dt
  end select
  
  do local_id = 1, grid%nlmax
    if (associated(patch%imat)) then
      if (patch%imat(grid%nL2G(local_id)) <= 0) cycle
    endif

    istart = 1 + (local_id-1)*option%nflowdof
    if (volume_p(local_id)>1.D0) r_p (istart:istart+2)=r_p(istart:istart+2)/volume_p(local_id)
    if (r_p(istart) >1E20 .or. r_p(istart) <-1E20) print *, r_p (istart:istart+2)
  enddo

! print *,'finished rp vol scale'
!  if (option%use_isothermal) then
#ifdef ISOTHERMAL
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)   ! corresponding ghost index
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
      istart = 3 + (local_id-1)*option%nflowdof
      r_p(istart) = 0.D0 ! xx_loc_p(2 + (ng-1)*option%nflowdof) - yy_p(p1-1)
    enddo
#endif
!  endif
  !call VecRestoreArrayF90(r, r_p, ierr)


  if (patch%aux%Immis%inactive_cells_exist) then
    do i=1,patch%aux%Immis%n_zero_rows
      r_p(patch%aux%Immis%zero_rows_local(i)) = 0.d0
    enddo
  endif

  call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)
! call VecRestoreArrayF90(field%flow_yy, yy_p, ierr)
  call VecRestoreArrayF90(field%flow_xx_loc, xx_loc_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
!geh refactor  call VecRestoreArrayF90(field%porosity_loc, porosity_loc_p, ierr)
!geh refactor  call VecRestoreArrayF90(field%tortuosity_loc, tortuosity_loc_p, ierr)
!geh refactor  call VecRestoreArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
!geh refactor  call VecRestoreArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
!geh refactor  call VecRestoreArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
!geh refactor  call VecRestoreArrayF90(field%volume, volume_p, ierr)
  call VecRestoreArrayF90(field%ithrm_loc, ithrm_loc_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%icap_loc, icap_loc_p, ierr);CHKERRQ(ierr)

  if (realization%debug%vecview_residual) then
    string = 'Iresidual'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call VecView(r,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
  if (realization%debug%vecview_solution) then
    string = 'Ixx'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call VecView(xx,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
end subroutine ImmisResidualPatch

! ************************************************************************** !

subroutine ImmisJacobian(snes,xx,A,B,realization,ierr)
  ! 
  ! Computes the Jacobian
  ! 
  ! Author: Chuan Lu
  ! Date: 10/10/08
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module

  implicit none

  SNES :: snes
  Vec :: xx
  Mat :: A, B, J
  type(realization_subsurface_type) :: realization
  PetscErrorCode :: ierr
  
  type(patch_type), pointer :: cur_patch
  type(grid_type),  pointer :: grid
  
  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call ImmisJacobianPatch(snes,xx,A,B,realization,ierr)
    cur_patch => cur_patch%next
  enddo

end subroutine ImmisJacobian

! ************************************************************************** !

subroutine ImmisJacobianPatch(snes,xx,A,B,realization,ierr)
  ! 
  ! Computes the Jacobian
  ! 
  ! Author: Chuan Lu
  ! Date: 10/13/08
  ! 

  use Connection_module
  use Option_module
  use Grid_module
  use Realization_Subsurface_class
  use Patch_module
  use Coupler_module
  use Field_module
  use Debug_module
  
  implicit none

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  type(realization_subsurface_type) :: realization

  PetscErrorCode :: ierr
  PetscInt :: nvar,neq,nr
  PetscInt :: ithrm_up, ithrm_dn, i
  PetscInt :: ip1, ip2 

  PetscReal, pointer :: porosity_loc_p(:), volume_p(:), &
                          xx_loc_p(:), tortuosity_loc_p(:),&
                          perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)
  PetscReal, pointer :: icap_loc_p(:), ithrm_loc_p(:)
  PetscInt :: icap,icap_up,icap_dn
  PetscInt :: ii, jj
  PetscReal :: dw_kg,dw_mol,enth_src_co2,enth_src_h2o,rho
  PetscReal :: tsrc1,qsrc1,csrc1,hsrc1
  PetscReal :: dd_up, dd_dn, dd, f_up, f_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: dw_dp,dw_dt,hw_dp,hw_dt,dresT_dp,dresT_dt
  PetscReal :: D_up, D_dn  ! "Diffusion" constants upstream and downstream of a face.
  PetscReal :: zero, norm
  PetscReal :: upweight
  PetscReal :: max_dev  
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  PetscInt ::  natural_id_up,natural_id_dn
  
  PetscReal :: Jup(1:realization%option%nflowdof,1:realization%option%nflowdof), &
            Jdn(1:realization%option%nflowdof,1:realization%option%nflowdof)
  
  PetscInt :: istart, iend
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscBool :: enthalpy_flag
  PetscInt :: iconn, idof
  PetscInt :: sum_connection  
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity
  PetscReal :: Res(realization%option%nflowdof) 
  PetscReal :: xxbc(1:realization%option%nflowdof), delxbc(1:realization%option%nflowdof)
  PetscReal :: ResInc(realization%patch%grid%nlmax,realization%option%nflowdof,&
           realization%option%nflowdof)
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option 
  type(field_type), pointer :: field 
  type(Immis_parameter_type), pointer :: immis_parameter
  type(Immis_auxvar_type), pointer :: auxvars(:), auxvars_bc(:)
  
  PetscReal :: vv_darcy(realization%option%nphase), voltemp
  PetscReal :: ra(1:realization%option%nflowdof,1:realization%option%nflowdof*2) 
  PetscReal, pointer :: msrc(:)
! PetscReal :: msrc(1:realization%option%nflowspec)
  PetscReal :: psrc(1:realization%option%nphase), ss_flow(1:realization%option%nphase)
  PetscReal :: dddt, dddp, fg, dfgdp, dfgdt, eng, dhdt, dhdp, visc, dvdt,&
               dvdp, xphi
  PetscInt :: nsrcpara  
  
  PetscViewer :: viewer
  Vec :: debug_vec
  character(len=MAXSTRINGLENGTH) :: string

!-----------------------------------------------------------------------
! R stand for residual
!  ra       1              2              3              4          5              6            7      8
! 1: p     dR/dpi         dR/dTi          dR/dci        dR/dsi   dR/dpim        dR/dTim
! 2: T
! 3: c
! 4  s         
!-----------------------------------------------------------------------

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  immis_parameter => patch%aux%Immis%immis_parameter
  auxvars => patch%aux%Immis%auxvars
  auxvars_bc => patch%aux%Immis%auxvars_bc
  
! dropped derivatives:
!   1.D0 gas phase viscocity to all p,t,c,s
!   2. Average molecular weights to p,t,s

#if 0
!  call ImmisNumericalJacobianTest(xx,realization)
#endif

 ! print *,'*********** In Jacobian ********************** '
  call MatZeroEntries(A,ierr);CHKERRQ(ierr)

  call VecGetArrayF90(field%flow_xx_loc, xx_loc_p, ierr);CHKERRQ(ierr)
!geh refactor  call VecGetArrayF90(field%porosity_loc, porosity_loc_p, ierr)
!geh refactor  call VecGetArrayF90(field%tortuosity_loc, tortuosity_loc_p, ierr)
!geh refactor  call VecGetArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
!geh refactor  call VecGetArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
!geh refactor  call VecGetArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
!geh refactor  call VecGetArrayF90(field%volume, volume_p, ierr)

  call VecGetArrayF90(field%ithrm_loc, ithrm_loc_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%icap_loc, icap_loc_p, ierr);CHKERRQ(ierr)

 ResInc = 0.D0
#if 1
  ! Accumulation terms ------------------------------------
  do local_id = 1, grid%nlmax  ! For each local node do...
     ghosted_id = grid%nL2G(local_id)
     !geh - Ignore inactive cells with inactive materials
     if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
     endif
     iend = local_id*option%nflowdof
     istart = iend-option%nflowdof+1
     icap = int(icap_loc_p(ghosted_id))
     
     do nvar =1, option%nflowdof
        call ImmisAccumulation(auxvars(ghosted_id)%auxvar_elem(nvar), &
             porosity_loc_p(ghosted_id), &
             volume_p(local_id), &
             immis_parameter%dencpr(int(ithrm_loc_p(ghosted_id))), &
             option,ONE_INTEGER, res) 
        ResInc( local_id,:,nvar) =  ResInc(local_id,:,nvar) + Res(:)
     enddo
     
  enddo
#endif
#if 1
  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sink_list%first 
  do 
    if (.not.associated(source_sink)) exit
    
    ! check whether enthalpy dof is included
  !  if (source_sink%flow_condition%num_sub_conditions > 3) then
      enthalpy_flag = PETSC_TRUE
   ! else
   !   enthalpy_flag = PETSC_FALSE
   ! endif

    if (associated(source_sink%flow_condition%pressure)) then
      psrc(:) = source_sink%flow_condition%pressure%dataset%rarray(:)
    endif
    tsrc1 = source_sink%flow_condition%temperature%dataset%rarray(1)
    csrc1 = source_sink%flow_condition%concentration%dataset%rarray(1)
 !   hsrc1=0.D0
    if (enthalpy_flag) hsrc1 = source_sink%flow_condition%enthalpy%dataset%rarray(1)

   ! qsrc1 = qsrc1 / FMWH2O ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
   ! csrc1 = csrc1 / FMWCO2
!     msrc(:)= psrc(:)
!     msrc(1) =  msrc(1) / FMWH2O
!     msrc(2) =  msrc(2) / FMWCO2

!clu add
    select case(source_sink%flow_condition%itype(1))
      case(MASS_RATE_SS)
        msrc => source_sink%flow_condition%rate%dataset%rarray
        nsrcpara= 2
      case(WELL_SS)
        msrc => source_sink%flow_condition%well%dataset%rarray
        nsrcpara = 7 + option%nflowspec 
      case default
        print *, 'ims mode does not support source/sink type: ', source_sink%flow_condition%itype(1)
        stop  
    end select
 
    cur_connection_set => source_sink%connection_set
 
    do iconn = 1, cur_connection_set%num_connections      
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif
!     if (enthalpy_flag) then
!       r_p(local_id*option%nflowdof) = r_p(local_id*option%nflowdof) - hsrc1 * option%flow_dt   
!     endif         
      do nvar =1, option%nflowdof
        call ImmisSourceSink(msrc,nsrcpara,psrc,tsrc1,hsrc1,auxvars(ghosted_id)%auxvar_elem(nvar),&
        source_sink%flow_condition%itype(1),Res,ss_flow,enthalpy_flag, option)
      
        ResInc(local_id,jh2o,nvar)=  ResInc(local_id,jh2o,nvar) - Res(jh2o)
        ResInc(local_id,jco2,nvar)=  ResInc(local_id,jco2,nvar) - Res(jco2)
        if (enthalpy_flag) & 
           ResInc(local_id,option%nflowdof,nvar)=&
           ResInc(local_id,option%nflowdof,nvar)- Res(option%nflowdof) 

      enddo 
    enddo
    source_sink => source_sink%next
  enddo
#endif
! Boundary conditions
#if 1
  ! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection_set
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif

      if (ghosted_id<=0) then
        print *, "Wrong boundary node index... STOP!!!"
        stop
      endif

      ithrm_dn = int(ithrm_loc_p(ghosted_id))
      D_dn = immis_parameter%ckwet(ithrm_dn)

      ! for now, just assume diagonal tensor
      perm_dn = perm_xx_loc_p(ghosted_id)*abs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id)*abs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id)*abs(cur_connection_set%dist(3,iconn))
      ! dist(0,iconn) = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3,iconn) = vector(3) - unit vector
      distance_gravity = cur_connection_set%dist(0,iconn) * &
                         dot_product(option%gravity, &
                                     cur_connection_set%dist(1:3,iconn))
      icap_dn = int(icap_loc_p(ghosted_id))

! Then need fill up increments for BCs
      delxbc=0.D0;
      do idof =1, option%nflowdof
        select case(boundary_condition%flow_condition%itype(idof))
        case(DIRICHLET_BC)
          xxbc(idof) = boundary_condition%flow_aux_real_var(idof,iconn)
          delxbc(idof) = 0.D0
        case(HYDROSTATIC_BC,SEEPAGE_BC)
          xxbc(MPH_PRESSURE_DOF) = boundary_condition%flow_aux_real_var(MPH_PRESSURE_DOF,iconn)
          if (idof >= MPH_TEMPERATURE_DOF) then
            xxbc(idof) = xx_loc_p((ghosted_id-1)*option%nflowdof+idof)
            delxbc(idof) = patch%aux%Immis%delx(idof,ghosted_id)
          endif
        case(NEUMANN_BC, ZERO_GRADIENT_BC)
          ! solve for pb from Darcy's law given qb /= 0
          xxbc(idof) = xx_loc_p((ghosted_id-1)*option%nflowdof+idof)
          delxbc(idof) = patch%aux%Immis%delx(idof,ghosted_id)
        end select
      enddo
    !print *,'BC:',boundary_condition%flow_condition%itype, xxbc, delxbc

 
      call ImmisAuxVarCompute_Ninc(xxbc,auxvars_bc(sum_connection)%auxvar_elem(0),&
         patch%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr,&
         realization%fluid_properties, option)
      call ImmisAuxVarCompute_Winc(xxbc,delxbc,&
         auxvars_bc(sum_connection)%auxvar_elem(1:option%nflowdof),&
         patch%saturation_function_array(int(icap_loc_p(ghosted_id)))%ptr,&
         realization%fluid_properties,option)
    
      do nvar=1,option%nflowdof
        call ImmisBCFlux(boundary_condition%flow_condition%itype, &
          boundary_condition%flow_aux_real_var(:,iconn), &
          auxvars_bc(sum_connection)%auxvar_elem(nvar), &
          auxvars(ghosted_id)%auxvar_elem(nvar), &
          porosity_loc_p(ghosted_id), &
          tortuosity_loc_p(ghosted_id), &
          immis_parameter%sir(:,icap_dn), &
          cur_connection_set%dist(0,iconn),perm_dn,D_dn, &
          cur_connection_set%area(iconn), &
          distance_gravity,option, &
          vv_darcy,Res)
        ResInc(local_id,1:option%nflowdof,nvar) = ResInc(local_id,1:option%nflowdof,nvar) - Res(1:option%nflowdof)
      enddo
    enddo
    boundary_condition => boundary_condition%next
   enddo
#endif
! Set matrix values related to single node terms: Accumulation, Source/Sink, BC
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif

    ra=0.D0
    max_dev=0.D0
    do neq=1, option%nflowdof
      do nvar=1, option%nflowdof
        ra(neq,nvar)=(ResInc(local_id,neq,nvar)-patch%aux%Immis%res_old_AR(local_id,neq))/patch%aux%Immis%delx(nvar,ghosted_id)
        if (max_dev < dabs(ra(3,nvar))) max_dev = dabs(ra(3,nvar))
      enddo
    enddo
   
    select case(option%idt_switch)
      case(1) 
        ra(1:option%nflowdof,1:option%nflowdof) =ra(1:option%nflowdof,1:option%nflowdof) /option%flow_dt
      case(-1)
        if (option%flow_dt>1) ra(1:option%nflowdof,1:option%nflowdof) =ra(1:option%nflowdof,1:option%nflowdof) /option%flow_dt
    end select

    Jup=ra(1:option%nflowdof,1:option%nflowdof)
    if (volume_p(local_id)>1.D0 ) Jup=Jup / volume_p(local_id)
   
     ! if (n==1) print *,  blkmat11, volume_p(n), ra
    call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup,ADD_VALUES, &
                                  ierr);CHKERRQ(ierr)
  end do

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    string = 'jacobian_srcsink'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
#if 1
  ! Interior Flux Terms -----------------------------------  
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0    
  ResInc = 0.D0
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id_up) <= 0 .or. &
            patch%imat(ghosted_id_dn) <= 0) cycle
      endif

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping   
     ! natural_id_up = grid%nG2N(ghosted_id_up)
     ! natural_id_dn = grid%nG2N(ghosted_id_dn)
   
      fraction_upwind = cur_connection_set%dist(-1,iconn)
      distance = cur_connection_set%dist(0,iconn)
      ! distance = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3,iconn) = vector(3) - unit vector
      distance_gravity = distance * &
                         dot_product(option%gravity, &
                                     cur_connection_set%dist(1:3,iconn))
      dd_up = distance*fraction_upwind
      dd_dn = distance-dd_up ! should avoid truncation error
      ! upweight could be calculated as 1.d0-fraction_upwind
      ! however, this introduces ever so slight error causing pflow-overhaul not
      ! to match pflow-orig.  This can be changed to 1.d0-fraction_upwind
      upweight = dd_dn/(dd_up+dd_dn)
    
      ! for now, just assume diagonal tensor
      perm_up = perm_xx_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id_up)*abs(cur_connection_set%dist(3,iconn))

      perm_dn = perm_xx_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(1,iconn))+ &
                perm_yy_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(2,iconn))+ &
                perm_zz_loc_p(ghosted_id_dn)*abs(cur_connection_set%dist(3,iconn))
    
      ithrm_up = int(ithrm_loc_p(ghosted_id_up))
      ithrm_dn = int(ithrm_loc_p(ghosted_id_dn))
      D_up = immis_parameter%ckwet(ithrm_up)
      D_dn = immis_parameter%ckwet(ithrm_dn)
    
      icap_up = int(icap_loc_p(ghosted_id_up))
      icap_dn = int(icap_loc_p(ghosted_id_dn))
      
      do nvar = 1, option%nflowdof 
         call ImmisFlux(auxvars(ghosted_id_up)%auxvar_elem(nvar),porosity_loc_p(ghosted_id_up), &
                          tortuosity_loc_p(ghosted_id_up),immis_parameter%sir(:,icap_up), &
                          dd_up,perm_up,D_up, &
                          auxvars(ghosted_id_dn)%auxvar_elem(0),porosity_loc_p(ghosted_id_dn), &
                          tortuosity_loc_p(ghosted_id_dn),immis_parameter%sir(:,icap_dn), &
                          dd_dn,perm_dn,D_dn, &
                          cur_connection_set%area(iconn),distance_gravity, &
                          upweight, option, vv_darcy, Res)
            ra(:,nvar)= (Res(:)-patch%aux%Immis%res_old_FL(iconn,:))/patch%aux%Immis%delx(nvar,ghosted_id_up)

         call ImmisFlux(auxvars(ghosted_id_up)%auxvar_elem(0),porosity_loc_p(ghosted_id_up), &
                          tortuosity_loc_p(ghosted_id_up),immis_parameter%sir(:,icap_up), &
                          dd_up,perm_up,D_up, &
                          auxvars(ghosted_id_dn)%auxvar_elem(nvar),porosity_loc_p(ghosted_id_dn),&
                          tortuosity_loc_p(ghosted_id_dn),immis_parameter%sir(:,icap_dn), &
                          dd_dn,perm_dn,D_dn, &
                          cur_connection_set%area(iconn),distance_gravity, &
                          upweight, option, vv_darcy, Res)
         ra(:,nvar+option%nflowdof)= (Res(:)-patch%aux%Immis%res_old_FL(iconn,:))/patch%aux%Immis%delx(nvar,ghosted_id_dn)
    enddo

    select case(option%idt_switch)
    case(1)
       ra =ra / option%flow_dt
    case(-1)  
       if (option%flow_dt>1)  ra =ra / option%flow_dt
    end select
    
    if (local_id_up > 0) then
       voltemp=1.D0
       if (volume_p(local_id_up)>1.D0)then
         voltemp = 1.D0/volume_p(local_id_up)
       endif
       Jup(:,1:option%nflowdof)= ra(:,1:option%nflowdof)*voltemp !11
       jdn(:,1:option%nflowdof)= ra(:, 1 + option%nflowdof:2 * option%nflowdof)*voltemp !12

       call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_up-1, &
            Jup,ADD_VALUES,ierr);CHKERRQ(ierr)
       call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_dn-1, &
            Jdn,ADD_VALUES,ierr);CHKERRQ(ierr)
    endif
    if (local_id_dn > 0) then
       voltemp=1.D0
       if (volume_p(local_id_dn)>1.D0)then
         voltemp=1.D0/volume_p(local_id_dn)
       endif
       Jup(:,1:option%nflowdof)= -ra(:,1:option%nflowdof)*voltemp !21
       jdn(:,1:option%nflowdof)= -ra(:, 1 + option%nflowdof:2 * option%nflowdof)*voltemp !22

 
       call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_dn-1, &
            Jdn,ADD_VALUES,ierr);CHKERRQ(ierr)
       call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_up-1, &
            Jup,ADD_VALUES,ierr);CHKERRQ(ierr)
    endif
 enddo
    cur_connection_set => cur_connection_set%next
  enddo
#endif
  if (realization%debug%matview_Jacobian_detailed) then
 ! print *,'end inter flux'
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    string = 'jacobian_flux'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
#if 0
  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call PetscViewerASCIIOpen(option%mycomm,'jacobian_bcflux.out',viewer, &
                              ierr);CHKERRQ(ierr)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
#endif
  
  call VecRestoreArrayF90(field%flow_xx_loc, xx_loc_p, ierr);CHKERRQ(ierr)
!geh refactor  call VecRestoreArrayF90(field%porosity_loc, porosity_loc_p, ierr)
!geh refactor  call VecRestoreArrayF90(field%tortuosity_loc, tortuosity_loc_p, ierr)
!geh refactor  call VecRestoreArrayF90(field%perm_xx_loc, perm_xx_loc_p, ierr)
!geh refactor  call VecRestoreArrayF90(field%perm_yy_loc, perm_yy_loc_p, ierr)
!geh refactor  call VecRestoreArrayF90(field%perm_zz_loc, perm_zz_loc_p, ierr)
!geh refactor  call VecRestoreArrayF90(field%volume, volume_p, ierr)

   
  call VecRestoreArrayF90(field%ithrm_loc, ithrm_loc_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%icap_loc, icap_loc_p, ierr);CHKERRQ(ierr)
! print *,'end jac'
  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
 ! call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
#if 0
! zero out isothermal and inactive cells
#ifdef ISOTHERMAL
  zero = 0.d0
  call MatZeroRowsLocal(A,n_zero_rows,zero_rows_local_ghosted,zero, &
                        PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                        ierr);CHKERRQ(ierr)
  do i=1, n_zero_rows
    ii = mod(zero_rows_local(i),option%nflowdof)
    ip1 = zero_rows_local_ghosted(i)
    if (ii == 0) then
      ip2 = ip1-1
    elseif (ii == option%nflowdof-1) then
      ip2 = ip1+1
    else
      ip2 = ip1
    endif
    call MatSetValuesLocal(A,1,ip1,1,ip2,1.d0,INSERT_VALUES, &
                           ierr);CHKERRQ(ierr)
  enddo

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
#else
#endif
#endif

  if (patch%aux%Immis%inactive_cells_exist) then
    f_up = 1.d0
    call MatZeroRowsLocal(A,patch%aux%Immis%n_zero_rows, &
                          patch%aux%Immis%zero_rows_local_ghosted,f_up, &
                          PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
                          ierr);CHKERRQ(ierr)
  endif

  if (realization%debug%matview_Jacobian) then
    string = 'Ijacobian'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
  if (realization%debug%norm_Jacobian) then
    call MatNorm(A,NORM_1,norm,ierr);CHKERRQ(ierr)
    write(option%io_buffer,'("1 norm: ",es11.4)') norm
    call printMsg(option)
    call MatNorm(A,NORM_FROBENIUS,norm,ierr);CHKERRQ(ierr)
    write(option%io_buffer,'("2 norm: ",es11.4)') norm
    call printMsg(option)
    call MatNorm(A,NORM_INFINITY,norm,ierr);CHKERRQ(ierr)
    write(option%io_buffer,'("inf norm: ",es11.4)') norm
    call printMsg(option)
!    call GridCreateVector(grid,ONEDOF,debug_vec,GLOBAL)
!    call MatGetRowMaxAbs(A,debug_vec,PETSC_NULL_INTEGER,ierr)
!    call VecMax(debug_vec,i,norm,ierr)
!    call VecDestroy(debug_vec,ierr)
  endif
end subroutine ImmisJacobianPatch

! ************************************************************************** !

subroutine ImmisMaxChange(realization,dpmax,dtmpmax,dsmax)
  ! 
  ! Computes the maximum change in the solution vector
  ! 
  ! Author: Chuan Lu
  ! Date: 01/15/08
  ! 

  use Realization_Subsurface_class
  use Field_module
  use Option_module
  use Field_module

  implicit none
  
  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(field_type), pointer :: field
  PetscReal :: dpmax, dtmpmax, dsmax
  PetscErrorCode :: ierr 

  option => realization%option
  field => realization%field

  dpmax=0.D0
  dtmpmax=0.D0 
  dsmax=0.D0

  call VecWAXPY(field%flow_dxx,-1.d0,field%flow_xx,field%flow_yy, &
                ierr);CHKERRQ(ierr)
  call VecStrideNorm(field%flow_dxx,ZERO_INTEGER,NORM_INFINITY,dpmax, &
                     ierr);CHKERRQ(ierr)
  call VecStrideNorm(field%flow_dxx,ONE_INTEGER,NORM_INFINITY,dtmpmax, &
                     ierr);CHKERRQ(ierr)
  call VecStrideNorm(field%flow_dxx,TWO_INTEGER,NORM_INFINITY,dsmax, &
                     ierr);CHKERRQ(ierr)

  !print *, 'Max changes=', option%dpmax,option%dtmpmax, option%dcmax,option%dsmax
end subroutine ImmisMaxChange

! ************************************************************************** !

function ImmisGetTecplotHeader(realization, icolumn)
  ! 
  ! Returns Richards contribution to
  ! Tecplot file header
  ! 
  ! Author: Chuan Lu
  ! Date: 10/13/08
  ! 

  use Realization_Subsurface_class
  use Option_module
  use Field_module

  implicit none
  
  character(len=MAXSTRINGLENGTH) :: ImmisGetTecplotHeader
  type(realization_subsurface_type) :: realization
  PetscInt :: icolumn
  
  character(len=MAXSTRINGLENGTH) :: string, string2
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  PetscInt :: i
  
  option => realization%option
  field => realization%field
  
  string = ''

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-T [C]"'')') icolumn
  else
    write(string2,'('',"T [C]"'')')
  endif
  string = trim(string) // trim(string2)
  
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-P [Pa]"'')') icolumn
  else
    write(string2,'('',"P [Pa]"'')')
  endif
  string = trim(string) // trim(string2)
  
! if (icolumn > -1) then
!   icolumn = icolumn + 1
!   write(string2,'('',"'',i2,''-PHASE"'')') icolumn
! else
!   write(string2,'('',"PHASE"'')')
! endif
! string = trim(string) // trim(string2)
  
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-S(l)"'')') icolumn
  else
    write(string2,'('',"S(l)"'')')
  endif
  string = trim(string) // trim(string2)

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-S(g)"'')') icolumn
  else
    write(string2,'('',"S(g)"'')')
  endif
  string = trim(string) // trim(string2)
    
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-D(l)"'')') icolumn
  else
    write(string2,'('',"D(l)"'')')
  endif
  string = trim(string) // trim(string2)

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-D(g)"'')') icolumn
  else
    write(string2,'('',"D(g)"'')')
  endif
  string = trim(string) // trim(string2)
    
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-U(l)"'')') icolumn
  else
    write(string2,'('',"U(l)"'')')
  endif
  string = trim(string) // trim(string2)

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-U(g)"'')') icolumn
  else
    write(string2,'('',"U(g)"'')')
  endif
  string = trim(string) // trim(string2)
    
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-Vis(l)"'')') icolumn
  else
    write(string2,'('',"Vis(l)"'')')
  endif
  string = trim(string) // trim(string2)

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-Vis(g)"'')') icolumn
  else
    write(string2,'('',"Vis(g)"'')')
  endif
  string = trim(string) // trim(string2)
    
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-Mob(l)"'')') icolumn
  else
    write(string2,'('',"Mob(l)"'')')
  endif
  string = trim(string) // trim(string2)

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-Mob(g)"'')') icolumn
  else
    write(string2,'('',"Mob(g)"'')')
  endif
  string = trim(string) // trim(string2)
  
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-PHASE"'')') icolumn
  else
    write(string2,'('',"PHASE"'')')
  endif
  string = trim(string) // trim(string2)

  ImmisGetTecplotHeader = string

end function ImmisGetTecplotHeader

! ************************************************************************** !

subroutine ImmisSetPlotVariables(list)
  ! 
  ! Adds variables to be printed to list
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/12
  ! 
  
  use Output_Aux_module
  use Variables_module

  implicit none

  type(output_variable_list_type), pointer :: list

  type(output_variable_type) :: output_variable
  character(len=MAXWORDLENGTH) :: name, units
  
  if (associated(list%first)) then
    return
  endif

  name = 'Temperature'
  units = 'C'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               TEMPERATURE)
  
  name = 'Liquid Pressure'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_PRESSURE,units, &
                               LIQUID_PRESSURE)

  name = 'Gas Pressure'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_PRESSURE,units, &
                               GAS_PRESSURE)

  name = 'Liquid Saturation'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
                               LIQUID_SATURATION)

  name = 'Gas Saturation'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
                               GAS_SATURATION)

  name = 'Liquid Density'
  units = 'kg/m^3'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               LIQUID_DENSITY)

  name = 'Gas Density'
  units = 'kg/m^3'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               GAS_DENSITY)

  name = 'Liquid Energy'
  units = 'kJ/mol'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               LIQUID_ENERGY)

  name = 'Gas Energy'
  units = 'kJ/mol'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               GAS_ENERGY)

  name = 'Liquid Viscosity'
  units = 'Pa.s'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               LIQUID_VISCOSITY)

  name = 'Gas Viscosity'
  units = 'Pa.s'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               GAS_VISCOSITY)

  name = 'Liquid Mobility'
  units = '1/Pa.s'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               LIQUID_MOBILITY)

  name = 'Gas Mobility'
  units = '1/Pa.s'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               GAS_MOBILITY)

! name = 'Phase'
! units = ''
! output_variable%iformat = 1 ! integer
! call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
!                              PHASE)

end subroutine ImmisSetPlotVariables

! ************************************************************************** !

subroutine ImmisDestroy(realization)
  ! 
  ! Deallocates variables associated with Immis
  ! 
  ! Author: Chuan Lu
  ! Date: 10/14/08
  ! 

  use Realization_Subsurface_class

  implicit none
  
  type(realization_subsurface_type) :: realization
  
  ! need to free array in aux vars
  !call ImmisAuxDestroy(patch%aux%Immis)

end subroutine ImmisDestroy


#if 0

! ************************************************************************** !

subroutine ImmisCheckpointWrite(discretization, viewer)
  ! 
  ! Writes vecs to checkpoint file
  ! Author: Chuan Lu
  ! 

  use Discretization_module

  implicit none
  
  type(discretization_type) :: discretization
  PetscViewer :: viewer
  
  Vec :: global_var
  PetscErrorCode :: ierr
  
  call VecView(global_var,viewer,ierr);CHKERRQ(ierr)
  call VecDestroy(global_var,ierr);CHKERRQ(ierr)
  
  
end subroutine ImmisCheckpointWrite

! ************************************************************************** !

subroutine ImmisCheckpointRead(discretization,viewer)
  ! 
  ! Reads vecs from checkpoint file
  ! Author: Chuan Lu
  ! 

  use Discretization_module

  implicit none
  
  type(discretization_type) :: discretization
  PetscViewer :: viewer
  
  Vec :: global_var
  PetscErrorCode :: ierr
  
  call VecLoad(global_var, viewer, ierr);CHKERRQ(ierr)
  call VecDestroy(global_var,ierr);CHKERRQ(ierr)
  
end subroutine ImmisCheckpointRead

#endif

end module Immis_module
