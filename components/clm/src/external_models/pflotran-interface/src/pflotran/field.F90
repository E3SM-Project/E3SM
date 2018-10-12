module Field_module

! IMPORTANT NOTE: This module can have no dependencies on other modules!!!
 
#include "petsc/finclude/petscvec.h"
  use petscvec
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: field_type 
    
    !get material id
    ! 1 degree of freedom
    Vec :: porosity0
    Vec :: porosity_base_store
    Vec :: porosity_t
    Vec :: porosity_tpdt
    Vec :: porosity_geomech_store
    Vec :: tortuosity0
    Vec :: ithrm_loc
    Vec :: icap_loc
    Vec :: iphas_loc, iphas_old_loc

    Vec :: perm0_xx, perm0_yy, perm0_zz
    !geh: required for higher order, but not supported at this time.
!    Vec :: perm0_xz, perm0_xy, perm0_yz
    
    Vec :: work, work_loc

    Vec :: volume0
    Vec :: compressibility0
    
    ! residual vectors
    Vec :: flow_r          
    Vec :: tran_r
    
    ! Solution vectors (yy = previous solution, xx = current iterate)
    Vec :: flow_xx, flow_xx_loc, flow_dxx, flow_yy, flow_accum, flow_accum2
    Vec :: tran_xx, tran_xx_loc, tran_dxx, tran_yy, tran_accum

    ! vectors for operator splitting
    Vec :: tran_rhs
    Vec :: tran_rhs_coef
    
    Vec :: tran_log_xx, tran_work_loc
    
    ! mass transfer
    Vec :: flow_mass_transfer
    Vec :: tran_mass_transfer
    
    Vec :: flow_ts_mass_balance, flow_total_mass_balance
    Vec :: tran_ts_mass_balance, tran_total_mass_balance

    ! vector that holds the second layer of ghost cells for tvd
    Vec :: tvd_ghosts

    ! vectors to save temporally average quantities
    Vec, pointer :: avg_vars_vec(:)
    PetscInt :: nvars

    ! vectors to save temporally average flowrates
    Vec :: flowrate_inst
    Vec :: flowrate_aveg

    ! vectors to save velocity at face
    Vec :: vx_face_inst
    Vec :: vy_face_inst
    Vec :: vz_face_inst

    Vec, pointer :: max_change_vecs(:)

  end type field_type

  public :: FieldCreate, &
            FieldDestroy

contains

! ************************************************************************** !

function FieldCreate()
  ! 
  ! Allocates and initializes a new Field object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  ! 

  implicit none
  
  type(field_type), pointer :: FieldCreate
  
  type(field_type), pointer :: field
  
  allocate(field)
  
  ! nullify PetscVecs
  field%porosity0 = PETSC_NULL_VEC
  field%porosity_base_store = PETSC_NULL_VEC
  field%porosity_geomech_store = PETSC_NULL_VEC
  field%porosity_t = PETSC_NULL_VEC
  field%porosity_tpdt = PETSC_NULL_VEC
  field%tortuosity0 = PETSC_NULL_VEC
  field%ithrm_loc = PETSC_NULL_VEC
  field%icap_loc = PETSC_NULL_VEC
  field%iphas_loc = PETSC_NULL_VEC
  field%iphas_old_loc = PETSC_NULL_VEC

  field%perm0_xx = PETSC_NULL_VEC
  field%perm0_yy = PETSC_NULL_VEC
  field%perm0_zz = PETSC_NULL_VEC
  
  field%work = PETSC_NULL_VEC
  field%work_loc = PETSC_NULL_VEC

  field%volume0 = PETSC_NULL_VEC
  field%compressibility0 = PETSC_NULL_VEC
  
  field%flow_r = PETSC_NULL_VEC
  field%flow_xx = PETSC_NULL_VEC
  field%flow_xx_loc = PETSC_NULL_VEC
  field%flow_dxx = PETSC_NULL_VEC
  field%flow_yy = PETSC_NULL_VEC
  field%flow_accum = PETSC_NULL_VEC
  field%flow_accum2 = PETSC_NULL_VEC

  field%tran_r = PETSC_NULL_VEC
  field%tran_log_xx = PETSC_NULL_VEC
  field%tran_xx = PETSC_NULL_VEC
  field%tran_xx_loc = PETSC_NULL_VEC
  field%tran_dxx = PETSC_NULL_VEC
  field%tran_yy = PETSC_NULL_VEC
  field%tran_accum = PETSC_NULL_VEC
  field%tran_work_loc = PETSC_NULL_VEC

  field%tvd_ghosts = PETSC_NULL_VEC

  field%tran_rhs = PETSC_NULL_VEC
  field%tran_rhs_coef = PETSC_NULL_VEC
  
  field%flow_mass_transfer = PETSC_NULL_VEC
  field%tran_mass_transfer = PETSC_NULL_VEC
  
  field%flow_ts_mass_balance = PETSC_NULL_VEC
  field%flow_total_mass_balance = PETSC_NULL_VEC
  field%tran_ts_mass_balance = PETSC_NULL_VEC
  field%tran_total_mass_balance = PETSC_NULL_VEC

  nullify(field%avg_vars_vec)
  field%nvars = 0

  field%flowrate_inst = PETSC_NULL_VEC
  field%flowrate_aveg = PETSC_NULL_VEC

  field%vx_face_inst = PETSC_NULL_VEC
  field%vy_face_inst = PETSC_NULL_VEC
  field%vz_face_inst = PETSC_NULL_VEC

  nullify(field%max_change_vecs)

  FieldCreate => field
  
end function FieldCreate

! ************************************************************************** !

subroutine FieldDestroy(field)
  ! 
  ! Deallocates a field object
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/15/07
  ! 

  implicit none
  
  type(field_type), pointer :: field
  
  PetscErrorCode :: ierr
  PetscInt :: ivar
  PetscInt :: num_vecs

  if (.not.associated(field)) return

  ! Destroy PetscVecs
  if (field%porosity0 /= PETSC_NULL_VEC) then
    call VecDestroy(field%porosity0,ierr);CHKERRQ(ierr)
  endif
  if (field%porosity_base_store /= PETSC_NULL_VEC) then
    call VecDestroy(field%porosity_base_store,ierr);CHKERRQ(ierr)
  endif
  if (field%porosity_geomech_store /= PETSC_NULL_VEC) then
    call VecDestroy(field%porosity_geomech_store,ierr);CHKERRQ(ierr)
  endif
  if (field%porosity_t /= PETSC_NULL_VEC) then
    call VecDestroy(field%porosity_t,ierr);CHKERRQ(ierr)
  endif
  if (field%porosity_tpdt /= PETSC_NULL_VEC) then
    call VecDestroy(field%porosity_tpdt,ierr);CHKERRQ(ierr)
  endif
  if (field%tortuosity0 /= PETSC_NULL_VEC) then
    call VecDestroy(field%tortuosity0,ierr);CHKERRQ(ierr)
  endif
  if (field%ithrm_loc /= PETSC_NULL_VEC) then
    call VecDestroy(field%ithrm_loc,ierr);CHKERRQ(ierr)
  endif
  if (field%icap_loc /= PETSC_NULL_VEC) then
    call VecDestroy(field%icap_loc,ierr);CHKERRQ(ierr)
  endif
  if (field%iphas_loc /= PETSC_NULL_VEC) then
    call VecDestroy(field%iphas_loc,ierr);CHKERRQ(ierr)
  endif
  if (field%iphas_old_loc /= PETSC_NULL_VEC) then
    call VecDestroy(field%iphas_old_loc,ierr);CHKERRQ(ierr)
  endif

  if (field%perm0_xx /= PETSC_NULL_VEC) then
    call VecDestroy(field%perm0_xx,ierr);CHKERRQ(ierr)
  endif
  if (field%perm0_yy /= PETSC_NULL_VEC) then
    call VecDestroy(field%perm0_yy,ierr);CHKERRQ(ierr)
  endif
  if (field%perm0_zz /= PETSC_NULL_VEC) then
    call VecDestroy(field%perm0_zz,ierr);CHKERRQ(ierr)
  endif
  
  if (field%work /= PETSC_NULL_VEC) then
    call VecDestroy(field%work,ierr);CHKERRQ(ierr)
  endif
  if (field%work_loc /= PETSC_NULL_VEC) then
    call VecDestroy(field%work_loc,ierr);CHKERRQ(ierr)
  endif

  if (field%volume0 /= PETSC_NULL_VEC) then
    call VecDestroy(field%volume0,ierr);CHKERRQ(ierr)
  endif

  if (field%compressibility0 /= PETSC_NULL_VEC) then
    call VecDestroy(field%compressibility0,ierr);CHKERRQ(ierr)
  endif
 
  if (field%flow_r /= PETSC_NULL_VEC) then
    call VecDestroy(field%flow_r,ierr);CHKERRQ(ierr)
  endif
  if (field%flow_xx /= PETSC_NULL_VEC) then
    call VecDestroy(field%flow_xx,ierr);CHKERRQ(ierr)
  endif
  if (field%flow_xx_loc /= PETSC_NULL_VEC) then
    call VecDestroy(field%flow_xx_loc,ierr);CHKERRQ(ierr)
  endif
  if (field%flow_dxx /= PETSC_NULL_VEC) then
    call VecDestroy(field%flow_dxx,ierr);CHKERRQ(ierr)
  endif
  if (field%flow_yy /= PETSC_NULL_VEC) then
    call VecDestroy(field%flow_yy,ierr);CHKERRQ(ierr)
  endif
  if (field%flow_accum /= PETSC_NULL_VEC) then
    call VecDestroy(field%flow_accum,ierr);CHKERRQ(ierr)
  endif
  if (field%flow_accum2 /= PETSC_NULL_VEC) then
    call VecDestroy(field%flow_accum2,ierr);CHKERRQ(ierr)
  endif
  
  if (field%tran_r /= PETSC_NULL_VEC) then
    call VecDestroy(field%tran_r,ierr);CHKERRQ(ierr)
  endif
  if (field%tran_log_xx /= PETSC_NULL_VEC) then
    call VecDestroy(field%tran_log_xx,ierr);CHKERRQ(ierr)
  endif
  if (field%tran_xx /= PETSC_NULL_VEC) then
    call VecDestroy(field%tran_xx,ierr);CHKERRQ(ierr)
  endif
  if (field%tran_xx_loc /= PETSC_NULL_VEC) then
    call VecDestroy(field%tran_xx_loc,ierr);CHKERRQ(ierr)
  endif
  if (field%tran_dxx /= PETSC_NULL_VEC) then
    call VecDestroy(field%tran_dxx,ierr);CHKERRQ(ierr)
  endif
  if (field%tran_yy /= PETSC_NULL_VEC) then
    call VecDestroy(field%tran_yy,ierr);CHKERRQ(ierr)
  endif
  if (field%tran_accum /= PETSC_NULL_VEC) then
    call VecDestroy(field%tran_accum,ierr);CHKERRQ(ierr)
  endif
  if (field%tran_work_loc /= PETSC_NULL_VEC) then
    call VecDestroy(field%tran_work_loc,ierr);CHKERRQ(ierr)
  endif
  
  if (field%tran_rhs /= PETSC_NULL_VEC) then
    call VecDestroy(field%tran_rhs,ierr);CHKERRQ(ierr)
  endif
  if (field%tran_rhs_coef /= PETSC_NULL_VEC) then
    call VecDestroy(field%tran_rhs_coef,ierr);CHKERRQ(ierr)
  endif

  if (field%flow_mass_transfer /= PETSC_NULL_VEC) then
    call VecDestroy(field%flow_mass_transfer,ierr);CHKERRQ(ierr)
  endif
  if (field%tran_mass_transfer /= PETSC_NULL_VEC) then
    call VecDestroy(field%tran_mass_transfer,ierr);CHKERRQ(ierr)
  endif

  if (field%flow_ts_mass_balance /= PETSC_NULL_VEC) then
    call VecDestroy(field%flow_ts_mass_balance,ierr);CHKERRQ(ierr)
  endif
  if (field%flow_total_mass_balance /= PETSC_NULL_VEC) then
    call VecDestroy(field%flow_total_mass_balance,ierr);CHKERRQ(ierr)
  endif
  if (field%tran_ts_mass_balance /= PETSC_NULL_VEC) then
    call VecDestroy(field%tran_ts_mass_balance,ierr);CHKERRQ(ierr)
  endif
  if (field%tran_total_mass_balance /= PETSC_NULL_VEC) then
    call VecDestroy(field%tran_total_mass_balance,ierr);CHKERRQ(ierr)
  endif
    
  if (field%tvd_ghosts /= PETSC_NULL_VEC) then
    call VecDestroy(field%tvd_ghosts,ierr);CHKERRQ(ierr)
  endif

  do ivar = 1,field%nvars
    call VecDestroy(field%avg_vars_vec(ivar),ierr);CHKERRQ(ierr)
  enddo

  if (field%flowrate_inst /= PETSC_NULL_VEC) then
    call VecDestroy(field%flowrate_inst,ierr);CHKERRQ(ierr)
  endif
  if (field%flowrate_aveg /= PETSC_NULL_VEC) then
    call VecDestroy(field%flowrate_aveg,ierr);CHKERRQ(ierr)
  endif

  if (field%vx_face_inst /= PETSC_NULL_VEC) then
    call VecDestroy(field%vx_face_inst,ierr);CHKERRQ(ierr)
  endif
  if (field%vy_face_inst /= PETSC_NULL_VEC) then
    call VecDestroy(field%vy_face_inst,ierr);CHKERRQ(ierr)
  endif
  if (field%vz_face_inst /= PETSC_NULL_VEC) then
    call VecDestroy(field%vz_face_inst,ierr);CHKERRQ(ierr)
  endif

  if (associated(field%max_change_vecs)) then
    !geh: kludge as the compiler returns i4 in 64-bit
    num_vecs = size(field%max_change_vecs)
    call VecDestroyVecsF90(num_vecs,field%max_change_vecs,ierr);CHKERRQ(ierr)
  endif

  deallocate(field)
  nullify(field)
  
end subroutine FieldDestroy

end module Field_module
