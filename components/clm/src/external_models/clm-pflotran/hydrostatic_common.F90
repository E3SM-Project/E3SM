module Hydrostatic_Common_module
 
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  ! commented out, use those defined in pflotran_constants
  !LIQ = water to be consistent with the remainder fo the code
  !PetscInt, parameter, public :: HYDRO_LIQ_PHASE = 1  
  !PetscInt, parameter, public :: HYDRO_GAS_PHASE = 2
  !PetscInt, parameter, public :: HYDRO_OIL_PHASE = 3 


  type, public :: one_dim_grid_type
    PetscReal :: delta_z
    PetscReal :: min_z
    PetscReal :: max_z
    PetscReal, pointer :: z(:)
    PetscInt :: idatum
  contains 
    procedure :: ElevationIdLoc
    procedure, public :: InterpFromWellConnTo1DGrid
  end type one_dim_grid_type


  public :: GetCellOnPhaseConact, &
            PhaseHydrostaticPressure, &
            PhaseDensity, &
            PressInterp, &
            PressGrad, &  
            CompVertTempProfile, &
            CreateOneDimGrid, &
            DestroyOneDimGrid
                    
 
contains

! ************************************************************************** !

function GetCellOnPhaseConact(connection_set,grid,imat,z_phase_contact,option)

  ! Returns the ghosted_id of the closest cell to phase contact  
  ! Where more cells have the same minimum distance from z_phase_contact,
  ! the first cell in the list will be selected as for the minloc function 
  !
  ! Author: Paolo Orsini
  ! Date: 12/31/15
  ! 

  use Grid_module
  !use Coupler_module
  use Connection_module
  use Utility_module
  use Option_module

  implicit none

  !type(coupler_type), intent(in) :: coupler
  type(connection_set_type), intent(in) :: connection_set
  type(grid_type), intent(in) :: grid
  PetscInt, pointer, intent(in) :: imat(:)
  PetscReal, intent(in) :: z_phase_contact
  type(option_type) :: option

  PetscInt :: GetCellOnPhaseConact

  PetscInt :: iconn, local_id, ghosted_id
  PetscInt :: iconn_min(1)
  PetscReal, pointer :: phase_contact_dist(:) 

  allocate(phase_contact_dist(connection_set%num_connections))
  phase_contact_dist = 1.d20
  
  do iconn=1,connection_set%num_connections
    local_id = connection_set%id_dn(iconn)
    ghosted_id = grid%nL2G(local_id)
    if (imat(ghosted_id) <= 0) cycle
    phase_contact_dist(iconn) = dabs(grid%z(ghosted_id) - z_phase_contact)
  end do

  iconn_min = minloc(phase_contact_dist(:))

  local_id = connection_set%id_dn(iconn_min(1))

  !write(*,*) "in get ghost -rank=", option%myrank, "icon_min=", iconn_min(1), &
  !           " ghost =", grid%nL2G(local_id), "size dist = ", size(phase_contact_dist(:))

  GetCellOnPhaseConact = grid%nL2G(local_id) 

  call DeallocateArray(phase_contact_dist)
 

end function GetCellOnPhaseConact

! ************************************************************************** !

subroutine PhaseHydrostaticPressure(one_d_grid,gravity,iphase,press_start, &
                                    id_start,xm_nacl,temp,press,den_kg)
  ! 
  ! Compute an hydrostatic pressure profile for a given phase and 1D 
  ! 1D discretisation
  ! The "starting point" for the computation of the pressure profile
  ! is the elevation where a reference pressure is given, it can be the datum,
  ! or one of the phase contact interface (e.g. OWC, OGC, etc)
  !
  ! Author: Paolo Orsini
  ! Date: 12/21/15
  ! 

  implicit none

  class(one_dim_grid_type) :: one_d_grid
  PetscReal, intent(in) :: gravity(:) ! this is option%gravity
  PetscInt, intent(in) :: iphase
  PetscReal, intent(in) :: press_start
  PetscInt, intent(in) :: id_start
  PetscReal, intent(in) :: temp(:)
  PetscReal, intent(in) :: xm_nacl
  PetscReal, intent(out) :: press(:)
  PetscReal, intent(out) :: den_kg(:)

  PetscReal :: pressure, pressure0, rho, rho_one, rho_kg, rho_zero
  PetscInt :: ipressure, num_iteration
  
  if(iphase == GAS_PHASE ) then
    print *, "PhaseHydrostaticPressure does not support gas"
    stop
  end if

  rho_kg = PhaseDensity(iphase,press_start,temp(id_start),xm_nacl)

  ! fill properties for reference pressure
  den_kg(id_start) = rho_kg
  press(id_start) = press_start 

  pressure0 = press_start 
  rho_zero = rho_kg 
  do ipressure=id_start+1,size(one_d_grid%z(:))
    rho_kg = PhaseDensity(iphase,pressure0,temp(ipressure),xm_nacl)
    num_iteration = 0
    do 
      pressure = pressure0 + 0.5d0*(rho_kg+rho_zero) * &
                 gravity(Z_DIRECTION) * one_d_grid%delta_z
      rho_one = PhaseDensity(iphase,pressure,temp(ipressure),xm_nacl)
      !check convergence on density
      if (dabs(rho_kg-rho_one) < 1.d-10) exit
      rho_kg = rho_one
      num_iteration = num_iteration + 1
      if (num_iteration > 100) then
        print *,'Phase-Hydrostatic iteration failed to converge', &
                 num_iteration,rho_one,rho_kg
        !print *, condition%name, idatum
        !print *, pressure_array
        stop
      endif
    enddo
    rho_zero = rho_kg
    press(ipressure) = pressure
    den_kg(ipressure) = rho_kg
    pressure0 = pressure
  enddo

  ! compute pressures below one_d_grid%z(id_start), if any
  pressure0 = press(id_start)
  rho_zero = den_kg(id_start)
  do ipressure=id_start-1,1,-1
    rho_kg = PhaseDensity(iphase,pressure0,temp(ipressure),xm_nacl)
    num_iteration = 0
    do                   ! notice the negative sign (-) here
      pressure = pressure0 - 0.5d0*(rho_kg+rho_zero) * &
                 gravity(Z_DIRECTION) * one_d_grid%delta_z
      rho_one = PhaseDensity(iphase,pressure,temp(ipressure),xm_nacl)
      !check convergence on density
      if (dabs(rho_kg-rho_one) < 1.d-10) exit
      rho_kg = rho_one
      num_iteration = num_iteration + 1
      if (num_iteration > 100) then
        print *,'Phase-Hydrostatic iteration failed to converge', &
                 num_iteration,rho_one,rho_kg
        !print *, condition%name, idatum
        !print *, pressure_array
        stop
      endif
    enddo
    rho_zero = rho_kg
    press(ipressure) = pressure
    den_kg(ipressure) = rho_kg
    pressure0 = pressure
  enddo

  !write(*,*) "den fine bottom", den_kg(1) 
  !write(*,*) "den fine top", den_kg(size(one_d_grid%z(:)))

end subroutine PhaseHydrostaticPressure

! ************************************************************************** !
function PhaseDensity(iphase,p,t,xm_nacl) 
  !
  ! computes phase density given the specified phase
  !
  ! Author: Paolo Orsini
  ! Date: 12/21/15
  ! 

  use EOS_Oil_module
  use EOS_Water_module
  use EOS_Gas_module ! when gas is considered

  implicit none

  PetscInt, intent(in) :: iphase
  PetscReal, intent(in) :: p
  PetscReal, intent(in) :: t
  PetscReal, intent(in) :: xm_nacl
 
  PetscReal :: PhaseDensity ! kg/m3 
  PetscReal :: dw_mol
  PetscReal :: aux(1)

  PetscErrorCode :: ierr

  select case(iphase)
    case(LIQUID_PHASE)
      aux(1) = xm_nacl
      call EOSWaterDensityExt(t,p,aux,PhaseDensity,dw_mol,ierr)
    case(GAS_PHASE)
      !call EOSGasDensityNoDerive(t,p,PhaseDensity,ierr)
      ! rho_kg = rho * GAS_FMW (to get gas FMW currenlty mode specific)
      ! gas_fmw should be defined in gas_eos 
    case(OIL_PHASE)
      call EOSOilDensity(t,p,PhaseDensity,ierr)
      PhaseDensity = PhaseDensity * EOSOilGetFMW() 
  end select

end function PhaseDensity 

! ************************************************************************** !

function PressInterp(ipressure,dist_x,dist_y,dist_z_for_pressure,gravity, &
                     pressure_array,density_array,pressure_gradient)

  ! Computes hydrostatic pressure over dz from a reference pressure and
  ! a density, combining possible horizontal gradients  
  !
  ! Author: Paolo Orsini
  ! Date: 12/21/15


  implicit none

  PetscInt, intent(in) :: ipressure
  PetscReal, intent(in) :: dist_x
  PetscReal, intent(in) :: dist_y
  PetscReal, intent(in) :: dist_z_for_pressure
  PetscReal, intent(in) :: gravity(:)
  PetscReal, intent(in) :: pressure_array(:)
  PetscReal, intent(in) :: density_array(:)
  PetscReal, intent(in) :: pressure_gradient(:)

  PetscReal :: PressInterp

  PressInterp = pressure_array(ipressure) + &
                density_array(ipressure) * gravity(Z_DIRECTION) * &
                dist_z_for_pressure + &
                pressure_gradient(X_DIRECTION) * dist_x + & ! gradient in Pa/m
                pressure_gradient(Y_DIRECTION) * dist_y

end function PressInterp

! ************************************************************************** !

!function HydrostaticPressOverDz(dz,gravity,press_ref,density_ref)
!
!  implicit none
!
!  PetscReal, intent(in) :: dz
!  PetscReal, intent(in) :: gravity(:)
!  PetscReal, intent(in) :: press_ref
!  PetscReal, intent(in) :: density_ref
!
!  PetscReal :: HydrostaticPressOverDz
!
!  HydrostaticPressOverDz = press_ref + density_ref * gravity(Z_DIRECTION) * dz
!
!end function HydrostaticPressOverDz

! ************************************************************************** !

function PressGrad(dist_x,dist_y,dist_z,press_ref,pressure_gradient)

  ! Computes pressure from a reference pressure, distances in the three 
  ! directions, and a 3d local pressure gradient  
  !
  ! Author: Paolo Orsini
  ! Date: 12/21/15


  implicit none

  PetscReal, intent(in) :: dist_x
  PetscReal, intent(in) :: dist_y
  PetscReal, intent(in) :: dist_z
  PetscReal, intent(in) :: press_ref
  PetscReal, intent(in) :: pressure_gradient(:)

  PetscReal :: PressGrad

  PressGrad = press_ref + &
              pressure_gradient(X_DIRECTION)*dist_x + & ! gradient in Pa/m
              pressure_gradient(Y_DIRECTION)*dist_y + &
              pressure_gradient(Z_DIRECTION)*dist_z 

end function PressGrad

! ************************************************************************** !


subroutine CompVertTempProfile(one_d_grid,temp_grad,temp_at_datum, &
                               temp_profile)
  ! 
  ! Computes the temperature vertical profile on the 1D grid use for the 
  ! hydrostatic pressure calculation
  !
  ! Author: Paolo Orsini
  ! Date: 12/21/15
  ! 

  implicit none

  class(one_dim_grid_type) :: one_d_grid
  PetscReal, intent(in) :: temp_grad(:)
  PetscReal, intent(in) :: temp_at_datum
  PetscReal, intent(out) :: temp_profile(:)
 
  PetscInt :: i_z

  temp_profile(one_d_grid%idatum) = temp_at_datum

  do i_z=one_d_grid%idatum+1,size(one_d_grid%z(:))
    temp_profile(i_z) = temp_profile(i_z-1) + &
                        temp_grad(Z_DIRECTION)*one_d_grid%delta_z   
  end do

  do i_z=one_d_grid%idatum-1,1,-1
    temp_profile(i_z) = temp_profile(i_z+1) - & ! note the (-) sign
                        temp_grad(Z_DIRECTION)*one_d_grid%delta_z   
  end do

end subroutine CompVertTempProfile
! ************************************************************************** !

function CreateOneDimGrid(min_z,max_z,datum)
  ! 
  ! Computes 1D grid for interpolation needed in hydrostatic pressure
  ! computation
  !
  ! Author: Paolo Orsini
  ! Date: 12/21/15
  ! 

  implicit none

  class(one_dim_grid_type), pointer :: one_d_grid
  PetscReal,intent(in) :: datum(3)
  PetscReal, intent(in) :: min_z, max_z

  class(one_dim_grid_type), pointer :: CreateOneDimGrid

  PetscInt :: num_z, i_z
  PetscReal :: dist_z

  allocate(one_d_grid)

  one_d_grid%min_z = min_z
  one_d_grid%max_z = max_z

  one_d_grid%delta_z = min((max_z-min_z)/500.d0,1.d0)
  ! if zero, assign 1.d0 to avoid divide by zero below. essentially the grid
  ! is flat.
  if (one_d_grid%delta_z < 1.d-40) one_d_grid%delta_z = 1.d0

  num_z = int((max_z-min_z)/one_d_grid%delta_z) + 1

  allocate(one_d_grid%z(num_z))

  one_d_grid%idatum = int((datum(Z_DIRECTION)-min_z)/(max_z-min_z) * &
                        dble(num_z))+1  

  one_d_grid%z(one_d_grid%idatum) = datum(Z_DIRECTION)  

  dist_z = 0.d0
  do i_z=one_d_grid%idatum+1,num_z
    dist_z = dist_z + one_d_grid%delta_z
    one_d_grid%z(i_z) = one_d_grid%z(one_d_grid%idatum) + dist_z
  enddo

  dist_z = 0.d0
  do i_z = one_d_grid%idatum-1,1,-1
    dist_z = dist_z + one_d_grid%delta_z
    one_d_grid%z(i_z) = one_d_grid%z(one_d_grid%idatum) - dist_z
  enddo

  CreateOneDimGrid => one_d_grid

end function CreateOneDimGrid

! ************************************************************************** !
function ElevationIdLoc(this,elevation)
  ! 
  ! Detect id location in one_d_grid given an absolute elevation
  ! Note: absolute elevation, not relative to the datum 
  !
  ! Author: Paolo Orsini
  ! Date: 12/21/15
  ! 

  implicit none

  class(one_dim_grid_type) :: this
  PetscReal, intent(in) :: elevation !this is the global elevation
  
  PetscInt :: ElevationIdLoc
   
  ElevationIdLoc = int( (elevation-this%min_z)/(this%max_z-this%min_z) * &
                        dble(size(this%z(:))) ) + 1 

  !  idatum = int((datum(Z_DIRECTION)-min_z)/(max_z-min_z) * &
  !               dble(num_pressures))+1

end function ElevationIdLoc

! ************************************************************************** !

!subroutine InterpFromWellConnTo1DGrid(this,grid,connection_set, w_conn_z &
!                                      well_conn_val,interp_val)
subroutine InterpFromWellConnTo1DGrid(this,w_conn_z,ord, &
                                      well_conn_val,interp_val)

  !   
  ! Interpolate well connection varìaiable values to the vertical fine grid
  ! Note that the vertical fine grid is ordered for ascending z
  ! This function works for deviated wells -
  ! Needs to be extended for horizontal wells  
  !
  ! Author: Paolo Orsini
  ! Date: 06/17/16
  ! 

  implicit none

  !class(one_dim_grid_type) :: one_d_grid
  class(one_dim_grid_type) :: this
  !type(connection_set_type), intent(in) :: connection_set
  !type(grid_type), intent(in) :: grid
  PetscReal, intent(in) :: w_conn_z(:)
  PetscReal, intent(in) :: well_conn_val(:) ! same size of w_conn_z but not z-ordered
  PetscInt, intent(in) :: ord(:)           ! same size of w_conn_z -> well_num_conns
  PetscReal, intent(out) :: interp_val(:)  ! same size of this%z(:)

  PetscInt :: inode,iconn  
  PetscReal, parameter :: tol_z = 1.0d-6 ! m - to spot nearly overlapping nodes

  !well_conn_val(ord(iconn)) -> gives conn_val in z ascending order 

  !for horizontal and mixed deviated/horizontal wells
  ! - loop over w_conn_z and identify aligned connections using tol_z
  ! - for aligned well connections perform an average of well_conn_val
  !   over the number of aligned connections 
  

  do inode=1,size(this%z(:))
    if ( this%z(inode) < (w_conn_z(ONE_INTEGER)-tol_z) ) then
      interp_val(inode) = well_conn_val(ord(ONE_INTEGER)) 
    else if ( this%z(inode) > (w_conn_z(size(ord(:)))+tol_z) ) then
      interp_val(inode) = well_conn_val(ord(size(ord(:))))
    else 
      do iconn=1,size(ord(:))-1 !loop to well_num_conns-1
        !nearly overlapping node
        if ( (this%z(inode) > (w_conn_z(iconn)-tol_z) ).and. &
             (this%z(inode) < (w_conn_z(iconn)+tol_z) ) &
           ) then
          interp_val(inode) = well_conn_val(ord(iconn)) ! assign well conn value
        else if ( ( this%z(inode) > (w_conn_z(iconn)+tol_z) ).and. &
                  ( this%z(inode) < (w_conn_z(iconn+1)-tol_z) ) &
                ) then
          !perform linear inteprolation
          interp_val(inode) = well_conn_val(ord(iconn)) + &
                              ( well_conn_val(ord(iconn+1)) - &
                                well_conn_val(ord(iconn)) ) / &
                              ( w_conn_z(iconn+1) - w_conn_z(iconn) ) * &
                              ( this%z(inode) - w_conn_z(iconn) ) 
        end if
      end do !end loop over well connections  
    end if !end if extrapolation 
  end do !end loop over vertical fine grid nodes 

end subroutine InterpFromWellConnTo1DGrid

! ************************************************************************** !

subroutine DestroyOneDimGrid(one_d_grid)
  !
  ! destroys OneDimGrid
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/30/15
  ! 

  use Utility_module 

  implicit none

  class(one_dim_grid_type), pointer :: one_d_grid

  if (.not.associated(one_d_grid) ) return

  call DeallocateArray(one_d_grid%z)

  deallocate(one_d_grid)
  nullify(one_d_grid)

end subroutine DestroyOneDimGrid

! ************************************************************************** !

end module Hydrostatic_Common_module

