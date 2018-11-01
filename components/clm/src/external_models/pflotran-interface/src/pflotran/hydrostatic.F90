module Hydrostatic_module
 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private


  public :: HydrostaticUpdateCoupler, &
            HydrostaticTest
 
contains

! ************************************************************************** !

subroutine HydrostaticUpdateCoupler(coupler,option,grid)
  ! 
  ! Computes the hydrostatic initial/boundary
  ! condition
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/28/07
  ! 

  use EOS_Water_module

  use Option_module
  use Grid_module
  use Coupler_module
  use Condition_module
  use Connection_module
  use Region_module
  use Grid_Structured_module
  use Utility_module, only : DotProduct
  use Dataset_Gridded_HDF5_class
  use Dataset_Common_HDF5_class
  use Dataset_Ascii_class
  

  implicit none

  type(coupler_type) :: coupler
  type(option_type) :: option
  type(grid_type) :: grid
  
  PetscInt :: local_id, ghosted_id, iconn
  PetscInt :: num_iteration, ipressure, idatum, num_pressures
  PetscReal :: dist_x, dist_y, dist_z, delta_z, dist_z_for_pressure
  PetscReal :: dx_conn, dy_conn, dz_conn
  PetscReal :: rho_kg, rho_one, rho_zero, pressure, pressure0, pressure_at_datum
  PetscReal :: temperature_at_datum, temperature
  PetscReal :: concentration_at_datum
  PetscReal :: gas_pressure
  PetscReal :: xm_nacl
  PetscReal :: max_z, min_z, temp_real
  PetscInt :: num_faces, face_id_ghosted, conn_id, num_regions
  type(connection_set_type), pointer :: conn_set_ptr
  PetscReal, pointer :: pressure_array(:)
  PetscReal, allocatable :: density_array(:), z(:)
  PetscReal :: pressure_gradient(3), piezometric_head_gradient(3), datum(3)
  PetscReal :: temperature_gradient(3), concentration_gradient(3)
  PetscReal :: gravity_magnitude
  PetscReal :: z_offset
  PetscReal :: aux(1), dummy
  PetscErrorCode :: ierr
  
  class(dataset_gridded_hdf5_type), pointer :: datum_dataset
  PetscReal :: datum_dataset_rmax
  PetscReal :: datum_dataset_rmin
  
  type(flow_condition_type), pointer :: condition
  
  type(connection_set_type), pointer :: cur_connection_set
  
  condition => coupler%flow_condition
  
  xm_nacl = option%m_nacl * FMWNACL
  xm_nacl = xm_nacl /(1.d3 + xm_nacl)
  aux(1) = xm_nacl
  
  nullify(pressure_array)
  nullify(datum_dataset)
  
  delta_z = min((grid%z_max_global-grid%z_min_global)/500.d0,1.d0)
  ! if zero, assign 1.d0 to avoid divide by zero below. essentially the grid
  ! is flat.
  if (delta_z < 1.d-40) delta_z = 1.d0
  temperature_at_datum = option%reference_temperature
  concentration_at_datum = 0.d0
  datum = 0.d0

  gas_pressure = 0.d0
  pressure_gradient = 0.d0
  temperature_gradient = 0.d0
  piezometric_head_gradient = 0.d0
  concentration_gradient = 0.d0
  
  select case(option%iflowmode)

    case default
      ! for now, just set it; in future need to account for a different 
      ! temperature datum
      !geh: this is a trick to determine if the dataset is hdf5 type.
      if (associated(condition%temperature)) then
        if (associated(DatasetCommonHDF5Cast(condition%&
                                             temperature%dataset))) then
          option%io_buffer = 'HDF5-type datasets for temperature are not &
            &supported for hydrostatic, seepage, or conductance boundary &
            &conditions.'
          call printErrMsg(option)
        endif
        if (condition%temperature%itype == DIRICHLET_BC) then
#ifndef THDIRICHLET_TEMP_BC_HACK
          temperature_at_datum = &
            condition%temperature%dataset%rarray(1)
#else
          if (associated(condition%temperature%dataset%rarray)) then
            temperature_at_datum = &
              condition%temperature%dataset%rarray(1)
          else
            temperature_at_datum = option%reference_temperature
          endif
#endif
          if (associated(condition%temperature%gradient)) then
            temperature_gradient(1:3) = &
              condition%temperature%gradient%rarray(1:3)
          endif
        endif
      endif
      if (associated(condition%concentration)) then
        if (condition%concentration%itype == DIRICHLET_BC .and. &
            associated(condition%concentration%dataset)) then
            concentration_at_datum = &
              condition%concentration%dataset%rarray(1)
            if (associated(condition%concentration%gradient)) then
              concentration_gradient(1:3) = &
                condition%concentration%gradient%rarray(1:3)
            endif
        else
          concentration_at_datum = UNINITIALIZED_DOUBLE
          concentration_gradient = 0.d0
        endif
      endif

      pressure_at_datum = &
        condition%pressure%dataset%rarray(1)
      ! gradient is in m/m; needs conversion to Pa/m
      if (associated(condition%pressure%gradient)) then
        piezometric_head_gradient(1:3) = &
          condition%pressure%gradient%rarray(1:3)
      endif
  end select      

  if (associated(condition%datum)) then
    nullify(datum_dataset)
    select type(dataset=>condition%datum)
      class is (dataset_ascii_type)
        datum = dataset%rarray(1:3)
      class is (dataset_gridded_hdf5_type)
        datum_dataset => dataset
        ! set datum here equal to estimated mid value of dataset
        datum(1:3) = UNINITIALIZED_DOUBLE
        datum_dataset_rmax = maxval(datum_dataset%rarray)
        datum_dataset_rmin = minval(datum_dataset%rarray)
        datum(3) = 0.5d0*(datum_dataset_rmax+datum_dataset_rmin)
        ! round the number to the nearest whole number
        datum(3) = int(datum(3))
      class default
        option%io_buffer = &
          'Incorrect dataset type in HydrostaticUpdateCoupler. Dataset "' // &
          trim(condition%datum%name) // '" in file "' // &
          trim(condition%datum%filename) // '".'
        call printErrMsg(option)
    end select
  endif
      
  call EOSWaterDensityExt(temperature_at_datum,pressure_at_datum, &
                          aux,rho_kg,dummy,ierr)
  if (ierr /= 0) then
    call printMsgByCell(option,-1, &
                        'Error in HydrostaticUpdateCoupler->EOSWaterDensity')
  endif
  
  gravity_magnitude = sqrt(DotProduct(option%gravity,option%gravity))
  
  if (dabs(gravity_magnitude-EARTH_GRAVITY) > 0.1d0) then
    option%io_buffer = 'Magnitude of gravity vector is not near 9.81.'
    call printErrMsg(option)
  endif
  
  ! if a pressure gradient is prescribed in Z (at this point it will be a
  ! piezometric head gradient), the units of the pressure gradient are
  ! Pa/m and the pressure gradient does not need conversion
  if (dabs(piezometric_head_gradient(Z_DIRECTION)) < 1.d-40) then
    pressure_gradient(1:3) = piezometric_head_gradient(1:3)* &
                             rho_kg*gravity_magnitude
  else
    pressure_gradient(1:3) = piezometric_head_gradient(1:3)
  endif

  if (dabs(pressure_gradient(Z_DIRECTION)) < 1.d-40) then
    ! compute the vertical gradient based on a 1 meter vertical spacing and
    ! interpolate the values from that array
    if (associated(datum_dataset)) then
      temp_real = grid%z_max_global - grid%z_min_global
      max_z = max(grid%z_max_global,datum_dataset_rmax+temp_real)+1.d0
      min_z = min(grid%z_min_global,datum_dataset_rmin-temp_real)-1.d0
    else
      max_z = max(grid%z_max_global,datum(Z_DIRECTION))+1.d0 ! add 1m buffer
      min_z = min(grid%z_min_global,datum(Z_DIRECTION))-1.d0
    endif
    
    num_pressures = int((max_z-min_z)/delta_z) + 1
    allocate(pressure_array(num_pressures))
    allocate(density_array(num_pressures))
    allocate(z(num_pressures))
    pressure_array = 0.d0
    density_array = 0.d0
    z = 0.d0
    ! place this pressure in the array
    idatum = int((datum(Z_DIRECTION)-min_z)/(max_z-min_z) * &
                 dble(num_pressures))+1
    pressure_array(idatum) = pressure_at_datum
    call EOSWaterDensityExt(temperature_at_datum,pressure_at_datum, &
                            aux,rho_kg,dummy,ierr)
    if (ierr /= 0) then
      call printMsgByCell(option,-2, &
                        'Error in HydrostaticUpdateCoupler->EOSWaterDensity')
    endif
    temperature = temperature_at_datum
    pressure0 = pressure_at_datum
    density_array(idatum) = rho_kg
    z(idatum) = datum(Z_DIRECTION)
    ! compute pressures above datum, if any
    dist_z = 0.d0
    rho_zero = rho_kg
    do ipressure=idatum+1,num_pressures
      dist_z = dist_z + delta_z
      select case(option%iflowmode)
        case(TH_MODE)
          temperature = temperature + temperature_gradient(Z_DIRECTION)*delta_z
      end select
      call EOSWaterDensityExt(temperature,pressure0, &
                              aux,rho_kg,dummy,ierr)
      if (ierr /= 0) then
        call printMsgByCell(option,-3, &
                      'Error in HydrostaticUpdateCoupler->EOSWaterDensity')
      endif
      num_iteration = 0
      do 
        pressure = pressure0 + 0.5d0*(rho_kg+rho_zero) * &
                   option%gravity(Z_DIRECTION) * delta_z
        call EOSWaterDensityExt(temperature,pressure,aux,rho_one,dummy,ierr)
        if (ierr /= 0) then
          call printMsgByCell(option,-4, &
                        'Error in HydrostaticUpdateCoupler->EOSWaterDensity')
        endif
!geh        call EOSWaterDensityNaCl(temperature,pressure,xm_nacl,rho_one) 
        if (dabs(rho_kg-rho_one) < 1.d-10) exit
        rho_kg = rho_one
        num_iteration = num_iteration + 1
        if (num_iteration > 100) then
          print *,'Hydrostatic iteration failed to converge',num_iteration, &
                   rho_one,rho_kg
          print *, condition%name, idatum
          print *, pressure_array
          stop
        endif
      enddo
      rho_zero = rho_kg
      pressure_array(ipressure) = pressure
      density_array(ipressure) = rho_kg
      z(ipressure) = z(idatum)+dist_z
      pressure0 = pressure
    enddo

    ! compute pressures below datum, if any
    pressure0 = pressure_array(idatum)
    select case(option%iflowmode)
      case(TH_MODE)
        temperature = temperature_at_datum
    end select
    dist_z = 0.d0
    rho_zero = density_array(idatum)
    do ipressure=idatum-1,1,-1
      dist_z = dist_z + delta_z
      select case(option%iflowmode)
        case(TH_MODE)
          temperature = temperature - temperature_gradient(Z_DIRECTION)*delta_z
      end select
      call EOSWaterDensityExt(temperature,pressure0,aux,rho_kg,dummy,ierr)
      if (ierr /= 0) then
        call printMsgByCell(option,-5, &
                      'Error in HydrostaticUpdateCoupler->EOSWaterDensity')
      endif
      num_iteration = 0
      do                   ! notice the negative sign (-) here
        pressure = pressure0 - 0.5d0*(rho_kg+rho_zero) * &
                   option%gravity(Z_DIRECTION) * delta_z
        call EOSWaterDensityExt(temperature,pressure,aux,rho_one,dummy,ierr)
        if (ierr /= 0) then
          call printMsgByCell(option,-6, &
                        'Error in HydrostaticUpdateCoupler->EOSWaterDensity')
        endif
        if (dabs(rho_kg-rho_one) < 1.d-10) exit
        rho_kg = rho_one
        num_iteration = num_iteration + 1
        if (num_iteration > 100) then
          print *,'Hydrostatic iteration failed to converge',num_iteration, &
                  rho_one,rho_kg
          print *, condition%name, idatum
          print *, pressure_array
          stop
        endif
      enddo
      rho_zero = rho_kg
      pressure_array(ipressure) = pressure
      density_array(ipressure) = rho_kg
      z(ipressure) = z(idatum)-dist_z
      pressure0 = pressure
    enddo
  endif

  dx_conn = 0.d0
  dy_conn = 0.d0
  dz_conn = 0.d0

  num_faces = coupler%connection_set%num_connections

  do iconn=1, num_faces !geh: this should really be num_faces!
    local_id = coupler%connection_set%id_dn(iconn)
    ghosted_id = grid%nL2G(local_id)
  
    ! geh: note that this is a boundary connection, thus the entire 
    !      distance is between the face and cell center
    if (associated(coupler%connection_set%dist)) then
      dx_conn = coupler%connection_set%dist(0,iconn)* &
        coupler%connection_set%dist(1,iconn)
      dy_conn = coupler%connection_set%dist(0,iconn)* &
        coupler%connection_set%dist(2,iconn)
      dz_conn = coupler%connection_set%dist(0,iconn)* &
        coupler%connection_set%dist(3,iconn)
    endif
    if (associated(datum_dataset)) then
      ! correct datum based on dataset value
      ! if we interpolate in x and y, then we can use grid%x/y - dx/y_conn 
      ! for x and y then we set dist_x and dist_y = 0.
      dist_x = 0.d0
      dist_y = 0.d0
      call DatasetGriddedHDF5InterpolateReal(datum_dataset, &
                                          grid%x(ghosted_id)-dx_conn, &
                                          grid%y(ghosted_id)-dy_conn, &
                                          0.d0,temp_real,option)
      ! temp_real is now the real datum
      dist_z = grid%z(ghosted_id)-dz_conn-temp_real
      z_offset = temp_real-datum(Z_DIRECTION)
    else
      ! note the negative (-) d?_conn is required due to the offset of 
      ! the boundary face
      dist_x = grid%x(ghosted_id)-dx_conn-datum(X_DIRECTION)
      dist_y = grid%y(ghosted_id)-dy_conn-datum(Y_DIRECTION)
      dist_z = grid%z(ghosted_id)-dz_conn-datum(Z_DIRECTION)
      z_offset = 0.d0
    endif

    if (associated(pressure_array)) then
      ipressure = idatum+int(dist_z/delta_z)
      dist_z_for_pressure = grid%z(ghosted_id)-dz_conn-(z(ipressure) + z_offset)
      pressure = pressure_array(ipressure) + &
                 density_array(ipressure)*option%gravity(Z_DIRECTION) * &
                 dist_z_for_pressure + &
!                 (grid%z(ghosted_id)-z(ipressure)) + &
                 pressure_gradient(X_DIRECTION)*dist_x + & ! gradient in Pa/m
                 pressure_gradient(Y_DIRECTION)*dist_y
    else
      pressure = pressure_at_datum + &
                 pressure_gradient(X_DIRECTION)*dist_x + & ! gradient in Pa/m
                 pressure_gradient(Y_DIRECTION)*dist_y + &
                 pressure_gradient(Z_DIRECTION)*dist_z 
    endif
   
 
    if (pressure < option%minimum_hydrostatic_pressure) &
      pressure = option%minimum_hydrostatic_pressure

    ! assign pressure
    select case(option%iflowmode)
      case default
        if (condition%pressure%itype == SEEPAGE_BC) then
          coupler%flow_aux_real_var(1,iconn) = &
            max(pressure,option%reference_pressure)
        else if (condition%pressure%itype == CONDUCTANCE_BC) then
          coupler%flow_aux_real_var(1,iconn) = &
            max(pressure,option%reference_pressure)
          select case(option%iflowmode)
            case(TH_MODE)
              coupler%flow_aux_real_var(TH_CONDUCTANCE_DOF,iconn) = &
                condition%pressure%aux_real(1)
          end select
        else
          coupler%flow_aux_real_var(1,iconn) = pressure
        endif
    end select

    ! assign other dofs
    select case(option%iflowmode)

      case(TH_MODE)
        temperature = temperature_at_datum + &
                    ! gradient in K/m
                    temperature_gradient(X_DIRECTION)*dist_x + & 
                    temperature_gradient(Y_DIRECTION)*dist_y + &
                    temperature_gradient(Z_DIRECTION)*dist_z
        coupler%flow_aux_real_var(TH_TEMPERATURE_DOF,iconn) = temperature
        coupler%flow_aux_int_var(TH_PRESSURE_DOF,iconn) = condition%iphase

      case default
        coupler%flow_aux_int_var(COUPLER_IPHASE_INDEX,iconn) = 1
    end select

  enddo

  if (associated(pressure_array)) deallocate(pressure_array)
  nullify(pressure_array)
  if (allocated(z)) deallocate(z)
  if (allocated(density_array)) deallocate(density_array)
  !geh: Do not deallocate datum_dataset as it is soleley a pointer to an
  !     external dataset.
  nullify(datum_dataset)

end subroutine HydrostaticUpdateCoupler

! ************************************************************************** !

subroutine HydrostaticTest()
  ! 
  ! Computes the hydrostatic initial/boundary
  ! condition (more accurately than before0
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/28/07
  ! 

  use EOS_Water_module
  
  implicit none
  
  PetscInt :: iz, i, i_increment, num_increment
  PetscInt :: max_num_pressures, i_up, i_dn, num_iteration
  PetscReal :: rho_kg, rho_one, rho_zero, pressure0, pressure, temperature
  PetscReal :: increment(4)
  PetscReal :: xm_nacl, dist_z, dist
  PetscReal :: aux(1)
  PetscReal :: dummy
  PetscErrorCode :: ierr

  PetscReal, pointer :: density_array(:,:), pressure_array(:,:)
  
  increment(1) = 1.d-1
  increment(2) = 1.d-0
  increment(3) = 1.d+1
  increment(4) = 1.d+2
  num_increment = size(increment)
  
  temperature = 25.d0

  xm_nacl = 0.d0
  aux(1) = xm_nacl
  
  max_num_pressures = int(1000.d0/increment(1)+0.5d0)+1
  
  allocate(density_array(max_num_pressures,num_increment))
  allocate(pressure_array(max_num_pressures,num_increment)) 
  density_array = 0.d0
  pressure_array = 0.d0
  
  do i_increment = 1, num_increment
    pressure = 101325.d0
    call EOSWaterDensityExt(temperature,pressure,aux,rho_kg,dummy,ierr)
    dist_z = 0.d0
    pressure_array(1,i_increment) = pressure
    density_array(1,i_increment) = rho_kg
    do iz=1,int(1000.d0/increment(i_increment)+0.5d0)
      dist_z = dist_z + increment(i_increment)
      pressure0 = pressure
      num_iteration = 0
      do
        pressure = pressure0 + rho_kg * EARTH_GRAVITY * increment(i_increment)
        call EOSWaterDensityExt(temperature,pressure,aux,rho_one,dummy,ierr)
        if (dabs(rho_kg-rho_one) < 1.d-10) exit
        rho_kg = rho_one
        num_iteration = num_iteration + 1
        if (num_iteration > 100) then
          print *,'HydrostaticInitCondition failed to converge', &
                  num_iteration,rho_one,rho_kg
          stop
        endif
      enddo
      i = int(dist_z/increment(1)+0.5d0)+1
      pressure_array(i,i_increment) = pressure
      density_array(i,i_increment) = rho_kg
      pressure0 = pressure
      rho_zero = rho_kg  
    enddo
  enddo

  do i_increment=2,num_increment
    dist_z = 0.d0
    i = 1
    i_up = 1
    do iz = 1,int(1000.d0/increment(i_increment)+0.5d0)
      i_dn = i_up + int(increment(i_increment)/increment(1)+1.d-6)
      dist = increment(1)
      do 
        i = i + 1
        if (dist >= 0.9999d0*increment(i_increment)) exit
        pressure_array(i,i_increment) = pressure_array(i_up,i_increment)* &
                                        (1.d0-dist/increment(i_increment))+ &
                                        pressure_array(i_dn,i_increment)* &
                                        dist/increment(i_increment)
        density_array(i,i_increment) = density_array(i_up,i_increment)* &
                                       (1.d0-dist/increment(i_increment))+ &
                                       density_array(i_dn,i_increment)* &
                                       dist/increment(i_increment)
        dist = dist + increment(1)
      enddo
      dist_z = dist_z + increment(i_increment)
      i_up = i_dn
    enddo
  enddo

  
  open(unit=86,file='pressures.dat')
  dist_z = 0.d0
  do iz = 1,max_num_pressures
    write(86,'(100(es16.10,x))') dist_z, &
                (density_array(iz,i_increment),i_increment=1,num_increment), &
                (pressure_array(iz,i_increment),i_increment=1,num_increment)
    dist_z = dist_z + increment(1)
  enddo
  close(86)
  
  deallocate(pressure_array)
  deallocate(density_array)
    
end subroutine HydrostaticTest

end module Hydrostatic_module

