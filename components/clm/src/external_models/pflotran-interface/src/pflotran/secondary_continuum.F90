! added by S. Karra 07/11/12

module Secondary_Continuum_module
  
#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use Secondary_Continuum_Aux_module

  use PFLOTRAN_Constants_module
  use Utility_module, only : Equal
  
  implicit none

  private

  ! secondary continuum cell type
  PetscInt, parameter, public :: SLAB = 0
  PetscInt, parameter, public :: NESTED_CUBE = 1
  PetscInt, parameter, public :: NESTED_SPHERE = 2

  PetscReal, parameter :: perturbation_tolerance = 1.d-5

  public :: SecondaryContinuumType, &
            SecondaryContinuumSetProperties, &
            SecondaryRTAuxVarInit, &
            SecondaryRTResJacMulti, &
            SecondaryRTAuxVarComputeMulti, &
            THCSecHeatAuxVarCompute, &
            THSecHeatAuxVarCompute, &
            SecondaryRTUpdateIterate, &
            SecondaryRTUpdateEquilState, &
            SecondaryRTUpdateKineticState, &
            SecondaryRTTimeCut

contains

! ************************************************************************** !

subroutine SecondaryContinuumType(sec_continuum,nmat,aream, &
            volm,dm1,dm2,aperture,epsilon,log_spacing,outer_spacing, &
            interfacial_area,option)
  ! 
  ! The area, volume, grid sizes for secondary continuum
  ! are calculated based on the input dimensions and geometry
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 07/11/12
  ! 

  use Option_module

  implicit none
  
  type(sec_continuum_type) :: sec_continuum

  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  PetscInt :: igeom, nmat, m
  PetscReal :: aream(nmat), volm(nmat), dm1(nmat), dm2(nmat)
  PetscReal :: dy, r0, r1, aream0, am0, vm0, interfacial_area
  PetscReal :: num_density, aperture, epsilon, fracture_spacing
  PetscReal :: outer_spacing, matrix_block_size
  PetscReal :: grid_spacing(nmat)
  PetscBool :: log_spacing
  PetscReal :: sum

  PetscInt, save :: icall

  data icall/0/

  igeom = sec_continuum%itype
  option%nsec_cells = nmat
    
  select case (igeom)      
    case(SLAB)
    
      dy = sec_continuum%slab%length/nmat
      aream0 = sec_continuum%slab%area
      do m = 1, nmat
        volm(m) = dy*aream0
      enddo
      am0 = aream0
      vm0 = nmat*dy*aream0
      interfacial_area = am0/vm0
     
      do m = 1, nmat
        aream(m) = aream0
        dm1(m) = 0.5d0*dy
        dm2(m) = 0.5d0*dy
      enddo

      if (icall == 0 .and. OptionPrintToFile(option)) then
        icall = 1
        string = 'DCDM Multiple Continuum Model'
        write(option%fid_out,'(/,2x,a,/)') trim(string)
        string = 'Slab'
        write(option%fid_out,'(2x,a,/)') trim(string)
        num_density = (1.d0-epsilon)/vm0
        write(option%fid_out,'(2x,"number density: ",11x,1pe12.4," m^(-3)")') num_density
        write(option%fid_out,'(2x,"matrix block size: ",8x,1pe12.4," m")') sec_continuum%slab%length
        write(option%fid_out,'(2x,"epsilon: ",18x,1pe12.4)') epsilon
        write(option%fid_out,'(2x,"specific interfacial area: ",1pe12.4," m^(-1)")') interfacial_area
        do m = 1, nmat
          if (m == 1) write(option%fid_out,'(/,2x,"node matrix volume fraction")') 
          write(option%fid_out,'(2x,i3,3x,1pe12.4)') m,volm(m)/vm0 !*(1.d0 - epsilon)
        enddo
!       aperture = r0*(1.d0/(1.d0-epsilon)**(1.d0/3.d0)-1.d0)
!       write(option%fid_out,'(2x,"aperture: ",17x,1pe12.4," m")') aperture
      endif
      
      ! Store the distances
      sec_continuum%distance(1) = dm1(1)
      do m = 2, nmat
        sec_continuum%distance(m) = sec_continuum%distance(m-1) + &
                                      dm2(m-1) + dm1(m)
      enddo
          
    case(NESTED_CUBE)

      if (sec_continuum%nested_cube%fracture_spacing > 0.d0) then

        fracture_spacing = sec_continuum%nested_cube%fracture_spacing
!        override epsilon if aperture defined
        if (aperture > 0.d0) then
          r0 = fracture_spacing - aperture
          epsilon = 1.d0 - (1.d0 + aperture/r0)**(-3.d0)
        else if (epsilon > 0.d0) then
          r0 = fracture_spacing*(1.d0-epsilon)**(1.d0/3.d0)
          aperture = r0*((1.d0-epsilon)**(-1.d0/3.d0)-1.d0)
        endif
                                            
      else if (sec_continuum%nested_cube%matrix_block_size > 0.d0) then

        r0 = sec_continuum%nested_cube%matrix_block_size

!        override epsilon if aperture defined
        if (aperture > 0.d0) then
          fracture_spacing = r0 + aperture
          epsilon = 1.d0 - (1.d0 + aperture/r0)**(-3.d0)
        else if (epsilon > 0.d0) then
          fracture_spacing = r0*(1.d0-epsilon)**(-1.d0/3.d0)
          aperture = fracture_spacing - r0
        endif
      endif
      
      if (log_spacing) then 
        
        matrix_block_size = r0
        call SecondaryContinuumCalcLogSpacing(matrix_block_size,outer_spacing, &
                                              nmat,grid_spacing,option)
        
        r0 = 2.d0*grid_spacing(1)
        dm1(1) = 0.5d0*grid_spacing(1)
        dm2(1) = 0.5d0*grid_spacing(1)
        volm(1) = r0**3.d0
        aream(1) = 6.d0*r0**2.d0
        do m = 2, nmat
          dm1(m) = 0.5d0*grid_spacing(m)
          dm2(m) = 0.5d0*grid_spacing(m)
          r1 = r0 + 2.d0*(dm1(m) + dm2(m))
          volm(m) = r1**3.d0 - r0**3.d0
          aream(m) = 6.d0*r1**2.d0
          r0 = r1
        enddo
        r0 = matrix_block_size
        am0 = 6.d0*r0**2.d0
        vm0 = r0**3.d0
        interfacial_area = am0/vm0

      else
        dy = r0/nmat/2.d0
     
        r0 = 2.d0*dy
        volm(1) = r0**3.d0
        do m = 2, nmat
          r1 = r0 + 2.d0*dy
          volm(m) = r1**3.d0 - r0**3.d0
          r0 = r1
        enddo
      
        r0 = 2.d0*dy
        aream(1) = 6.d0*r0**2.d0
        dm1(1) = 0.5d0*dy
        dm2(1) = 0.5d0*dy
        do m = 2, nmat
          dm1(m) = 0.5d0*dy
          dm2(m) = 0.5d0*dy
          r0 = r0 + 2.d0*dy
          aream(m) = 6.d0*r0**2.d0
        enddo
        r0 = real(2.d0*nmat)*dy
        am0 = 6.d0*r0**2.d0
        vm0 = r0**3.d0
        interfacial_area = am0/vm0
      endif

      if (icall == 0 .and. OptionPrintToFile(option)) then
        icall = 1
        string = 'DCDM Multiple Continuum Model'
        write(option%fid_out,'(/,2x,a,/)') trim(string)
        string = 'Nested Cubes'
        write(option%fid_out,'(2x,a,/)') trim(string)
        num_density = (1.d0-epsilon)/vm0
        write(option%fid_out,'(2x,"number density: ",11x,1pe12.4," m^(-3)")') num_density
        write(option%fid_out,'(2x,"matrix block size: ",8x,1pe12.4," m")') r0
        write(option%fid_out,'(2x,"epsilon: ",18x,1pe12.4)') epsilon
        write(option%fid_out,'(2x,"specific interfacial area: ",1pe12.4," m^(-1)")') interfacial_area
        write(option%fid_out,'(2x,"fracture aperture: ",8x,1pe12.4," m")') aperture
        write(option%fid_out,'(2x,"fracture spacing: ",9x,1pe12.4," m")') fracture_spacing
        write(option%fid_out,'(/,2x,"node  vol. frac.      dm1         dm2         aream       dy          y")')
        r0 = 0.d0
        do m = 1, nmat
          if (m == 1) then
            r0 = r0 + dm1(m)
          else
            r0 = r0 + dm2(m-1)+dm1(m)
          endif
          write(option%fid_out,'(2x,i3,3x,1p6e12.4)') m,volm(m)/vm0,dm1(m),dm2(m),aream(m), &
          dm1(m)+dm2(m),r0
        enddo
      endif

      ! Store the distances
      sec_continuum%distance(1) = dm1(1)
      do m = 2, nmat
        sec_continuum%distance(m) = sec_continuum%distance(m-1) + &
                                      dm2(m-1) + dm1(m)
      enddo     

    case(NESTED_SPHERE)
    
      dy = sec_continuum%nested_sphere%radius/nmat
      r0 = dy

      volm(1) = 4.d0/3.d0*pi*r0**3.d0
      do m = 2, nmat
        r1 = r0 + dy
        volm(m) = 4.d0/3.d0*pi*(r1**3.d0 - r0**3.d0)
        r0 = r1
      enddo
      
      r0 = dy
      aream(1) = 4.d0*pi*r0**2.d0
      dm1(1) = 0.5d0*dy
      dm2(1) = 0.5d0*dy
      do m = 2, nmat
        r0 = r0 + dy
        dm1(m) = 0.5d0*dy
        dm2(m) = 0.5d0*dy
        aream(m) = 4.d0*pi*r0**2.d0
      enddo
      r0 = 0.5d0*real(2.d0*nmat)*dy
      am0 = 4.d0*pi*r0**2.d0
      vm0 = am0*r0/3.d0
      interfacial_area = am0/vm0

      if (icall == 0 .and. OptionPrintToFile(option)) then
        icall = 1
        string = 'DCDM Multiple Continuum Model'
        write(option%fid_out,'(/,2x,a,/)') trim(string)
        string = 'Nested Spheres'
        write(option%fid_out,'(2x,a,/)') trim(string)
        num_density = (1.d0-epsilon)/vm0
        write(option%fid_out,'(2x,"number density: ",11x,1pe12.4," m^(-3)")') num_density
        write(option%fid_out,'(2x,"sphere radius: ",8x,1pe12.4," m")') sec_continuum%nested_sphere%radius
        write(option%fid_out,'(2x,"epsilon: ",18x,1pe12.4)') epsilon
        write(option%fid_out,'(2x,"specific interfacial area: ",1pe12.4," m^(-1)")') interfacial_area
        do m = 1, nmat
          if (m == 1) write(option%fid_out,'(/,2x,"node matrix volume fraction")') 
          write(option%fid_out,'(2x,i3,3x,1pe12.4)') m,volm(m)/vm0*(1.d0 - epsilon)
        enddo

!       aperture = r0*(1.d0/(1.d0-epsilon)**(1.d0/3.d0)-1.d0)
!       write(option%fid_out,'(2x,"aperture: ",17x,1pe12.4," m")') aperture
      endif
      
      ! Store the distances
      sec_continuum%distance(1) = dm1(1)
      do m = 2, nmat
        sec_continuum%distance(m) = sec_continuum%distance(m-1) + &
                                      dm2(m-1) + dm1(m)
      enddo
                        
  end select
  
  
  sum = 0.d0
  do m = 1,nmat
    if (volm(m)/vm0 > 1.d0) then
      print *, 'Error: volume fraction for cell', m, 'is greater than 1.'
      stop
    else 
      sum = sum + volm(m)/vm0
    endif
  enddo
  
  if (icall /= 2 .and. OptionPrintToFile(Option)) then
    icall = 2
    write(option%fid_out,'(/,"sum of volume fractions:",1x,1pe12.4)') sum
  endif
  
  if (abs(sum - 1.d0) > 1.d-6) then
    option%io_buffer = 'Error: Sum of the volume fractions of the' // &
                       ' secondary cells is not equal to 1.'
    call printErrMsg(option)
  endif
  
end subroutine SecondaryContinuumType

! ************************************************************************** !

subroutine SecondaryContinuumSetProperties(sec_continuum, &
                                           sec_continuum_name, & 
                                           sec_continuum_length, &
                                           sec_continuum_matrix_block_size, &
                                           sec_continuum_fracture_spacing, &
                                           sec_continuum_radius, &
                                           sec_continuum_area, &
                                           option)
  ! 
  ! The type, dimensions of the secondary
  ! continuum are set
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 07/17/12
  ! 
                                    
  use Option_module
  use String_module
  
  implicit none
  
  type(sec_continuum_type) :: sec_continuum
  type(option_type) :: option
  PetscReal :: sec_continuum_matrix_block_size
  PetscReal :: sec_continuum_fracture_spacing
  PetscReal :: sec_continuum_length
  PetscReal :: sec_continuum_area
  PetscReal :: sec_continuum_radius
  character(len=MAXWORDLENGTH) :: sec_continuum_name

  call StringToUpper(sec_continuum_name)
  
  select case(trim(sec_continuum_name))
    case("SLAB")
      sec_continuum%itype = SLAB
      sec_continuum%slab%length = sec_continuum_length
      if (Equal(sec_continuum_area,0.d0)) then
        option%io_buffer = 'Keyword "AREA" not specified for SLAB type ' // &
                           'under SECONDARY_CONTINUUM'
        call printErrMsg(option)
      endif
      sec_continuum%slab%area = sec_continuum_area
    case("NESTED_CUBES")
      sec_continuum%itype = NESTED_CUBE
      sec_continuum%nested_cube%matrix_block_size = sec_continuum_matrix_block_size
      sec_continuum%nested_cube%fracture_spacing = sec_continuum_fracture_spacing
    case("NESTED_SPHERES")
      sec_continuum%itype = NESTED_SPHERE
      sec_continuum%nested_sphere%radius = sec_continuum_radius
    case default
      option%io_buffer = 'Keyword "' // trim(sec_continuum_name) // '" not ' // &
                         'recognized in SecondaryContinuumSetProperties()'
      call printErrMsg(option)  
  end select
      
end subroutine SecondaryContinuumSetProperties  

! ************************************************************************** !

subroutine SecondaryContinuumCalcLogSpacing(matrix_size,outer_grid_size, &
                                            sec_num_cells,grid_spacing,option)
  ! 
  ! Given the matrix block size and the
  ! grid spacing of the outer most secondary continuum cell, a geometric
  ! series is assumed and the grid spacing of the rest of the cells is
  ! calculated
  ! 
  ! Equation:
  ! \frac{1 - \rho}{1 - \rho_M}*\rho*(M-1) = \frac{2\Delta\xi_m}{l_M}
  !
  ! where
  !   \Delta\xi_m: Grid spacing of the outer most continuum cell (INPUT)
  !   l_M        : Matrix block size (INPUT)
  !   M          : Number of secondary continuum cells (INPUT)
  !   \rho       : Logarithmic grid spacing factor (COMPUTED)
  !
  ! Author: Satish Karra, LANL
  ! Date: 07/17/12
  ! 
                                              
  use Option_module
  
  implicit none
  
  type(option_type) :: option
  PetscReal :: matrix_size, outer_grid_size
  PetscInt :: sec_num_cells
  PetscReal :: grid_spacing(sec_num_cells)
  PetscReal :: delta, delta_new, inner_grid_size
  PetscReal :: F, dF
  PetscReal, parameter :: tol = 1.d-12
  PetscInt, parameter :: maxit = 50
  PetscInt :: i 
  
  
  if (mod(sec_num_cells,2) /= 0) then
     option%io_buffer = 'NUM_CELLS under SECONDARY_CONTINUUM has to be' // &
                        ' even for logarithmic grid spacing'
      call printErrMsg(option)
  endif
  
  delta = 0.99d0
  
  do i = 1, maxit
    F = (1.d0 - delta)/(1.d0 - delta**sec_num_cells)*delta**(sec_num_cells - 1.d0) - &
        2.d0*outer_grid_size/matrix_size
    dF = (1.d0 + sec_num_cells*(delta - 1.d0) - delta**sec_num_cells)/ &
         (delta**sec_num_cells - 1.d0)**2.d0*delta**(sec_num_cells - 2.d0)
    delta_new = delta + F/dF
    if ((abs(F) < tol)) exit
    delta = delta_new
    if (delta < 0.d0) delta = 0.5d0
!   if (delta > 1.d0) delta = 0.9d0
  enddo
  
  if (i == maxit) then
     option%io_buffer = 'Log Grid spacing solution has not converged' // &
                        ' with given fracture values.'
     call printErrMsg(option)    
  endif

  inner_grid_size = outer_grid_size/delta**(sec_num_cells - 1)
  
  do i = 1, sec_num_cells
    grid_spacing(i) = inner_grid_size*delta**(i-1)
  enddo

!  write(option%fid_out,'("  Logarithmic grid spacing: delta = ",1pe12.4)') delta
    
end subroutine SecondaryContinuumCalcLogSpacing

! ************************************************************************** !

subroutine SecondaryRTTimeCut(realization)
  ! 
  ! Resets secondary concentrations to previous time
  ! step when there is a time cut
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 05/29/13
  ! 

  use Realization_Subsurface_class
  use Grid_module
  use Reaction_Aux_module
  
  implicit none
  class(realization_subsurface_type) :: realization
  type(reaction_type), pointer :: reaction
  type(sec_transport_type), pointer :: rt_sec_transport_vars(:)
  type(grid_type), pointer :: grid
  
  PetscInt :: local_id, ghosted_id
  PetscInt :: ngcells, ncomp
  PetscInt :: cell, comp

  reaction => realization%reaction
  rt_sec_transport_vars => realization%patch%aux%SC_RT%sec_transport_vars
  grid => realization%patch%grid  
  
  ncomp = reaction%naqcomp
  
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (realization%patch%imat(ghosted_id) <= 0) cycle
    do comp = 1, ncomp
      ngcells = rt_sec_transport_vars(local_id)%ncells
      do cell = 1, ngcells
        rt_sec_transport_vars(local_id)%updated_conc(comp,cell) = &
          rt_sec_transport_vars(local_id)%sec_rt_auxvar(cell)%pri_molal(comp)
      enddo
    enddo
  enddo
 
end subroutine SecondaryRTTimeCut

! ************************************************************************** !

subroutine SecondaryRTAuxVarInit(ptr,rt_sec_transport_vars,reaction, &
                                 initial_condition,constraint,option)
  ! 
  ! Initializes all the secondary continuum reactive
  ! transport variables
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 02/05/13
  ! 
  
  use Coupler_module
  use Transport_Constraint_module
  use Condition_module
  use Global_Aux_module
  use Material_module
  use Option_module
  use Reaction_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Material_Aux_class
  
  use EOS_Water_module
  
  implicit none 
  
  type(sec_transport_type) :: rt_sec_transport_vars
  type(material_property_type), pointer :: ptr
  type(reaction_type), pointer :: reaction
  type(coupler_type), pointer :: initial_condition
  type(option_type), pointer :: option
  type(reactive_transport_auxvar_type), pointer :: rt_auxvar
  type(global_auxvar_type), pointer :: global_auxvar
  class(material_auxvar_type), allocatable :: material_auxvar
  type(tran_constraint_type), pointer :: constraint
  type(flow_condition_type), pointer :: initial_flow_condition
  

  PetscReal :: equil_conc(reaction%mineral%nmnrl)
  PetscInt :: i, cell
  PetscReal :: area_per_vol
  PetscReal :: dum1
  PetscInt :: num_iterations
  PetscErrorCode :: ierr
  
  num_iterations = 0

  allocate(material_auxvar)
  call MaterialAuxVarInit(material_auxvar,option)
  material_auxvar%porosity = option%reference_porosity

  call SecondaryContinuumSetProperties( &
        rt_sec_transport_vars%sec_continuum, &
        ptr%secondary_continuum_name, &
        ptr%secondary_continuum_length, &
        ptr%secondary_continuum_matrix_block_size, &
        ptr%secondary_continuum_fracture_spacing, &
        ptr%secondary_continuum_radius, &
        ptr%secondary_continuum_area, &
        option)
        
  rt_sec_transport_vars%ncells = ptr%secondary_continuum_ncells
  rt_sec_transport_vars%aperture = ptr%secondary_continuum_aperture
  rt_sec_transport_vars%epsilon = ptr%secondary_continuum_epsilon 
  rt_sec_transport_vars%log_spacing = ptr%secondary_continuum_log_spacing
  rt_sec_transport_vars%outer_spacing = ptr%secondary_continuum_outer_spacing    
        
  allocate(rt_sec_transport_vars%area(rt_sec_transport_vars%ncells))
  allocate(rt_sec_transport_vars%vol(rt_sec_transport_vars%ncells))
  allocate(rt_sec_transport_vars%dm_minus(rt_sec_transport_vars%ncells))
  allocate(rt_sec_transport_vars%dm_plus(rt_sec_transport_vars%ncells))
  allocate(rt_sec_transport_vars%sec_continuum% &
             distance(rt_sec_transport_vars%ncells))
    
  call SecondaryContinuumType(rt_sec_transport_vars%sec_continuum, &
                              rt_sec_transport_vars%ncells, &
                              rt_sec_transport_vars%area, &
                              rt_sec_transport_vars%vol, &
                              rt_sec_transport_vars%dm_minus, &
                              rt_sec_transport_vars%dm_plus, &
                              rt_sec_transport_vars%aperture, &
                              rt_sec_transport_vars%epsilon, &
                              rt_sec_transport_vars%log_spacing, &
                              rt_sec_transport_vars%outer_spacing, &
                              area_per_vol,option)                                
  rt_sec_transport_vars%interfacial_area = area_per_vol* &
         (1.d0 - rt_sec_transport_vars%epsilon)*ptr% &
         secondary_continuum_area_scaling
  
  ! Initializing the secondary RT auxvars
  allocate(rt_sec_transport_vars%sec_rt_auxvar(rt_sec_transport_vars%ncells))
  do cell = 1, rt_sec_transport_vars%ncells
    call RTAuxVarInit(rt_sec_transport_vars%sec_rt_auxvar(cell),reaction,option)
  enddo

  allocate(rt_sec_transport_vars%sec_jac(reaction%naqcomp,reaction%naqcomp))    
           
  ! Allocate diagonal terms
  allocate(rt_sec_transport_vars%cxm(reaction%naqcomp,reaction%naqcomp,&
           rt_sec_transport_vars%ncells)) 
  allocate(rt_sec_transport_vars%cxp(reaction%naqcomp,reaction%naqcomp,&
           rt_sec_transport_vars%ncells))  
  allocate(rt_sec_transport_vars%cdl(reaction%naqcomp,reaction%naqcomp,&
           rt_sec_transport_vars%ncells)) 
  allocate(rt_sec_transport_vars% &
           r(reaction%naqcomp*rt_sec_transport_vars%ncells))
  allocate(rt_sec_transport_vars% &
           updated_conc(reaction%naqcomp,rt_sec_transport_vars%ncells))
           
  
  initial_flow_condition => initial_condition%flow_condition
  do cell = 1, rt_sec_transport_vars%ncells
    global_auxvar => initial_condition%tran_condition% &
                       constraint_coupler_list%global_auxvar
    rt_auxvar => rt_sec_transport_vars%sec_rt_auxvar(cell)
    if (associated(initial_flow_condition)) then
      if (associated(initial_flow_condition%pressure)) then
        if (associated(initial_flow_condition%pressure%dataset)) then
          global_auxvar%pres = &
            initial_flow_condition%pressure%dataset%rarray(1)
        else
          global_auxvar%pres = option%reference_pressure
        endif
      else 
        global_auxvar%pres = option%reference_pressure
      endif
      if (associated(initial_flow_condition%temperature)) then
        if (associated(initial_flow_condition%temperature%dataset)) then
          global_auxvar%temp  = &
            initial_flow_condition%temperature%dataset%rarray(1)
        else
          global_auxvar%temp = option%reference_temperature
        endif
      else
        global_auxvar%temp = option%reference_temperature
      endif
        
      call EOSWaterDensity(global_auxvar%temp, &
                           global_auxvar%pres(1), &
                           global_auxvar%den_kg(1), &
                           dum1,ierr)
    else
      global_auxvar%pres = option%reference_pressure
      global_auxvar%temp = option%reference_temperature
      global_auxvar%den_kg(option%liquid_phase) = &
        option%reference_density(option%liquid_phase)

    endif
    global_auxvar%sat = option%reference_saturation

    if (option%transport%nphase > option%nphase) then
      ! gas phase not considered explicitly on flow side
      global_auxvar%den_kg(option%gas_phase) = &
        option%reference_density(option%gas_phase)
      global_auxvar%sat(option%gas_phase) = &
        1.d0 - global_auxvar%sat(option%liquid_phase)
    endif
                      
    call ReactionEquilibrateConstraint(rt_auxvar,global_auxvar, &
                          material_auxvar, &
                          reaction,constraint%name, &
                          constraint%aqueous_species, &
                          constraint%free_ion_guess, &
                          constraint%minerals, &
                          constraint%colloids, &
                          constraint%immobile_species, &
                          num_iterations, &
                          PETSC_FALSE,option)   
                          
    rt_sec_transport_vars%updated_conc(:,cell) =  rt_auxvar%pri_molal   
       
  enddo                                    
  
  call MaterialAuxVarStrip(material_auxvar)
  deallocate(material_auxvar)
  
  rt_sec_transport_vars%sec_jac_update = PETSC_FALSE
  rt_sec_transport_vars%sec_jac = 0.d0
  rt_sec_transport_vars%cxm = 0.d0
  rt_sec_transport_vars%cxp = 0.d0
  rt_sec_transport_vars%cdl = 0.d0
  rt_sec_transport_vars%r = 0.d0
      
end subroutine SecondaryRTAuxVarInit  

! ************************************************************************** !

subroutine SecondaryRTResJacMulti(sec_transport_vars,auxvar, &
                                  global_auxvar,prim_vol, &
                                  reaction,diffusion_coefficient, &
                                  porosity,option,res_transport)
  ! 
  ! RTSecondaryTransportMulti:  Calculates the source term contribution due to
  ! secondary continuum in the primary continuum residual for multicomponent
  ! system assuming only aqueous reaction
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 1/31/13
  ! 
                               
                            
  use Option_module 
  use Global_Aux_module
  use Block_Solve_module
  use Block_Tridiag_module
  use Utility_module
  use Reaction_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Material_Aux_class

  implicit none
  
  type(sec_transport_type) :: sec_transport_vars
  type(reactive_transport_auxvar_type) :: auxvar
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(reaction_type), pointer :: reaction
  type(option_type) :: option
  PetscReal :: coeff_left(reaction%naqcomp,reaction%naqcomp, &
                          sec_transport_vars%ncells)
  PetscReal :: coeff_diag(reaction%naqcomp,reaction%naqcomp, &
                          sec_transport_vars%ncells)
  PetscReal :: coeff_right(reaction%naqcomp,reaction%naqcomp, &
                           sec_transport_vars%ncells)
  PetscReal :: res(sec_transport_vars%ncells*reaction%naqcomp)
  PetscReal :: rhs(sec_transport_vars%ncells*reaction%naqcomp)
  PetscReal :: D_M(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: identity(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: b_M(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: sec_jac(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: inv_D_M(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: conc_upd(reaction%naqcomp,sec_transport_vars%ncells) 
  PetscReal :: total_upd(reaction%naqcomp,sec_transport_vars%ncells)
  PetscReal :: total_prev(reaction%naqcomp,sec_transport_vars%ncells)
  PetscReal :: conc_current_M(reaction%naqcomp)
  PetscReal :: total_current_M(reaction%naqcomp)
  PetscReal :: res_transport(reaction%naqcomp)
  PetscReal :: total_primary_node(reaction%naqcomp)
  PetscReal :: area(sec_transport_vars%ncells)
  PetscReal :: vol(sec_transport_vars%ncells)
  PetscReal :: dm_plus(sec_transport_vars%ncells)
  PetscReal :: dm_minus(sec_transport_vars%ncells)
  PetscReal :: res_react(reaction%naqcomp)
  PetscReal :: jac_react(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: dtotal(reaction%naqcomp,reaction%naqcomp,sec_transport_vars%ncells)
  PetscReal :: dtotal_prim(reaction%naqcomp,reaction%naqcomp)
  PetscInt :: i, j, k, n, l
  PetscInt :: ngcells, ncomp
  PetscReal :: area_fm
  PetscReal :: diffusion_coefficient
  PetscReal :: porosity
  PetscReal :: arrhenius_factor
  PetscReal :: pordt, pordiff
  PetscReal :: prim_vol ! volume of primary grid cell
  PetscReal :: dCsec_dCprim(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: dPsisec_dCprim(reaction%naqcomp,reaction%naqcomp)
  PetscInt :: jcomp, lcomp, kcomp, icplx, ncompeq
  PetscReal :: sec_sec_molal_M(reaction%neqcplx)   ! secondary species molality of secondary continuum
  
  PetscInt :: pivot(reaction%naqcomp,sec_transport_vars%ncells)
  PetscInt :: indx(reaction%naqcomp)
  PetscInt :: d, ier
  PetscReal :: m
  
  ! Quantities for numerical jacobian
  PetscReal :: conc_prim(reaction%naqcomp)
  PetscReal :: conc_prim_pert(reaction%naqcomp)
  PetscReal :: sec_jac_num(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: conc_current_M_pert(reaction%naqcomp)
  PetscReal :: total_current_M_pert(reaction%naqcomp)
  PetscReal :: res_transport_pert(reaction%naqcomp)
  PetscReal :: total_primary_node_pert(reaction%naqcomp)
  PetscReal :: dtotal_prim_num(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: dPsisec_dCprim_num(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: pert
  PetscReal :: coeff_diag_dm(reaction%naqcomp,reaction%naqcomp, &
                          sec_transport_vars%ncells)
  PetscReal :: coeff_left_dm(reaction%naqcomp,reaction%naqcomp, &
                          sec_transport_vars%ncells)
  PetscReal :: coeff_right_dm(reaction%naqcomp,reaction%naqcomp, &
                          sec_transport_vars%ncells)
  PetscReal :: coeff_left_pert(reaction%naqcomp,reaction%naqcomp, &
                          sec_transport_vars%ncells)
  PetscReal :: coeff_diag_pert(reaction%naqcomp,reaction%naqcomp, &
                          sec_transport_vars%ncells)
  PetscReal :: coeff_right_pert(reaction%naqcomp,reaction%naqcomp, &
                           sec_transport_vars%ncells)
  PetscReal :: coeff_left_copy(reaction%naqcomp,reaction%naqcomp, &
                          sec_transport_vars%ncells)
  PetscReal :: coeff_diag_copy(reaction%naqcomp,reaction%naqcomp, &
                          sec_transport_vars%ncells)
  PetscReal :: coeff_right_copy(reaction%naqcomp,reaction%naqcomp, &
                           sec_transport_vars%ncells)

  PetscReal :: total_sorb_upd(reaction%naqcomp,sec_transport_vars%ncells) 
  PetscReal :: total_sorb_prev(reaction%naqcomp,sec_transport_vars%ncells)
  PetscReal :: dtotal_sorb_upd(reaction%naqcomp,reaction%naqcomp,sec_transport_vars%ncells)

  class(material_auxvar_type), allocatable :: material_auxvar
  
  ngcells = sec_transport_vars%ncells
  area = sec_transport_vars%area
  vol = sec_transport_vars%vol          
  dm_plus = sec_transport_vars%dm_plus
  dm_minus = sec_transport_vars%dm_minus
  area_fm = sec_transport_vars%interfacial_area
  ncomp = reaction%naqcomp
  
  do j = 1, ncomp
    do i = 1, ngcells
      total_prev(j,i) = sec_transport_vars%sec_rt_auxvar(i)%total(j,1)
      if (reaction%neqsorb > 0) then
        total_sorb_prev(j,i) = sec_transport_vars%sec_rt_auxvar(i)%total_sorb_eq(j)
      endif
    enddo
  enddo
  conc_upd = sec_transport_vars%updated_conc
    
  ! Note that sec_transport_vars%sec_rt_auxvar(i)%pri_molal(j) units are in mol/kg
  ! Need to convert to mol/L since the units of total. in the Thomas 
  ! algorithm are in mol/L 
  
  coeff_left = 0.d0
  coeff_diag = 0.d0
  coeff_right = 0.d0
  res = 0.d0
  rhs = 0.d0
  D_M = 0.d0
  identity = 0.d0
  b_M = 0.d0
  inv_D_M = 0.d0
  total_current_M = 0.d0
  dPsisec_dCprim = 0.d0
  dCsec_dCprim = 0.d0
  
  total_primary_node = auxvar%total(:,1)                         ! in mol/L 
  dtotal_prim = auxvar%aqueous%dtotal(:,:,1)
  pordt = porosity/option%tran_dt
  pordiff = porosity*diffusion_coefficient

  call RTAuxVarInit(rt_auxvar,reaction,option)
  do i = 1, ngcells
    call RTAuxVarCopy(rt_auxvar,sec_transport_vars%sec_rt_auxvar(i),option)
    rt_auxvar%pri_molal = conc_upd(:,i)
    call RTotalAqueous(rt_auxvar,global_auxvar,reaction,option)
    if (reaction%neqsorb > 0) then
      call SecondaryRTotalSorb(rt_auxvar,global_auxvar,material_auxvar,reaction,option)
    endif
    total_upd(:,i) = rt_auxvar%total(:,1)
    dtotal(:,:,i) = rt_auxvar%aqueous%dtotal(:,:,1)
    if (reaction%neqsorb > 0) then 
      total_sorb_upd(:,i) = rt_auxvar%total_sorb_eq(:)
      dtotal_sorb_upd(:,:,i) = rt_auxvar%dtotal_sorb_eq(:,:)
    endif
  enddo 
                          
!================ Calculate the secondary residual =============================        

  do j = 1, ncomp
      
    ! Accumulation
    do i = 1, ngcells
      n = j + (i-1)*ncomp
      res(n) = pordt*(total_upd(j,i) - total_prev(j,i))*vol(i)    ! in mol/L*m3/s
      if (reaction%neqsorb > 0) then 
        res(n) = res(n) + vol(i)/option%tran_dt*(total_sorb_upd(j,i) - total_sorb_prev(j,i))
      endif      
    enddo
  
    ! Flux terms
    do i = 2, ngcells - 1
      n = j + (i-1)*ncomp
      res(n) = res(n) - pordiff*area(i)/(dm_minus(i+1) + dm_plus(i))* &
                        (total_upd(j,i+1) - total_upd(j,i))
      res(n) = res(n) + pordiff*area(i-1)/(dm_minus(i) + dm_plus(i-1))* &
                        (total_upd(j,i) - total_upd(j,i-1))                      
    enddo
         
              
    ! Apply boundary conditions
    ! Inner boundary
    res(j) = res(j) - pordiff*area(1)/(dm_minus(2) + dm_plus(1))* &
                      (total_upd(j,2) - total_upd(j,1))
                                      
    ! Outer boundary
    res(j+(ngcells-1)*ncomp) = res(j+(ngcells-1)*ncomp) - &
                               pordiff*area(ngcells)/dm_plus(ngcells)* &
                               (total_primary_node(j) - total_upd(j,ngcells))
    res(j+(ngcells-1)*ncomp) = res(j+(ngcells-1)*ncomp) + &
                               pordiff*area(ngcells-1)/(dm_minus(ngcells) &
                               + dm_plus(ngcells-1))*(total_upd(j,ngcells) - &
                               total_upd(j,ngcells-1))  
                               
  enddo
                         
  res = res*1.d3 ! Convert mol/L*m3/s to mol/s                                                    
                                                                                                          
!================ Calculate the secondary jacobian =============================        


  do j = 1, ncomp
    do k = 1, ncomp  
        ! Accumulation
        do i = 1, ngcells 
          coeff_diag(j,k,i) = coeff_diag(j,k,i) + pordt*vol(i)
          if (reaction%neqsorb > 0) then
            coeff_diag(j,k,i) = coeff_diag(j,k,i) + vol(i)/option%tran_dt*(dtotal_sorb_upd(j,k,i))
          endif
        enddo
  
        ! Flux terms
        do i = 2, ngcells-1
          coeff_diag(j,k,i) = coeff_diag(j,k,i) + &
                              pordiff*area(i)/(dm_minus(i+1) + dm_plus(i)) + &
                              pordiff*area(i-1)/(dm_minus(i) + dm_plus(i-1))
          coeff_left(j,k,i) = coeff_left(j,k,i) - &
                              pordiff*area(i-1)/(dm_minus(i) + dm_plus(i-1))
          coeff_right(j,k,i) = coeff_right(j,k,i) - &
                               pordiff*area(i)/(dm_minus(i+1) + dm_plus(i))
        enddo
  
  
        ! Apply boundary conditions
        ! Inner boundary
        coeff_diag(j,k,1) = coeff_diag(j,k,1) + &
                            pordiff*area(1)/(dm_minus(2) + dm_plus(1))
                   
        coeff_right(j,k,1) = coeff_right(j,k,1) - &
                             pordiff*area(1)/(dm_minus(2) + dm_plus(1))
  
        ! Outer boundary -- closest to primary node
        coeff_diag(j,k,ngcells) = coeff_diag(j,k,ngcells) + &
                                  pordiff*area(ngcells-1)/(dm_minus(ngcells) &
                                  + dm_plus(ngcells-1)) + &
                                  pordiff*area(ngcells)/dm_plus(ngcells)
        coeff_left(j,k,ngcells) = coeff_left(j,k,ngcells) - &
                                  pordiff*area(ngcells-1)/(dm_minus(ngcells) + &
                                  dm_plus(ngcells-1)) 

    enddo    
  enddo

!============================= Include dtotal ==================================        
  
  ! Include dtotal (units of kg water/ L water)
  i = 1
  do j = 1, ncomp
    do k = 1, ncomp
      coeff_diag(j,k,i) = coeff_diag(j,k,i)*dtotal(j,k,i) ! m3/s*kg/L
      coeff_right(j,k,i) = coeff_right(j,k,i)*dtotal(j,k,i+1)
    enddo
  enddo
  do i = 2, ngcells-1
    do j = 1, ncomp
      do k = 1, ncomp
        coeff_diag(j,k,i) = coeff_diag(j,k,i)*dtotal(j,k,i) ! m3/s*kg/L
        coeff_left(j,k,i) = coeff_left(j,k,i)*dtotal(j,k,i-1)
        coeff_right(j,k,i) = coeff_right(j,k,i)*dtotal(j,k,i+1)
      enddo
    enddo
  enddo
  i = ngcells
  do j = 1, ncomp
    do k = 1, ncomp
      coeff_diag(j,k,i) = coeff_diag(j,k,i)*dtotal(j,k,i) ! m3/s*kg/L
      coeff_left(j,k,i) = coeff_left(j,k,i)*dtotal(j,k,i-1)
    enddo
  enddo
  
  ! Sorption
  do j = 1, ncomp
    do k = 1, ncomp  
      ! Accumulation
      do i = 1, ngcells 
        if (reaction%neqsorb > 0) then
          coeff_diag(j,k,i) = coeff_diag(j,k,i) + vol(i)/option%tran_dt*(dtotal_sorb_upd(j,k,i))
        endif
      enddo
    enddo
  enddo
    
  
  ! Convert m3/s*kg/L to kg water/s
  coeff_right = coeff_right*1.d3
  coeff_left = coeff_left*1.d3
  coeff_diag = coeff_diag*1.d3
  
!====================== Add reaction contributions =============================        
  
  ! Reaction 
  allocate(material_auxvar)
  call MaterialAuxVarInit(material_auxvar,option)
  do i = 1, ngcells
    res_react = 0.d0
    jac_react = 0.d0
    call RTAuxVarCopy(rt_auxvar,sec_transport_vars%sec_rt_auxvar(i), &
                      option)
    rt_auxvar%pri_molal = conc_upd(:,i) ! in mol/kg
    call RTotalAqueous(rt_auxvar,global_auxvar,reaction,option)
    material_auxvar%porosity = porosity
    material_auxvar%volume = vol(i)
    call RReaction(res_react,jac_react,PETSC_TRUE, &
                   rt_auxvar,global_auxvar,material_auxvar,reaction,option)                     
    do j = 1, ncomp
      res(j+(i-1)*ncomp) = res(j+(i-1)*ncomp) + res_react(j) 
    enddo
    coeff_diag(:,:,i) = coeff_diag(:,:,i) + jac_react  ! in kg water/s
  enddo  
  call MaterialAuxVarStrip(material_auxvar)
  deallocate(material_auxvar)
         
!============================== Forward solve ==================================        
                        
  rhs = -res   
  
  if (reaction%use_log_formulation) then
  ! scale the jacobian by concentrations
    i = 1
    do k = 1, ncomp
      coeff_diag(:,k,i) = coeff_diag(:,k,i)*conc_upd(k,i) ! m3/s*kg/L
      coeff_right(:,k,i) = coeff_right(:,k,i)*conc_upd(k,i+1)
    enddo
    do i = 2, ngcells-1
      do k = 1, ncomp
        coeff_diag(:,k,i) = coeff_diag(:,k,i)*conc_upd(k,i) ! m3/s*kg/L
        coeff_left(:,k,i) = coeff_left(:,k,i)*conc_upd(k,i-1)
        coeff_right(:,k,i) = coeff_right(:,k,i)*conc_upd(k,i+1)
      enddo
    enddo
    i = ngcells
      do k = 1, ncomp
        coeff_diag(:,k,i) = coeff_diag(:,k,i)*conc_upd(k,i) ! m3/s*kg/L
        coeff_left(:,k,i) = coeff_left(:,k,i)*conc_upd(k,i-1)
      enddo
  endif 
  
  ! First do an LU decomposition for calculating D_M matrix
  coeff_diag_dm = coeff_diag
  coeff_left_dm = coeff_left
  coeff_right_dm = coeff_right
  
  select case (option%secondary_continuum_solver)
    case(1) 
      do i = 2, ngcells
        coeff_left_dm(:,:,i-1) = coeff_left_dm(:,:,i)
      enddo
      coeff_left_dm(:,:,ngcells) = 0.d0
      call bl3dfac(ngcells,ncomp,coeff_right_dm,coeff_diag_dm,coeff_left_dm,pivot)  
    case(2)
      call decbt(ncomp,ngcells,ncomp,coeff_diag_dm,coeff_right_dm,coeff_left_dm,pivot,ier)
      if (ier /= 0) then
        print *,'error in matrix decbt: ier = ',ier
        stop
      endif
    case(3)
      ! Thomas algorithm for tridiagonal system
      ! Forward elimination
      if (ncomp /= 1) then
        option%io_buffer = 'THOMAS algorithm can be used only with single '// &
                           'component chemistry'
        call printErrMsg(option)
      endif
      do i = 2, ngcells
        m = coeff_left_dm(ncomp,ncomp,i)/coeff_diag_dm(ncomp,ncomp,i-1)
        coeff_diag_dm(ncomp,ncomp,i) = coeff_diag_dm(ncomp,ncomp,i) - &
                                    m*coeff_right_dm(ncomp,ncomp,i-1)
      enddo        
    case default
      option%io_buffer = 'SECONDARY_CONTINUUM_SOLVER can be only ' // &
                         'HINDMARSH or KEARST. For single component'// &
                         'chemistry THOMAS can be used.'
      call printErrMsg(option)  
  end select
  
  ! Set the values of D_M matrix and create identity matrix of size ncomp x ncomp  
  do i = 1, ncomp
    do j = 1, ncomp
      D_M(i,j) = coeff_diag_dm(i,j,ngcells)
      if (j == i) then
        identity(i,j) = 1.d0
      else
        identity(i,j) = 0.d0
      endif
    enddo
  enddo
  
  ! Find the inverse of D_M
  call ludcmp(D_M,ncomp,indx,d) 
  do j = 1, ncomp
    call lubksb(D_M,ncomp,indx,identity(1,j))
  enddo  
  inv_D_M = identity      
  
  if (option%numerical_derivatives_multi_coupling) then  
    ! Store the coeffs for numerical jacobian
    coeff_diag_copy = coeff_diag
    coeff_left_copy = coeff_left
    coeff_right_copy = coeff_right
  endif

  select case (option%secondary_continuum_solver)
    case(1) 
      do i = 2, ngcells
        coeff_left(:,:,i-1) = coeff_left(:,:,i)
      enddo
      coeff_left(:,:,ngcells) = 0.d0
      call bl3dfac(ngcells,ncomp,coeff_right,coeff_diag,coeff_left,pivot)  
      call bl3dsolf(ngcells,ncomp,coeff_right,coeff_diag,coeff_left,pivot, &
                    ONE_INTEGER,rhs)
    case(2)
      call decbt(ncomp,ngcells,ncomp,coeff_diag,coeff_right,coeff_left, &
                 pivot,ier)
      if (ier /= 0) then
        print *,'error in matrix decbt: ier = ',ier
        stop
      endif
      call solbtf(ncomp,ngcells,ncomp,coeff_diag,coeff_right,coeff_left, &
                  pivot,rhs)
    case(3)
      ! Thomas algorithm for tridiagonal system
      ! Forward elimination
      if (ncomp /= 1) then
        option%io_buffer = 'THOMAS algorithm can be used only with single '// &
                           'component chemistry'
        call printErrMsg(option)
      endif
      do i = 2, ngcells
        m = coeff_left(ncomp,ncomp,i)/coeff_diag(ncomp,ncomp,i-1)
        coeff_diag(ncomp,ncomp,i) = coeff_diag(ncomp,ncomp,i) - &
                                    m*coeff_right(ncomp,ncomp,i-1)
        rhs(i) = rhs(i) - m*rhs(i-1)
      enddo        
      rhs(ngcells) = rhs(ngcells)/coeff_diag(ncomp,ncomp,ngcells)
    case default
      option%io_buffer = 'SECONDARY_CONTINUUM_SOLVER can be only ' // &
                         'HINDMARSH or KEARST. For single component'// &
                         'chemistry THOMAS can be used.'
      call printErrMsg(option)  
  end select
    
  ! Update the secondary concentrations
  do i = 1, ncomp
    if (reaction%use_log_formulation) then
      ! convert log concentration to concentration
      rhs(i+(ngcells-1)*ncomp) = dsign(1.d0,rhs(i+(ngcells-1)*ncomp))* &
        min(dabs(rhs(i+(ngcells-1)*ncomp)),reaction%max_dlnC)
      conc_current_M(i) = conc_upd(i,ngcells)*exp(rhs(i+(ngcells-1)*ncomp))
    else
      conc_current_M(i) = conc_upd(i,ngcells) + rhs(i+(ngcells-1)*ncomp)
    endif
  enddo

  ! Update the secondary continuum totals at the outer matrix node
  call RTAuxVarCopy(rt_auxvar,sec_transport_vars%sec_rt_auxvar(ngcells), &
                    option)
  rt_auxvar%pri_molal = conc_current_M ! in mol/kg
  call RTotalAqueous(rt_auxvar,global_auxvar,reaction,option)
  total_current_M = rt_auxvar%total(:,1)
  if (reaction%neqcplx > 0) sec_sec_molal_M = rt_auxvar%sec_molal
  call RTAuxVarStrip(rt_auxvar)
  

  b_m = pordiff/dm_plus(ngcells)*area(ngcells)*inv_D_M ! in m3/kg
  b_m = b_m*1.d3 ! in L/kg  For log formulation, L/mol
  
  dCsec_dCprim = b_m*dtotal_prim
      
  ! Calculate the dervative of outer matrix node total with respect to the 
  ! primary node concentration
  
  if (reaction%use_log_formulation) then ! log formulation
    do j = 1, ncomp
      do l = 1, ncomp
        dPsisec_dCprim(j,l) = dCsec_dCprim(j,l)*conc_current_M(j)
      enddo
    enddo
   
    if (reaction%neqcplx > 0) then
      do icplx = 1, reaction%neqcplx
        ncompeq = reaction%eqcplxspecid(0,icplx)
        do j = 1, ncompeq
          jcomp = reaction%eqcplxspecid(j,icplx)
          do l = 1, ncompeq
            lcomp = reaction%eqcplxspecid(l,icplx)
            do k = 1, ncompeq
              kcomp = reaction%eqcplxspecid(k,icplx)
              dPsisec_dCprim(jcomp,lcomp) = dPsisec_dCprim(jcomp,lcomp) + &
                                            reaction%eqcplxstoich(j,icplx)* &
                                            reaction%eqcplxstoich(k,icplx)* &
                                            dCsec_dCprim(kcomp,lcomp)* &
                                            sec_sec_molal_M(icplx)
            enddo
          enddo      
        enddo
      enddo
    endif
   
  else   ! linear case  
  
    dPsisec_dCprim = dCsec_dCprim       ! dimensionless
  
    if (reaction%neqcplx > 0) then
      do icplx = 1, reaction%neqcplx
        ncompeq = reaction%eqcplxspecid(0,icplx)
        do j = 1, ncompeq
          jcomp = reaction%eqcplxspecid(j,icplx)
          do l = 1, ncompeq
            lcomp = reaction%eqcplxspecid(l,icplx)
            do k = 1, ncompeq
              kcomp = reaction%eqcplxspecid(k,icplx)
              dPsisec_dCprim(jcomp,lcomp) = dPsisec_dCprim(jcomp,lcomp) + &
                                            reaction%eqcplxstoich(j,icplx)* &
                                            reaction%eqcplxstoich(k,icplx)* &
                                            dCsec_dCprim(kcomp,lcomp)* &
                                            sec_sec_molal_M(icplx)/ &
                                            conc_current_M(kcomp)
            enddo
          enddo      
        enddo
      enddo
    endif
  
  endif
  
  dPsisec_dCprim = dPsisec_dCprim*global_auxvar%den_kg(1)*1.d-3 ! in kg/L
            
  ! Calculate the coupling term
  res_transport = pordiff/dm_plus(ngcells)*area_fm* &
                  (total_current_M - total_primary_node)*prim_vol*1.d3 ! in mol/s
                                           
  ! Calculate the jacobian contribution due to coupling term
  sec_jac = area_fm*pordiff/dm_plus(ngcells)*(dPsisec_dCprim - dtotal_prim)* &
            prim_vol*1.d3 ! in kg water/s
      
  ! Store the contribution to the primary jacobian term
  sec_transport_vars%sec_jac = sec_jac 
  sec_transport_vars%sec_jac_update = PETSC_TRUE
  
  ! Store the coefficients from LU decomposition of the block tridiagonal
  ! sytem. These will be called later to perform backsolve to the get the
  ! updated secondary continuum concentrations at the end of the timestep
  sec_transport_vars%cxm = coeff_left
  sec_transport_vars%cxp = coeff_right
  sec_transport_vars%cdl = coeff_diag
  
  ! Store the solution of the forward solve
  sec_transport_vars%r = rhs
  
!============== Numerical jacobian for coupling term ===========================


  if (option%numerical_derivatives_multi_coupling) then

    call RTAuxVarInit(rt_auxvar,reaction,option)
    conc_prim = auxvar%pri_molal
    conc_prim_pert = conc_prim
  
    do l = 1, ncomp
      
      conc_prim_pert = conc_prim
      pert = conc_prim(l)*perturbation_tolerance
      conc_prim_pert(l) = conc_prim_pert(l) + pert
  
      res = 0.d0
      rhs = 0.d0
    
      coeff_diag_pert = coeff_diag_copy
      coeff_left_pert = coeff_left_copy
      coeff_right_pert = coeff_right_copy

      call RTAuxVarCopy(rt_auxvar,auxvar,option)
      rt_auxvar%pri_molal = conc_prim_pert ! in mol/kg
      call RTotalAqueous(rt_auxvar,global_auxvar,reaction,option)
      total_primary_node_pert = rt_auxvar%total(:,1)
                          
!================ Calculate the secondary residual =============================        

      do j = 1, ncomp
      
        ! Accumulation
        do i = 1, ngcells
          n = j + (i-1)*ncomp
          res(n) = pordt*(total_upd(j,i) - total_prev(j,i))*vol(i)    ! in mol/L*m3/s
        enddo
  
        ! Flux terms
        do i = 2, ngcells - 1
          n = j + (i-1)*ncomp
          res(n) = res(n) - pordiff*area(i)/(dm_minus(i+1) + dm_plus(i))* &
                            (total_upd(j,i+1) - total_upd(j,i))
          res(n) = res(n) + pordiff*area(i-1)/(dm_minus(i) + dm_plus(i-1))* &
                            (total_upd(j,i) - total_upd(j,i-1))                      
        enddo
         
              
        ! Apply boundary conditions
        ! Inner boundary
        res(j) = res(j) - pordiff*area(1)/(dm_minus(2) + dm_plus(1))* &
                          (total_upd(j,2) - total_upd(j,1))
                                      
        ! Outer boundary
        res(j+(ngcells-1)*ncomp) = res(j+(ngcells-1)*ncomp) - &
                                   pordiff*area(ngcells)/dm_plus(ngcells)* &
                                   (total_primary_node_pert(j) -  &
                                   total_upd(j,ngcells))
        res(j+(ngcells-1)*ncomp) = res(j+(ngcells-1)*ncomp) + &
                                   pordiff*area(ngcells-1)/(dm_minus(ngcells) &
                                   + dm_plus(ngcells-1))*(total_upd(j,ngcells) - &
                                   total_upd(j,ngcells-1))  
                               
      enddo
                         
      res = res*1.d3 ! Convert mol/L*m3/s to mol/s                                                    
                                                                                                                  
!============================== Forward solve ==================================        
                        
      rhs = -res   
           
    select case (option%secondary_continuum_solver)
      case(1) 
        call bl3dfac(ngcells,ncomp,coeff_right_pert,coeff_diag_pert, &
                      coeff_left_pert,pivot)  
        call bl3dsolf(ngcells,ncomp,coeff_right_pert,coeff_diag_pert, &
                       coeff_left_pert,pivot,ONE_INTEGER,rhs)
      case(2)
        call decbt(ncomp,ngcells,ncomp,coeff_diag_pert,coeff_right_pert, &
                    coeff_left_pert,pivot,ier)
        if (ier /= 0) then
          print *,'error in matrix decbt: ier = ',ier
          stop
        endif
        call solbtf(ncomp,ngcells,ncomp,coeff_diag_pert,coeff_right_pert, &
                     coeff_left_pert,pivot,rhs)
      case(3)
        ! Thomas algorithm for tridiagonal system
        ! Forward elimination
        if (ncomp /= 1) then
          option%io_buffer = 'THOMAS algorithm can be used only with '// &
                             'single component chemistry'
          call printErrMsg(option)
        endif
        do i = 2, ngcells
          m = coeff_left_pert(ncomp,ncomp,i)/coeff_diag_pert(ncomp,ncomp,i-1)
          coeff_diag_pert(ncomp,ncomp,i) = coeff_diag_pert(ncomp,ncomp,i) - &
                                      m*coeff_right_pert(ncomp,ncomp,i-1)
          rhs(i) = rhs(i) - m*rhs(i-1)
        enddo        
        rhs(ngcells) = rhs(ngcells)/coeff_diag(ncomp,ncomp,ngcells)
      case default
        option%io_buffer = 'SECONDARY_CONTINUUM_SOLVER can be only ' // &
                           'HINDMARSH or KEARST. For single component'// &
                           'chemistry THOMAS can be used.'
        call printErrMsg(option)  
      end select      
    
      ! Update the secondary concentrations
      do i = 1, ncomp
        if (reaction%use_log_formulation) then
          ! convert log concentration to concentration
          rhs(i+(ngcells-1)*ncomp) = dsign(1.d0,rhs(i+(ngcells-1)*ncomp))* &
            min(dabs(rhs(i+(ngcells-1)*ncomp)),reaction%max_dlnC)
          conc_current_M_pert(i) = conc_upd(i,ngcells)* &
                                     exp(rhs(i+(ngcells-1)*ncomp))
        else
          conc_current_M_pert(i) = conc_upd(i,ngcells) + &
                                     rhs(i+(ngcells-1)*ncomp)
        endif
      enddo

      ! Update the secondary continuum totals at the outer matrix node
      call RTAuxVarCopy(rt_auxvar,sec_transport_vars%sec_rt_auxvar(ngcells), &
                        option)
      rt_auxvar%pri_molal = conc_current_M_pert ! in mol/kg
      call RTotalAqueous(rt_auxvar,global_auxvar,reaction,option)
      total_current_M_pert = rt_auxvar%total(:,1)
             
      ! Calculate the coupling term
      res_transport_pert = pordiff/dm_plus(ngcells)*area_fm* &
                           (total_current_M_pert - total_primary_node_pert)* &
                           prim_vol*1.d3 ! in mol/s
  
      dtotal_prim_num(:,l) = (total_primary_node_pert(:) - &
                               total_primary_node(:))/pert
  
      dPsisec_dCprim_num(:,l) = (total_current_M_pert(:) - &
                                  total_current_M(:))/pert
  
      sec_jac_num(:,l) = (res_transport_pert(:) - res_transport(:))/pert
    
    enddo    

    call RTAuxVarStrip(rt_auxvar)
    sec_transport_vars%sec_jac = sec_jac_num 

  endif
  

end subroutine SecondaryRTResJacMulti

! ************************************************************************** !

subroutine SecondaryRTUpdateIterate(line_search,P0,dP,P1,dX_changed, &
                                    X1_changed,realization,ierr)
  ! 
  ! Checks update after the update is done
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 02/22/13
  ! 

  use Realization_Subsurface_class
  use Option_module
  use Grid_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
 
  implicit none
  
  SNESLineSearch :: line_search
  Vec :: P0
  Vec :: dP
  Vec :: P1
  class(realization_subsurface_type) :: realization
  ! ignore changed flag for now.
  PetscBool :: dX_changed
  PetscBool :: X1_changed
  
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  type(sec_transport_type), pointer :: rt_sec_transport_vars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(reaction_type), pointer :: reaction
  PetscInt :: local_id, ghosted_id
  PetscReal :: sec_diffusion_coefficient
  PetscReal :: sec_porosity
  PetscErrorCode :: ierr
  PetscReal :: inf_norm_sec
  PetscReal :: max_inf_norm_sec
  
  option => realization%option
  grid => realization%patch%grid
  rt_auxvars => realization%patch%aux%RT%auxvars
  global_auxvars => realization%patch%aux%Global%auxvars
  reaction => realization%reaction
  if (option%use_mc) then
    rt_sec_transport_vars => realization%patch%aux%SC_RT%sec_transport_vars
  endif  
  
  dX_changed = PETSC_FALSE
  X1_changed = PETSC_FALSE
  
  max_inf_norm_sec = 0.d0
  
  if (option%use_mc) then
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (realization%patch%imat(ghosted_id) <= 0) cycle
      sec_diffusion_coefficient = realization%patch% &
                                  material_property_array(1)%ptr% &
                                  secondary_continuum_diff_coeff
      sec_porosity = realization%patch%material_property_array(1)%ptr% &
                    secondary_continuum_porosity

      call SecondaryRTAuxVarComputeMulti(&
                                    rt_sec_transport_vars(local_id), &
                                    reaction, &
                                    option)              
 
      call SecondaryRTCheckResidual(rt_sec_transport_vars(local_id), &
                                    rt_auxvars(ghosted_id), &
                                    global_auxvars(ghosted_id), &
                                    reaction,sec_diffusion_coefficient, &
                                    sec_porosity,option,inf_norm_sec)
                                      
      max_inf_norm_sec = max(max_inf_norm_sec,inf_norm_sec)                                                                   
    enddo 
    call MPI_Allreduce(max_inf_norm_sec,option%infnorm_res_sec,ONE_INTEGER_MPI, &
                       MPI_DOUBLE_PRECISION, &
                       MPI_MAX,option%mycomm,ierr)
  endif
  
      
end subroutine SecondaryRTUpdateIterate

! ************************************************************************** !

subroutine SecondaryRTUpdateEquilState(sec_transport_vars,global_auxvars, &
                                       reaction,option) 
  ! 
  ! Updates the equilibrium secondary continuum
  ! variables
  ! at the end of time step
  ! 
  ! Author: Satish Karra, LANL; Glenn Hammond (modification)
  ! Date: 02/22/13; 06/27/13
  ! 
                                     

  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Reaction_module
  use Global_Aux_module
 
  implicit none
  

  type(option_type), pointer :: option
  type(sec_transport_type) :: sec_transport_vars
  type(global_auxvar_type) :: global_auxvars
  type(reaction_type), pointer :: reaction
  PetscInt :: ngcells,ncomp
  PetscInt :: i,j
  
  ngcells = sec_transport_vars%ncells
  ncomp = reaction%naqcomp
                   
  do j = 1, ncomp
    do i = 1, ngcells
      sec_transport_vars%sec_rt_auxvar(i)%pri_molal(j) = sec_transport_vars%&
        updated_conc(j,i)
    enddo
  enddo
    
  do i = 1, ngcells
    call RTotalAqueous(sec_transport_vars%sec_rt_auxvar(i),global_auxvars, &
                reaction,option)
  enddo
  
end subroutine SecondaryRTUpdateEquilState

! ************************************************************************** !

subroutine SecondaryRTUpdateKineticState(sec_transport_vars,global_auxvars, &
                                         reaction,porosity,option) 
  ! 
  ! Updates the kinetic secondary continuum
  ! variables at the end of time step
  ! 
  ! Author: Satish Karra, LANL; Glenn Hammond (modification)
  ! Date: 02/22/13; 06/27/13
  ! 
                                     

  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Reaction_module
  use Global_Aux_module
  use Material_Aux_class
 
  implicit none
  

  type(option_type), pointer :: option
  type(sec_transport_type) :: sec_transport_vars
  type(global_auxvar_type) :: global_auxvars
  type(reaction_type), pointer :: reaction
  PetscReal :: porosity
  PetscInt :: ngcells
  PetscReal :: vol(sec_transport_vars%ncells)
  PetscReal :: res_react(reaction%naqcomp)
  PetscReal :: jac_react(reaction%naqcomp,reaction%naqcomp)
  PetscInt :: i,j
  class(material_auxvar_type), allocatable :: material_auxvar
  
  ngcells = sec_transport_vars%ncells
  vol = sec_transport_vars%vol     
                   
  res_react = 0.d0
  jac_react = 0.d0 ! These are not used anyway
  allocate(material_auxvar)
  call MaterialAuxVarInit(material_auxvar,option)
  do i = 1, ngcells
    material_auxvar%porosity = porosity
    material_auxvar%volume = vol(i)
    call RReaction(res_react,jac_react,PETSC_FALSE, &
                   sec_transport_vars%sec_rt_auxvar(i), &
                   global_auxvars,material_auxvar,reaction,option)
  enddo
  call MaterialAuxVarStrip(material_auxvar)
  deallocate(material_auxvar)
  
  if (reaction%mineral%nkinmnrl > 0) then
    do i = 1, ngcells
      do j = 1, reaction%mineral%nkinmnrl
        sec_transport_vars%sec_rt_auxvar(i)%mnrl_volfrac(j) = &
          sec_transport_vars%sec_rt_auxvar(i)%mnrl_volfrac(j) + &
          sec_transport_vars%sec_rt_auxvar(i)%mnrl_rate(j)* &
          reaction%mineral%kinmnrl_molar_vol(j)* &
          option%tran_dt
          if (sec_transport_vars%sec_rt_auxvar(i)%mnrl_volfrac(j) < 0.d0) &
            sec_transport_vars%sec_rt_auxvar(i)%mnrl_volfrac(j) = 0.d0
      enddo
    enddo
  endif

  
end subroutine SecondaryRTUpdateKineticState

! ************************************************************************** !

subroutine SecondaryRTCheckResidual(sec_transport_vars,auxvar, &
                                    global_auxvar, &
                                    reaction,diffusion_coefficient, &
                                    porosity,option,inf_norm_sec)
  ! 
  ! The residual of the secondary domain are checked
  ! to ensure convergence
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 1/31/13
  ! 
                                    
  use Option_module 
  use Global_Aux_module
  use Block_Solve_module
  use Block_Tridiag_module
  use Utility_module
  use Reaction_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Material_Aux_class

  implicit none
  
  type(sec_transport_type) :: sec_transport_vars
  type(reactive_transport_auxvar_type) :: auxvar
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(reaction_type), pointer :: reaction
  type(option_type) :: option
  
  PetscReal :: res(sec_transport_vars%ncells*reaction%naqcomp)
  PetscReal :: conc_upd(reaction%naqcomp,sec_transport_vars%ncells) 
  PetscReal :: total_upd(reaction%naqcomp,sec_transport_vars%ncells)
  PetscReal :: total_prev(reaction%naqcomp,sec_transport_vars%ncells)
  PetscReal :: total_primary_node(reaction%naqcomp)
  PetscReal :: area(sec_transport_vars%ncells)
  PetscReal :: vol(sec_transport_vars%ncells)
  PetscReal :: dm_plus(sec_transport_vars%ncells)
  PetscReal :: dm_minus(sec_transport_vars%ncells)
  PetscReal :: res_react(reaction%naqcomp)
  PetscReal :: jac_react(reaction%naqcomp,reaction%naqcomp)
  PetscInt :: i, j, k, n
  PetscInt :: ngcells, ncomp
  PetscReal :: area_fm
  PetscReal :: diffusion_coefficient
  PetscReal :: porosity
  PetscReal :: arrhenius_factor
  PetscReal :: pordt, pordiff
  PetscReal :: inf_norm_sec
  class(material_auxvar_type), allocatable :: material_auxvar

  PetscReal :: total_sorb_upd(reaction%naqcomp,sec_transport_vars%ncells) 
  PetscReal :: total_sorb_prev(reaction%naqcomp,sec_transport_vars%ncells)
  
  ngcells = sec_transport_vars%ncells
  area = sec_transport_vars%area
  vol = sec_transport_vars%vol          
  dm_plus = sec_transport_vars%dm_plus
  dm_minus = sec_transport_vars%dm_minus
  area_fm = sec_transport_vars%interfacial_area
  ncomp = reaction%naqcomp

  do j = 1, ncomp
    do i = 1, ngcells
      total_prev(j,i) = sec_transport_vars%sec_rt_auxvar(i)%total(j,1)
      if (reaction%neqsorb > 0) then
        total_sorb_prev(j,i) = sec_transport_vars%sec_rt_auxvar(i)%total_sorb_eq(j)
      endif
    enddo
  enddo
  conc_upd = sec_transport_vars%updated_conc
    
  ! Note that sec_transport_vars%sec_rt_auxvar(i)%pri_molal(j) units are in mol/kg
  ! Need to convert to mol/L since the units of total. in the Thomas 
  ! algorithm are in mol/L 
  
  res = 0.d0
  
  total_primary_node = auxvar%total(:,1)                         ! in mol/L 
  pordt = porosity/option%tran_dt
  pordiff = porosity*diffusion_coefficient

  call RTAuxVarInit(rt_auxvar,reaction,option)
  do i = 1, ngcells
    call RTAuxVarCopy(rt_auxvar,sec_transport_vars%sec_rt_auxvar(i),option)
    rt_auxvar%pri_molal = conc_upd(:,i)
    call RTotalAqueous(rt_auxvar,global_auxvar,reaction,option)
    if (reaction%neqsorb > 0) then
      call SecondaryRTotalSorb(rt_auxvar,global_auxvar,material_auxvar,reaction,option)
    endif
    total_upd(:,i) = rt_auxvar%total(:,1)
    if (reaction%neqsorb > 0) then 
      total_sorb_upd(:,i) = rt_auxvar%total_sorb_eq(:)
    endif
  enddo
                                    
!================ Calculate the secondary residual =============================        

  do j = 1, ncomp
      
    ! Accumulation
    do i = 1, ngcells
      n = j + (i-1)*ncomp
      res(n) = pordt*(total_upd(j,i) - total_prev(j,i))*vol(i)    ! in mol/L*m3/s
      if (reaction%neqsorb > 0) then 
        res(n) = res(n) + vol(i)/option%tran_dt*(total_sorb_upd(j,i) - total_sorb_prev(j,i))
      endif 
    enddo
  
    ! Flux terms
    do i = 2, ngcells - 1
      n = j + (i-1)*ncomp
      res(n) = res(n) - pordiff*area(i)/(dm_minus(i+1) + dm_plus(i))* &
                        (total_upd(j,i+1) - total_upd(j,i))
      res(n) = res(n) + pordiff*area(i-1)/(dm_minus(i) + dm_plus(i-1))* &
                        (total_upd(j,i) - total_upd(j,i-1))                      
    enddo
         
              
    ! Apply boundary conditions
    ! Inner boundary
    res(j) = res(j) - pordiff*area(1)/(dm_minus(2) + dm_plus(1))* &
                      (total_upd(j,2) - total_upd(j,1))
                                      
    ! Outer boundary
    res(j+(ngcells-1)*ncomp) = res(j+(ngcells-1)*ncomp) - &
                               pordiff*area(ngcells)/dm_plus(ngcells)* &
                               (total_primary_node(j) - total_upd(j,ngcells))
    res(j+(ngcells-1)*ncomp) = res(j+(ngcells-1)*ncomp) + &
                               pordiff*area(ngcells-1)/(dm_minus(ngcells) &
                               + dm_plus(ngcells-1))*(total_upd(j,ngcells) - &
                               total_upd(j,ngcells-1))  
                               
  enddo
                         
  res = res*1.d3 ! Convert mol/L*m3/s to mol/s                                             
                                    
                                    
!====================== Add reaction contributions =============================        
  
  ! Reaction 
  allocate(material_auxvar)
  call MaterialAuxVarInit(material_auxvar,option)
  do i = 1, ngcells
    res_react = 0.d0
    jac_react = 0.d0
    call RTAuxVarCopy(rt_auxvar,sec_transport_vars%sec_rt_auxvar(i), &
                      option)
    rt_auxvar%pri_molal = conc_upd(:,i) ! in mol/kg
    call RTotalAqueous(rt_auxvar,global_auxvar,reaction,option)
    material_auxvar%porosity = porosity
    material_auxvar%volume = vol(i)
    call RReaction(res_react,jac_react,PETSC_FALSE, &
                   rt_auxvar,global_auxvar,material_auxvar,reaction,option)                     
    do j = 1, ncomp
      res(j+(i-1)*ncomp) = res(j+(i-1)*ncomp) + res_react(j) 
    enddo
  enddo  
  call MaterialAuxVarStrip(material_auxvar)         
  deallocate(material_auxvar)
  
 ! Need to decide how to scale the residual with volumes
  do i = 1, ngcells
    do j = 1, ncomp
      res(j+(i-1)*ncomp) = res(j+(i-1)*ncomp)/vol(i)
    enddo
  enddo
    
  inf_norm_sec = maxval(abs(res))  
  call RTAuxVarStrip(rt_auxvar)  
                                                                      
end subroutine SecondaryRTCheckResidual                                    

! ************************************************************************** !

subroutine SecondaryRTAuxVarComputeMulti(sec_transport_vars,reaction, &
                                         option)
  ! 
  ! Updates the secondary continuum
  ! concentrations at end of each time step for multicomponent system
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 2/1/13
  ! 
                               
                            
  use Option_module 
  use Reaction_Aux_module
  use Reaction_module
  use Reactive_Transport_Aux_module
  use Block_Solve_module
  use Block_Tridiag_module
  use Utility_module
  

  implicit none
  
  type(sec_transport_type) :: sec_transport_vars
  type(reaction_type), pointer :: reaction
  type(option_type) :: option
  PetscReal :: coeff_left(reaction%naqcomp,reaction%naqcomp, &
                 sec_transport_vars%ncells)
  PetscReal :: coeff_diag(reaction%naqcomp,reaction%naqcomp, &
                 sec_transport_vars%ncells)
  PetscReal :: coeff_right(reaction%naqcomp,reaction%naqcomp, &
                 sec_transport_vars%ncells)
  PetscReal :: rhs(sec_transport_vars%ncells*reaction%naqcomp)
  PetscReal :: conc_upd(reaction%naqcomp,sec_transport_vars%ncells) 
  PetscInt :: i, j, n
  PetscInt :: ngcells, ncomp
  PetscInt :: pivot(reaction%naqcomp,sec_transport_vars%ncells)
  PetscInt :: indx(reaction%naqcomp)
  PetscInt :: d
    
  ngcells = sec_transport_vars%ncells
  ncomp = reaction%naqcomp
  ! Note that sec_transport_vars%sec_conc units are in mol/kg
  ! Need to convert to mol/L since the units of conc. in the Thomas 
  ! algorithm are in mol/L 
  
  coeff_left = 0.d0
  coeff_diag = 0.d0
  coeff_right = 0.d0
  rhs = 0.d0

  conc_upd = sec_transport_vars%updated_conc
           
  ! Use the stored coefficient matrices from LU decomposition of the
  ! block triagonal sytem
  coeff_left = sec_transport_vars%cxm
  coeff_right = sec_transport_vars%cxp
  coeff_diag = sec_transport_vars%cdl
  rhs = sec_transport_vars%r
    
  select case (option%secondary_continuum_solver)
    case(1) 
      call bl3dsolb(ngcells,ncomp,coeff_right,coeff_diag,coeff_left,pivot, &
                    ONE_INTEGER,rhs)
    case(2)
      call solbtb(ncomp,ngcells,ncomp,coeff_diag,coeff_right,coeff_left, &
                  pivot,rhs)
    case(3)
      do i = ngcells-1, 1, -1
        rhs(i) = (rhs(i) - coeff_right(ncomp,ncomp,i)*rhs(i+1))/ &
                             coeff_diag(ncomp,ncomp,i)
      enddo
    case default
      option%io_buffer = 'SECONDARY_CONTINUUM_SOLVER can be only ' // &
                         'HINDMARSH or KEARST. For single component'// &
                         'chemistry THOMAS can be used.'
      call printErrMsg(option)  
  end select  
  
  do j = 1, ncomp
    do i = 1, ngcells
      n = j + (i - 1)*ncomp
      if (reaction%use_log_formulation) then 
        ! convert log concentration to concentration
        rhs(n) = dsign(1.d0,rhs(n))*min(dabs(rhs(n)),reaction%max_dlnC) 
        conc_upd(j,i) = exp(rhs(n))*conc_upd(j,i)
      else
        conc_upd(j,i) = rhs(n) + conc_upd(j,i)
      endif
      if (conc_upd(j,i) < 0.d0) conc_upd(j,i) = 1.d-8
    enddo
  enddo
  
  sec_transport_vars%updated_conc = conc_upd
    
end subroutine SecondaryRTAuxVarComputeMulti

! ************************************************************************** !

subroutine THCSecHeatAuxVarCompute(sec_heat_vars,global_auxvar, &
                                   therm_conductivity,dencpr, &
                                   option)
  ! 
  ! Computes secondary auxillary variables for each
  ! grid cell for heat transfer only
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/5/12
  ! 

  use Option_module 
  use Global_Aux_module
  
  implicit none
  
  type(sec_heat_type) :: sec_heat_vars
  type(global_auxvar_type) :: global_auxvar
  type(option_type) :: option
  PetscReal :: coeff_left(sec_heat_vars%ncells)
  PetscReal :: coeff_diag(sec_heat_vars%ncells)
  PetscReal :: coeff_right(sec_heat_vars%ncells)
  PetscReal :: rhs(sec_heat_vars%ncells)
  PetscReal :: sec_temp(sec_heat_vars%ncells)
  PetscReal :: area(sec_heat_vars%ncells)
  PetscReal :: vol(sec_heat_vars%ncells)
  PetscReal :: dm_plus(sec_heat_vars%ncells)
  PetscReal :: dm_minus(sec_heat_vars%ncells)
  PetscInt :: i, ngcells
  PetscReal :: area_fm
  PetscReal :: alpha, therm_conductivity, dencpr
  PetscReal :: temp_primary_node
  PetscReal :: m
  
  ngcells = sec_heat_vars%ncells
  area = sec_heat_vars%area
  vol = sec_heat_vars%vol
  dm_plus = sec_heat_vars%dm_plus
  dm_minus = sec_heat_vars%dm_minus
  area_fm = sec_heat_vars%interfacial_area
  temp_primary_node = global_auxvar%temp
  
  coeff_left = 0.d0
  coeff_diag = 0.d0
  coeff_right = 0.d0
  rhs = 0.d0
  sec_temp = 0.d0
  
  alpha = option%flow_dt*therm_conductivity/dencpr

  
  ! Setting the coefficients
  do i = 2, ngcells-1
    coeff_left(i) = -alpha*area(i-1)/((dm_minus(i) + dm_plus(i-1))*vol(i))
    coeff_diag(i) = alpha*area(i-1)/((dm_minus(i) + dm_plus(i-1))*vol(i)) + &
                    alpha*area(i)/((dm_minus(i+1) + dm_plus(i))*vol(i)) + 1.d0
    coeff_right(i) = -alpha*area(i)/((dm_minus(i+1) + dm_plus(i))*vol(i))
  enddo
  
  coeff_diag(1) = alpha*area(1)/((dm_minus(2) + dm_plus(1))*vol(1)) + 1.d0
  coeff_right(1) = -alpha*area(1)/((dm_minus(2) + dm_plus(1))*vol(1))
  
  coeff_left(ngcells) = -alpha*area(ngcells-1)/ &
                       ((dm_minus(ngcells) + dm_plus(ngcells-1))*vol(ngcells))
  coeff_diag(ngcells) = alpha*area(ngcells-1)/ &
                       ((dm_minus(ngcells) + dm_plus(ngcells-1))*vol(ngcells)) &
                       + alpha*area(ngcells)/(dm_plus(ngcells)*vol(ngcells)) &
                       + 1.d0

                        
  rhs = sec_heat_vars%sec_temp  ! secondary continuum values from previous time step
  rhs(ngcells) = rhs(ngcells) + & 
                 alpha*area(ngcells)/(dm_plus(ngcells)*vol(ngcells))* &
                 temp_primary_node
                
  ! Thomas algorithm for tridiagonal system
  ! Forward elimination
  do i = 2, ngcells
    m = coeff_left(i)/coeff_diag(i-1)
    coeff_diag(i) = coeff_diag(i) - m*coeff_right(i-1)
    rhs(i) = rhs(i) - m*rhs(i-1)
  enddo

  ! Back substitution
  ! Calculate temperature in the secondary continuum
  sec_temp(ngcells) = rhs(ngcells)/coeff_diag(ngcells)
  do i = ngcells-1, 1, -1
    sec_temp(i) = (rhs(i) - coeff_right(i)*sec_temp(i+1))/coeff_diag(i)
  enddo
  
  sec_heat_vars%sec_temp = sec_temp
            
end subroutine THCSecHeatAuxVarCompute

! ************************************************************************** !

subroutine THSecHeatAuxVarCompute(sec_heat_vars,global_auxvar, &
                                   therm_conductivity,dencpr, &
                                   option)
  ! 
  ! Computes secondary auxillary variables for each
  ! grid cell for heat transfer only
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/5/12
  ! 

  use Option_module 
  use Global_Aux_module
  
  implicit none
  
  type(sec_heat_type) :: sec_heat_vars
  type(global_auxvar_type) :: global_auxvar
  type(option_type) :: option
  PetscReal :: coeff_left(sec_heat_vars%ncells)
  PetscReal :: coeff_diag(sec_heat_vars%ncells)
  PetscReal :: coeff_right(sec_heat_vars%ncells)
  PetscReal :: rhs(sec_heat_vars%ncells)
  PetscReal :: sec_temp(sec_heat_vars%ncells)
  PetscReal :: area(sec_heat_vars%ncells)
  PetscReal :: vol(sec_heat_vars%ncells)
  PetscReal :: dm_plus(sec_heat_vars%ncells)
  PetscReal :: dm_minus(sec_heat_vars%ncells)
  PetscInt :: i, ngcells
  PetscReal :: area_fm
  PetscReal :: alpha, therm_conductivity, dencpr
  PetscReal :: temp_primary_node
  PetscReal :: m
  
  ngcells = sec_heat_vars%ncells
  area = sec_heat_vars%area
  vol = sec_heat_vars%vol
  dm_plus = sec_heat_vars%dm_plus
  dm_minus = sec_heat_vars%dm_minus
  area_fm = sec_heat_vars%interfacial_area
  temp_primary_node = global_auxvar%temp
  
  coeff_left = 0.d0
  coeff_diag = 0.d0
  coeff_right = 0.d0
  rhs = 0.d0
  sec_temp = 0.d0
  
  alpha = option%flow_dt*therm_conductivity/dencpr

  
  ! Setting the coefficients
  do i = 2, ngcells-1
    coeff_left(i) = -alpha*area(i-1)/((dm_minus(i) + dm_plus(i-1))*vol(i))
    coeff_diag(i) = alpha*area(i-1)/((dm_minus(i) + dm_plus(i-1))*vol(i)) + &
                    alpha*area(i)/((dm_minus(i+1) + dm_plus(i))*vol(i)) + 1.d0
    coeff_right(i) = -alpha*area(i)/((dm_minus(i+1) + dm_plus(i))*vol(i))
  enddo
  
  coeff_diag(1) = alpha*area(1)/((dm_minus(2) + dm_plus(1))*vol(1)) + 1.d0
  coeff_right(1) = -alpha*area(1)/((dm_minus(2) + dm_plus(1))*vol(1))
  
  coeff_left(ngcells) = -alpha*area(ngcells-1)/ &
                       ((dm_minus(ngcells) + dm_plus(ngcells-1))*vol(ngcells))
  coeff_diag(ngcells) = alpha*area(ngcells-1)/ &
                       ((dm_minus(ngcells) + dm_plus(ngcells-1))*vol(ngcells)) &
                       + alpha*area(ngcells)/(dm_plus(ngcells)*vol(ngcells)) &
                       + 1.d0

                        
  rhs = sec_heat_vars%sec_temp  ! secondary continuum values from previous time step
  rhs(ngcells) = rhs(ngcells) + & 
                 alpha*area(ngcells)/(dm_plus(ngcells)*vol(ngcells))* &
                 temp_primary_node
                
  ! Thomas algorithm for tridiagonal system
  ! Forward elimination
  do i = 2, ngcells
    m = coeff_left(i)/coeff_diag(i-1)
    coeff_diag(i) = coeff_diag(i) - m*coeff_right(i-1)
    rhs(i) = rhs(i) - m*rhs(i-1)
  enddo

  ! Back substitution
  ! Calculate temperature in the secondary continuum
  sec_temp(ngcells) = rhs(ngcells)/coeff_diag(ngcells)
  do i = ngcells-1, 1, -1
    sec_temp(i) = (rhs(i) - coeff_right(i)*sec_temp(i+1))/coeff_diag(i)
  enddo
  
  sec_heat_vars%sec_temp = sec_temp
            
end subroutine THSecHeatAuxVarCompute

! ************************************************************************** !

subroutine SecondaryRTotalSorb(rt_auxvar,global_auxvar,material_auxvar,reaction, &
                               option)
  ! 
  ! Computes the secondary total sorbed component concentrations and
  ! derivative with respect to free-ion
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 02/20/2014
  ! 

  use Option_module
  use Global_Aux_module
  use Reaction_Aux_module
  use Reaction_module
  use Reactive_Transport_Aux_module
  use Material_Aux_class
  
  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  call RZeroSorb(rt_auxvar)
  
  if (reaction%neqkdrxn > 0) then
    call SecondaryRTotalSorbKD(rt_auxvar,global_auxvar,material_auxvar, &
                               reaction,option)
  endif
  
end subroutine SecondaryRTotalSorb

! ************************************************************************** !

subroutine SecondaryRTotalSorbKD(rt_auxvar,global_auxvar,material_auxvar,reaction, &
                        option)
  ! 
  ! Computes the total sorbed component concentrations and
  ! derivative with respect to free-ion for the linear
  ! K_D model
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 02/20/2014
  ! 

  use Option_module
  use Reaction_Aux_module
  use Reaction_module
  use Reactive_Transport_Aux_module
  use Material_Aux_class
  use Global_Aux_module

  implicit none

  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  PetscInt :: irxn
  PetscInt :: icomp
  PetscReal :: res
  PetscReal :: dres_dc
  PetscReal :: activity
  PetscReal :: molality
  PetscReal :: tempreal
  PetscReal :: one_over_n
  PetscReal :: activity_one_over_n

  ! Surface Complexation
  do irxn = 1, reaction%neqkdrxn
    icomp = reaction%eqkdspecid(irxn)
    molality = rt_auxvar%pri_molal(icomp)
    activity = molality*rt_auxvar%pri_act_coef(icomp) ! Activity coefficient needs?
    select case(reaction%sec_cont_eqkdtype(irxn))
      case(SORPTION_LINEAR)
        ! Csorb = Kd*Caq
        res = reaction%sec_cont_eqkddistcoef(irxn)*activity
        dres_dc = res/molality
      case(SORPTION_LANGMUIR)
        ! Csorb = K*Caq*b/(1+K*Caq)
        tempreal = reaction%sec_cont_eqkddistcoef(irxn)*activity
        res = tempreal*reaction%sec_cont_eqkdlangmuirb(irxn) / (1.d0 + tempreal)
        dres_dc = res/molality - &
                  res / (1.d0 + tempreal) * tempreal / molality
      case(SORPTION_FREUNDLICH)
        ! Csorb = Kd*Caq**(1/n)
        one_over_n = 1.d0/reaction%sec_cont_eqkdfreundlichn(irxn)
        activity_one_over_n = activity**one_over_n
        res = reaction%sec_cont_eqkddistcoef(irxn)* &
                activity**one_over_n
        dres_dc = res/molality*one_over_n
      case default
        res = 0.d0
        dres_dc = 0.d0
    end select
    rt_auxvar%total_sorb_eq(icomp) = rt_auxvar%total_sorb_eq(icomp) + res
    rt_auxvar%dtotal_sorb_eq(icomp,icomp) = &
      rt_auxvar%dtotal_sorb_eq(icomp,icomp) + dres_dc 
  enddo

end subroutine SecondaryRTotalSorbKD


end module Secondary_Continuum_module
            
