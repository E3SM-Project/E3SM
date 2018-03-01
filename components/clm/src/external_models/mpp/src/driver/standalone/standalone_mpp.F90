! ----------------------------------------------------------------------
function remove_filename_extension(filename)
  !
  ! Remove the extension from a file name to get a base filename.
  !
  ! We start at the end of the filename and assume that the extension
  ! is marked by a period.
  !    
  implicit none
    
  character(len=256), intent(in) :: filename
  character(len=256) :: remove_filename_extension
  integer :: ext_index
    
  ext_index = scan(filename, '.', .true.)
  if (ext_index == 0) then
     ! no period marking an extension...
     ext_index = len(trim(filename)) + 1
     
  end if
  remove_filename_extension = filename(1:ext_index-1)

end function remove_filename_extension

! ----------------------------------------------------------------------
program standalone_mpp
  !
#include <petsc/finclude/petsc.h>
  !
  use mass_and_heat_model_problem , only : run_mass_and_heat_model_problem
  use mass_and_heat_model_problem , only : output_regression_mass_and_heat_model_problem
  use heat_transport_1D_problem   , only : run_heat_transport_1D_problem
  use heat_transport_1D_problem   , only : output_regression_heat_transport_1D_problem
  use vsfm_celia1990_problem      , only : run_vsfm_celia1990_problem
  use vsfm_celia1990_problem      , only : output_regression_vsfm_celia1990_problem
  use vsfm_vchannel_problem       , only : run_vsfm_vchannel_problem
  use vsfm_vchannel_problem       , only : output_regression_vsfm_vchannel_problem
  use vsfm_spac_problem           , only : run_vsfm_spac_problem
  use vsfm_spac_problem           , only : output_regression_vsfm_spac_problem
  use vsfm_spac_campbell_problem  , only : run_vsfm_spac_campbell_problem
  use vsfm_spac_campbell_problem  , only : output_regression_vsfm_spac_campbell_problem
  use petscsys
  !
  implicit none
  !
  !
  PetscErrorCode               :: ierr
  character(len=256)           :: namelist_filename
  character(len=256)           :: filename_base
  character(len=256), external :: remove_filename_extension
  character(len=256)           :: problem_type
  character(len=256)           :: ioerror_msg
  character(len=2560)          :: namelist_buffer
  logical                      :: write_regression_output
  PetscInt                     :: num_cells
  PetscBool                    :: flg
  integer                      :: nml_unit, nml_error

  namelist / mpp_driver / problem_type
  namelist / regression_test / write_regression_output, num_cells
  nml_unit = 16

  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

  PETSC_COMM_WORLD     = MPI_COMM_WORLD
  PETSC_COMM_SELF      = MPI_COMM_SELF
  namelist_filename    = ''
  write_regression_output = .false.
  num_cells            = 0

  call PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-namelist', &
       namelist_filename, flg, ierr)

  if (len(trim(adjustl(namelist_filename))) ==  0) then
     write(*,*)'ERROR: -namelist <filename> was not defined. Bailing out.'
     call exit(-1)
  endif

  filename_base = remove_filename_extension(namelist_filename)

  !
  ! Read namelist in buffer
  !
  open(unit=nml_unit, file=trim(namelist_filename), action='read', access='stream', &
       form='unformatted', iostat=nml_error)
  if (nml_error /= 0) then
     write(*,*)'ERROR: Unable to open namelist file: ',trim(namelist_filename)
     call exit(-1)
  endif

  read(unit=nml_unit, iostat=nml_error, iomsg=ioerror_msg) namelist_buffer
  if (.not. is_iostat_end(nml_error)) then
     write(*,*)"ERROR: Unable to read '",trim(namelist_filename),"' till EOF"
     call exit(-1)
  else
     write(*, '(a, a, a)') "Read '", trim(namelist_filename), "' until EOF."
  endif
     
  !
  ! Read namelist: mpp_driver
  !
  read(namelist_buffer, nml=mpp_driver, iostat=nml_error, iomsg=ioerror_msg)
  if (nml_error /= 0) then
     write(*,*)'ERROR: Unable to read "mpp_driver" in namelist file '
     call exit(-1)
  endif

  !
  ! Read namelist: regression_test
  !
  read(namelist_buffer, nml=regression_test, iostat=nml_error, iomsg=ioerror_msg)

  if (trim(problem_type) == 'mass_and_heat') then
     call run_mass_and_heat_model_problem()
     if (write_regression_output) then
        call output_regression_mass_and_heat_model_problem(filename_base, num_cells)
     endif

  else if (trim(problem_type) == 'heat_transport_1D') then
     call run_heat_transport_1D_problem()        
     if (write_regression_output) then
        call output_regression_heat_transport_1D_problem(filename_base, num_cells)
     endif

  else if(trim(problem_type) == 'vsfm_celia1990') then
     call run_vsfm_celia1990_problem()

     if (write_regression_output) then
        call output_regression_vsfm_celia1990_problem(filename_base, num_cells)
     endif

  else if(trim(problem_type) == 'vsfm_vchannel') then
     call run_vsfm_vchannel_problem()

     if (write_regression_output) then
        call output_regression_vsfm_vchannel_problem(filename_base, num_cells)
     endif

  else if(trim(problem_type) == 'vsfm_spac') then
     call run_vsfm_spac_problem()

     if (write_regression_output) then
        call output_regression_vsfm_spac_problem(filename_base, num_cells)
     endif

  else if(trim(problem_type) == 'vsfm_spac_campbell') then
     call run_vsfm_spac_campbell_problem()

     if (write_regression_output) then
        call output_regression_vsfm_spac_campbell_problem(filename_base, num_cells)
     endif

  else
     write(*,*)"problem_type = '", trim(problem_type), "' is unsupported."
  endif

  close(nml_unit)

  call PetscFinalize(ierr)

  call exit(0)
  
end program standalone_mpp
