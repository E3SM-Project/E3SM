module globals
   integer, parameter :: naer=10 ! number of aerosol species

   integer :: ncidi  = -1        ! input necdf file id
   integer :: londimidi = -1     ! longitude dimension id input file
   integer :: latdimidi = -1     ! latitude dimension id input file
   integer :: levdimidi = -1     ! level dimension id input file
   integer :: timedimidi = -1    ! time dimension id input file

!   integer :: nxi = -1           ! x-dimension size input file
!   integer :: nyi = -1           ! y-dimension size input file

   integer :: ncido  = -1        ! output necdf file id
   integer :: londimido = -1     ! longitude dimension id output file
   integer :: latdimido = -1     ! latitude dimension id output file
   integer :: ncoldimido = -1    ! ncol dimension id output file only for 1d horizontal files
   integer :: levdimido = -1     ! level dimension id output file
   integer :: ilevdimido = -1    ! interface dimension id output file
   integer :: timedimido = -1    ! time dimension id output file

!   integer :: nxo = -1           ! x-dimension size output file
!   integer :: nyo = -1           ! y-dimension size output file

!   integer :: nz = -1            ! z-dimension size input and output files
!   integer :: ntime = -1         ! time-dimension size input and output files
end module globals

