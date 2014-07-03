module fv_control_mod
  use shr_kind_mod, only : r8=> shr_kind_r8

! !PUBLIC MEMBER FUNCTIONS:
  public dyn_readnl        ! Set and/or get all runtime options


! !PUBLIC DATA MEMBERS:

  real(r8) :: tmass0
  real(r8) :: zgsint
  ! f-v dynamics specific
  ! _ord = 1: first order upwind
  ! _ord = 2: 2nd order van Leer (Lin et al 1994)
  ! _ord = 3: standard PPM 
  ! _ord = 4: enhanced PPM (default)

  integer :: nsplit = 0                  ! Lagrangian time splits (Lin-Rood only)
  integer :: nspltrac = 0                ! Tracer time splits (Lin-Rood only)
  integer :: nspltvrm = 0                ! Vertical re-mapping time splits (Lin-Rood only)
  integer :: iord = 4                    ! scheme to be used in E-W direction
  integer :: jord = 4                    ! scheme to be used in N-S direction
  integer :: kord = 4                    ! scheme to be used for vertical mapping
  logical :: dyn_conservative = .false.  ! Flag indicating whether the dynamics is conservative
  integer :: filtcw = 0                  ! flag for filtering c-grid winds
  integer :: ct_overlap = 0              ! nonzero for overlap of cd_core and trac2d, 0 otherwise
  integer :: trac_decomp = 1             ! size of tracer domain decomposition for trac2d
  integer :: fft_flt = 1                 ! 0 => FFT/algebraic filter; 1 => FFT filter
  integer :: div24del2flag = 2           ! 2 for 2nd order div damping, 4 for 4th order div damping,
                                         ! 42 for 4th order div damping plus 2nd order velocity damping
  real(r8):: del2coef = 3.e5_r8          ! strength of 2nd order velocity damping

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!-------------------------------------------------------------------------
!BOP
! !ROUTINE:  dyn_readnl --- Read dynamics-specific namelist variables
!
! !INTERFACE:
subroutine dyn_readnl( nlfilename )
! !USES:
   use units,           only: getunit, freeunit
   use namelist_utils,  only: find_group_name
   use spmd_utils,      only: masterproc
   use abortutils,      only: endrun
   use cam_logfile,     only: iulog
#if defined(SPMD)
   use mpishorthand
#endif

!
! !PARAMETERS:
   character(len=*), intent(in) :: nlfilename

! !DESCRIPTION: Read in FV-specific namelist variables.  Must be 
!               performed before dyn\_init
!
! !REVISION HISTORY:
!   2010.05.15   Sawyer  Creation
!
!EOP
!=========================================================================
!BOC
! Local variables
   integer :: ierr           ! error code
   integer :: unitn          ! namelist unit number
   character(len=*), parameter ::  subname = "dyn_readnl"

   namelist /dyn_fv_inparm/ nsplit, nspltrac, nspltvrm, filtcw, fft_flt,   &
                            iord, jord, kord, dyn_conservative,            &
                            div24del2flag, del2coef,ct_overlap, trac_decomp
   if (masterproc) then
      write(iulog,*) 'Read in dyn_fv_inparm namelist from: ', trim(nlfilename)
      unitn = getunit()
      open( unitn, file=trim(nlfilename), status='old' )

      ! Look for dyn_fv_inparm group name in the input file.  If found, leave the
      ! file positioned at that namelist group.
      call find_group_name(unitn, 'dyn_fv_inparm', status=ierr)
      if (ierr == 0) then  ! found dyn_fv_inparm
         read(unitn, dyn_fv_inparm, iostat=ierr)  ! read the dyn_fv_inparm namelist group
         if (ierr /= 0) then
            call endrun( subname//':: namelist read returns an'// &
                          ' error condition for dyn_fv_inparm' )
         end if
      else
         call endrun(subname // ':: can''t find dyn_fv_inparm in file ' // trim(nlfilename))
      end if
      close( unitn )
      call freeunit( unitn )
   endif

#if defined(SPMD)
   ! f-v dynamics specific
   call mpibcast (ct_overlap  ,1,mpiint,0,mpicom)
   call mpibcast (trac_decomp ,1,mpiint,0,mpicom)
   call mpibcast (nsplit  ,1,mpiint,0,mpicom)
   call mpibcast (nspltrac,1,mpiint,0,mpicom)
   call mpibcast (nspltvrm,1,mpiint,0,mpicom)
   call mpibcast (iord    ,1,mpiint,0,mpicom)
   call mpibcast (jord    ,1,mpiint,0,mpicom)
   call mpibcast (kord    ,1,mpiint,0,mpicom)
   call mpibcast (dyn_conservative,1,mpilog,0,mpicom)
   call mpibcast (filtcw  ,1,mpiint,0,mpicom)
   call mpibcast (fft_flt ,1,mpiint,0,mpicom)
   call mpibcast (div24del2flag ,1,mpiint,0,mpicom)
   call mpibcast (del2coef ,1,mpir8,0,mpicom)

#endif

!
! Finite volume code only: Set Lagrangian time splits.  A default of zero indicates the number
! should be automatically computed unless the user enters something.
!

!EOC
end subroutine dyn_readnl
!-------------------------------------------------------------------------

end module fv_control_mod
