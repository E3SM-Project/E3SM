module history_scam
!----------------------------------------------------------------------- 
! 
! Purpose: SCAM specific history code.
!
! Public functions/subroutines:
!   bldfld, h_default
! 
! Author: anonymous from code in cam_history.F90
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8, r4 => shr_kind_r4
   use constituents, only: pcnst, cnst_name, cnst_longname, hadvnam, vadvnam, &
	                   tendnam, tottnam, fixcnam
   use ppgrid,       only: pver, pverp
   use pmgrid,       only: plev
   use cam_history,  only: phys_decomp, dyn_decomp,  addfld, outfld, add_default

   use scamMod, only :divq3d,divt3d,wfld,divq,divt,divu,divv

   implicit none

PRIVATE

   public :: scm_intht

!#######################################################################
CONTAINS
   subroutine scm_intht()
!----------------------------------------------------------------------- 
! 
! Purpose: 
!
! add master list fields to scm
! 
! Method: Call a subroutine to add each field
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
      use ppgrid, only: pver, pverp
      use constituents
      use cam_history, only: add_default
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
! Local variables
!
      integer m,j        ! Indices
      real(r8) dummy
!
! Call addfld to add each field to the Master Field List.
!
      call addfld ('TDIFF   ','K      ',plev,    'A','difference from observed temp', dyn_decomp)

      call addfld ('TOBS    ','K      ',plev,    'A','observed temp', phys_decomp)
      call addfld ('QDIFF   ','kg/kg   ',plev,    'A','difference from observed water',dyn_decomp)

      call addfld ('QOBS    ','kg/kg   ',plev,    'A','observed water',phys_decomp)
      call addfld ('PRECOBS','mm/day',plev,    'A','Total (convective and large-scale) precipitation rate', phys_decomp)
      call addfld ('DIVQ    ','kg/kg/s ',plev,    'A','Q advection tendency (horizontal)', phys_decomp)
      call addfld ('DIVQ3D  ','kg/kg/s ',pver,    'A','Q advection tendency (horiz/vert combined)', dyn_decomp)
      call addfld ('DIVV  ','m/s2    ',plev,    'A','V advection tendency (horizontal)', phys_decomp)
      call addfld ('DIVU  ','m/s2    ',plev,    'A','U advection tendency (horizontal)', phys_decomp)
      call addfld ('DIVT   ','K/s     ',plev,    'A','T advection tendency (horizontal)', phys_decomp)
      call addfld ('DIVT3D ','K/s     ',pver,    'A','T advection tendency (horiz/vert combined)', dyn_decomp)

      call addfld ('SHFLXOBS','W/m2    ',1,    'A','Obs Surface sensible heat flux',phys_decomp)
      call addfld ('LHFLXOBS','W/m2    ',1,    'A','Obs Surface latent heat flux',phys_decomp)
      call addfld ('TRELAX  ','K      ',plev,    'A','t relaxation amount', dyn_decomp)
      call addfld ('QRELAX  ','kg/kg  ',plev,    'A','q relaxation amount', dyn_decomp)
      call addfld ('TAURELAX','seconds',plev,    'A','relaxation time constant', dyn_decomp)
      call add_default ('TDIFF     ', 1, ' ')
      call add_default ('QDIFF     ', 1, ' ')
   end subroutine scm_intht

!#######################################################################
 end module history_scam
