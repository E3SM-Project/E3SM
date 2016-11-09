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
   use cam_history,  only:    addfld, horiz_only, outfld, add_default

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
      call addfld ('TDIFF',(/ 'lev' /),    'A','K','difference from observed temp', gridname='gauss_grid')

      call addfld ('TOBS',(/ 'lev' /),    'A','K','observed temp')
      call addfld ('QDIFF',(/ 'lev' /),    'A','kg/kg','difference from observed water',gridname='gauss_grid')

      call addfld ('QOBS',(/ 'lev' /),    'A','kg/kg','observed water')
      call addfld ('PRECOBS',(/ 'lev' /),    'A','mm/day','Total (convective and large-scale) precipitation rate')
      call addfld ('DIVQ',(/ 'lev' /),    'A','kg/kg/s','Q advection tendency (horizontal)')
      call addfld ('DIVQ3D',(/ 'lev' /),    'A','kg/kg/s','Q advection tendency (horiz/vert combined)', gridname='gauss_grid')
      call addfld ('DIVV',(/ 'lev' /),    'A','m/s2','V advection tendency (horizontal)')
      call addfld ('DIVU',(/ 'lev' /),    'A','m/s2','U advection tendency (horizontal)')
      call addfld ('DIVT',(/ 'lev' /),    'A','K/s','T advection tendency (horizontal)')
      call addfld ('DIVT3D',(/ 'lev' /),    'A','K/s','T advection tendency (horiz/vert combined)', gridname='gauss_grid')

      call addfld ('SHFLXOBS',horiz_only,    'A','W/m2','Obs Surface sensible heat flux')
      call addfld ('LHFLXOBS',horiz_only,    'A','W/m2','Obs Surface latent heat flux')
      call addfld ('TRELAX',(/ 'lev' /),    'A','K','t relaxation amount', gridname='gauss_grid')
      call addfld ('QRELAX',(/ 'lev' /),    'A','kg/kg','q relaxation amount', gridname='gauss_grid')
      call addfld ('TAURELAX',(/ 'lev' /),    'A','seconds','relaxation time constant', gridname='gauss_grid')
      call add_default ('TDIFF     ', 1, ' ')
      call add_default ('QDIFF     ', 1, ' ')
   end subroutine scm_intht

!#######################################################################
 end module history_scam
