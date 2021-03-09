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
      use dycore, only: dycore_is
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
! Local variables
!
      integer m,j        ! Indices
      real(r8) dummy
      character(len=100) dyngrid

      if (dycore_is('EUL')) then
        dyngrid = 'gauss_grid'
      else
        dyngrid = 'physgrid'
      endif
!
! Call addfld to add each field to the Master Field List.
!
      !+ Make this have backwards compatibility with Eulerian core
      call addfld ('TDIFF',(/ 'lev' /),    'A','K','difference from observed temp',gridname=trim(dyngrid))

      call addfld ('TOBS',(/ 'lev' /),    'A','K','observed temp',gridname=trim(dyngrid))
      call addfld ('QDIFF',(/ 'lev' /),    'A','kg/kg','difference from observed water',gridname=trim(dyngrid))

      call addfld ('QOBS',(/ 'lev' /),    'A','kg/kg','observed water',gridname=trim(dyngrid))
      call addfld ('PRECOBS',(/ 'lev' /),    'A','mm/day','Total (convective and large-scale) precipitation rate',gridname=trim(dyngrid))
      call addfld ('DIVQ',(/ 'lev' /),    'A','kg/kg/s','Q advection tendency (horizontal)',gridname=trim(dyngrid))
      call addfld ('DIVQ3D',(/ 'lev' /),    'A','kg/kg/s','Q advection tendency (horiz/vert combined)',gridname=trim(dyngrid))
      call addfld ('DIVV',(/ 'lev' /),    'A','m/s2','V advection tendency (horizontal)',gridname=trim(dyngrid))
      call addfld ('DIVU',(/ 'lev' /),    'A','m/s2','U advection tendency (horizontal)',gridname=trim(dyngrid))
      call addfld ('DIVT',(/ 'lev' /),    'A','K/s','T advection tendency (horizontal)',gridname=trim(dyngrid))
      call addfld ('DIVT3D',(/ 'lev' /),    'A','K/s','T advection tendency (horiz/vert combined)',gridname=trim(dyngrid))

      call addfld ('SHFLXOBS',horiz_only,    'A','W/m2','Obs Surface sensible heat flux',gridname=trim(dyngrid))
      call addfld ('LHFLXOBS',horiz_only,    'A','W/m2','Obs Surface latent heat flux',gridname=trim(dyngrid))
      call addfld ('TRELAX',(/ 'lev' /),    'A','K','t relaxation amount',gridname=trim(dyngrid))
      call addfld ('QRELAX',(/ 'lev' /),    'A','kg/kg','q relaxation amount',gridname=trim(dyngrid))
      call addfld ('TAURELAX',(/ 'lev' /),    'A','seconds','relaxation time constant',gridname=trim(dyngrid))
      call add_default ('TDIFF     ', 1, ' ')
      call add_default ('QDIFF     ', 1, ' ')
   end subroutine scm_intht   

!#######################################################################
 end module history_scam
