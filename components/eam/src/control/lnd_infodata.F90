module lnd_infodata

  !Purpose: This module contains info that the
  ! atmosphere model needs from the land model
  ! during the model initialization such as
  ! namelist variables and flags

use shr_kind_mod, only: cs => shr_kind_cs

implicit none

!keep everything private
private

!Precipitation downscaling method used in the land model (current possible options: ERMM (default), FNM)
!(Initialize to UNSET to avoid uninitialized behavior)
character(len = cs), public :: precip_downscaling_method = 'UNSET'

end module lnd_infodata
