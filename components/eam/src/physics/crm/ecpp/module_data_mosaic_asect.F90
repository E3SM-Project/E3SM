module module_data_mosaic_asect

  implicit none
  !-----------------------------------------------------------------------------
  ! The variables in this module provide a means of organizing and accessing
  ! aerosol species in the "chem" array by their chemical component, 
  ! size bin (or mode), "type", and "phase"
  !
  ! Their purpose is to allow flexible coding of process modules, 
  ! compared to "hard-coding" using the chem array p_xxx indices
  ! (e.g., p_so4_a01, p_so4_a02, ...; p_num_a01, ...)
  !
  !-----------------------------------------------------------------------------
  ! The aerosol type will allow treatment of an externally mixed 
  ! aerosol.  The current MOSAIC code has only 1 type, with the implicit
  ! assumption of internal mixing.  Eventually, multiple types 
  ! could treat fresh primary BC/OC, fresh SO4 from nucleation, 
  ! aged BC/OC/SO4/... mixture, soil dust, sea salt, ... 
  !
  ! [Note:  the value of "xx_phase" will be between 1 and nphase_aer 
  ! for phases that are active in a simulation.  The others
  ! will have non-positive values.]
  !
  !-----------------------------------------------------------------------------
  ! [Note:  dens_aer(c,t) == dens_mastercomp_aer(mastercompptr_aer(c,t))
  ! The dens_mastercomp_aer is used in some initialization routines.
  ! The dens_aer is used in most other places because of convenience.]
  !-----------------------------------------------------------------------------

  

end module module_data_mosaic_asect
