module shr_nuopc_multiinst_mod

  implicit none
  public

  ! NOTE: NUM_COMP_INST_XXX are cpp variables set in buildlib.csm_share

  integer, parameter :: num_inst_atm = NUM_COMP_INST_ATM
  integer, parameter :: num_inst_lnd = NUM_COMP_INST_LND
  integer, parameter :: num_inst_ocn = NUM_COMP_INST_OCN
  integer, parameter :: num_inst_ice = NUM_COMP_INST_ICE
  integer, parameter :: num_inst_glc = NUM_COMP_INST_GLC
  integer, parameter :: num_inst_wav = NUM_COMP_INST_WAV
  integer, parameter :: num_inst_rof = NUM_COMP_INST_ROF
  integer, parameter :: num_inst_esp = NUM_COMP_INST_ESP
  integer, parameter :: num_inst_total = &
       num_inst_atm + num_inst_lnd + num_inst_ocn + num_inst_ice + &
       num_inst_glc + num_inst_wav + num_inst_rof + num_inst_esp + 1

  integer :: num_inst_min, num_inst_max

end module shr_nuopc_multiinst_mod
