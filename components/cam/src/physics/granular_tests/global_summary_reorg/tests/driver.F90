program main

  use test_core, only: initialize, finalize, test_glb_verif_smry

  call initialize
  call test_glb_verif_smry
  call finalize

end program
