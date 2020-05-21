#include "share/scream_session.hpp"

extern "C"
{

  void eam_init_pio(const MPI_Fint& mpicomm, const int& atm_id, const int& numdim, const int& numvar);
  

/* ========================================= */
  int eam_init_pio_c(const MPI_Fint& mpicomm, const int& atm_id, const int& numdim, const int& numvar) {

  return 0;

  } //eam_init_pio_c
/* ========================================= */
} // extern "C"
