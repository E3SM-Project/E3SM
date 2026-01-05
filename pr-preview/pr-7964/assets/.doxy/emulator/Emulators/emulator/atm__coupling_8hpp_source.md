

# File atm\_coupling.hpp

[**File List**](files.md) **>** [**components**](dir_409f97388efe006bc3438b95e9edef48.md) **>** [**emulator\_comps**](dir_cd6ef227c082afa5b90fe3621cc9f093.md) **>** [**eatm**](dir_54689134e1a693092e83f56806593839.md) **>** [**src**](dir_1c3b735e18de9b9534f50214e18facf2.md) **>** [**impl**](dir_6975f7b28201ba7a9e865ff30c48a340.md) **>** [**atm\_coupling.hpp**](atm__coupling_8hpp.md)

[Go to the documentation of this file](atm__coupling_8hpp.md)


```C++


#ifndef ATM_COUPLING_HPP
#define ATM_COUPLING_HPP

#include "../../../common/src/coupling_fields.hpp"
#include "atm_field_manager.hpp"

namespace emulator {
namespace impl {

struct AtmCouplingIndices {
  // =========================================================================
  // Export indices (a2x)
  // =========================================================================
  int Sa_z = -1;       
  int Sa_u = -1;       
  int Sa_v = -1;       
  int Sa_tbot = -1;    
  int Sa_ptem = -1;    
  int Sa_shum = -1;    
  int Sa_dens = -1;    
  int Sa_pbot = -1;    
  int Sa_pslv = -1;    
  int Faxa_lwdn = -1;  
  int Faxa_rainc = -1; 
  int Faxa_rainl = -1; 
  int Faxa_snowc = -1; 
  int Faxa_snowl = -1; 
  int Faxa_swndr = -1; 
  int Faxa_swvdr = -1; 
  int Faxa_swndf = -1; 
  int Faxa_swvdf = -1; 
  int Faxa_swnet = -1; 

  // =========================================================================
  // Import indices (x2a)
  // =========================================================================
  int Sx_t = -1;      
  int So_t = -1;      
  int Faxx_sen = -1;  
  int Faxx_evap = -1; 
  int Faxx_lat = -1;  
  int Faxx_taux = -1; 
  int Faxx_tauy = -1; 
  int Faxx_lwup = -1; 
  int Sx_avsdr = -1;  
  int Sx_anidr = -1;  
  int Sx_avsdf = -1;  
  int Sx_anidf = -1;  
  int Sl_snowh = -1;  
  int Si_snowh = -1;  
  int Sx_tref = -1;   
  int Sx_qref = -1;   
  int Sx_u10 = -1;    
  int Sf_ifrac = -1;  
  int Sf_ofrac = -1;  
  int Sf_lfrac = -1;  

  void initialize(CouplingFieldsBase &fields);
};

void import_atm_fields(const double *import_data, int ncols, int nfields,
                       const AtmCouplingIndices &idx, AtmFieldManager &fields);

void export_atm_fields(double *export_data, int ncols, int nfields,
                       const AtmCouplingIndices &idx,
                       const AtmFieldManager &fields);

} // namespace impl
} // namespace emulator

#endif // ATM_COUPLING_HPP
```


