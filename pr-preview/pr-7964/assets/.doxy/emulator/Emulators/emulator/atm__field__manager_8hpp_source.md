

# File atm\_field\_manager.hpp

[**File List**](files.md) **>** [**components**](dir_409f97388efe006bc3438b95e9edef48.md) **>** [**emulator\_comps**](dir_cd6ef227c082afa5b90fe3621cc9f093.md) **>** [**eatm**](dir_54689134e1a693092e83f56806593839.md) **>** [**src**](dir_1c3b735e18de9b9534f50214e18facf2.md) **>** [**impl**](dir_6975f7b28201ba7a9e865ff30c48a340.md) **>** [**atm\_field\_manager.hpp**](atm__field__manager_8hpp.md)

[Go to the documentation of this file](atm__field__manager_8hpp.md)


```C++


#ifndef ATM_FIELD_MANAGER_HPP
#define ATM_FIELD_MANAGER_HPP

#include <map>
#include <string>
#include <vector>

namespace emulator {
namespace impl {

class AtmFieldManager {
public:
  // =========================================================================
  // Constants
  // =========================================================================

  static constexpr int N_INPUT_CHANNELS = 39;

  static constexpr int N_OUTPUT_CHANNELS = 44;

  AtmFieldManager() = default;
  ~AtmFieldManager() = default;

  // =========================================================================
  // Allocation
  // =========================================================================

  void allocate(int ncols);

  void deallocate();

  void set_defaults(int ncols);

  bool is_allocated() const { return m_allocated; }

  // =========================================================================
  // Generic Field Access
  // =========================================================================

  std::vector<double> *get_field_ptr(const std::string &name);

  void register_dynamic_field(const std::string &name);

  // =========================================================================
  // AI Model Tensors
  // =========================================================================

  std::vector<double> net_inputs;  
  std::vector<double> net_outputs; 

  // =========================================================================
  // Imported Fields (x2a - from coupler)
  // =========================================================================

  std::vector<double> shf;          
  std::vector<double> cflx;         
  std::vector<double> lhf;          
  std::vector<double> wsx;          
  std::vector<double> wsy;          
  std::vector<double> lwup;         
  std::vector<double> asdir;        
  std::vector<double> aldir;        
  std::vector<double> asdif;        
  std::vector<double> aldif;        
  std::vector<double> ts;           
  std::vector<double> sst;          
  std::vector<double> snowhland;    
  std::vector<double> snowhice;     
  std::vector<double> tref;         
  std::vector<double> qref;         
  std::vector<double> u10;          
  std::vector<double> u10withgusts; 
  std::vector<double> icefrac;      
  std::vector<double> ocnfrac;      
  std::vector<double> lndfrac;      

  // =========================================================================
  // Exported Fields (a2x - to coupler)
  // =========================================================================

  std::vector<double> zbot;  
  std::vector<double> ubot;  
  std::vector<double> vbot;  
  std::vector<double> tbot;  
  std::vector<double> ptem;  
  std::vector<double> shum;  
  std::vector<double> dens;  
  std::vector<double> pbot;  
  std::vector<double> pslv;  
  std::vector<double> lwdn;  
  std::vector<double> rainc; 
  std::vector<double> rainl; 
  std::vector<double> snowc; 
  std::vector<double> snowl; 
  std::vector<double> swndr; 
  std::vector<double> swvdr; 
  std::vector<double> swndf; 
  std::vector<double> swvdf; 
  std::vector<double> swnet; 

  // =========================================================================
  // Dynamic Fields
  // =========================================================================

  std::map<std::string, std::vector<double>> dynamic_fields;

private:
  bool m_allocated = false; 
  int m_ncols = 0;          

  void init_field_map();
  std::map<std::string, std::vector<double> *> m_field_map;
};

} // namespace impl
} // namespace emulator

#endif // ATM_FIELD_MANAGER_HPP
```


