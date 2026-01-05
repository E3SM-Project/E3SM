

# File emulator\_comp.hpp

[**File List**](files.md) **>** [**common**](dir_03a82009d675450cd7b26f7887b718e6.md) **>** [**src**](dir_cf65ed25a3cb4ff7a3950496111c5548.md) **>** [**emulator\_comp.hpp**](emulator__comp_8hpp.md)

[Go to the documentation of this file](emulator__comp_8hpp.md)


```C++


#ifndef EMULATOR_COMP_HPP
#define EMULATOR_COMP_HPP

#include "emulator_logger.hpp"
#include <mpi.h>
#include <string>
#include <vector>

namespace emulator {

enum class CompType {
  ATM = 0, 
  OCN = 1, 
  ICE = 2, 
  LND = 3  
};

class EmulatorComp {
public:
  explicit EmulatorComp(CompType type);

  virtual ~EmulatorComp() = default;

  // =========================================================================
  // Public Interface (called from Fortran via C interface)
  // =========================================================================

  void create_instance(MPI_Comm comm, int comp_id, const char *input_file,
                       int run_type, int start_ymd, int start_tod);

  void set_grid_data(int nx, int ny, int num_local_cols, int num_global_cols,
                     const int *col_gids, const double *lat, const double *lon,
                     const double *area);

  void setup_coupling(double *import_data, double *export_data, int num_imports,
                      int num_exports, int field_size);

  void initialize();

  void run(int dt);

  void finalize();

  // =========================================================================
  // Grid Accessors
  // =========================================================================

  int get_num_local_cols() const { return m_num_local_cols; }

  int get_num_global_cols() const { return m_num_global_cols; }

  int get_nx() const { return m_nx; }

  int get_ny() const { return m_ny; }

  void get_local_col_gids(int *gids) const;

  void get_cols_latlon(double *lat, double *lon) const;

  void get_cols_area(double *area) const;

  // =========================================================================
  // Component Accessors
  // =========================================================================

  CompType type() const { return m_type; }

  MPI_Comm comm() const { return m_comm; }

  int comp_id() const { return m_comp_id; }

  int rank() const { return m_rank; }

  bool is_root() const { return m_rank == 0; }

protected:
  // =========================================================================
  // Virtual Methods (implement in derived classes)
  // =========================================================================

  virtual void init_impl() = 0;

  virtual void run_impl(int dt) = 0;

  virtual void final_impl() = 0;

  virtual void run_inference(const std::vector<double> &inputs,
                             std::vector<double> &outputs) = 0;

  virtual void import_from_coupler() {}

  virtual void export_to_coupler() {}

  // =========================================================================
  // Protected Data
  // =========================================================================

  MPI_Comm m_comm; 
  int m_comp_id;   
  int m_rank;      
  int m_nprocs;    
  CompType m_type; 
  int m_run_type;  
  Logger m_logger; 

  // Time and Steps
  int m_current_ymd = 0; 
  int m_current_tod = 0; 
  int m_step_count = 0;  

  // Grid data
  int m_num_local_cols = 0;    
  int m_num_global_cols = 0;   
  int m_nx = 0;                
  int m_ny = 0;                
  std::vector<int> m_col_gids; 
  std::vector<double> m_lat;   
  std::vector<double> m_lon;   
  std::vector<double> m_area;  

  // Coupling data
  double *m_import_data = nullptr; 
  double *m_export_data = nullptr; 
  int m_num_imports = 0;           
  int m_num_exports = 0;           
  int m_field_size = 0;            

  std::string m_input_file; 

private:
  void read_grid_file(const std::string &config_file);
  void distribute_grid_data(const std::vector<double> &lon_global,
                            const std::vector<double> &lat_global,
                            const std::vector<double> &area_global);
  void setup_default_grid();
  void advance_time(int dt);

  bool m_initialized = false;
};

// ===========================================================================
// Utility Functions
// ===========================================================================

inline std::string get_import_prefix(CompType type) {
  switch (type) {
  case CompType::ATM:
    return "x2a";
  case CompType::OCN:
    return "x2o";
  case CompType::ICE:
    return "x2i";
  case CompType::LND:
    return "x2l";
  default:
    return "x2x";
  }
}

inline std::string get_export_prefix(CompType type) {
  switch (type) {
  case CompType::ATM:
    return "a2x";
  case CompType::OCN:
    return "o2x";
  case CompType::ICE:
    return "i2x";
  case CompType::LND:
    return "l2x";
  default:
    return "x2x";
  }
}

inline std::string get_comp_name(CompType type) {
  switch (type) {
  case CompType::ATM:
    return "atm";
  case CompType::OCN:
    return "ocn";
  case CompType::ICE:
    return "ice";
  case CompType::LND:
    return "lnd";
  default:
    return "xxx";
  }
}

} // namespace emulator

#endif // EMULATOR_COMP_HPP
```


