#ifndef SCREAM_PRESCRIBED_AEROSOL_HPP
#define SCREAM_PRESCRIBED_AEROSOL_HPP

#include "physics/spa/spa_functions.hpp"
#include "share/atm_process/atmosphere_process.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/grid/remap/abstract_remapper.hpp"
#include <ekat/ekat_parameter_list.hpp>

#include <string>

namespace scream
{

/*
 * The class responsible to handle the calculation of the subgrid cloud fractions
 *
 * The AD should store exactly ONE instance of this class stored
 * in its list of subcomponents (the AD should make sure of this).
*/

class SPA : public AtmosphereProcess
{
public:
  using SPAFunc  = spa::SPAFunctions<Real, DefaultDevice>;
  using Spack    = SPAFunc::Spack;
  using KT       = ekat::KokkosTypes<DefaultDevice>;

  using view_1d = typename SPAFunc::view_1d<Spack>;
  using view_2d = typename SPAFunc::view_2d<Spack>;

  template<typename ScalarT>
  using uview_2d = Unmanaged<typename KT::template view_2d<ScalarT>>;

  // Constructors
  SPA (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The type of subcomponent
  AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }

  // The name of the subcomponent
  std::string name () const { return "spa"; }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

  // Structure for storing local variables initialized using the ATMBufferManager
  struct Buffer {
    // Used to store temporary data during spa_main
    SPAFunc::SPAInput spa_temp;

    // Temporary to use
    uview_2d<Spack> p_mid_src;
  };
protected:

  // The three main overrides for the subcomponent
  void initialize_impl (const RunType run_type);
  void initialize_spa_impl ();
  void run_impl        (const double dt);
  void finalize_impl   ();

  // Computes total number of bytes needed for local variables
  size_t requested_buffer_size_in_bytes() const;

  // Set local variables using memory provided by
  // the ATMBufferManager
  void init_buffers(const ATMBufferManager &buffer_manager);

  // Keep track of field dimensions and the iteration count
  int m_num_cols;
  int m_num_levs;
  int m_num_src_levs;
  int m_nswbands = 14;
  int m_nlwbands = 16;

  // Struct which contains temporary variables used during spa_main
  Buffer m_buffer;

  // IO structure to read in data for standard grids (keep it around to avoid re-creating PIO decomps)
  std::shared_ptr<AtmosphereInput>   SPADataReader;
  // Similar to above, but stores info to read data for IOP grid
  std::shared_ptr<SPAFunc::IOPReader>  SPAIOPDataReader;

  // Structures to store the data used for interpolation
  std::shared_ptr<AbstractRemapper>  SPAHorizInterp;

  SPAFunc::SPATimeState     SPATimeState;
  SPAFunc::SPAInput         SPAData_start;
  SPAFunc::SPAInput         SPAData_end;
  SPAFunc::SPAOutput        SPAData_out;

  std::shared_ptr<const AbstractGrid>   m_grid;
}; // class SPA

} // namespace scream

#endif // SCREAM_SPA_HPP
