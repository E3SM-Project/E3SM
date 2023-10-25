#ifndef SCREAM_IMPORTER_HPP
#define SCREAM_IMPORTER_HPP

#include "share/atm_process/atmosphere_process.hpp"
#include "ekat/ekat_parameter_list.hpp"
#include "share/atm_process/SCDataManager.hpp"

#include "surface_coupling_utils.hpp"

#include <string>

namespace scream
{

/*
 * The class responsible to handle the calculation of the subgrid cloud fractions
 *
 * The AD should store exactly ONE instance of this class stored
 * in its list of subcomponents (the AD should make sure of this).
*/

class SurfaceCouplingImporter : public AtmosphereProcess
{
public:

  template<typename DevT, typename DataT>
  using view_1d = typename KokkosTypes<DevT>::template view_1d<DataT>;
  template<typename DevT, typename DataT>
  using view_2d = typename KokkosTypes<DevT>::template view_2d<DataT>;

  template<typename DevT, typename ScalarT>
  using uview_1d = Unmanaged<view_1d<DevT, ScalarT>>;
  template<typename DevT, typename ScalarT>
  using uview_2d = Unmanaged<view_2d<DevT, ScalarT>>;

   using name_t = char[32];

  // Constructors
  SurfaceCouplingImporter (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The type of subcomponent
  AtmosphereProcessType type () const {
    return AtmosphereProcessType::SurfaceCouplingImporter;
  }

  // The name of the subcomponent
  std::string name () const { return "SurfaceCouplingImporter"; }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

  // Function which performes the import to scream fields.
  // If calling in initialize_impl(), set
  // called_during_initialization=true to avoid importing to fields
  // which are not needed during initialization.
  void do_import(const bool called_during_initialization=false);

  // Take and store data from SCDataManager
  void setup_surface_coupling_data(const SCDataManager &sc_data_manager);

protected:

  // The three main overrides for the subcomponent
  void initialize_impl (const RunType run_type);
  void run_impl        (const double dt);
  void finalize_impl   ();

  // Keep track of field dimensions
  Int m_num_cols; 

  // Number of fields in cpl data
  Int m_num_cpl_imports;

  // Number of imports to SCREAM
  Int m_num_scream_imports;

  // Views storing a 2d array with dims (num_cols,num_fields) for import data.
  // The field idx strides faster, since that's what mct does (so we can "view" the
  // pointer to the whole x2a array from Fortran)
  view_2d <DefaultDevice, Real> m_cpl_imports_view_d;
  uview_2d<HostDevice,    Real> m_cpl_imports_view_h;

  // Array storing the field names for imports
  name_t* m_import_field_names;

  // Views storing information for each import
  uview_1d<HostDevice, int>  m_cpl_indices_view;
  uview_1d<HostDevice, int>  m_vector_components_view;
  uview_1d<HostDevice, Real> m_constant_multiple_view;
  uview_1d<HostDevice, bool> m_do_import_during_init_view;

  // Column info used during import
  view_1d<DefaultDevice, SurfaceCouplingColumnInfo> m_column_info_d;
  decltype(m_column_info_d)::HostMirror             m_column_info_h;


  // The grid is needed for property checks
  std::shared_ptr<const AbstractGrid> m_grid;

}; // class SurfaceCouplingImporter

} // namespace scream

#endif // SCREAM_CLD_FRACTION_HPP
