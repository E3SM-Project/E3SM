#ifndef SCREAM_EXPORTER_HPP
#define SCREAM_EXPORTER_HPP

#include "surface_coupling_utils.hpp"

#include "share/atm_process/atmosphere_process.hpp"
#include "share/util/eamxx_common_physics_functions.hpp"
#include "share/util/eamxx_time_interpolation.hpp"
#include "share/atm_process/ATMBufferManager.hpp"
#include "share/atm_process/SCDataManager.hpp"

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

// enum to track how exported fields will be set.
enum ExportType:int {
  FROM_MODEL =  0,  // Variable will be derived from atmosphere model state
  FROM_FILE  =  1,  // Variable will be set given data from a file
  CONSTANT   =  2   // Set variable to a constant value
};


class SurfaceCouplingExporter : public AtmosphereProcess
{
public:

  using PF      = scream::PhysicsFunctions<DefaultDevice>;
  using KT      = ekat::KokkosTypes<DefaultDevice>;
  using Spack   = ekat::Pack<Real,SCREAM_SMALL_PACK_SIZE>;

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
  SurfaceCouplingExporter (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The type of subcomponent
  AtmosphereProcessType type () const {
    return AtmosphereProcessType::SurfaceCouplingExporter;
  }

  // The name of the subcomponent
  std::string name () const { return "SurfaceCouplingExporter"; }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

  // Structure for storing local variables initialized using the ATMBufferManager
  struct Buffer {
    static constexpr int num_2d_vector_mid = 2;
    static constexpr int num_2d_vector_int = 1;

    uview_2d<DefaultDevice, Spack> dz;
    uview_2d<DefaultDevice, Spack> z_mid;
    uview_2d<DefaultDevice, Spack> z_int;
  };

  // Function which performes the export from scream fields,
  // If calling in initialize_impl(), set
  // called_during_initialization=true to avoid exporting fields
  // which do not have valid entries.
  void do_export(const double dt, const bool called_during_initialization=false);             // Main export routine
  void compute_eamxx_exports(const double dt, const bool called_during_initialization=false); // Export vars are derived from eamxx state
  void set_constant_exports();                                                                // Export vars are set to a constant
  void set_from_file_exports();                                                               // Export vars are set by interpolation of data from files
  void do_export_to_cpl(const bool called_during_initialization=false);                       // Finish export by copying data to cpl structures.

  // Take and store data from SCDataManager
  void setup_surface_coupling_data(const SCDataManager &sc_data_manager);
protected:

  // The three main overrides for the subcomponent
  void initialize_impl (const RunType run_type);
  void run_impl        (const double dt);
  void finalize_impl   ();

  // Creates an helper field, not to be shared with the AD's FieldManager
  void create_helper_field (const std::string& name,
                            const FieldLayout& layout,
                            const std::string& grid_name);

  // Query if a local field exists
  bool has_helper_field (const std::string& name) const { return m_helper_fields.find(name)!=m_helper_fields.end(); }

  // Computes total number of bytes needed for local variables
  size_t requested_buffer_size_in_bytes() const;

  // Set local variables using memory provided by
  // the ATMBufferManager
  void init_buffers(const ATMBufferManager &buffer_manager);

  std::shared_ptr<const AbstractGrid> m_grid;

  // Keep track of field dimensions and the iteration count
  Int m_num_cols; 
  Int m_num_levs;

  // Some helper fields.
  std::map<std::string,Field> m_helper_fields;

  // Struct which contains local variables
  Buffer m_buffer;

  // Number of fields in cpl data
  Int m_num_cpl_exports;

  // Number of exports from EAMxx and how they will be handled
  int                               m_num_scream_exports;
  view_1d<DefaultDevice,ExportType> m_export_source;
  view_1d<HostDevice,ExportType>    m_export_source_h;
  std::map<std::string,Real>        m_export_constants;
  int                               m_num_from_model_exports=0;
  int                               m_num_const_exports=0;
  // For exporting from file
  int                               m_num_from_file_exports=0;
  util::TimeInterpolation           m_time_interp;
  std::vector<std::string>          m_export_from_file_field_names;

  // Views storing a 2d array with dims (num_cols,num_fields) for cpl export data.
  // The field idx strides faster, since that's what mct does (so we can "view" the
  // pointer to the whole a2x array from Fortran)
  view_2d <DefaultDevice, Real> m_cpl_exports_view_d;
  uview_2d<HostDevice,    Real> m_cpl_exports_view_h;

#ifdef HAVE_MOAB
  // Views storing a 2d array with dims (num_fields, num_cols) for moab cpl export data.
  // The field cols strides faster, since that's what moab does (so we can "view" the
  // pointer to the whole a2x_am(:,:) array from Fortran)
  view_2d <DefaultDevice, Real> m_moab_cpl_exports_view_d;
  uview_2d<HostDevice,    Real> m_moab_cpl_exports_view_h;
#endif
  // Array storing the field names for exports
  name_t*                   m_export_field_names;
  std::vector<std::string>  m_export_field_names_vector;

  // Views storing information for each export
  uview_1d<HostDevice, int>  m_cpl_indices_view;
  uview_1d<HostDevice, int>  m_vector_components_view;
  uview_1d<HostDevice, Real> m_constant_multiple_view;
  uview_1d<HostDevice, bool> m_do_export_during_init_view;

  // Column info used during export
  view_1d<DefaultDevice, SurfaceCouplingColumnInfo> m_column_info_d;
  decltype(m_column_info_d)::HostMirror             m_column_info_h;

}; // class SurfaceCouplingExporter

} // namespace scream

#endif // SCREAM_CLD_FRACTION_HPP
