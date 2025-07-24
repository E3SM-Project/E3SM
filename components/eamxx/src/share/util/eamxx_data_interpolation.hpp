#ifndef EAMXX_DATA_INTERPOLATION_HPP
#define EAMXX_DATA_INTERPOLATION_HPP

#include "share/grid/abstract_grid.hpp"
#include "share/grid/remap/abstract_remapper.hpp"
#include "share/util/eamxx_time_stamp.hpp"
#include "share/field/field.hpp"

namespace scream
{

class AtmosphereInput;

class DataInterpolation
{
public:
  using strvec_t = std::vector<std::string>;
  // This enum helps in two parts of DataInterpolation lifetime:
  //  - when building a VerticalRemapper, to correctly initialize it
  //  - at runtime, to load source data vertical profile coordinates
  // The first only applies if the user does not provide a pre-built vertical remapper
  enum VRemapType {
    None,         // This signales "no vert remap"
    Static1D,     // Uses a constant 1d (vertical) src pressure from input data
    Dynamic3D,    // Uses a time-dep 3d src pressure from input data
    Dynamic3DRef, // Reconstructs a reference 3d src pressure from time-dep PS in input data
    Custom        // The user will provide a vert remapper, and there won't be any need
                  // to set up the src pressure profile (the user will take care of it)
  };

  struct VertRemapData {
    VertRemapData() = default;

    VRemapType vr_type = None;
    std::string extrap_top = "P0";
    std::string extrap_bot = "P0";
    Real mask_value = std::numeric_limits<Real>::quiet_NaN(); // Unused for P0 extrapolation
    std::string pname; // What we need to load from nc file
    Field pmid, pint;  // The model pmid/pint
    std::shared_ptr<AbstractRemapper> custom_remapper; // Use this custom remapper
  };

  // Constructor(s) & Destructor
  DataInterpolation (const std::shared_ptr<const AbstractGrid>& model_grid,
                     const std::vector<Field>& fields);

  ~DataInterpolation () = default;

  void toggle_debug_output (bool enable_dbg_output) { m_dbg_output = enable_dbg_output; }

  void setup_time_database (const strvec_t& input_files,
                            const util::TimeLine timeline,
                            const util::TimeStamp& ref_ts = util::TimeStamp());

  // In case the input files store col/lev dims with exhotic names, the user can provide them here
  void set_input_files_dimname (const FieldTag t, const std::string& name) { m_input_files_dimnames[t] = name; }

  void create_horiz_remappers (const std::string& map_file = "");
  void create_horiz_remappers (const Real iop_lat, const Real iop_lon);
  void create_vert_remapper ();
  void create_vert_remapper (const VertRemapData& data);

  void register_fields_in_remappers ();

  void init_data_interval (const util::TimeStamp& t0);

  void run (const util::TimeStamp& ts);

  std::shared_ptr<AbstractGrid> get_grid_after_hremap () const { return m_grid_after_hremap; }

protected:

  void shift_data_interval ();
  void update_end_fields ();

  int get_input_files_dimlen (const std::string& dimname) const;

  // ----------- Internal data types ---------- //

  struct DataSlice {
    util::TimeStamp time;
    std::string     filename;
    int             time_idx; // slice index within the input file
  };

  struct TimeDatabase {
    strvec_t                files;
    std::vector<DataSlice>  slices;
    util::TimeLine          timeline;

    int size () const { return slices.size(); }
    int get_next_idx (int prev_idx) const;

    // Find interval containing t
    int find_interval (const util::TimeStamp& t) const;
  };

  // --------------- Internal data ------------- //

  std::shared_ptr<AtmosphereInput> m_reader;

  std::shared_ptr<const AbstractGrid> m_model_grid;
  std::shared_ptr<AbstractGrid>       m_grid_after_hremap; // nonconst b/c we may need to set some geo data

  std::vector<Field>                  m_fields;

  // Use two horiz remappers, so we only set them up once (it may be costly)
  std::shared_ptr<AbstractRemapper> m_horiz_remapper_beg;
  std::shared_ptr<AbstractRemapper> m_horiz_remapper_end;
  std::shared_ptr<AbstractRemapper> m_vert_remapper;

  // These are inited as the usual "ncol" and "lev" at construction, but the user
  // can reset them in case the input files store funky dimensions
  std::map<FieldTag,std::string>    m_input_files_dimnames;

  // If vertical remap happens, at runtime we may need to access some
  // versions of certain perssure fields. Store them here for convenient access
  std::map<std::string,Field>     m_helper_pressure_fields;

  VRemapType            m_vr_type;
  int                   m_nfields;

  util::TimeInterval    m_data_interval;
  std::pair<int,int>    m_curr_interval_idx;

  TimeDatabase          m_time_database;

  ekat::Comm            m_comm;
  ekat::ParameterList   m_params;

  bool                  m_time_db_created   = false;
  bool                  m_data_initialized  = false;

  bool m_dbg_output = false;
};

} // namespace scream

#endif // EAMXX_DATA_INTERPOLATION_HPP
