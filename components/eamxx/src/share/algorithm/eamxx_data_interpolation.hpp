#ifndef EAMXX_DATA_INTERPOLATION_HPP
#define EAMXX_DATA_INTERPOLATION_HPP

#include "share/grid/abstract_grid.hpp"
#include "share/remap/abstract_remapper.hpp"
#include "share/util/eamxx_time_stamp.hpp"
#include "share/field/field.hpp"
#include "share/field/field_reader.hpp"

#include <ekat_logger.hpp>

namespace scream
{

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

  enum TimeInterpType {
    Linear,
    Nearest
  };

  struct VertRemapData {
    VertRemapData() = default;

    VRemapType vr_type = None;
    std::string extrap_top = "P0";
    std::string extrap_bot = "P0";
    std::string pname; // What we need to load from nc file
    Field pmid, pint;  // The model pmid/pint
    std::shared_ptr<AbstractRemapper> custom_remapper; // Use this custom remapper
  };

  // Constructor(s) & Destructor
  DataInterpolation (const std::shared_ptr<const AbstractGrid>& model_grid,
                     const std::vector<Field>& fields);

  ~DataInterpolation () = default;

  void set_logger (const std::shared_ptr<ekat::logger::LoggerBase>& logger);

  // Setup time database for LinearHistory time-dependent data interpolation
  void setup_linear_time_database (const strvec_t& input_files,
                                   const TimeInterpType interp_type = Linear,
                                   const util::TimeStamp& ref_ts = util::TimeStamp());

  // Setup time database for YearlyPeriodic time-dependent data interpolation
  // ref_ts shifts raw file time values; start_ts selects the first logical
  // slice (the slice at or after start_ts becomes index 0 via rotation).
  void setup_periodic_time_database (const strvec_t& input_files,
                                     const TimeInterpType interp_type = Linear,
                                     const util::TimeStamp& start_ts = util::TimeStamp(),
                                     const util::TimeStamp& ref_ts = util::TimeStamp());

  // In case the input files store dims with exhotic names, the user can provide them here
  void set_input_files_dimname (const std::string& name, const std::string& nc_name);

  void create_horiz_remappers (const std::string& map_file = "");
  void create_horiz_remappers (const Real iop_lat, const Real iop_lon);
  void create_vert_remapper ();
  void create_vert_remapper (const VertRemapData& data);

  void register_fields_in_remappers ();

  void init_data_interval (const util::TimeStamp& t0);

  void run (const util::TimeStamp& ts);

  std::shared_ptr<AbstractGrid> get_grid_after_hremap () const { return m_grid_after_hremap; }

  void set_name (const std::string& name) { m_name = name; }

protected:

  void shift_data_interval ();
  void update_end_fields ();

  // Common body shared by setup_linear/periodic_time_database.
  // Reads timestamps from the given files and builds m_time_database.slices/files.
  void build_time_database_slices (const strvec_t& input_files,
                                   util::TimeLine timeline,
                                   const util::TimeStamp& ref_ts);
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

  std::shared_ptr<FieldReader> m_reader;

  std::shared_ptr<const AbstractGrid> m_model_grid;

  // The following are nonconst, since we may need to modify a few things
  std::shared_ptr<AbstractGrid>       m_data_grid;
  std::shared_ptr<AbstractGrid>       m_grid_after_hremap;

  std::vector<Field>                  m_fields;

  // Use two horiz remappers, so we only set them up once (it may be costly)
  std::shared_ptr<AbstractRemapper> m_horiz_remapper_beg;
  std::shared_ptr<AbstractRemapper> m_horiz_remapper_end;
  std::shared_ptr<AbstractRemapper> m_vert_remapper;

  // These are inited as the usual "ncol" and "lev" at construction, but the user
  // can reset them in case the input files store funky dimensions
  std::map<std::string,std::string>    m_input_files_dimnames;

  // If vertical remap happens, at runtime we may need to access some
  // versions of certain perssure fields. Store them here for convenient access
  std::map<std::string,Field>     m_helper_pressure_fields;

  VRemapType            m_vr_type;
  int                   m_nfields;

  util::TimeInterval    m_data_interval;
  std::pair<int,int>    m_curr_interval_idx;

  TimeDatabase          m_time_database;
  TimeInterpType        m_time_interp_type;

  ekat::Comm            m_comm;

  bool                  m_time_db_created   = false;
  bool                  m_data_initialized  = false;

  bool                  m_fields_have_col_dim = false;
  bool                  m_fields_have_lev_dim = false;
  bool                  m_fields_have_ilev_dim = false;

  std::string           m_name = "DataInterp";

  std::shared_ptr<ekat::logger::LoggerBase> m_logger;
};

} // namespace scream

#endif // EAMXX_DATA_INTERPOLATION_HPP
