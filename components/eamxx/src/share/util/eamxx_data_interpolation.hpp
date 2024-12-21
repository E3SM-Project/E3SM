#ifndef EAMXX_DATA_INTERPOLATION_HPP
#define EAMXX_DATA_INTERPOLATION_HPP

#include "share/grid/abstract_grid.hpp"
#include "share/grid/remap/abstract_remapper.hpp"
#include "share/util/scream_time_stamp.hpp"
#include "share/field/field.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/util/scream_time_stamp.hpp"

namespace scream{

class VerticalRemapper;

class DataInterpolation
{
public:
  using strvec_t = std::vector<std::string>;
  enum VRemapType {
    None,
    Static1D,
    Dynamic3D
  };

  // Constructor(s) & Destructor
  DataInterpolation (const std::shared_ptr<const AbstractGrid>& model_grid,
                     const std::vector<Field>& fields);

  ~DataInterpolation () = default;

  void toggle_debug_output (bool enable_dbg_output) { m_dbg_output = enable_dbg_output; }

  void setup_time_database (const strvec_t& input_files, const util::TimeLine timeline);

  void setup_remappers (const std::string& hremap_filename,
                        const VRemapType vr_type,
                        const std::string& data_pname,
                        const Field& model_pmid,
                        const Field& model_pint);

  void setup_remappers (const std::string& hremap_filename,
                        const VRemapType vr_type,
                        const std::string& extrap_type_top,
                        const std::string& extrap_type_bot,
                        const Real mask_value,
                        const std::string& data_pname,
                        const Field& model_pmid,
                        const Field& model_pint);

  void init_data_interval (const util::TimeStamp& t0);

  void run (const util::TimeStamp& ts);

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

  std::vector<Field>                  m_fields;

  // Use two horiz remappers, so we only set them up once (it may be costly)
  std::shared_ptr<AbstractRemapper> m_horiz_remapper_beg;
  std::shared_ptr<AbstractRemapper> m_horiz_remapper_end;
  std::shared_ptr<AbstractRemapper> m_vert_remapper;
  std::shared_ptr<VerticalRemapper> m_vremap;

  VRemapType            m_vr_type;
  int                   m_nfields;

  util::TimeInterval    m_data_interval;
  std::pair<int,int>    m_curr_interval_idx;

  TimeDatabase          m_time_database;

  ekat::Comm            m_comm;
  ekat::ParameterList   m_params;

  bool                  m_time_db_created   = false;
  bool                  m_remappers_created = false;
  bool                  m_data_initialized  = false;

  bool m_dbg_output = false;
};

} // namespace scream

#endif // EAMXX_DATA_INTERPOLATION_HPP
