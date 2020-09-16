#ifndef SCREAM_SCORPIO_HPP
#define SCREAM_SCORPIO_HPP

#include "scream_config.h"
#include "ekat/ekat_parameter_list.hpp"
#include "ekat/std_meta/ekat_std_enable_shared_from_this.hpp"
#include "ekat/mpi/ekat_comm.hpp"
#include "ekat/std_meta/ekat_std_utils.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"

#include "share/io/scream_scorpio_interface.hpp"

#include "share/field/field_repository.hpp"
#include "share/field/field.hpp"
#include "share/field/field_identifier.hpp"

/* There are a number of things that need to be cleaned up - in no particular order -
 *  1) Gather dimension data from grid manager or field repo.  We may need to use field repo since grid manager
 *     is unaware of say the number of components, for example.  --DONE--
 *  2) Clean-up init step, perhaps with subfunctions.  --DONE--
 *  3) Gather degree's of freedom from grid manager. Need to know the global-id index for
 *     each value in a view, in order if flattened.
 *  4) Handle MPI communication to allow for PIO_STRIDE (i.e. fewer pio ranks than total ranks).
 *  5) Make it so that all dimensions that are registered are automatically registered as variables and the values written (only once, not with time).  Probably need this information from the grid manager(?)
 *  6) Write Restart output.  Extra-savy would be to use DAG from AD to determine which fields were essential and just write those.
 *  7) Create an AtmosphereInput class built on these same concepts.  Can we make a master SCORPIO class with output and input as sub-classes?
 *  8) Create average, min, max and instantaneous options for output.
 *  9) Assure that IO is compatible with PACKing of variables.
 * 10) Add control of netCDF header information to the YAML file and this class.
 * 11) Better naming convention for output files?  Currently just keeps a counter, but in EAM the date/time range is used.
 */

namespace scream
{

class AtmosphereOutput 
{
public:
  using device_type = DefaultDevice; // may need to template class on this
  using Int = scream::Int;
  using Real = scream::Real;

  virtual ~AtmosphereOutput () = default;

  // Constructor
  AtmosphereOutput(const ekat::Comm& comm, const ekat::ParameterList& params)
  {
    m_comm   = comm;
    m_params = params;
  }

  // Main Functions
  void init(const FieldRepository<Real, device_type>& field_repo);
  void run(const FieldRepository<Real, device_type>& field_repo, const Real time);
  void finalize();

  // Helper Functions
  void check_status();
  std::map<std::string,Int> get_status() const { return m_status; }

protected:
  // Internal functions
  void add_identifier(const FieldRepository<Real, device_type>& field_repo, const std::string name);
  void register_variables();
  void set_degrees_of_freedom();
  // Internal variables
  ekat::ParameterList m_params;
  ekat::Comm          m_comm;
  
  std::string m_filename;
  std::string m_avg_type;
  std::string m_grid_name;

  Int m_out_max_steps;
  Int m_out_frequency;
  std::string m_out_units;

  std::map<std::string,FieldIdentifier>  m_fields;
  std::map<std::string,std::vector<Int>> m_dofs;
  std::map<std::string,Int>              m_dims;

  bool m_is_init = false;

  std::map<std::string,Int> m_status = {
                                  {"Init",    0},
                                  {"Run",     0},
                                  {"Finalize",0},
                                       }; 
private:

}; // Class AtmosphereOutput

// ====================== IMPLEMENTATION ===================== //
inline void AtmosphereOutput::init(const FieldRepository<Real, device_type>& field_repo /* , grid_manager */) 
{
  using namespace scream;
  using namespace scream::scorpio;

  m_status["Init"] += 1;
  m_filename = m_params.get<std::string>("FILENAME")+"_"+std::to_string(m_status["Init"])+".nc";  //TODO: The filename should be treated as a prefix to enable multiple files for the same control.  Like in the case of monthly output with 1 snap/file.
  EKAT_REQUIRE_MSG(!m_is_init,"Error! File " + m_filename + " has already been initialized.\n");

  // Parse the yaml file that controls this output instance.
  m_avg_type        = m_params.get<std::string>("AVERAGING TYPE");
  m_grid_name       = m_params.get<std::string>("GRID");
  auto& freq_params = m_params.sublist("FREQUENCY");
  m_out_max_steps   = freq_params.get<Int>("OUT_MAX_STEPS");
  m_out_frequency   = freq_params.get<Int>("OUT_N");
  m_out_units       = freq_params.get<std::string>("OUT_OPTION");

  // Create map of fields in this output with the field_identifier in the field repository.
  auto& var_params = m_params.sublist("FIELDS");
  for (int var_i=0; var_i<var_params.get<Int>("Number of Fields");++var_i)
  {
    /* Determine the variable name */
    std::string var_name = var_params.get<std::string>(ekat::util::strint("field",var_i+1));
    /* Find the <FieldIdentifier,Field> pair that corresponds with this variable name on the appropriate grid */
    add_identifier(field_repo,var_name);
  }

  // Register new netCDF file for output.
  register_outfile(m_filename);

  // Register dimensions with netCDF file.
  for (auto it : m_dims)
  {
    register_dimension(m_filename,it.first,it.first,it.second);
  }
  register_dimension(m_filename,"time","time",0);  // Note that time has an unknown length, setting the "length" to 0 tells the interface to set this dimension as having an unlimited length, thus allowing us to write as many timesnaps to file as we desire.
  
  // Register variables with netCDF file.  Must come after dimensions are registered.
  register_variables();
  set_degrees_of_freedom();

  // Finish the definition phase for this file.
  eam_pio_enddef  (m_filename); 
  m_is_init = true;

} // init
/* ---------------------------------------------------------- */
inline void AtmosphereOutput::run(const FieldRepository<Real, device_type>& field_repo, const Real time) 
{
  using namespace scream;
  using namespace scream::scorpio;
  using Int = scream::Int;
  using Real = scream::Real;


  if (!m_is_init) {init(field_repo);} // If m_is_init is false than the previous step was at max_n_steps.  Need a new file.

  pio_update_time(m_filename,time);
  // Cycle through all fields in this output file, grab the view and write to file.
  for (auto const& f_map : m_fields)
  {
    grid_write_data_array(m_filename,f_map.first,m_dofs.at(f_map.first)[0],field_repo.get_field(f_map.second).get_view().data());
  } 
  sync_outfile(m_filename);

  // Check if the maximum number of steps has been reached.  If so, close this file and flag for a new one if needed.
  m_status["Run"] += 1;
  if (m_status["Run"]% m_out_max_steps==0)
  {
    eam_pio_closefile(m_filename);
    m_is_init = false;
  }
} // run
/* ---------------------------------------------------------- */
inline void AtmosphereOutput::finalize() 
{
  using namespace scream;
  using namespace scream::scorpio;
  if (m_is_init) {eam_pio_closefile(m_filename);} // if m_is_init is true then the file was not closed during run step, close it now.
  m_status["Finalize"] += 1;
} // finalize
/* ---------------------------------------------------------- */
inline void AtmosphereOutput::check_status()
{
  printf("IO Status for Rank %5d, File - %.40s: (Init: %2d), (Run: %5d), (Finalize: %2d)\n",m_comm.rank(),m_filename.c_str(),m_status["Init"],m_status["Run"],m_status["Finalize"]);
} // check_status
/* ---------------------------------------------------------- */
inline void AtmosphereOutput::add_identifier(const FieldRepository<Real, device_type> &field_repo, const std::string name)
{
  for (auto myalias=field_repo.begin();myalias!=field_repo.end();++myalias)
  {
    if (myalias->first == name) {
      auto& map_it = myalias->second;
      // map_it is now a map<FieldIdentifier, Field>. Loop over it, and pick all fields defined on the output grid
      // map_it<Field::FieldHeader::FieldIdentifier,Field>
      for (const auto& it : map_it) {
        auto& field_id = it.first;
        if (it.first.get_grid_name()!=m_grid_name) {
          continue;
        }
        // ok, this field is on the correct mesh. add it to pio output
        m_fields.emplace(name,field_id);
        // check to see if all the dims for this field are already set to be registered.
        std::map<std::string,Int>::iterator tag_loc;
        for (int ii=0; ii<field_id.get_layout().rank(); ++ii)
        {
          // check tag against m_dims map.  If not in there, then add it.
          auto& tag = field_id.get_layout().tags()[ii];
          tag_loc = m_dims.find(tag2string(tag));
          if (tag_loc == m_dims.end()) 
          { 
            m_dims.emplace(std::make_pair(tag2string(tag),field_id.get_layout().dim(ii)));
          } else {  
            EKAT_REQUIRE_MSG(m_dims.at(tag2string(tag))==field_id.get_layout().dim(ii),
              "Error! Dimension " + tag2string(tag) + " on field " + name + " has conflicting lengths");
          }
        }
        return;
      }
      // Got this far means we found the field but not on the correct grid.
      EKAT_REQUIRE_MSG(true,"Error! Field " + name + " found in repo, but not on grid " + m_grid_name + ".\n");
    }
  }
  // Got this far means that the field was never found in the repo.
  EKAT_REQUIRE_MSG(true,"Error! Field " + name + " not found in repo.\n");
} // add_identifier
/* ---------------------------------------------------------- */
inline void AtmosphereOutput::register_variables()
{
  using namespace scream;
  using namespace scream::scorpio;
  // Cycle through all fields and register.
  for (auto it : m_fields)
  {
    auto  name = it.first;
    auto& fid  = it.second;
    // Determine the IO-decomp and construct a vector of dimension ids for this variable:
    std::string io_decomp_tag = "Real";  // Note, for now we only assume REAL variables.  This may change in the future.
    std::vector<std::string> vec_of_dims;
    for (int ii=0;ii<fid.get_layout().rank();++ii) {
      io_decomp_tag += "-" + tag2string(fid.get_layout().tag(ii)); // Concatenate the dimension string to the io-decomp string
      vec_of_dims.push_back(tag2string(fid.get_layout().tag(ii))); // Add dimensions string to vector of dims.
    }
    io_decomp_tag += "-time";  // TODO: Do we expect all vars to have a time dimension?  If not then how to trigger?  Should we register dimension variables (such as ncol and lat/lon) elsewhere in the dimension registration?  These won't have time.
    std::reverse(vec_of_dims.begin(),vec_of_dims.end()); // TODO: Reverse order of dimensions to match flip between C++ -> F90 -> PIO, may need to delete this line when switching to fully C++/C implementation.
    vec_of_dims.push_back("time");  //TODO: See the above comment on time.
    register_variable(m_filename, name, name, vec_of_dims.size(), vec_of_dims, PIO_REAL, io_decomp_tag);  // TODO  Need to change dtype to allow for other variables.  Currently the field_repo only stores Real variables so it is not an issue, but in the future if non-Real variables are added we will want to accomodate that.
  }
  // Finish by registering time as a variable.  TODO: Should this really be something registered during the reg. dimensions step? 
  register_variable(m_filename,"time","time",1,{"time"},  PIO_REAL,"t");
} // register_variables
/* ---------------------------------------------------------- */
inline void AtmosphereOutput::set_degrees_of_freedom()
{
  using namespace scream;
  using namespace scream::scorpio;
  // Cycle through all fields and set dof.
  for (auto it : m_fields)
  {
    auto  name = it.first;
    auto& fid  = it.second;
    Int dof_len, dof_start, dof_stop;
    dof_len = fid.get_layout().size()/m_comm.size();
    Int extra_procs = fid.get_layout().size() % dof_len;
    if (extra_procs>0) { dof_len += 1; }
    dof_start = m_comm.rank()*dof_len;
    if (m_comm.rank() == m_comm.size()-1) { dof_len = std::max((Int) fid.get_layout().size()-dof_start,0); }
    dof_stop = std::min(dof_start + dof_len-1,fid.get_layout().size()-1);
    dof_len = std::max(dof_stop-dof_start + 1,0);
    // determine which dof's this rank is responsible for.  TODO: This should really all be done via the grid manager which maps to global-ids.
    //                                                      note: global-ids from grid manager are for columns, but PIO thinks of them in 
    //                                                      reference to a flattened array of N-dimensions.  We will need to handle the other dimensions properly. 
    Int var_dof[dof_len];
    for (Int ii=0;ii<dof_len;++ii) {var_dof[ii] = dof_start+ii;}
    set_dof(m_filename,name,dof_len,var_dof);
    std::vector<Int> dof_vec = {dof_len,dof_start,dof_stop};
    m_dofs.emplace(std::make_pair(name,dof_vec));
  }
  // Set degree of freedom for "time"
  set_dof(m_filename,"time",0,0);
  /* TODO: 
   * Adjust DOF to accomodate packing for fields 
   * Gather DOF info directly from grid manager
  */
} // set_degrees_of_freedom
/* ---------------------------------------------------------- */
} //namespace scream
#endif // SCREAM_SCORPIO_HPP
