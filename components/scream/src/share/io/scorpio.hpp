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
 * 1) Gather dimension data from grid manager or field repo.  We may need to use field repo since grid manager
 *    is unaware of say the number of components, for example.
 * 2) Clean-up init step, perhaps with subfunctions.
 * 3) Gather degree's of freedom from grid manager. Need to know the global-id index for
 *    each value in a view, in order if flattened.
 * 4) Handle MPI communication to allow for PIO_STRIDE (i.e. fewer pio ranks than total ranks).
 * 5) Make it so that all dimensions that are registered are automatically registered as variables and the values written (only once, not with time).
 * 6) Write Restart output.  Extra-savy would be to use DAG from AD to determine which fields were essential and just write those.
 * 7) Create an AtmosphereInput class built on these same concepts.  Can we make a master SCORPIO class with output and input as sub-classes?
 * 8) Create average, min, max and instantaneous options for output.
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

  void init(const FieldRepository<Real, device_type>& field_repo);
  void run(const FieldRepository<Real, device_type>& field_repo, const Real time);
  void finalize();

  // Helper Functions
  void check_status();
  std::map<std::string,Int> get_status() const { return m_status; }

  AtmosphereOutput(const ekat::Comm& comm, const ekat::ParameterList& params)
  {
    m_comm   = comm;
    m_params = params;
  }
protected:
  ekat::ParameterList m_params;
  ekat::Comm          m_comm;
  
  std::string m_filename;
  std::string m_avg_type;
  std::string m_grid_name;

  Int m_out_max_steps;
  Int m_out_frequency;
  std::string m_out_units;

  std::map<std::string,FieldIdentifier> m_fields;
  std::map<std::string,std::vector<Int>>          m_dofs;

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

  m_filename = m_params.get<std::string>("FILENAME")+"_"+std::to_string(m_status["Run"])+".nc";  //TODO: The filename should be treated as a prefix to enable multiple files for the same control.  Like in the case of monthly output with 1 snap/file.
  EKAT_REQUIRE_MSG(m_status.at("Init")==0,"Error! File " + m_filename + " has already been initialized.\n");
  m_status["Init"] += 1;
  printf("Init for file %s\n",m_filename.c_str());

  // Parse the yaml file that controls this output instance.
  m_avg_type        = m_params.get<std::string>("AVERAGING TYPE");
  m_grid_name       = m_params.get<std::string>("GRID");
  auto& freq_params = m_params.sublist("FREQUENCY");
  m_out_max_steps   = freq_params.get<Int>("OUT_MAX_STEPS");
  m_out_frequency   = freq_params.get<Int>("OUT_N");
  m_out_units       = freq_params.get<std::string>("OUT_OPTION");

  // Register new netCDF file for output.
  register_outfile(m_filename);

  // Register dimensions with netCDF file.
  //   ->  register_dimension(filename,dim_shortname,dim_longname,dim_length);
  /* Use grids manager to get num_dofs, use param list to get ylen and zlen as "NLVL" and "NTRACERS" */
  int xlen=10, ylen=5, zlen=2;  //TODO: Switch this to work with a "grid"
  register_dimension(m_filename,"COL","horizontal distance",xlen);
  register_dimension(m_filename,"VL","vertical distance",ylen);
  register_dimension(m_filename,"CMP","height",zlen);
  register_dimension(m_filename,"time","time",0);  // Note that time has an unknown length, setting the "length" to 0 tells the interface to set this dimension as having an unlimited length, thus allowing us to write as many timesnaps to file as we desire.
  
  // Register variables with netCDF file.  Must come after dimensions are registered.
  // Additionally, determine the degrees-of-freedom for this variable that this MPI 
  // rank is responsible for writing.
  //    ->  register_variable(filename, var_shortname, var_longname, var_num_of_dims (rank of field), vector_of_dimensions, var_dtype, io_decomp_tag);
  //    field_utils.hpp -> get_layout_type
  /* Obtain the list of fields for this output stream from parameter list */
  auto& var_params = m_params.sublist("FIELDS");
  /* Cycle through the fields and register the variables */
  for (int var_i=0; var_i<var_params.get<Int>("Number of Fields");++var_i)
  {
    /* Determine the variable name */
    std::string name = var_params.get<std::string>(ekat::util::strint("field",var_i+1));
    /* Search field repository for this variable, if not found, throw an error */
    bool found_field = false;
    for (auto myalias=field_repo.begin();myalias!=field_repo.end();++myalias)
    {
      if (myalias->first == name) {
        auto& map_it = myalias->second;
        // map_it is now a map<FieldIdentifier, Field>. Loop over it, and pick all fields defined on the output grid
        // map_it<Field::FieldHeader::FieldIdentifier,Field>
        bool found_one = false;
        for (const auto& it : map_it) {
          auto& field_id = it.first;
          if (it.first.get_grid_name()!=m_grid_name) {
            continue;
          }
          // ok, this field is on the correct mesh. add it to pio output
          m_fields.emplace(name,field_id);
          std::string io_decomp_tag = "Real";
          std::vector<std::string> vec_of_dims;
          for (int kk=0;kk<field_id.get_layout().rank();++kk) {
            auto& l_layout = field_id.get_layout();
            const auto& l_tag = l_layout.tag(kk);
            std::string l_dimname = tag2string(l_tag);
            io_decomp_tag += "-" + tag2string(field_id.get_layout().tag(kk));
            vec_of_dims.push_back(tag2string(field_id.get_layout().tag(kk)));
          }
          io_decomp_tag += "-time";
          std::reverse(vec_of_dims.begin(),vec_of_dims.end()); // TODO: Reverse order of dimensions to match flip between C++ -> F90 -> PIO, may need to delete this line when switching to fully C++/C implementation.
          vec_of_dims.push_back("time");
          found_field = true;
          found_one = true;
          register_variable(m_filename, name, name, field_id.get_layout().rank()+1, vec_of_dims, PIO_REAL, io_decomp_tag);  // TODO  Need to change dtype to allow for other variables.  Currently the field_repo only stores Real variables so it is not an issue, but in the future if non-Real variables are added we will want to accomodate that.
          // Get and then set the DOF's for this rank.
          // note: m_comm.rank() and m_comm.size() = my mpi rank and the total number of mpi ranks, respectively.
          // TODO: This DOF information should be available from the grid manager given the grid and field layout.  So fix this simple implementation.  
          // Note: When we begin to allow for # PIO procs < Total Procs then we will need to pass DOF information from multiple procs to the writer
          //       procs.
          Int dof_len, dof_start, dof_stop;
          dof_len = field_id.get_layout().size()/m_comm.size();
          Int extra_procs = field_id.get_layout().size() % dof_len;
          if (extra_procs>0) { dof_len += 1; }
          dof_start = m_comm.rank()*dof_len;
          if (m_comm.rank() == m_comm.size()-1) { dof_len = std::max((Int) field_id.get_layout().size()-dof_start,0); }
          dof_stop = std::min(dof_start + dof_len-1,field_id.get_layout().size()-1);
          dof_len = std::max(dof_stop-dof_start + 1,0);
          // determine which dof's this rank is responsible for
          Int var_dof[dof_len];
          for (Int ii=0;ii<dof_len;++ii) {var_dof[ii] = dof_start+ii;}
          set_dof(m_filename,name,dof_len,var_dof);
          std::vector<Int> dof_vec = {dof_len,dof_start,dof_stop};
          m_dofs.emplace(std::make_pair(name,dof_vec));
        }
        EKAT_REQUIRE_MSG(found_one,"Error! Field " + name + " found in repo, but not on grid " + m_grid_name + ".\n");
      } 
    }
    EKAT_REQUIRE_MSG(found_field,"Error! Field " + name + " not found in repo.\n");
  } // output_names
  // Register time as a variable TODO: Is there a better way then hard-code?
  register_variable(m_filename,"time","time",1,{"time"},  PIO_REAL,"t");
  set_dof(m_filename,"time",0,0);
  /* TODO: 
 * Adjust DOF to accomodate packing for fields 
 * Gather DOF info directly from grid manager
 * Clean-up this class, maybe add some sub-functions so that the init function is leaner.
  */

  // Finish the definition phase for this file.
  eam_pio_enddef  (m_filename); 

} // init
/* ---------------------------------------------------------- */
inline void AtmosphereOutput::run(const FieldRepository<Real, device_type>& field_repo, const Real time) 
{
  using namespace scream;
  using namespace scream::scorpio;
  using Int = scream::Int;
  using Real = scream::Real;

  m_status["Run"] += 1;
  // Check if the maximum number of steps has been reached.  If so, close this file and create a new one.
  if (m_status["Run"]% m_out_max_steps==1)
  {
    printf("ASD - we reached max steps\n");
    m_status["Init"] = 0;
    m_status["Run"]  = 1;
    eam_pio_closefile(m_filename);
    init(field_repo);
  }

  pio_update_time(m_filename,time);
  // Cycle through all fields in this output file, grab the view and write to file.
  for (auto const& f_map : m_fields)
  {
    grid_write_data_array(m_filename,f_map.first,m_dofs.at(f_map.first)[0],field_repo.get_field(f_map.second).get_view().data());
  } 
  sync_outfile(m_filename);
} // run
/* ---------------------------------------------------------- */
inline void AtmosphereOutput::finalize() 
{
  using namespace scream;
  using namespace scream::scorpio;
  eam_pio_closefile(m_filename);
  m_status["Finalize"] += 1;
} // finalize
/* ---------------------------------------------------------- */
inline void AtmosphereOutput::check_status()
{
  printf("IO Status for Rank %5d: %.40s - (%1d, %5d, %2d)\n",m_comm.rank(),m_filename.c_str(),m_status["Init"],m_status["Run"],m_status["Finalize"]);
} // check_status
} //namespace scream
#endif // SCREAM_SCORPIO_HPP
