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
  void run();
  void finalize();

  AtmosphereOutput(const ekat::Comm& comm, const ekat::ParameterList& params)
  {
    m_comm   = comm;
    m_params = params;
  }
protected:
  ekat::ParameterList m_params;
  ekat::Comm          m_comm;
  
  std::string m_avg_type;
  std::string m_grid_name;

  Int m_out_max_steps;
  Int m_out_frequency;
  std::string m_out_units;

private:
  std::string m_filename;

#if 0
  std::string filename;
  std::string output_type;
  std::string grid_type;
  Int Out_N;
  std::string Out_Freq;
  Int Out_Max_N;

  std::vector<std::string> List_of_Fields;
  std::vector<std::string> List_of_Dimensions;

  void set_field_dof( const std::string name, const Int dof_len, const Int dof_start, const Int dof_stop );
  void set_comm( const Int comm );
  void set_decomp(const std::string name, const std::string decomp);
  void set_dtype(const std::string name, const Int);

  std::string get_decomp(const std::string name);
  Int get_dtype(const std::string name);
  std::vector<Int> get_field_dof( const std::string name);

  void set_dimensions (const std::vector<std::string> dim_array_in)
  {
    for (std::string dim : dim_array_in)
    {
      List_of_Dimensions.push_back(dim);
    }
  }

  /* AtmosphereOutput constructor takes an output YAML file as input to set all local values */
  AtmosphereOutput(const char* input_yaml_file) 
  {
  
    ekat::ParameterList file_params("Output File Parameters");
    parse_yaml_file (input_yaml_file, file_params); // Parse the YAML file associated with that name

    filename = file_params.get<std::string>("FILENAME");
    output_type = file_params.get<std::string>("AVERAGING TYPE");
    grid_type = file_params.get<std::string>("GRID");
    
    /* frequency control */
    auto& freq_params = file_params.sublist("FREQUENCY");
    Out_Freq  = freq_params.get<std::string>("OUT_OPTION");
    Out_N     = freq_params.get<Int>("OUT_N");
    Out_Max_N = freq_params.get<Int>("OUT_MAX_STEPS"); 
 
    /* Variables for output */  //TODO: switch to string array for field entries in YAML and loop over size of array.
    auto& field_params = file_params.sublist("FIELDS");
    const auto& num_of_fields = field_params.get<Int>("Number of Fields"); 
    for (Int it = 1;it <= num_of_fields;it++) {
      auto& fieldname = field_params.get<std::string>("f"+std::to_string(it));
      List_of_Fields.push_back(fieldname);
      ekat::ParameterList loc_param(fieldname);
      field_outputs[fieldname] = loc_param;
    }
  }

private:
  Int myrank;
  Int numranks;

  std::map<std::string,std::string> var_iodecomp;
  std::map<std::string,Int> var_dtype;

  std::map<std::string,std::vector<Int>> var_dof;
  std::map<std::string,std::vector<Int>> dof_vec;
  std::map<std::string,ekat::ParameterList> field_outputs;
#endif
}; // Class AtmosphereOutput

// ====================== IMPLEMENTATION ===================== //
inline void AtmosphereOutput::init(const FieldRepository<Real, device_type>& field_repo /* grids_manager reference */) 
{
  using namespace scream;
  using namespace scream::scorpio;

  // Parse the yaml file that controls this output instance.
  m_params.print();

  m_avg_type        = m_params.get<std::string>("AVERAGING TYPE");
  m_grid_name       = m_params.get<std::string>("GRID");
  auto& freq_params = m_params.sublist("FREQUENCY");
  m_out_max_steps   = freq_params.get<Int>("OUT_MAX_STEPS");
  m_out_frequency   = freq_params.get<Int>("OUT_N");
  m_out_units       = freq_params.get<std::string>("OUT_OPTION");


  // Register new netCDF file for output.
  m_filename = m_params.get<std::string>("FILENAME");  //TODO: The filename should be treated as a prefix to enable multiple files for the same control.  Like in the case of monthly output with 1 snap/file.
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
  auto& var_params = m_params.sublist("FIELDS");
//  for (const auto& name : m_output_names) 
//  {
    std::string name = "data_3d";
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
          printf("ASD grid = %s\n",it.first.get_grid_name().c_str());
          if (it.first.get_grid_name()!=m_grid_name) {
            continue;
          }
          // ok, this field is on the correct mesh. add it to pio output
          std::string io_decomp_tag = "Real";
          //const char* vector_of_dims[field_id.get_layout().rank()];
          std::vector<std::string> vec_of_dims;
          for (int kk=0;kk<field_id.get_layout().rank();++kk) {
            auto& l_layout = field_id.get_layout();
            const auto& l_tag = l_layout.tag(kk);
            std::string l_dimname = tag2string(l_tag);
            io_decomp_tag += "-" + tag2string(field_id.get_layout().tag(kk));
          //  vector_of_dims[kk] = tag2string(field_id.get_layout().tag(kk)).c_str();
            vec_of_dims.push_back(tag2string(field_id.get_layout().tag(kk)));
            printf(" --- %s = %d\n",field_id.get_id_string().c_str(),l_layout.dims()[kk]);
          }
          for (int kk = 0;kk<field_id.get_layout().rank();++kk) { printf("ASD - var_dims first: %d -> %s\n",kk,vec_of_dims[kk].c_str()); }
          printf("ASD - Found field %s -> decomp = %s\n",name.c_str(),io_decomp_tag.c_str());
          found_field = true;
          found_one = true;
          // AaronDonahue TODO: What you are doing here is now that you have the field identifier, gather the info that you need to register the variable with IO.
          register_variable(m_filename, name, name, field_id.get_layout().rank(), vec_of_dims, PIO_REAL, io_decomp_tag);  // TODO  Need to change dtype to allow for other variables
          // Get and then set the DOF's for this rank.
          // note: m_comm.rank() and m_comm.size() = my mpi rank and the total number of mpi ranks, respectively.
          Int dof_len, dof_start, dof_stop;
          dof_len = field_id.get_layout().rank()/m_comm.size();
          Int extra_procs = field_id.get_layout().rank() % dof_len;
          if (extra_procs>0) { dof_len += 1; }
          dof_start = m_comm.rank()*dof_len;
          if (m_comm.rank() == m_comm.size()-1) { dof_len = std::max((Int) field_id.get_layout().rank()-dof_start,0); }
          dof_stop = std::min(dof_start + dof_len-1,field_id.get_layout().rank()-1);
          dof_len = std::max(dof_stop-dof_start + 1,0);
          Int var_dof[dof_len];
          for (Int ii=0;ii<dof_len;++ii) {var_dof[ii] = dof_start+ii;}
          set_dof(m_filename,name,dof_len,var_dof);
        }
        //throw error if found_one is false maybe?
      } 
    }
   /*   auto& map_it = std::find(repo_start,repo_end,name); */
      EKAT_REQUIRE_MSG(found_field,"Error! Field not found in repo.\n");
//  } // output_names

  
  /* TODO: Adjust this to accomodate packing for fields */
  //get_dof(length_of_var-1d, myrank, numranks, OUT length_of_dof, OUT dof_start, OUT dof_stop);
  //set_dof(filename,var_shortname,length_of_dof, array_of_dof_indices);
//  FieldIdentifier& field_repo.get_field("x");

  // Finish the definition phase for this file.
  eam_pio_enddef  (m_filename); 
}
/* ---------------------------------------------------------- */
inline void AtmosphereOutput::run() 
{
  // Do nothing
}
/* ---------------------------------------------------------- */
inline void AtmosphereOutput::finalize() 
{
  // Do nothing
}
#if 0
inline void AtmosphereOutput::set_comm( const Int comm )
{
  MPI_Comm_rank(comm, &myrank);  // Store rank and total number of ranks for determining which chunk of global array this rank is responsible for reading.
  MPI_Comm_size(comm, &numranks);
}

//inline void AtmosphereOutput::set_field( )
//{
//
//}

inline void AtmosphereOutput::set_field_dof( const std::string name, const Int dof_len, const Int dof_start, const Int dof_stop )
{
  var_dof[name] = {dof_len,dof_start};
  for (Int ii=dof_start;ii<=dof_stop;ii++) { dof_vec[name].push_back(ii+1); }
}
inline void AtmosphereOutput::set_decomp(const std::string name, const std::string decomp)
{
  var_iodecomp[name] = decomp;
}
inline void AtmosphereOutput::set_dtype(const std::string name, const Int dtype)
{
  var_dtype[name] = dtype;
}

inline std::vector<Int> AtmosphereOutput::get_field_dof( const std::string name )
{
  return var_dof[name];
}
inline std::string AtmosphereOutput::get_decomp(const std::string name)
{
  return var_iodecomp[name];
}
inline Int AtmosphereOutput::get_dtype(const std::string name)
{
  return var_dtype[name];
}
#endif
} //namespace scream
#endif // SCREAM_SCORPIO_HPP
