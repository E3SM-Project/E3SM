#ifndef SCREAM_IO_HPP
#define SCREAM_IO_HPP

#include "ekat/ekat_parameter_list.hpp"
#include "ekat/std_meta/ekat_std_enable_shared_from_this.hpp"
#include "ekat/mpi/ekat_comm.hpp"
#include "ekat/std_meta/ekat_std_utils.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"

class AtmosphereOutput 
{

public:
  using Int = scream::Int;

  virtual ~AtmosphereOutput () = default;

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
//    file_params.print();  // Print what is in this file: TODO Delete eventually once we know things are working.

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
//    for (Int it=0;it<List_of_Fields.size();it++) { printf("-- %s\n",List_of_Fields[it].c_str()); }
  }

private:
  Int myrank;
  Int numranks;

  std::map<std::string,std::string> var_iodecomp;
  std::map<std::string,Int> var_dtype;

  std::map<std::string,std::vector<Int>> var_dof;
  std::map<std::string,std::vector<Int>> dof_vec;
  std::map<std::string,ekat::ParameterList> field_outputs;
}; // Class AtmosphereOutput

// ====================== IMPLEMENTATION ===================== //
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
#endif // SCREAM_IO_HPP
