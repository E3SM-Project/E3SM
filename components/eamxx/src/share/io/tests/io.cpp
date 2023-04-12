#include <catch2/catch.hpp>
#include <memory>

#include "share/io/scream_output_manager.hpp"
#include "share/io/scorpio_output.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/io/scream_scorpio_interface.hpp"

#include "diagnostics/register_diagnostics.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/grid/point_grid.hpp"

#include "share/field/field_identifier.hpp"
#include "share/field/field_header.hpp"
#include "share/field/field.hpp"
#include "share/field/field_manager.hpp"

#include "share/util/scream_time_stamp.hpp"
#include "share/scream_types.hpp"
#include "scream_config.h"

#include "ekat/util/ekat_units.hpp"
#include "ekat/ekat_parameter_list.hpp"
#include "ekat/ekat_scalar_traits.hpp"
#include "ekat/ekat_assert.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"
#include "ekat/util/ekat_string_utils.hpp"

namespace {

using namespace scream;
using namespace ekat::units;
// Make sure packsize isn't bigger than the packsize for this machine, but not so big that we end up with only 1 pack.
const int packsize = SCREAM_SMALL_PACK_SIZE;
using Pack         = ekat::Pack<Real,packsize>;

using KT = KokkosTypes<DefaultDevice>;
template <typename S>
using view_2d = typename KT::template view_2d<S>;
using FID = FieldIdentifier;
using FL  = FieldLayout;

std::shared_ptr<FieldManager>
get_test_fm(std::shared_ptr<const AbstractGrid> grid);

std::shared_ptr<GridsManager>
get_test_gm(const ekat::Comm& io_comm, const Int num_gcols, const Int num_levs);

ekat::ParameterList get_in_params(const std::string& type,
                                  const std::string& freq,
                                  const ekat::Comm& comm,
                                  const std::string& t_first_write);

view_2d<Real>::HostMirror
get_diagnostic_input(const ekat::Comm& comm, const std::shared_ptr<GridsManager>& gm,
                     const int time_index, const std::string& filename);

int get_current_t(const int tt, const int dt, const int freq,  const std::string& frequency_units);

Real generate_data_xy(const Int time, const Int i, const Int j);
Real check_data_xy   (const Int time, const Int dt, const Int i, const Int j, const std::string& avg_Type);

class DiagTest : public AtmosphereDiagnostic
{
public:
  DiagTest (const ekat::Comm& comm, const ekat::ParameterList&params)
    : AtmosphereDiagnostic(comm,params)
  {
    //Do nothing
  }

  std::string name() const { return "DiagnosticTest"; }

  void set_grids (const std::shared_ptr<const GridsManager> gm) {
    using namespace ekat::units;
    using namespace ShortFieldTagsNames;
    using FL = FieldLayout;

    const auto grid = gm->get_grid("Physics");
    const auto& grid_name = grid->name();
    m_num_cols  = grid->get_num_local_dofs(); // Number of columns on this rank
    m_num_levs  = grid->get_num_vertical_levels();  // Number of levels per column

    std::vector<FieldTag> tag_2d = {COL,LEV};
    std::vector<Int>     dims_2d = {m_num_cols,m_num_levs};
    FL lt( tag_2d, dims_2d );

    constexpr int ps = Pack::n;
    add_field<Required>("field_packed",lt,kg/m,grid_name,ps);

    // We have to initialize the m_diagnostic_output:
    FieldIdentifier fid (name(), lt, K, grid_name);
    m_diagnostic_output = Field(fid);
    auto& C_ap = m_diagnostic_output.get_header().get_alloc_properties();
    C_ap.request_allocation(ps);
    m_diagnostic_output.allocate_view();
  }

protected:

  void compute_diagnostic_impl () {
    const auto& v_A  = get_field_in("field_packed").get_view<const Real**,Host>();
    auto v_me = m_diagnostic_output.get_view<Real**,Host>();
    // Have the dummy diagnostic just manipulate the field_packed data arithmetically.  In this case (x*3 + 2)
    for (int i=0; i<m_num_cols; ++i) {
      for (int k=0; k<m_num_levs; ++k) {
        v_me(i,k) = v_A(i,k)*3.0 + 2.0;
      }
    }
    m_diagnostic_output.sync_to_dev();
  }

  void initialize_impl (const RunType /* run_type */ ) {
    m_diagnostic_output.get_header().get_tracking().update_time_stamp(timestamp());
  }

  // Clean up
  void finalize_impl ( /* inputs */ ) {}

  // Internal variables
  int m_num_cols, m_num_levs;

};

/*===================================================================================================*/
void run_multisnap(const std::string& output_freq_units) {
  const std::string output_type = "multisnap";
  ekat::Comm io_comm(MPI_COMM_WORLD);  // MPI communicator group used for I/O set as ekat object.
  Int num_gcols = 2*io_comm.size();
  Int num_levs = 2 + SCREAM_SMALL_PACK_SIZE;

  // Initialize the pio_subsystem for this test:
  MPI_Fint fcomm = MPI_Comm_c2f(io_comm.mpi_comm());  // MPI communicator group used for I/O.  In our simple test we use MPI_COMM_WORLD, however a subset could be used.
  scorpio::eam_init_pio_subsystem(fcomm);   // Gather the initial PIO subsystem data creater by component coupler

  // First set up a field manager and grids manager to interact with the output functions
  auto gm = get_test_gm(io_comm,num_gcols,num_levs);
  auto grid = gm->get_grid("Point Grid");
  int num_lcols = grid->get_num_local_dofs();

  // Need to register the Test Diagnostic so that IO can access it
  auto& diag_factory = AtmosphereDiagnosticFactory::instance();
  diag_factory.register_product("DummyDiagnostic",&create_atmosphere_diagnostic<DiagTest>);
  diag_factory.register_product("FieldAtLevel",&create_atmosphere_diagnostic<DiagTest>);

  // Construct a timestamp
  util::TimeStamp t0 ({2000,1,1},{0,0,0});
  IOControl io_control; // Needed for testing input.
  io_control.timestamp_of_last_write = t0;
  io_control.nsamples_since_last_write = 0;
  io_control.frequency_units         = output_freq_units;
  std::vector<std::string> output_stamps; 

  // Gather testing data based on the output frequency units
  ekat::ParameterList control_params;
  ekat::parse_yaml_file("io_test_control.yaml",control_params);
  EKAT_REQUIRE_MSG(control_params.isSublist(output_freq_units),"ERROR! output frequency units " + output_freq_units + " does not have an entry in io_test_control.yaml");
  const auto& freq_params = control_params.sublist(output_freq_units);
  const Int max_steps = freq_params.get<int>("Max Steps");
  const Int dt        = freq_params.get<int>("dt");
  {
    util::TimeStamp time  = t0;
    util::TimeStamp time0(t0);

    // Re-create the fm anew, so the fields are re-inited for each output type
    auto field_manager = get_test_fm(grid);
    field_manager->init_fields_time_stamp(t0);

    // Set up parameter list control for output
    ekat::ParameterList params;
    ekat::parse_yaml_file("io_test_" + output_type + ".yaml",params);
    params.set<std::string>("Floating Point Precision","real");
    auto& params_sub = params.sublist("output_control");
    params_sub.set<std::string>("frequency_units",output_freq_units);
    io_control.frequency = params_sub.get<int>("Frequency");

    // Set up output manager.
    OutputManager om;
    om.setup(io_comm,params,field_manager,gm,t0,t0,false);
    io_comm.barrier();

    const auto& out_fields = field_manager->get_groups_info().at("output");
    using namespace ShortFieldTagsNames;
    // Time loop
    for (Int ii=0;ii<max_steps;++ii) {
      time += dt;
      Int time_in_sec = time.seconds_from(time0);
      // Update the fields
      for (const auto& fname : out_fields->m_fields_names) {
        auto f  = field_manager->get_field(fname);
        f.sync_to_host();
        auto fl = f.get_header().get_identifier().get_layout();
        switch (fl.rank()) {
          case 1:
            {
              auto v = f.get_view<Real*,Host>();
              // If field tag is columns use generate_x, levels use generate_y
              for (int i=0; i<fl.dim(0); ++i) {
                if (fl.has_tag(COL)) {
                  v(i) = generate_data_xy(time_in_sec,i,0);
                } else {
                  v(i) = generate_data_xy(time_in_sec,0,i);
                }
              }
            }
            break;
          case 2:
            {
              auto v = f.get_view<Real**,Host>();
              for (int i=0; i<fl.dim(0); ++i) {
                for (int j=0; j<fl.dim(1); ++j) {
                  v(i,j) = generate_data_xy(time_in_sec,i,j);
                }
              }
            }
            break;
          default:
            EKAT_ERROR_MSG ("Error! Unexpected field rank.\n");
        }
        f.sync_to_dev();
      }

      // Run the output manager for this time step
      om.run(time);
      if (io_control.is_write_step(time)) {
        output_stamps.push_back(time.to_string());
        io_control.nsamples_since_last_write = 0;
        io_control.timestamp_of_last_write = time;
      }
    }

    // Finalize the output manager (close files)
    om.finalize();
  }

  // Get a fresh new field manager
  auto field_manager = get_test_fm(grid);

  // We can use the produced output files to simultaneously check output quality and the
  // ability to read input.
  // NOTE: single-snap file start outputing at the final time only,
  //       while multi-snap starts to output at the first time step
  REQUIRE(output_stamps.size()>0);
  auto input_params = get_in_params(output_type,output_freq_units, io_comm,output_stamps.front());
  const Real tol = std::numeric_limits<Real>::epsilon();

  // TODO: Create a small nc dummy file and a separate unit test which tests all input functions.
  // Test that pio_inq_dimlen is correct, using a file from one of the above parameter lists.
  {
    auto test_filename = input_params.get<std::string>("Filename");
    scorpio::register_file(test_filename,scorpio::Read);
    Int test_gcols_len = scorpio::get_dimlen_c2f(test_filename.c_str(),"ncol");
    REQUIRE(test_gcols_len==num_gcols);
    scorpio::eam_pio_closefile(test_filename);
  }

  auto f1 = field_manager->get_field("field_1");
  auto f2 = field_manager->get_field("field_2");
  auto f3 = field_manager->get_field("field_3");
  auto f4 = field_manager->get_field("field_packed");

  // Add field @ level. Note: cannot use subfield, since
  //  - we don't want to alias f3
  //  - we cannot use subfield along 2nd dim for rank-2 fields
  using namespace ShortFieldTagsNames;
  const auto& f3_lt = f3.get_header().get_identifier().get_layout();
  EKAT_REQUIRE_MSG (f3_lt.rank()>0,"WHAT?!?\n");
  const auto units = f3.get_header().get_identifier().get_units();
  auto f3_tom = Field(FID("field_3@tom",f3_lt.strip_dim(LEV),units,grid->name()));
  auto f3_bot = Field(FID("field_3@bot",f3_lt.strip_dim(LEV),units,grid->name()));
  auto f3_lev_2 = Field(FID("field_3@lev_2",f3_lt.strip_dim(LEV),units,grid->name()));
  f3_tom.allocate_view();
  f3_bot.allocate_view();
  f3_lev_2.allocate_view();
  field_manager->add_field(f3_tom);
  field_manager->add_field(f3_bot);
  field_manager->add_field(f3_lev_2);

  // Get host views
  auto f1_host = f1.get_view<Real*,Host>();
  auto f2_host = f2.get_view<Real*,Host>();
  auto f3_host = f3.get_view<Real**,Host>();
  auto f4_host = f4.get_view<Real**,Host>();
  auto f3_tom_host = f3_tom.get_view<Real*,Host>();
  auto f3_bot_host = f3_bot.get_view<Real*,Host>();
  auto f3_lev_2_host = f3_lev_2.get_view<Real*,Host>();

  // Check multisnap output; note, tt starts at 1 instead of 0 to follow netcdf time dimension indexing.
  AtmosphereInput multi_input(input_params,field_manager);
  for (int tt = 0; tt<std::min(max_steps,10); tt++) {
    multi_input.read_variables(tt);
    f1.sync_to_host();
    f2.sync_to_host();
    f3.sync_to_host();
    f4.sync_to_host();
    f3_tom.sync_to_host();
    f3_bot.sync_to_host();
    f3_lev_2.sync_to_host();

    int tt1 = tt + 1;
    // Here tt1 is the snap, we need to figure out what time that is in seconds given the frequency units
    const int current_t = get_current_t(tt1,dt,io_control.frequency,output_freq_units);
    // Check timestamp, note, with multisnap we also verify that we can grab the timestamp at location tt1
    {
      auto test_filename = input_params.get<std::string>("Filename");
      Real time_val = scorpio::read_time_at_index_c2f(test_filename.c_str(),tt1);
      Real time_in_days = current_t/86400.; // current_t is in seconds, need to convert to days.
      REQUIRE(time_val==time_in_days);
    }
    for (int ii=0;ii<num_lcols;++ii) {
      REQUIRE(std::abs(f1_host(ii)-check_data_xy(current_t,dt,ii,0,"instant"))<tol);
      for (int jj=0;jj<num_levs;++jj) {
        REQUIRE(std::abs(f3_host(ii,jj)-check_data_xy(current_t,dt,ii,jj,"instant"))<tol);
        REQUIRE(std::abs(f4_host(ii,jj)-check_data_xy(current_t,dt,ii,jj,"instant"))<tol);
      }

      REQUIRE (f3_tom_host(ii)==f3_host(ii,0));
      REQUIRE (f3_bot_host(ii)==f3_host(ii,num_levs-1));
      REQUIRE (f3_lev_2_host(ii)==f3_host(ii,2));
    }
    for (int jj=0;jj<num_levs;++jj) {
      REQUIRE(std::abs(f2_host(jj)-check_data_xy(current_t,dt,0,jj,"instant"))<tol);
    }
  }
  multi_input.finalize();
  // All Done 
  scorpio::eam_pio_finalize();
} // end function run_multisnap()
/*===================================================================================================*/
void run(const std::string& output_type,const std::string& output_freq_units) {
  ekat::Comm io_comm(MPI_COMM_WORLD);  // MPI communicator group used for I/O set as ekat object.
  Int num_gcols = 2*io_comm.size();
  Int num_levs = 2 + SCREAM_SMALL_PACK_SIZE;

  // Initialize the pio_subsystem for this test:
  MPI_Fint fcomm = MPI_Comm_c2f(io_comm.mpi_comm());  // MPI communicator group used for I/O.  In our simple test we use MPI_COMM_WORLD, however a subset could be used.
  scorpio::eam_init_pio_subsystem(fcomm);   // Gather the initial PIO subsystem data creater by component coupler

  // First set up a field manager and grids manager to interact with the output functions
  auto gm = get_test_gm(io_comm,num_gcols,num_levs);
  auto grid = gm->get_grid("Point Grid");
  int num_lcols = grid->get_num_local_dofs();

  // Need to register the Test Diagnostic so that IO can access it
  auto& diag_factory = AtmosphereDiagnosticFactory::instance();
  diag_factory.register_product("DummyDiagnostic",&create_atmosphere_diagnostic<DiagTest>);

  // Construct a timestamp
  util::TimeStamp t0 ({2000,1,1},{0,0,0});
  IOControl io_control; // Needed for testing input.
  io_control.timestamp_of_last_write = t0;
  io_control.nsamples_since_last_write = 0;
  io_control.frequency_units         = output_freq_units;
  std::vector<std::string> output_stamps; 

  // Gather testing data based on the output frequency units
  ekat::ParameterList control_params;
  ekat::parse_yaml_file("io_test_control.yaml",control_params);
  EKAT_REQUIRE_MSG(control_params.isSublist(output_freq_units),"ERROR! output frequency units " + output_freq_units + " does not have an entry in io_test_control.yaml");
  const auto& freq_params = control_params.sublist(output_freq_units);
  const Int max_steps = freq_params.get<int>("Max Steps");
  const Int dt        = freq_params.get<int>("dt");

  {
    util::TimeStamp time = t0;
    util::TimeStamp time0(t0);

    // Re-create the fm anew, so the fields are re-inited for each output type
    auto field_manager = get_test_fm(grid);
    field_manager->init_fields_time_stamp(t0);
    // Set up parameter list control for output
    ekat::ParameterList params;
    ekat::parse_yaml_file("io_test_template.yaml",params);
    params.set<std::string>("Averaging Type",output_type);
    auto& params_sub = params.sublist("output_control");
    params_sub.set<std::string>("frequency_units",output_freq_units);
    io_control.frequency = params_sub.get<int>("Frequency");
    // Set up output manager.
    OutputManager om;
    om.setup(io_comm,params,field_manager,gm,t0,t0,false);
    io_comm.barrier();

    const auto& out_fields = field_manager->get_groups_info().at("output");
    using namespace ShortFieldTagsNames;
    // Time loop
    for (Int ii=0;ii<max_steps;++ii) {
      time += dt;
      Int time_in_sec = time.seconds_from(time0);
      // Update the fields
      for (const auto& fname : out_fields->m_fields_names) {
        auto f  = field_manager->get_field(fname);
        f.sync_to_host();
        auto fl = f.get_header().get_identifier().get_layout();
        switch (fl.rank()) {
          case 1:
            {
              auto v = f.get_view<Real*,Host>();
              for (int i=0; i<fl.dim(0); ++i) {
                if (fl.has_tag(COL)) {
                  v(i) = generate_data_xy(time_in_sec,i,0);
                } else {
                  v(i) = generate_data_xy(time_in_sec,0,i);
                }
              }
            }
            break;
          case 2:
            {
              auto v = f.get_view<Real**,Host>();
              for (int i=0; i<fl.dim(0); ++i) {
                for (int j=0; j<fl.dim(1); ++j) {
                  v(i,j) = generate_data_xy(time_in_sec,i,j);
                }
              }
            }
            break;
          default:
            EKAT_ERROR_MSG ("Error! Unexpected field rank.\n");
        }
        f.sync_to_dev();
      }

      // Run the output manager for this time step
      om.run(time);
      if (io_control.is_write_step(time)) {
        output_stamps.push_back(time.to_string());
        io_control.nsamples_since_last_write = 0;
        io_control.timestamp_of_last_write = time;
      }
    }

    // Finalize the output manager (close files)
    om.finalize();
  }

  // Get a fresh new field manager
  auto field_manager = get_test_fm(grid);

  // We can use the produced output files to simultaneously check output quality and the
  // ability to read input.
  // NOTE: single-snap file start outputing at the final time only,
  //       while multi-snap starts to output at the first time step
  // First we re-construct the timestamp for the input file:
  REQUIRE(output_stamps.size()>0);
  auto input_params = get_in_params(output_type,output_freq_units, io_comm,output_stamps.front());
  const Real tol = 1000*std::numeric_limits<Real>::epsilon();

  // TODO: Create a small nc dummy file and a separate unit test which tests all input functions.
  // Test that pio_inq_dimlen is correct, using a file from one of the above parameter lists.
  {
    auto test_filename = input_params.get<std::string>("Filename");
    scorpio::register_file(test_filename,scorpio::Read);
    Int test_gcols_len = scorpio::get_dimlen_c2f(test_filename.c_str(),"ncol");
    REQUIRE(test_gcols_len==num_gcols);
    scorpio::eam_pio_closefile(test_filename);
  }

  auto f1 = field_manager->get_field("field_1");
  auto f2 = field_manager->get_field("field_2");
  auto f3 = field_manager->get_field("field_3");
  auto f4 = field_manager->get_field("field_packed");

  // Add field @ level. Note: cannot use subfield, since
  //  - we don't want to alias f3
  //  - we cannot use subfield along 2nd dim for rank-2 fields
  using namespace ShortFieldTagsNames;
  const auto units = f3.get_header().get_identifier().get_units();
  const auto& f3_lt = f3.get_header().get_identifier().get_layout();
  EKAT_REQUIRE_MSG (f3_lt.rank()>0,"WHAT?!?\n");
  auto f3_tom = Field(FID("field_3@tom",f3_lt.strip_dim(LEV),units,grid->name()));
  auto f3_bot = Field(FID("field_3@bot",f3_lt.strip_dim(LEV),units,grid->name()));
  auto f3_lev_2 = Field(FID("field_3@lev_2",f3_lt.strip_dim(LEV),units,grid->name()));
  f3_tom.allocate_view();
  f3_bot.allocate_view();
  f3_lev_2.allocate_view();
  field_manager->add_field(f3_tom);
  field_manager->add_field(f3_bot);
  field_manager->add_field(f3_lev_2);

  // Get host views
  auto f1_host = f1.get_view<Real*,Host>();
  auto f2_host = f2.get_view<Real*,Host>();
  auto f3_host = f3.get_view<Real**,Host>();
  auto f4_host = f4.get_view<Real**,Host>();
  auto f3_tom_host = f3_tom.get_view<Real*,Host>();
  auto f3_bot_host = f3_bot.get_view<Real*,Host>();
  auto f3_lev_2_host = f3_lev_2.get_view<Real*,Host>();

  // Read data
  AtmosphereInput test_input(input_params,field_manager);
  test_input.read_variables();
  test_input.finalize();
  // Check instant output
  f1.sync_to_host();
  f2.sync_to_host();
  f3.sync_to_host();
  f4.sync_to_host();
  f3_tom.sync_to_host();
  f3_bot.sync_to_host();
  f3_lev_2.sync_to_host();

  if (output_type == "instant") {
    // The diagnostic is not present in the field manager.  So we can't use the scorpio_input class
    // to read in the data.  Here we use raw IO routines to gather the data for testing.
    auto f_diag_ins_h = get_diagnostic_input(io_comm, gm, 0, input_params.get<std::string>("Filename"));
    for (int ii=0;ii<num_lcols;++ii) {
      for (int jj=0;jj<num_levs;++jj) {
        REQUIRE(std::abs(f_diag_ins_h(ii,jj)-(3.0*f4_host(ii,jj)+2.0))<tol);
      }
    }
  }
  Int current_t = max_steps*dt;
  // Check timestamp
  {
    auto test_filename = input_params.get<std::string>("Filename");
    scorpio::register_file(test_filename,scorpio::Read);
    Real time_val = scorpio::read_curr_time_c2f(test_filename.c_str());
    scorpio::eam_pio_closefile(test_filename);
    Real time_in_days = current_t/86400.; // current_t is in seconds, need to convert to days.
    REQUIRE(time_val==time_in_days);
  }
  // Check values
  for (int ii=0;ii<num_lcols;++ii) {
    REQUIRE(std::abs(f1_host(ii)-check_data_xy(current_t,dt,ii,0,output_type))<tol);
    for (int jj=0;jj<num_levs;++jj) {
      REQUIRE(std::abs(f3_host(ii,jj)-check_data_xy(current_t,dt,ii,jj,output_type))<tol);
      REQUIRE(std::abs(f4_host(ii,jj)-check_data_xy(current_t,dt,ii,jj,output_type))<tol);
    }

    REQUIRE (f3_tom_host(ii)==f3_host(ii,0));
    REQUIRE (f3_bot_host(ii)==f3_host(ii,num_levs-1));
    REQUIRE (f3_lev_2_host(ii)==f3_host(ii,2));
  }
  for (int jj=0;jj<num_levs;++jj) {
    REQUIRE(std::abs(f2_host(jj)-check_data_xy(current_t,dt,0,jj,output_type))<tol);
  }
  // All Done 
  scorpio::eam_pio_finalize();
} // end function run()
/*===================================================================================================*/
std::shared_ptr<FieldManager> get_test_fm(std::shared_ptr<const AbstractGrid> grid)
{
  using namespace ShortFieldTagsNames;
  using FL = FieldLayout;
  using FR = FieldRequest;

  // Create a fm
  auto fm = std::make_shared<FieldManager>(grid);

  const int num_lcols = grid->get_num_local_dofs();
  const int num_levs = grid->get_num_vertical_levels();

  // Create some fields for this fm
  std::vector<FieldTag> tag_h  = {COL};
  std::vector<FieldTag> tag_v  = {LEV};
  std::vector<FieldTag> tag_2d = {COL,LEV};

  std::vector<Int>     dims_h  = {num_lcols};
  std::vector<Int>     dims_v  = {num_levs};
  std::vector<Int>     dims_2d = {num_lcols,num_levs};

  const std::string& gn = grid->name();

  FieldIdentifier fid1("field_1",FL{tag_h,dims_h},m,gn);
  FieldIdentifier fid2("field_2",FL{tag_v,dims_v},kg,gn);
  FieldIdentifier fid3("field_3",FL{tag_2d,dims_2d},kg/m,gn);
  FieldIdentifier fid4("field_packed",FL{tag_2d,dims_2d},kg/m,gn);

  // Register fields with fm
  // Make sure packsize isn't bigger than the packsize for this machine, but not so big that we end up with only 1 pack.
  fm->registration_begins();
  fm->register_field(FR{fid1,"output"});
  fm->register_field(FR{fid2,"output"});
  fm->register_field(FR{fid3,"output"});
  fm->register_field(FR{fid4,"output",Pack::n}); // Register field as packed
  fm->registration_ends();

  // Make sure that field 4 is in fact a packed field
  auto field4 = fm->get_field(fid4);
  if (SCREAM_SMALL_PACK_SIZE > 1) {
    auto fid4_padding = field4.get_header().get_alloc_properties().get_padding();
    REQUIRE(fid4_padding > 0);
  }

  // Initialize these fields
  auto f1 = fm->get_field(fid1);
  auto f2 = fm->get_field(fid2);
  auto f3 = fm->get_field(fid3);
  auto f4 = fm->get_field(fid4);
  auto f1_host = f1.get_view<Real*,Host>();
  auto f2_host = f2.get_view<Real*,Host>();
  auto f3_host = f3.get_view<Real**,Host>();
  auto f4_host = f4.get_view<Pack**,Host>();

  for (int ii=0;ii<num_lcols;++ii) {
    f1_host(ii) = generate_data_xy(0,ii,0);
    for (int jj=0;jj<num_levs;++jj) {
      f2_host(jj) = generate_data_xy(0,0,jj);
      f3_host(ii,jj) = generate_data_xy(0,ii,jj);
      int ipack = jj / packsize;
      int ivec  = jj % packsize;
      f4_host(ii,ipack)[ivec] = generate_data_xy(0,ii,jj);
    }
  }
  // Update timestamp
  util::TimeStamp time ({2000,1,1},{0,0,0});
  fm->init_fields_time_stamp(time);
  // Sync back to device
  f1.sync_to_dev();
  f2.sync_to_dev();
  f3.sync_to_dev();
  f4.sync_to_dev();

  return fm;
}
/*==========================================================================================================*/
  Real generate_data_xy(const Int time, const Int i, const Int j) {
    return i + (j+1)*10.0 + time;
  }

  Real check_data_xy(const Int time, const Int dt, const Int i, const Int j, const std::string& avg_type)
  {
    Real avg_val;
    if (avg_type=="instant") {
      avg_val = generate_data_xy(time,i,j);
    } else if (avg_type=="average") {
      Int curr_time = dt;
      Int N = 0;
      avg_val = 0;
      while (curr_time<=time) {
        avg_val += generate_data_xy(curr_time,i,j);
        curr_time += dt;
        N += 1;
      }
      avg_val /= N;
    } else if (avg_type=="min") {
      avg_val = generate_data_xy(dt,i,j);
      Int curr_time = 2*dt;
      while (curr_time<=time) {
        Real tmp_data = generate_data_xy(curr_time,i,j);
        avg_val = std::min(avg_val,tmp_data);
        curr_time += dt;
      }
    } else if (avg_type=="max") {
      avg_val = generate_data_xy(0,i,j);
      Int curr_time = 2*dt;
      while (curr_time<=time) {
        Real tmp_data = generate_data_xy(curr_time,i,j);
        avg_val = std::max(avg_val,tmp_data);
        curr_time += dt;
      }
    } else {
      EKAT_REQUIRE_MSG(false,"Error! Incorrect type for check_data in io_test: " + avg_type);
    }
    return avg_val;
  }
/*==========================================================================================================*/
std::shared_ptr<GridsManager> get_test_gm(const ekat::Comm& io_comm, const Int num_gcols, const Int num_levs)
{
  ekat::ParameterList gm_params;
  gm_params.set("number_of_global_columns",num_gcols);
  gm_params.set("number_of_vertical_levels",num_levs);
  auto gm = create_mesh_free_grids_manager(io_comm,gm_params);
  gm->build_grids();
  return gm;
}
/*==================================================================================================*/
ekat::ParameterList get_in_params(const std::string& type,
                                  const std::string& freq,
                                  const ekat::Comm& comm,
                                  const std::string& t_first_write)
{
  using vos_type = std::vector<std::string>;
  ekat::ParameterList in_params("Input Parameters");
  bool multisnap = type=="multisnap";

  std::string filename =
        "io_" + std::string(multisnap ? "multisnap_test" : "output_test")
      + "." + (multisnap ? ekat::upper_case("Instant") : ekat::upper_case(type))
      + "." + freq + "_x1" + std::string(multisnap ? "" : "0")
      + ".np" + std::to_string(comm.size())
      + "." + t_first_write + ".nc";

  in_params.set<std::string>("Filename",filename);
  in_params.set<vos_type>("Field Names",
      {"field_1", "field_2", "field_3", "field_packed", "field_3@tom", "field_3@bot", "field_3@lev_2"});
  in_params.set<std::string>("Floating Point Precision","real");

  return in_params;
}
/*========================================================================================================*/
view_2d<Real>::HostMirror
get_diagnostic_input(const ekat::Comm& comm, const std::shared_ptr<GridsManager>& gm,
                     const int time_index, const std::string& filename)
{
  using namespace ShortFieldTagsNames;
  using view_1d = typename KT::template view_1d<Real>;

  auto grid = gm->get_grid("Point Grid");
  int ncols = grid->get_num_local_dofs();
  int nlevs = grid->get_num_vertical_levels();

  view_2d<Real> f_diag("DummyDiag",ncols,nlevs);
  auto f_diag_h = Kokkos::create_mirror_view(f_diag);

  std::vector<std::string> fnames = {"DummyDiagnostic"};
  std::map<std::string,view_1d::HostMirror> host_views;
  std::map<std::string,FieldLayout>  layouts;
  host_views["DummyDiagnostic"] = view_1d::HostMirror(f_diag_h.data(),f_diag_h.size());
  layouts.emplace("DummyDiagnostic",FieldLayout( {COL,LEV}, {ncols,nlevs} ) );

  ekat::ParameterList in_params;
  in_params.set("Field Names",fnames);
  in_params.set("Filename",filename);
  in_params.set<std::string>("Floating Point Precision","real");
  AtmosphereInput input(comm,in_params);
  input.init(grid,host_views,layouts);
  input.read_variables(time_index);
  input.finalize();

  return f_diag_h;
}
/*========================================================================================================*/
int get_current_t(const int tt, const int dt, const int freq, const std::string& frequency_units) {
      if (frequency_units == "nsteps") {
        // Just use dt 
        return tt*dt;
      // We will need to use timestamp information
      } else if (frequency_units == "nsecs") {
        return tt*freq;
      } else if (frequency_units == "nmins") {
        return tt*freq*60;
      } else if (frequency_units == "nhours") {
        return tt*freq*3600;
      } else if (frequency_units == "ndays") {
        return tt*freq*3600*24;
      }
  EKAT_REQUIRE_MSG(false,"Error: unknown frequency unit passed to get_current_t");
}
/*========================================================================================================*/
TEST_CASE("input_output_basic","io")
{
  ekat::Comm comm (MPI_COMM_WORLD);

  // Register all potential diagnostics
  register_diagnostics();

  std::vector<std::string> output_types =
    { "instant","average", "max", "min"};
  std::vector<std::string> output_freq =
    { "nsteps", "nsecs", "nmins", "nhours", "ndays"};
  for (const auto& of : output_freq) {
    if (comm.am_i_root()) {
      printf("Testing output freq %s\n",of.c_str());
    }
    for (const auto& ot : output_types) {
      if (comm.am_i_root()) {
        printf("  Testing output type %s...",ot.c_str());
      }
      run(ot,of);
      if (comm.am_i_root()) {
        printf("Done!\n");
      }
    } 
    if (comm.am_i_root()) {
      printf("  Testing output type multisnap...");
    }
    run_multisnap(of);
    if (comm.am_i_root()) {
      printf("Done!\n");
    }
  }
}
} // anonymous namespace
