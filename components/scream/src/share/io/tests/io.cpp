#include <catch2/catch.hpp>
#include <memory>

#include "ekat/ekat_parse_yaml_file.hpp"
#include "ekat/util/ekat_string_utils.hpp"
#include "share/io/scream_output_manager.hpp"
#include "share/io/scorpio_output.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/io/scream_scorpio_interface.hpp"

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

namespace {

using namespace scream;
using namespace ekat::units;
// Make sure packsize isn't bigger than the packsize for this machine, but not so big that we end up with only 1 pack.
const int packsize = SCREAM_SMALL_PACK_SIZE;
using Pack         = ekat::Pack<Real,packsize>;

using KT = KokkosTypes<DefaultDevice>;
template <typename S>
using view_2d = typename KT::template view_2d<S>;

std::shared_ptr<FieldManager>
get_test_fm(std::shared_ptr<const AbstractGrid> grid);

std::shared_ptr<GridsManager>
get_test_gm(const ekat::Comm& io_comm, const Int num_gcols, const Int num_levs);

ekat::ParameterList get_in_params(const std::string& type,
                                  const ekat::Comm& comm,
                                  const util::TimeStamp& t0,
                                  const int dt, const int max_steps);

view_2d<Real>::HostMirror
get_diagnostic_input(const ekat::Comm& comm, const std::shared_ptr<GridsManager>& gm,
                     const int time_index, const std::string& filename);

class DiagTest : public AtmosphereDiagnostic
{
public:
  DiagTest (const ekat::Comm& comm, const ekat::ParameterList&params)
    : AtmosphereDiagnostic(comm,params)
  {
    //Do nothing
  }

  std::string name() const { return "DiagnosticTest"; }

  // Get the required grid for the diagnostic
  std::set<std::string> get_required_grids () const {
    static std::set<std::string> s;
    s.insert(m_params.get<std::string>("Grid"));
    return s;
  }

  void set_grids (const std::shared_ptr<const GridsManager> gm) {
    using namespace ekat::units;
    using namespace ShortFieldTagsNames;
    using FL = FieldLayout;

    const auto& grid_name = m_params.get<std::string>("Grid");
    const auto grid = gm->get_grid(grid_name);
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

TEST_CASE("input_output_basic","io")
{

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

  std::vector<std::string> output_types =
    { "instant","average", "max", "min", "multisnap" };

  const Int max_steps = 10;
  const Int dt = 1;
  for (const auto& type : output_types) {
    util::TimeStamp time = t0;

    // Re-create the fm anew, so the fields are re-inited for each output type
    auto field_manager = get_test_fm(grid);
    field_manager->init_fields_time_stamp(t0);
    ekat::ParameterList params;
    ekat::parse_yaml_file("io_test_" + type + ".yaml",params);
    OutputManager om;
    om.setup(io_comm,params,field_manager,gm,t0,t0,false);
    io_comm.barrier();

    const auto& out_fields = field_manager->get_groups_info().at("output");
    // Time loop
    for (Int ii=0;ii<max_steps;++ii) {

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
                v(i) += dt;
              }
            }
            break;
          case 2:
            {
              auto v = f.get_view<Real**,Host>();
              for (int i=0; i<fl.dim(0); ++i) {
                for (int j=0; j<fl.dim(1); ++j) {
                  v(i,j) += dt;
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
      time += dt;
      om.run(time);
    }

    // Finalize the output manager (close files)
    om.finalize();
  }

  // Get a fresh new field manager
  auto field_manager = get_test_fm(grid);
  const auto& out_fields = field_manager->get_groups_info().at("output");
  auto reset_fields = [&] () {
    // Because the next test for input in this sequence relies on the same fields.  We set all field values
    // to nan to ensure no false-positive tests if a field is simply not read in as input.
    for (const auto& fname : out_fields->m_fields_names) {
      auto f = field_manager->get_field(fname);
      f.deep_copy(ekat::ScalarTraits<Real>::invalid());
    }
  };
  reset_fields();

  // At this point we should have 5 files output:
  // 1 file each for averaged, instantaneous, min and max data.
  // And 1 file with multiple time snaps of instantaneous data.
  // Cycle through each output and make sure it is correct.
  // We can use the produced output files to simultaneously check output quality and the
  // ability to read input.
  // NOTE: single-snap file start outputing at the final time only,
  //       while multi-snap starts to output at the first time step
  auto ins_params = get_in_params("instant",io_comm,t0,dt,max_steps);
  auto avg_params = get_in_params("average",io_comm,t0,dt,max_steps);
  auto min_params = get_in_params("min",io_comm,t0,dt,max_steps);
  auto max_params = get_in_params("max",io_comm,t0,dt,max_steps);
  auto multi_params = get_in_params("multisnap",io_comm,t0,dt,max_steps);
  const Real tol = 100*std::numeric_limits<Real>::epsilon();
  // TODO: Create a small nc dummy file and a separate unit test which tests all input functions.
  // Test that pio_inq_dimlen is correct, using a file from one of the above parameter lists.
  {
    auto test_filename = ins_params.get<std::string>("Filename");
    scorpio::register_file(test_filename,scorpio::Read);
    Int test_gcols_len = scorpio::get_dimlen_c2f(test_filename.c_str(),"ncol");
    REQUIRE(test_gcols_len==num_gcols);
    scorpio::eam_pio_closefile(test_filename);
  }

  auto f1 = field_manager->get_field("field_1");
  auto f2 = field_manager->get_field("field_2");
  auto f3 = field_manager->get_field("field_3");
  auto f4 = field_manager->get_field("field_packed");
  auto f1_host = f1.get_view<Real*,Host>();
  auto f2_host = f2.get_view<Real*,Host>();
  auto f3_host = f3.get_view<Real**,Host>();
  auto f4_host = f4.get_view<Real**,Host>();


  // Check instant output
  AtmosphereInput ins_input(ins_params,field_manager);
  ins_input.read_variables();
  f1.sync_to_host();
  f2.sync_to_host();
  f3.sync_to_host();
  f4.sync_to_host();

  // The diagnostic is not present in the field manager.  So we can't use the scorpio_input class
  // to read in the data.  Here we use raw IO routines to gather the data for testing.
  auto f_diag_ins_h = get_diagnostic_input(io_comm, gm, 1, ins_params.get<std::string>("Filename"));

  for (int ii=0;ii<num_lcols;++ii) {
    REQUIRE(std::abs(f1_host(ii)-(max_steps*dt+ii))<tol);
    for (int jj=0;jj<num_levs;++jj) {
      REQUIRE(std::abs(f3_host(ii,jj)-(ii+max_steps*dt + (jj+1)/10.))<tol);
      REQUIRE(std::abs(f4_host(ii,jj)-(ii+max_steps*dt + (jj+1)/10.))<tol);
      REQUIRE(std::abs(f_diag_ins_h(ii,jj)-(3.0*f4_host(ii,jj)+2.0))<tol);
    }
  }
  for (int jj=0;jj<num_levs;++jj) {
    REQUIRE(std::abs(f2_host(jj)-(max_steps*dt + (jj+1)/10.))<tol);
  }
  ins_input.finalize();
  reset_fields();

  // Check average output
  AtmosphereInput avg_input(avg_params,field_manager);
  avg_input.read_variables();
  f1.sync_to_host();
  f2.sync_to_host();
  f3.sync_to_host();
  Real avg_val;
  for (int ii=0;ii<num_lcols;++ii) {
    avg_val = (max_steps+1)/2.0*dt + ii; // Sum(x0+i*dt,i=1...N) = N*x0 + dt*N*(N+1)/2, AVG = Sum/N, note x0=ii in this case
    REQUIRE(std::abs(f1_host(ii)-avg_val)<tol);
    for (int jj=0;jj<num_levs;++jj) {
      avg_val = (max_steps+1)/2.0*dt + (jj+1)/10.+ii;  //note x0=(jj+1)/10+ii in this case.
      REQUIRE(std::abs(f3_host(ii,jj)-avg_val)<tol);
    }
  }
  for (int jj=0;jj<num_levs;++jj) {
    avg_val = (max_steps+1)/2.0*dt + (jj+1)/10.;  //note x0=(jj+1)/10 in this case.
    REQUIRE(std::abs(f2_host(jj)-avg_val)<tol);
  }
  avg_input.finalize();
  reset_fields();

  // Check max output
  // The max should be equivalent to the instantaneous because this function is monotonically increasing.
  AtmosphereInput max_input(max_params,field_manager);
  max_input.read_variables();
  f1.sync_to_host();
  f2.sync_to_host();
  f3.sync_to_host();
  for (int ii=0;ii<num_lcols;++ii) {
    REQUIRE(std::abs(f1_host(ii)-(max_steps*dt+ii))<tol);
    for (int jj=0;jj<num_levs;++jj) {
      REQUIRE(std::abs(f3_host(ii,jj)-(ii+max_steps*dt + (jj+1)/10.))<tol);
    }
  }
  for (int jj=0;jj<num_levs;++jj) {
    REQUIRE(std::abs(f2_host(jj)-(max_steps*dt + (jj+1)/10.))<tol);
  }
  max_input.finalize();
  // Check min output
  // The min should be equivalent to the first step because this function is monotonically increasing.
  AtmosphereInput min_input(min_params,field_manager);
  min_input.read_variables();
  f1.sync_to_host();
  f2.sync_to_host();
  f3.sync_to_host();
  for (int ii=0;ii<num_lcols;++ii) {
    REQUIRE(std::abs(f1_host(ii)-(dt+ii))<tol);
    for (int jj=0;jj<num_levs;++jj) {
      REQUIRE(std::abs(f3_host(ii,jj)-(dt+ii + (jj+1)/10.))<tol);
    }
  }
  for (int jj=0;jj<num_levs;++jj) {
    REQUIRE(std::abs(f2_host(jj)-(dt+(jj+1)/10.))<tol);
  }
  min_input.finalize();
  reset_fields();

  // Check multisnap output; note, tt starts at 1 instead of 0 to follow netcdf time dimension indexing.
  AtmosphereInput multi_input(multi_params,field_manager);
  for (int tt = 1; tt<=std::min(max_steps,10); tt++) {
    multi_input.read_variables(tt);
    f1.sync_to_host();
    f2.sync_to_host();
    f3.sync_to_host();
    f4.sync_to_host();

    for (int ii=0;ii<num_lcols;++ii) {
      REQUIRE(std::abs(f1_host(ii)-(tt*dt+ii))<tol);
      for (int jj=0;jj<num_levs;++jj) {
        REQUIRE(std::abs(f3_host(ii,jj)-(ii+tt*dt + (jj+1)/10.))<tol);
        REQUIRE(std::abs(f4_host(ii,jj)-(ii+tt*dt + (jj+1)/10.))<tol);
      }
    }
    for (int jj=0;jj<num_levs;++jj) {
      REQUIRE(std::abs(f2_host(jj)-(tt*dt + (jj+1)/10.))<tol);
    }
  }
  multi_input.finalize();

  // All Done 
  scorpio::eam_pio_finalize();
}

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
    f1_host(ii) = ii;
    for (int jj=0;jj<num_levs;++jj) {
      f2_host(jj) = (jj+1)/10.0;
      f3_host(ii,jj) = (ii) + (jj+1)/10.0;
      int ipack = jj / packsize;
      int ivec  = jj % packsize;
      f4_host(ii,ipack)[ivec] = (ii) + (jj+1)/10.0;
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
std::shared_ptr<GridsManager> get_test_gm(const ekat::Comm& io_comm, const Int num_gcols, const Int num_levs)
{
  ekat::ParameterList gm_params;
  gm_params.set("Number of Global Columns",num_gcols);
  gm_params.set("Number of Vertical Levels",num_levs);
  auto gm = create_mesh_free_grids_manager(io_comm,gm_params);
  gm->build_grids(std::set<std::string>{"Point Grid"});
  return gm;
}
/*==================================================================================================*/
ekat::ParameterList get_in_params(const std::string& type,
                                  const ekat::Comm& comm,
                                  const util::TimeStamp& t0,
                                  const int dt, const int max_steps)
{
  using vos_type = std::vector<std::string>;
  ekat::ParameterList in_params("Input Parameters");
  bool multisnap = type=="multisnap";

  auto t_first_write = t0 + (multisnap ? dt : dt*max_steps);

  std::string filename =
        "io_" + std::string(multisnap ? "multisnap_test" : "output_test")
      + "." + (multisnap ? ekat::upper_case("Instant") : ekat::upper_case(type))
      + ".Steps_x1" + std::string(multisnap ? "" : "0")
      + ".np" + std::to_string(comm.size())
      + "." + t_first_write.to_string() + ".nc";

  in_params.set<std::string>("Filename",filename);
  in_params.set<vos_type>("Field Names",{"field_1", "field_2", "field_3", "field_packed"});
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
  AtmosphereInput input(comm,in_params);
  input.init(grid,host_views,layouts);
  input.read_variables(time_index);
  input.finalize();

  return f_diag_h;
}

} // anonymous namespace
