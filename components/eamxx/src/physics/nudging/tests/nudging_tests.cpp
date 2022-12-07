#include "catch2/catch.hpp"
#include "physics/nudging/atmosphere_nudging.hpp"
#include "share/grid/mesh_free_grids_manager.hpp"

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



using namespace scream;

//namespace {

std::shared_ptr<GridsManager>
create_gm (const ekat::Comm& comm, const int ncols, const int nlevs) {

  const int num_global_cols = ncols*comm.size();

  ekat::ParameterList gm_params;
  gm_params.set<int>("number_of_global_columns", num_global_cols);
  gm_params.set<int>("number_of_vertical_levels", nlevs);

  auto gm = create_mesh_free_grids_manager(comm,gm_params);
  gm->build_grids();

  return gm;
}


TEST_CASE("nudging") {
  std::cout<<"Get into nudging test"<<std::endl;

  using Pack = ekat::Pack<Real,SCREAM_PACK_SIZE>;
  using KT = KokkosTypes<DefaultDevice>;
  using view_1d = typename KT::template view_1d<Pack>;
  using view_1d_const = typename KT::template view_1d<const Pack>;

  using namespace ekat::units;
  using namespace ShortFieldTagsNames;
  using FL = FieldLayout;

  constexpr int packsize = SCREAM_PACK_SIZE;
  std::cout<<"packsize: "<<packsize<<std::endl;


  // A time stamp
  //util::TimeStamp t0 ({2022,1,1},{0,0,0});
  util::TimeStamp t0 ({2000,1,1},{0,0,0});
  
  //Fill data to interpolate   
  //Create .nc file with just T_mid
  //Then read that in probably in nudging code

  //const std::string output_type = "output";
  ekat::Comm io_comm(MPI_COMM_WORLD);  // MPI communicator group used for I/O set as ekat object.
  //Int num_gcols = 2*io_comm.size();
  Int num_levs = 33;


  // Initialize the pio_subsystem for this test:
  MPI_Fint fcomm = MPI_Comm_c2f(io_comm.mpi_comm());  // MPI communicator group used for I/O.  In our simple test we use MPI_COMM_WORLD, however a subset could be used.
  scorpio::eam_init_pio_subsystem(fcomm);   // Gather the initial PIO subsystem data creater by component coupler

  
  
  // First set up a field manager and grids manager to interact with the output functions
  //auto gm2 = create_gm(io_comm,ncols,num_levs);
  auto gm2 = create_gm(io_comm,3,num_levs);
  auto grid2 = gm2->get_grid("Point Grid");
  int num_lcols = grid2->get_num_local_dofs();

  // Construct a timestamp
  //util::TimeStamp t0 ({2000,1,1},{0,0,0});
  IOControl io_control; // Needed for testing input.
  io_control.timestamp_of_last_write = t0;
  io_control.nsamples_since_last_write = 0;
  io_control.frequency_units         = "nsteps";
  std::vector<std::string> output_stamps; 

  const Int dt        = 1;
  const Int max_steps = 2;
  {
    util::TimeStamp time  = t0;
    util::TimeStamp time0(t0);

    // Re-create the fm anew, so the fields are re-inited for each output type
    //auto field_manager = get_test_fm(grid);

    using namespace ShortFieldTagsNames;
    using FL = FieldLayout;
    using FR = FieldRequest;

    // Create a fm
    auto fm = std::make_shared<FieldManager>(grid2);

    const int num_lcols = grid2->get_num_local_dofs();
    const int num_levs = grid2->get_num_vertical_levels();

    // Create some fields for this fm
    std::vector<FieldTag> tag_h  = {COL};
    std::vector<FieldTag> tag_v  = {LEV};
    std::vector<FieldTag> tag_2d = {COL,LEV};

    std::vector<Int>     dims_h  = {num_lcols};
    std::vector<Int>     dims_v  = {num_levs};
    std::vector<Int>     dims_2d = {num_lcols,num_levs};

    const std::string& gn = grid2->name();

    FieldIdentifier fid3("T_mid_r",FL{tag_2d,dims_2d},K,gn);

    // Register fields with fm
    // Make sure packsize isn't bigger than the packsize for this machine, but not so big   that we end up with only 1 pack.
    fm->registration_begins();
    fm->register_field(FR{fid3,"output"});
    fm->registration_ends();


    // Initialize these fields
    auto f3 = fm->get_field(fid3);
    auto f3_host = f3.get_view<Real**,Host>();

    for (int ii=0;ii<num_lcols;++ii) {
      for (int jj=0;jj<num_levs;++jj) {
        f3_host(ii,jj) = 1;
      }
    }
    // Update timestamp
    //util::TimeStamp time ({2000,1,1},{0,0,0});
    fm->init_fields_time_stamp(time);
    // Sync back to device
    f3.sync_to_dev();

    fm->init_fields_time_stamp(t0);

    // Set up parameter list control for output
    
    ekat::ParameterList params;
    //ekat::parse_yaml_file("io_test_" + output_type + ".yaml",params);
    //params.set<std::string>("Floating Point Precision","real");
    params.set<std::string>("Casename","io_output_test");
    params.set<std::string>("Averaging Type","Instant");
    params.set<int>("Max Snapshots Per File",1);
    //auto& params_sub_f = params.sublist("Field Names");
    std::vector<std::string> fnames = {"T_mid_r"};
    //params_sub_f.set<std::vector<std::string>>("Field Names",fnames);
    params.set<std::vector<std::string>>("Field Names",fnames);
    //params_sub.set<std::string>("frequency_units","T_mid_r");
    auto& params_sub = params.sublist("output_control");
    //params_sub.set<std::string>("frequency_units",output_freq_units);
    params_sub.set<std::string>("frequency_units","nsteps");
    params_sub.set<int>("Frequency",1);
    params_sub.set<bool>("MPI Ranks in Filename",true);
    io_control.frequency = params_sub.get<int>("Frequency");
    //io_control.frequency = 1;
    

    // Set up output manager.
    OutputManager om;
    om.setup(io_comm,params,fm,gm2,t0,t0,false);
    io_comm.barrier();

    const auto& out_fields = fm->get_groups_info().at("output");
    using namespace ShortFieldTagsNames;
    // Time loop
    for (Int ii=0;ii<max_steps;++ii) {
      time += dt;
      Int time_in_sec = time.seconds_from(time0);
      std::cout<<"time_in_sec: "<<time_in_sec<<std::endl;
      // Update the fields
      for (const auto& fname : out_fields->m_fields_names) {
        auto f  = fm->get_field(fname);
	std::cout<<"fname: "<<fname<<std::endl;
        f.sync_to_host();
        auto fl = f.get_header().get_identifier().get_layout();
        switch (fl.rank()) {
          case 2:
            {
              auto v = f.get_view<Real**,Host>();
              for (int i=0; i<fl.dim(0); ++i) {
                for (int j=0; j<fl.dim(1); ++j) {
                  v(i,j) = (-j)*100 + i + 1;
		  std::cout<<"v(i,j): "<<v(i,j)<<std::endl;
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
	std::cout<<"Get in io_control.is_write_step(time)"<<std::endl;
        output_stamps.push_back(time.to_string());
        io_control.nsamples_since_last_write = 0;
        io_control.timestamp_of_last_write = time;
      }
    }

    // Finalize the output manager (close files)
    om.finalize();
  }
  

      
  //ekat::Comm comm(MPI_COMM_WORLD);

  // Create a grids manager
  const int ncols = 3;
  const int nlevs = packsize*2 + 1;  // Note, we need at least 3 levels for the test to work
  auto gm = create_gm(io_comm,ncols,nlevs);
  auto grid = gm->get_grid("Physics");

  const auto units = ekat::units::Units::invalid();

  std::cout<<"Columns outside: "<<ncols<<std::endl;
  std::cout<<"Levs outside: "<<nlevs<<std::endl;

  FieldIdentifier fid_mid ("T_mid",FL({COL,LEV},{ncols,nlevs}),units,grid->name());
  Field f_mid (fid_mid);
  f_mid.get_header().get_alloc_properties().request_allocation(packsize);
  f_mid.allocate_view();
  f_mid.get_header().get_tracking().update_time_stamp(t0);

  ekat::ParameterList params_mid;
  auto nudging_mid = std::make_shared<NUDGING>(io_comm,params_mid);
  nudging_mid->set_grids(gm);
  nudging_mid->set_required_field(f_mid);

  //fill data
  auto f_mid_v_h   = f_mid.get_view<Real**, Host>();
  for (int ilev=0; ilev<nlevs; ilev++){ 
    for (int icol=0; icol<ncols; icol++){
      f_mid_v_h(icol,ilev) = (-ilev)*100 + icol;
    }
  }
  f_mid.sync_to_dev();

  //initialize
  nudging_mid->initialize(t0,RunType::Initial);

  //run
  nudging_mid->run(100);
  f_mid.sync_to_host();

  //Now check that I was able to nudge it by a value of 1
  for (int ilev=0; ilev<nlevs; ilev++){ 
    for (int icol=0; icol<ncols; icol++){
      std::cout<<"f_mid_v_h(icol,ilev): "<<f_mid_v_h(icol,ilev)<<std::endl;
      REQUIRE(f_mid_v_h(icol,ilev) == (-ilev)*100 + icol + 1);
    }
  }
 
  
}


//} // empty namespace
