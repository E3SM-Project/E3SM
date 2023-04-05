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

  //This test makes a netcdf file with 3 columns and 34 levels
  //with 4 fields, T_mid, p_mid, qv, and horiz_winds
  //The p_mid field is filled based on the following equation:
  //2j+1, where j is the level
  //The T_mid, qv, and horiz_winds fields are filled based on the equations:
  //(i-1)*10000+200*j+10*(dt/250.)*ii, where i is the column, j is the level,
  //and ii in the timestep
  //It then runs then nudging module with the netcdf file as input to nudge
  //And then checks that the output fields of the nudging module match
  //what should be in the netcdf file
  
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;
  using FL = FieldLayout;

  // A time stamp
  util::TimeStamp t0 ({2000,1,1},{0,0,0});

  ekat::Comm io_comm(MPI_COMM_WORLD);  // MPI communicator group used for I/O set as ekat object.
  Int num_levs = 34;

  // Initialize the pio_subsystem for this test:
  // MPI communicator group used for I/O.  
  // In our simple test we use MPI_COMM_WORLD, however a subset could be used.
  MPI_Fint fcomm = MPI_Comm_c2f(io_comm.mpi_comm());  
  // Gather the initial PIO subsystem data creater by component coupler
  scorpio::eam_init_pio_subsystem(fcomm);   
  
  // First set up a field manager and grids manager to interact 
  // with the output functions
  auto gm2 = create_gm(io_comm,3,num_levs);
  auto grid2 = gm2->get_grid("Point Grid");
  int num_lcols = grid2->get_num_local_dofs();

  IOControl io_control; // Needed for testing input.
  io_control.timestamp_of_last_write = t0;
  io_control.nsamples_since_last_write = 0;
  io_control.frequency_units         = "nsteps";
  std::vector<std::string> output_stamps; 

  const Int dt        = 250;
  const Int max_steps = 12;
  {
    util::TimeStamp time  = t0;
    util::TimeStamp time0(t0);

    // Re-create the fm anew, so the fields are re-inited for each output type
    //auto field_manager = get_test_fm(grid);

    using namespace ShortFieldTagsNames;
    using FR = FieldRequest;

    // Create a fm
    auto fm = std::make_shared<FieldManager>(grid2);

    const int num_levs = grid2->get_num_vertical_levels();

    // Create some fields for this fm
    std::vector<FieldTag> tag_h  = {COL};
    std::vector<FieldTag> tag_v  = {LEV};
    std::vector<FieldTag> tag_2d = {COL,LEV};

    std::vector<Int>     dims_h  = {num_lcols};
    std::vector<Int>     dims_v  = {num_levs};
    std::vector<Int>     dims_2d = {num_lcols,num_levs};

    const std::string& gn = grid2->name();

    FieldIdentifier fid1("p_mid",FL{tag_2d,dims_2d},Pa,gn);
    FieldIdentifier fid2("T_mid",FL{tag_2d,dims_2d},K,gn);
    FieldIdentifier fid3("qv",FL{tag_2d,dims_2d},kg/kg,gn);
    FieldIdentifier fid4("u",FL{tag_2d,dims_2d},m/s,gn);
    FieldIdentifier fid5("v",FL{tag_2d,dims_2d},m/s,gn);

    // Register fields with fm
    fm->registration_begins();
    fm->register_field(FR{fid1,"output"});
    fm->register_field(FR{fid2,"output"});
    fm->register_field(FR{fid3,"output"});
    fm->register_field(FR{fid4,"output"});
    fm->register_field(FR{fid5,"output"});
    fm->registration_ends();

    // Initialize these fields
    auto f1      = fm->get_field(fid1);
    auto f1_host = f1.get_view<Real**,Host>();

    auto f2      = fm->get_field(fid2);
    auto f2_host = f2.get_view<Real**,Host>();

    auto f3      = fm->get_field(fid3);
    auto f3_host = f3.get_view<Real**,Host>();

    auto f4      = fm->get_field(fid4);
    auto f4_host = f4.get_view<Real**,Host>();

    auto f5      = fm->get_field(fid5);
    auto f5_host = f5.get_view<Real**,Host>();
    
    for (int ii=0;ii<num_lcols;++ii) {
      for (int jj=0;jj<num_levs;++jj) {
	f1_host(ii,jj) = 2*jj+1;
	f2_host(ii,jj) = (ii-1)*10000+200*jj+10*(-1);
	f3_host(ii,jj) = (ii-1)*10000+200*jj+10*(-1);
	f4_host(ii,jj) = (ii-1)*10000+200*jj+10*(-1);
	f5_host(ii,jj) = (ii-1)*10000+200*jj+10*(-1);
      }
    }
    fm->init_fields_time_stamp(time);
    f1.sync_to_dev();
    f2.sync_to_dev();
    f3.sync_to_dev();
    f4.sync_to_dev();
    f5.sync_to_dev();

    // Set up parameter list control for output
    ekat::ParameterList params;
    params.set<std::string>("filename_prefix","io_output_test");
    params.set<std::string>("Averaging Type","Instant");
    params.set<int>("Max Snapshots Per File",15);
    std::vector<std::string> fnames = {"T_mid","p_mid","qv","u","v"};
    params.set<std::vector<std::string>>("Field Names",fnames);
    auto& params_sub = params.sublist("output_control");
    params_sub.set<std::string>("frequency_units","nsteps");
    params_sub.set<int>("Frequency",1);
    params_sub.set<bool>("MPI Ranks in Filename",true);
    io_control.frequency = params_sub.get<int>("Frequency");

    // Set up output manager.
    OutputManager om;
    om.setup(io_comm,params,fm,gm2,t0,t0,false);
    io_comm.barrier();

    const auto& out_fields = fm->get_groups_info().at("output");
    using namespace ShortFieldTagsNames;
    // Time loop
    for (Int ii=0;ii<max_steps;++ii) {
      time += dt;
      // Update the fields
      for (const auto& fname : out_fields->m_fields_names) {
        auto f  = fm->get_field(fname);
        f.sync_to_host();
        auto fl = f.get_header().get_identifier().get_layout();
        switch (fl.rank()) {
          case 2:
            {
              auto v = f.get_view<Real**,Host>();
              for (int i=0; i<fl.dim(0); ++i) {
                for (int j=0; j<fl.dim(1); ++j) {
                  if (fname != "p_mid"){
		      v(i,j) = (i-1)*10000+200*j+10*(dt/250.)*ii;
                  }
                  if (fname == "p_mid"){
                    v(i,j) = 2*j+1;
		  }
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

  // Create a grids manager
  const int ncols = 3;
  const int nlevs = 35;  
  auto gm = create_gm(io_comm,ncols,nlevs);
  auto grid = gm->get_grid("Physics");

  ekat::ParameterList params_mid;
  std::string nudging_f = "io_output_test.INSTANT.nsteps_x1."\
                          "np1.2000-01-01-00000.nc";
  params_mid.set<std::string>("Nudging_Filename",nudging_f);
  auto nudging_mid = std::make_shared<Nudging>(io_comm,params_mid);

  nudging_mid->set_grids(gm);

  std::map<std::string,Field> input_fields;
  std::map<std::string,Field> output_fields;
  for (const auto& req : nudging_mid->get_required_field_requests()) {
    Field f(req.fid);
    auto & f_ap = f.get_header().get_alloc_properties();
    f_ap.request_allocation(1);
    f.allocate_view();
    const auto name = f.name();
    f.get_header().get_tracking().update_time_stamp(t0);
    nudging_mid->set_required_field(f);
    input_fields.emplace(name,f);
    if (name != "p_mid"){
      nudging_mid->set_computed_field(f);
      output_fields.emplace(name,f);
    }
  }

  //initialize
  nudging_mid->initialize(t0,RunType::Initial);
  Field p_mid       = input_fields["p_mid"];
  Field T_mid       = input_fields["T_mid"];
  Field qv          = input_fields["qv"];
  Field u          = input_fields["u"];
  Field v          = input_fields["v"];
  Field T_mid_o = output_fields["T_mid"];
  Field qv_mid_o = output_fields["qv"];
  Field u_o = output_fields["u"];
  Field v_o = output_fields["v"];

  //fill data
  //Don't fill T,qv,u,v because they will be nudged anyways  
  auto p_mid_v_h   = p_mid.get_view<Real**, Host>();
  for (int icol=0; icol<ncols; icol++){
    for (int ilev=0; ilev<nlevs; ilev++){ 
      p_mid_v_h(icol,ilev) = 2*ilev;
    }
  }
  T_mid.sync_to_dev();
  qv.sync_to_dev();
  u.sync_to_dev();
  v.sync_to_dev();
  p_mid.sync_to_dev();

  //10 timesteps of 100 s
  for (int time_s = 1; time_s < 10; time_s++){
    T_mid_o.sync_to_dev();
    qv_mid_o.sync_to_dev();
    u_o.sync_to_dev();
    v_o.sync_to_dev();
    auto T_mid_v_h_o   = T_mid_o.get_view<Real**, Host>();
    auto qv_h_o   = qv_mid_o.get_view<Real**, Host>();
    auto u_h_o   = u_o.get_view<Real**, Host>();
    auto v_h_o   = v_o.get_view<Real**, Host>();
    nudging_mid->run(100);
    T_mid_o.sync_to_host();
    qv_mid_o.sync_to_host();
    u_o.sync_to_host();
    v_o.sync_to_host();

    for (int icol=0; icol<ncols; icol++){
      for (int ilev=0; ilev<nlevs; ilev++){ 
        const int time_index = time_s*100./250.;

	//First deal with cases where destination pressure levels are outside
	//the range of the pressure levels from the external file
        //If destination pressure is 0 than it will use pressure value of 1
	//from external file since that is the lowest value. A time interpolation
	//is performed but there is no interpolation between levels necessary
	if (ilev == 0){
	  double val_before  = 10000*(icol-1) + 10*int(time_index-1);
	  double val_after   = 10000*(icol-1) + 10*int(time_index);
	  double w_aft       = time_s*100.-time_index*250.;
	  double w_bef       = (time_index+1)*250-time_s*100.;
	  double val_tim_avg = (val_before*w_bef + val_after*w_aft) / 250.;
          REQUIRE(abs(T_mid_v_h_o(icol,ilev) - val_tim_avg)<0.001);
          REQUIRE(abs(qv_h_o(icol,ilev) - val_tim_avg)<0.001);
          REQUIRE(abs(u_h_o(icol,ilev) - val_tim_avg)<0.001);
          REQUIRE(abs(v_h_o(icol,ilev) - val_tim_avg)<0.001);
	  continue;
	}

        //If destination pressure is 68 than it will use highest pressure value
	//from external file. A time interpolation is performed but there is
	//no interpolation between levels necessary
	if (ilev == (nlevs-1)){
	  double val_before  = 10000*(icol-1) + 200*(ilev-1) + 10*int(time_index-1);
	  double val_after   = 10000*(icol-1) + 200*(ilev-1) + 10*int(time_index);
	  double w_aft       = time_s*100.-time_index*250.;
	  double w_bef       = (time_index+1)*250-time_s*100.;
	  double val_tim_avg = (val_before*w_bef + val_after*w_aft) / 250.;
          REQUIRE(abs(T_mid_v_h_o(icol,ilev) - val_tim_avg)<0.001);
          REQUIRE(abs(qv_h_o(icol,ilev) - val_tim_avg)<0.001);
          REQUIRE(abs(u_h_o(icol,ilev) - val_tim_avg)<0.001);
          REQUIRE(abs(v_h_o(icol,ilev) - val_tim_avg)<0.001);
	  continue;
	}
	
        double val_before = 10000*(icol-1) + 200*(ilev-1) + 10*int(time_index-1);
        double val_after = 10000*(icol-1) + 200*(ilev-1) + 10*int(time_index);
        double w_aft = time_s*100.-time_index*250.;
        double w_bef = (time_index+1)*250-time_s*100.;
	double val_tim_avg = (val_before*w_bef + val_after*w_aft) / 250.;

        double val_before_n = 10000*(icol-1) + 200*(ilev) + 10*int(time_index-1);
        double val_after_n = 10000*(icol-1) + 200*(ilev) + 10*int(time_index);
        double w_aft_n = time_s*100.-time_index*250.;
        double w_bef_n = (time_index+1)*250-time_s*100.;
	double val_tim_avg_next = (val_before_n*w_bef_n + val_after_n*w_aft_n) / 250.;
	double val_avg = (val_tim_avg_next + val_tim_avg)/2.;

	REQUIRE(abs(T_mid_v_h_o(icol,ilev) - val_avg)<0.001);
	REQUIRE(abs(qv_h_o(icol,ilev) - val_avg)<0.001);
	REQUIRE(abs(u_h_o(icol,ilev) - val_avg)<0.001);
	REQUIRE(abs(v_h_o(icol,ilev) - val_avg)<0.001);
      }
    }

}

nudging_mid->finalize();
  
}

