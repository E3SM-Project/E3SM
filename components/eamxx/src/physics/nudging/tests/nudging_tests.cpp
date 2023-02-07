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

  //using Pack = ekat::Pack<Real,SCREAM_PACK_SIZE>;
  /*
  using Pack = ekat::Pack<Real,1>;
  using KT = KokkosTypes<DefaultDevice>;
  using view_1d = typename KT::template view_1d<Pack>;
  using view_1d_const = typename KT::template view_1d<const Pack>;
  */

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
  //Int num_levs = 33;
  Int num_levs = 34;


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

  const Int dt        = 250;
  const Int max_steps = 12;
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
    std::vector<FieldTag> tag_3d = {COL,CMP,LEV};

    std::vector<Int>     dims_h  = {num_lcols};
    std::vector<Int>     dims_v  = {num_levs};
    std::vector<Int>     dims_2d = {num_lcols,num_levs};
    std::vector<Int>     dims_3d = {num_lcols,2,num_levs};

    const std::string& gn = grid2->name();

    FieldIdentifier fid2("p_mid",FL{tag_2d,dims_2d},Pa,gn);
    FieldIdentifier fid3("T_mid",FL{tag_2d,dims_2d},K,gn);
    FieldIdentifier fid4("qv",FL{tag_2d,dims_2d},kg/kg,gn);
    FieldIdentifier fid5("horiz_winds",FL{tag_3d,dims_3d},m/s,gn);

    // Register fields with fm
    // Make sure packsize isn't bigger than the packsize for this machine, but not so big   that we end up with only 1 pack.
    fm->registration_begins();
    fm->register_field(FR{fid2,"output"});
    fm->register_field(FR{fid3,"output"});
    fm->register_field(FR{fid4,"output"});
    fm->register_field(FR{fid5,"output"});
    fm->registration_ends();


    // Initialize these fields
    auto f3 = fm->get_field(fid3);
    auto f3_host = f3.get_view<Real**,Host>();

    auto f2 = fm->get_field(fid2);
    auto f2_host = f2.get_view<Real**,Host>();

    auto f4 = fm->get_field(fid4);
    auto f4_host = f4.get_view<Real**,Host>();

    auto f5 = fm->get_field(fid5);
    auto f5_host = f5.get_view<Real***,Host>();
    
    for (int ii=0;ii<num_lcols;++ii) {
      for (int jj=0;jj<num_levs;++jj) {
        f4_host(ii,jj) = 1;
	f3_host(ii,jj) = 1;
        f2_host(ii,jj) = 1;
	f5_host(ii,0,jj) = 1;
	f5_host(ii,1,jj) = 1;
      }
    }
    // Update timestamp
    //util::TimeStamp time ({2000,1,1},{0,0,0});
    fm->init_fields_time_stamp(time);
    // Sync back to device
    f2.sync_to_dev();
    f3.sync_to_dev();
    f4.sync_to_dev();
    f5.sync_to_dev();

    fm->init_fields_time_stamp(t0);

    // Set up parameter list control for output
    
    ekat::ParameterList params;
    params.set<std::string>("Casename","io_output_test");
    params.set<std::string>("Averaging Type","Instant");
    params.set<int>("Max Snapshots Per File",10);
    //auto& params_sub_f = params.sublist("Field Names");
    std::vector<std::string> fnames = {"T_mid","p_mid","qv","horiz_winds"};
    //params_sub_f.set<std::vector<std::string>>("Field Names",fnames);
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
                  if (fname != "p_mid"){
		      v(i,j) = (i-1)*10000+200*j+10*(dt/250.)*ii;
		      std::cout<<"v("<<i<<","<<j<<"): "<<v(i,j)<<std::endl;
                  }
                  if (fname == "p_mid"){
                    v(i,j) = 2*j+1;
		  }
                }
              }
            }
	    break;
	case 3:
            {
              auto v = f.get_view<Real***,Host>();
              for (int i=0; i<fl.dim(0); ++i) {
                for (int j=0; j<fl.dim(2); ++j) {
		  v(i,0,j) = (i-1)*10000+200*j+10*(dt/250.)*ii;
		  v(i,1,j) = (i-1)*10000+200*j+10*(dt/250.)*ii;
		  std::cout<<"horiz 0 v("<<i<<","<<j<<"): "<<v(i,0,j)<<std::endl;
		  std::cout<<"horiz 1 v("<<i<<","<<j<<"): "<<v(i,1,j)<<std::endl;
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
  std::cout<<"Get Here 1: "<<std::endl;
  //ekat::Comm comm(MPI_COMM_WORLD);

  // Create a grids manager
  const int ncols = 3;
  //const int nlevs = packsize*2 + 1;  // Note, we need at least 3 levels for the test to work
  //const int nlevs = 33;  // Note, we need at least 3 levels for the test to work
  const int nlevs = 35;  
  auto gm = create_gm(io_comm,ncols,nlevs);
  auto grid = gm->get_grid("Physics");

  const auto units = ekat::units::Units::invalid();

  std::cout<<"Columns outside: "<<ncols<<std::endl;
  std::cout<<"Levs outside: "<<nlevs<<std::endl;

  ekat::ParameterList params_mid;
  std::vector<std::string> fnames = {"T_mid","p_mid","qv","horiz_winds"};
  //std::vector<std::string> fnames = {"T_mid","p_mid","qv"};
  //std::vector<std::string> fnames = {"T_mid"};
  std::string nudging_f = "io_output_test.INSTANT.nsteps_x1."\
                          "np1.2000-01-01-00250.nc";
  params_mid.set<std::string>("Nudging_Filename",nudging_f);
  params_mid.set<Int>("Time_Step_File",250);
  //std::vector<std::string> fnames = {"T_mid_r"};
  params_mid.set<std::vector<std::string>>("Field_Names",fnames);
  auto nudging_mid = std::make_shared<NUDGING>(io_comm,params_mid);
  nudging_mid->set_grids(gm);
  std::cout<<"Get Here 2: "<<std::endl;

  std::map<std::string,Field> input_fields;
  std::map<std::string,Field> output_fields;
  for (const auto& req : nudging_mid->get_required_field_requests()) {
    Field f(req.fid);
    auto & f_ap = f.get_header().get_alloc_properties();
    //f_ap.request_allocation(packsize);
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
  std::cout<<"Get Here 4: "<<std::endl;
  //initialize
  nudging_mid->initialize(t0,RunType::Initial);
  
  Field p_mid       = input_fields["p_mid"];
  Field f_mid       = input_fields["T_mid"];
  Field qv          = input_fields["qv"];
  Field horiz_winds = input_fields["horiz_winds"];
  //Field p_mid_o = output_fields["p_mid"];
  Field f_mid_o = output_fields["T_mid"];
  Field qv_mid_o = output_fields["qv"];
  Field hw_mid_o = output_fields["horiz_winds"];

  //fill data
  
  auto f_mid_v_h   = f_mid.get_view<Real**, Host>();
  auto p_mid_v_h   = p_mid.get_view<Real**, Host>();
  //Don't fill Temperature because it is going to be nudged anyways  
  for (int icol=0; icol<ncols; icol++){
    for (int ilev=0; ilev<nlevs; ilev++){ 
      //f_mid_v_h(icol,ilev) = 100*(icol-1) + 2*ilev;
      //p_mid_v_h(icol,ilev) = 2*ilev+2;
      p_mid_v_h(icol,ilev) = 2*ilev;
    }
  }
  f_mid.sync_to_dev();
  horiz_winds.sync_to_dev();
  qv.sync_to_dev();
  p_mid.sync_to_dev();
  std::cout<<"Get Here 5: "<<std::endl;
  //10 timesteps of 100 s
  for (int time_s = 1; time_s < 10; time_s++){
    f_mid_o.sync_to_dev();
    qv_mid_o.sync_to_dev();
    hw_mid_o.sync_to_dev();
    p_mid.sync_to_dev();
    auto f_mid_v_h_o   = f_mid_o.get_view<Real**, Host>();
    auto qv_h_o   = qv_mid_o.get_view<Real**, Host>();
    auto hw_h_o   = hw_mid_o.get_view<Real***, Host>();
    //nudging_mid->run(time_s*100);
    nudging_mid->run(100);
    f_mid_o.sync_to_host();
    qv_mid_o.sync_to_host();
    hw_mid_o.sync_to_host();
    p_mid.sync_to_host();

    if(time_s<4){continue;}
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
          REQUIRE(abs(f_mid_v_h_o(icol,ilev) - val_tim_avg)<0.1);
          REQUIRE(abs(qv_h_o(icol,ilev) - val_tim_avg)<0.1);
          REQUIRE(abs(hw_h_o(icol,0,ilev) - val_tim_avg)<0.1);
          REQUIRE(abs(hw_h_o(icol,1,ilev) - val_tim_avg)<0.1);
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
          REQUIRE(abs(f_mid_v_h_o(icol,ilev) - val_tim_avg)<0.1);
          REQUIRE(abs(qv_h_o(icol,ilev) - val_tim_avg)<0.1);
          REQUIRE(abs(hw_h_o(icol,0,ilev) - val_tim_avg)<0.1);
          REQUIRE(abs(hw_h_o(icol,1,ilev) - val_tim_avg)<0.1);
	  continue;
	}

	
        double val_before = 10000*(icol-1) + 200*(ilev-1) + 10*int(time_index-1);
        double val_after = 10000*(icol-1) + 200*(ilev-1) + 10*int(time_index);
        double w_aft = time_s*100.-time_index*250.;
        double w_bef = (time_index+1)*250-time_s*100.;
	//std::cout<<"val_before: "<<val_before<<std::endl;
        //std::cout<<"val_after: "<<val_after<<std::endl;
	//std::cout<<"w_before: "<<w_bef<<std::endl;
        //std::cout<<"w_after: "<<w_aft<<std::endl;
	double val_tim_avg = (val_before*w_bef + val_after*w_aft) / 250.;

        double val_before_n = 10000*(icol-1) + 200*(ilev) + 10*int(time_index-1);
        double val_after_n = 10000*(icol-1) + 200*(ilev) + 10*int(time_index);
        double w_aft_n = time_s*100.-time_index*250.;
        double w_bef_n = (time_index+1)*250-time_s*100.;
	double val_tim_avg_next = (val_before_n*w_bef_n + val_after_n*w_aft_n) / 250.;
        double val_avg = (val_tim_avg_next + val_tim_avg)/2.;

	std::cout<<"f_mid_v_h_o("<<icol<<","<<ilev<<"): "<<f_mid_v_h_o(icol,ilev)<<std::endl;
	std::cout<<"qv_h_o("<<icol<<","<<ilev<<"): "<<qv_h_o(icol,ilev)<<std::endl;
        std::cout<<"hw_h_o0("<<icol<<","<<ilev<<"): "<<hw_h_o(icol,0,ilev)<<std::endl;
	std::cout<<"hw_h_o1("<<icol<<","<<ilev<<"): "<<hw_h_o(icol,1,ilev)<<std::endl;
	//REQUIRE(f_mid_v_h_o(icol,ilev) == val_avg);
	REQUIRE(abs(f_mid_v_h_o(icol,ilev) - val_avg)<0.1);
	REQUIRE(abs(qv_h_o(icol,ilev) - val_avg)<0.1);
	REQUIRE(abs(hw_h_o(icol,0,ilev) - val_avg)<0.1);
	REQUIRE(abs(hw_h_o(icol,1,ilev) - val_avg)<0.1);
      }
    }

}

nudging_mid->finalize();
  
}


//} // empty namespace
