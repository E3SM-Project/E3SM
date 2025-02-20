#include <catch2/catch.hpp>
#include <memory>

#include "diagnostics/register_diagnostics.hpp"

#include "share/io/eamxx_output_manager.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/io/eamxx_scorpio_interface.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"

namespace {

using namespace scream;

constexpr int packsize = SCREAM_SMALL_PACK_SIZE;
using         Pack     = ekat::Pack<Real,packsize>;
using stratts_t = std::map<std::string,std::string>;

Real set_pressure(const Real p_top, const Real p_bot, const int nlevs, const int level);

std::shared_ptr<GridsManager>
get_test_gm(const ekat::Comm& io_comm, const Int num_gcols, const Int num_levs);

std::shared_ptr<FieldManager>
get_test_fm(std::shared_ptr<const AbstractGrid> grid, const bool midonly, const int p_ref=-1);

Real calculate_output(const Real pressure, const int col, const int cmp);

ekat::ParameterList set_output_params(const std::string& name, const std::string& remap_filename, const int p_ref, const bool vert_remap, const bool horiz_remap);
ekat::ParameterList set_input_params(const std::string& name, ekat::Comm& comm, const std::string& tstamp, const int p_ref);

bool approx(const Real a, const Real b) {
  const Real tol = std::numeric_limits<Real>::epsilon()*100000;
  if (std::abs(a-b) >= tol) {
    printf("Error::approx - difference of |%e - %e| = %e is greater than the max tolerance of %e\n",a,b,std::abs(a-b),tol);
  }
  return std::abs(a-b) < tol;
}

void print (const std::string& msg, const ekat::Comm& comm) {
  if (comm.am_i_root()) {
    printf("%s",msg.c_str());
  }
}

TEST_CASE("io_remap_test","io_remap_test")
{
  using namespace scream;

  // Setup the global structure
  ekat::Comm io_comm(MPI_COMM_WORLD);  // MPI communicator group used for I/O set as ekat object.
  print ("Starting io_remap_test ...\n",io_comm);


  print (" -> Test Setup ...\n",io_comm);
  scorpio::init_subsystem(io_comm);
  const int ncols_src = 64*io_comm.size();
  const int nlevs_src = 2*packsize + 1;
  const int dt = 10;

  // Construct a timestamp
  util::TimeStamp t0 ({2000,1,1},{0,0,0});
  // Setup for target levels.
  const int nlevs_tgt = packsize + 1;
  const Real p_top = 0.0;
  const Real p_bot = (nlevs_tgt-1)*(nlevs_src-1);

  // Create the grid and field manager for test
  // First set up a field manager and grids manager to interact with the output functions
  auto gm = get_test_gm(io_comm,ncols_src,nlevs_src);
  auto grid = gm->get_grid("Point Grid");
  const int  ncols_src_l = grid->get_num_local_dofs();
  auto field_manager = get_test_fm(grid, false);
  field_manager->init_fields_time_stamp(t0);
  print (" -> Test Setup ... done\n",io_comm);

  // Create remap data for both vertical and horizontal remapping
  // The strategy for remapping will be to map every 2 subsequent
  // columns to a single column. So we assume that `ncols_src` is
  // a multiple of 2.
  print (" -> Create remap file ... \n",io_comm);
  const int ncols_tgt_l = ncols_src_l/2;
  const int ncols_tgt = ncols_src/2;
  std::vector<int> col, row;
  std::vector<Real> S;
  const Real wgt = 0.4;
  for (int ii=0; ii<ncols_tgt_l; ii++) {
    const int src_col = 2*ii + ncols_src_l*io_comm.rank();
    row.push_back(1+ii+ncols_tgt_l*io_comm.rank());
    row.push_back(1+ii+ncols_tgt_l*io_comm.rank());
    col.push_back(1+src_col);
    col.push_back(1+src_col+1);
    S.push_back(wgt);
    S.push_back(1.0-wgt);
  }
  // For vertical remapping we will prokject onto a set of equally
  // spaced pressure levels from p_top to b_bot that is nearly half
  // the number of source columns.
  std::vector<Real> p_tgt;
  for (int ii=0; ii<nlevs_tgt; ++ii) {
    p_tgt.push_back(set_pressure(p_top, p_bot, nlevs_tgt, ii));
  }

  // Write remap data to file
  const std::string remap_filename = "remap_weights_np"+std::to_string(io_comm.size())+".nc";
  scorpio::register_file(remap_filename, scorpio::FileMode::Write);

  scorpio::define_dim(remap_filename,"n_a",ncols_src);
  scorpio::define_dim(remap_filename,"n_b",ncols_tgt);
  scorpio::define_dim(remap_filename,"n_s",ncols_src);
  scorpio::define_dim(remap_filename,"lev",nlevs_tgt);

  scorpio::define_var(remap_filename,"col",   {"n_s"},"int");
  scorpio::define_var(remap_filename,"row",   {"n_s"},"int");
  scorpio::define_var(remap_filename,"S",     {"n_s"},"real");
  scorpio::define_var(remap_filename,"p_levs",{"lev"},"real");

  scorpio::set_dim_decomp(remap_filename,"n_s",io_comm.rank()*ncols_src_l,ncols_src_l);
  scorpio::enddef(remap_filename);

  scorpio::write_var(remap_filename,"row",   row.data());
  scorpio::write_var(remap_filename,"col",   col.data());
  scorpio::write_var(remap_filename,"S",     S.data());
  scorpio::write_var(remap_filename,"p_levs",p_tgt.data());

  scorpio::release_file(remap_filename);
  print (" -> Create remap file ... done\n",io_comm);

  /*
   * Construct source data to be used for remapped output.
   * We want to test cases where some values may be masked, to accomplish
   * this we let the pressure profile at each column be a linear progression
   * between p_top to p_surf where p_surf is defined as:
   *
   *            /  p_bot                        for |x| > 4
   *   p_surf = |  0.5*(p_top+p_bot)            for |x| < 2
   *            \  p_bot - (-sign(x)*m + b)     otherwise, where m = (p_top+p_bot)/4.0 and b = p_top+p_bot
   *
   *
   *                             ----
   *                            /    \
   *                           /      \
   *                   --------        ----------
   */
  print (" -> Create pressure data ... \n",io_comm);
  print ("    -> setup x_src ... \n",io_comm);
  std::vector<Real> x_src;
  const auto& p_surf_f = field_manager->get_field("p_surf");
  const auto& p_surf   = p_surf_f.get_view<Real*,Host>();
  const Real dp_x = 16.0/(ncols_src-1);
  for (int ii=0; ii<ncols_src; ii++) {
    x_src.push_back(-8.0 + dp_x*ii);
  }

  print ("    -> p_surf ... \n",io_comm);
  Real slope = (p_top + p_bot) / 4.0;
  Real yint  = (p_top + p_bot);
  for (int ii=0; ii<ncols_src_l; ii++) {
    const int  g_dof = ii + io_comm.rank()*ncols_src_l;
    const Real loc_x = x_src[g_dof];
    if (loc_x < -4 || loc_x > 4) {
      p_surf(ii) = p_bot;
    } else if (loc_x > -2 && loc_x < 2) {
      p_surf(ii) =  0.5*(p_bot+p_top);
    } else {
      Real sign_x = loc_x < 0 ? 1 : -1;
      p_surf(ii) = p_bot -(sign_x * slope * loc_x + yint);
    }
  }
  p_surf_f.sync_to_dev();

  // With p_surf set we can set the actual data, starting with p_mid and p_int:
  print ("    -> p_mid and p_int ... \n",io_comm);
  const auto& pm_f = field_manager->get_field("p_mid");
  const auto& pi_f = field_manager->get_field("p_int");
  const auto& pm_v = pm_f.get_view<Real**,Host>();
  const auto& pi_v = pi_f.get_view<Real**,Host>();
  for (int ii=0; ii<ncols_src_l; ii++) {
    pi_v(ii,0) = set_pressure(p_top, p_surf(ii), nlevs_src+1, 0);
    for (int jj=0; jj<nlevs_src; jj++) {
      pi_v(ii,jj+1) = set_pressure(p_top, p_surf(ii), nlevs_src+1, jj+1);
      pm_v(ii,jj)   = 0.5*(pi_v(ii,jj)+pi_v(ii,jj+1));
    }
  }
  pm_f.sync_to_dev();
  pi_f.sync_to_dev();
  print (" -> Create pressure data ... done\n",io_comm);

  // We will assign the test values using an expression which is linear
  // in pressure, level and component.  This will make it easier to check
  // that the work of the remappers.
  print (" -> Create source data ... \n",io_comm);
  const auto& Yf_f = field_manager->get_field("Y_flat");
  const auto& Ym_f = field_manager->get_field("Y_mid");
  const auto& Yi_f = field_manager->get_field("Y_int");
  const auto& Vm_f = field_manager->get_field("V_mid");
  const auto& Vi_f = field_manager->get_field("V_int");

  const auto& Yf_v = Yf_f.get_view<Real*,Host>();
  const auto& Ym_v = Ym_f.get_view<Real**,Host>();
  const auto& Yi_v = Yi_f.get_view<Real**,Host>();
  const auto& Vm_v = Vm_f.get_view<Real***,Host>();
  const auto& Vi_v = Vi_f.get_view<Real***,Host>();

  for (int ii=0; ii<ncols_src_l; ii++) {
    Yf_v(ii)   = calculate_output(pm_v(ii,0),ii,0);
    Yi_v(ii,0) = calculate_output(pi_v(ii,0),ii,0);
    for (int cc=0; cc<2; cc++) {
      Vi_v(ii,cc,0) = calculate_output(pi_v(ii,0),ii,cc+1);
    }
    for (int jj=0; jj<nlevs_src; jj++) {
      Ym_v(ii,jj)   = calculate_output(pm_v(ii,jj),  ii,0);
      Yi_v(ii,jj+1) = calculate_output(pi_v(ii,jj+1),ii,0);
      for (int cc=0; cc<2; cc++) {
        Vm_v(ii,cc,jj)   = calculate_output(pm_v(ii,jj),  ii,cc+1);
        Vi_v(ii,cc,jj+1) = calculate_output(pi_v(ii,jj+1),ii,cc+1);
      }
    }
  }

  Yf_f.sync_to_dev();
  Ym_f.sync_to_dev();
  Yi_f.sync_to_dev();
  Vm_f.sync_to_dev();
  Vi_f.sync_to_dev();
  print (" -> Create source data ... done\n",io_comm);

  // Setup remapped output streams and run them
  print (" -> Create output ... \n",io_comm);
  register_diagnostics();
  OutputManager om_source, om_vert, om_horiz, om_vert_horiz;
  const int p_ref = (int)set_pressure(p_top, p_bot, nlevs_src+1,nlevs_src-1);

  print ("    -> source data ... \n",io_comm);
  auto source_remap_control = set_output_params("remap_source",remap_filename,p_ref,false,false);
  om_source.initialize(io_comm,source_remap_control,t0,false);
  om_source.setup(field_manager,gm);
  io_comm.barrier();
  om_source.init_timestep(t0,dt);
  om_source.run(t0+dt);
  om_source.finalize();
  print ("    -> source data ... done\n",io_comm);

  print ("    -> vertical remap ... \n",io_comm);
  auto vert_remap_control = set_output_params("remap_vertical",remap_filename,p_ref,true,false);
  om_vert.initialize(io_comm,vert_remap_control,t0,false);
  om_vert.setup(field_manager,gm);
  io_comm.barrier();
  om_vert.init_timestep(t0,dt);
  om_vert.run(t0+dt);
  om_vert.finalize();
  print ("    -> vertical remap ... done\n",io_comm);

  print ("    -> horizontal remap ... \n",io_comm);
  auto horiz_remap_control = set_output_params("remap_horizontal",remap_filename,p_ref,false,true);
  om_horiz.initialize(io_comm,horiz_remap_control,t0,false);
  om_horiz.setup(field_manager,gm);
  io_comm.barrier();
  om_horiz.init_timestep(t0,dt);
  om_horiz.run(t0+dt);
  om_horiz.finalize();
  print ("    -> horizontal remap ... done\n",io_comm);

  print ("    -> vertical-horizontal remap ... \n",io_comm);
  auto vert_horiz_remap_control = set_output_params("remap_vertical_horizontal",remap_filename,p_ref,true,true);
  om_vert_horiz.initialize(io_comm,vert_horiz_remap_control,t0,false);
  om_vert_horiz.setup(field_manager,gm);
  io_comm.barrier();
  om_vert_horiz.init_timestep(t0,dt);
  om_vert_horiz.run(t0+dt);
  om_vert_horiz.finalize();
  print ("    -> vertical-horizontal remap ... done\n",io_comm);
  print (" -> Create output ... done\n",io_comm);


  // Confirm that remapped fields are correct.
  print (" -> Test Remapped Output ... \n",io_comm);
  std::vector<std::string> fnames = {"Y_flat","Y_mid","Y_int","V_mid","V_int"};

  // ------------------------------------------------------------------------------------------------------
  //                                    ---  Vertical Remapping ---
  {
    // Note, the vertical remapper defaults to a mask value of std numeric limits scaled by 0.1;
    const float mask_val = vert_remap_control.isParameter("Fill Value")
                         ? vert_remap_control.get<double>("Fill Value") : constants::DefaultFillValue<float>().value;
    print ("    -> vertical remap ... \n",io_comm);
    auto gm_vert   = get_test_gm(io_comm,ncols_src,nlevs_tgt);
    auto grid_vert = gm_vert->get_grid("Point Grid");
    auto fm_vert   = get_test_fm(grid_vert,true,p_ref);
    auto vert_in   = set_input_params("remap_vertical",io_comm,t0.to_string(),p_ref);
    AtmosphereInput test_input(vert_in,fm_vert);
    test_input.read_variables();

    // Check the "test" metadata, which should match the field name
    // Note: the FieldAtPressureLevel diag should get the attribute from its input field,
    //       so the valuf for "Y_int"_at_XPa should be "Y_int"
    std::string att_val;
    const auto& filename = vert_in.get<std::string>("Filename");
    for (auto& fname : fnames) {
      att_val = scorpio::get_attribute<std::string>(filename,fname,"test");
      REQUIRE (att_val==fname);
    }
    std::string f_at_lev_name = "Y_int_at_" + std::to_string(p_ref) + "Pa";
    att_val = scorpio::get_attribute<std::string>(filename,f_at_lev_name,"test");
    REQUIRE (att_val=="Y_int");

    test_input.finalize();

    // Test vertically remapped output.
    // The single flat variable, "Y_flat" should match the source value exactly.  No vertical interpolation.
    // The other 4 variables should interpolate onto the expected value given the target pressure.  Note, if
    // there isn't a source pressure pair that contains a target pressure value then the verticallyi interpolated
    // value is expected to be masked.
    //
    // NOTE: For scorpio_output.cpp the mask value for vertical remapping is std::numeric_limits<Real>::max()/10.0
    const auto& Yf_f_vert = fm_vert->get_field("Y_flat");
    const auto& Ys_f_vert = fm_vert->get_field("Y_int_at_"+std::to_string(p_ref)+"Pa");
    const auto& Ym_f_vert = fm_vert->get_field("Y_mid");
    const auto& Yi_f_vert = fm_vert->get_field("Y_int");
    const auto& Vm_f_vert = fm_vert->get_field("V_mid");
    const auto& Vi_f_vert = fm_vert->get_field("V_int");

    const auto& Yf_v_vert = Yf_f_vert.get_view<Real*,Host>();
    const auto& Ys_v_vert = Ys_f_vert.get_view<Real*,Host>();
    const auto& Ym_v_vert = Ym_f_vert.get_view<Real**,Host>();
    const auto& Yi_v_vert = Yi_f_vert.get_view<Real**,Host>();
    const auto& Vm_v_vert = Vm_f_vert.get_view<Real***,Host>();
    const auto& Vi_v_vert = Vi_f_vert.get_view<Real***,Host>();

    for (int ii=0; ii<ncols_src_l; ii++) {
      const bool ref_masked = (p_ref>pi_v(ii,nlevs_src) || p_ref<pi_v(ii,0));
      const Real test_val = ref_masked ? mask_val : calculate_output(p_ref,ii,0);
      REQUIRE(approx(Ys_v_vert(ii),test_val));

      REQUIRE(approx(Yf_v_vert(ii), Yf_v(ii)));
      for (int jj=0; jj<nlevs_tgt; jj++) {
        auto p_jj = p_tgt[jj];
        const bool mid_masked = (p_jj>pm_v(ii,nlevs_src-1) || p_jj<pm_v(ii,0));
        const bool int_masked = (p_jj>pi_v(ii,nlevs_src)   || p_jj<pi_v(ii,0));
        REQUIRE(approx(Ym_v_vert(ii,jj),(mid_masked ? mask_val : calculate_output(p_jj,ii,0))));
        REQUIRE(approx(Yi_v_vert(ii,jj),(int_masked ? mask_val : calculate_output(p_jj,ii,0))));
        for (int cc=0; cc<2; cc++) {
          REQUIRE(approx(Vm_v_vert(ii,cc,jj), (mid_masked ? mask_val : calculate_output(p_jj,ii,cc+1))));
          REQUIRE(approx(Vi_v_vert(ii,cc,jj), (int_masked ? mask_val : calculate_output(p_jj,ii,cc+1))));
        }
      }
    }
    print ("    -> vertical remap ... done\n",io_comm);
  }
  // ------------------------------------------------------------------------------------------------------
  //                                    ---  Horizontal Remapping ---
  {
    // Note, the vertical remapper defaults to a mask value of std numeric limits scaled by 0.1;
    const float mask_val = horiz_remap_control.isParameter("Fill Value")
                         ? horiz_remap_control.get<double>("Fill Value") : constants::DefaultFillValue<float>().value;
    print ("    -> horizontal remap ... \n",io_comm);
    auto gm_horiz   = get_test_gm(io_comm,ncols_tgt,nlevs_src);
    auto grid_horiz = gm_horiz->get_grid("Point Grid");
    auto fm_horiz   = get_test_fm(grid_horiz,false,p_ref);
    auto horiz_in   = set_input_params("remap_horizontal",io_comm,t0.to_string(),p_ref);
    AtmosphereInput test_input(horiz_in,fm_horiz);
    test_input.read_variables();

    // Check the "test" metadata, which should match the field name
    // Note: the FieldAtPressureLevel diag should get the attribute from its input field,
    //       so the valuf for "Y_int"_at_XPa should be "Y_int"
    std::string att_val;
    const auto& filename = horiz_in.get<std::string>("Filename");
    for (auto& fname : fnames) {
      att_val = scorpio::get_attribute<std::string>(filename,fname,"test");
      REQUIRE (att_val==fname);
    }
    std::string f_at_lev_name = "Y_int_at_" + std::to_string(p_ref) + "Pa";
    att_val = scorpio::get_attribute<std::string>(filename,f_at_lev_name,"test");
    REQUIRE (att_val=="Y_int");
    test_input.finalize();

    // Test horizontally remapped output.
    // The remap we are testing is rather simple, each pair of subsequent columns are remapped to a single
    // column using `wgt` and `1-wgt` respectively.
    //
    // Note: For horizontal remapping we added the variable Y_min_at_XPa to check that this diagnostic does
    // provide some masking, since it applies vertical remapping.
    const auto& Yf_f_horiz = fm_horiz->get_field("Y_flat");
    const auto& Ys_f_horiz = fm_horiz->get_field("Y_int_at_"+std::to_string(p_ref)+"Pa");
    const auto& Ym_f_horiz = fm_horiz->get_field("Y_mid");
    const auto& Yi_f_horiz = fm_horiz->get_field("Y_int");
    const auto& Vm_f_horiz = fm_horiz->get_field("V_mid");
    const auto& Vi_f_horiz = fm_horiz->get_field("V_int");

    const auto& Yf_v_horiz = Yf_f_horiz.get_view<Real*,Host>();
    const auto& Ys_v_horiz = Ys_f_horiz.get_view<Real*,Host>();
    const auto& Ym_v_horiz = Ym_f_horiz.get_view<Real**,Host>();
    const auto& Yi_v_horiz = Yi_f_horiz.get_view<Real**,Host>();
    const auto& Vm_v_horiz = Vm_f_horiz.get_view<Real***,Host>();
    const auto& Vi_v_horiz = Vi_f_horiz.get_view<Real***,Host>();

    for (int ii=0; ii<ncols_tgt_l; ii++) {
      const int col1 = 2*ii;
      const int col2 = 2*ii+1;
      REQUIRE(approx(Yf_v_horiz(ii),Yf_v(col1)*wgt + Yf_v(col2)*(1.0-wgt)));
      REQUIRE(approx(Yi_v_horiz(ii,0), Yi_v(col1,0)*wgt + Yi_v(col2,0)*(1.0-wgt)));
      for (int cc=0; cc<2; cc++) {
        REQUIRE(approx(Vi_v_horiz(ii,cc,0), Vi_v(col1,cc,0)*wgt + Vi_v(col2,cc,0)*(1.0-wgt)));
      }
      for (int jj=0; jj<nlevs_src; jj++) {
        REQUIRE(approx(Ym_v_horiz(ii,jj),   Ym_v(col1,jj)*wgt   + Ym_v(col2,jj)*(1.0-wgt)));
        REQUIRE(approx(Yi_v_horiz(ii,jj+1), Yi_v(col1,jj+1)*wgt + Yi_v(col2,jj+1)*(1.0-wgt)));
        for (int cc=0; cc<2; cc++) {
          REQUIRE(approx(Vm_v_horiz(ii,cc,jj),   Vm_v(col1,cc,jj)*wgt   + Vm_v(col2,cc,jj)*(1.0-wgt)));
          REQUIRE(approx(Vi_v_horiz(ii,cc,jj+1), Vi_v(col1,cc,jj+1)*wgt + Vi_v(col2,cc,jj+1)*(1.0-wgt)));
        }
      }
      // For the pressured sliced variable we expect some masking which needs to be checked.
      Real Ys_exp = 0.0;
      Real Ys_wgt = 0.0;
      bool found  = false;
      if (p_ref<=pi_v(col1,nlevs_src) && p_ref>=pi_v(col1,0)) {
        found = true;
        Ys_exp += calculate_output(p_ref,col1,0)*wgt;
        Ys_wgt += wgt;
      }
      if (p_ref<=pi_v(col2,nlevs_src) && p_ref>=pi_v(col2,0)) {
        found = true;
        Ys_exp += calculate_output(p_ref,col2,0)*(1.0-wgt);
        Ys_wgt += (1.0 - wgt);
      }
      if (found) {
        Ys_exp /= Ys_wgt;
      } else {
        Ys_exp = mask_val;
      }
      REQUIRE(approx(Ys_v_horiz(ii), Ys_exp));
    }
    print ("    -> horizontal remap ... done\n",io_comm);
  }
  // ------------------------------------------------------------------------------------------------------
  //                                ---  Vertical + Horizontal Remapping ---
  {
    const float mask_val = vert_horiz_remap_control.isParameter("Fill Value")
                         ? vert_horiz_remap_control.get<double>("Fill Value") : constants::DefaultFillValue<float>().value;
    print ("    -> vertical + horizontal remap ... \n",io_comm);
    auto gm_vh   = get_test_gm(io_comm,ncols_tgt,nlevs_tgt);
    auto grid_vh = gm_vh->get_grid("Point Grid");
    auto fm_vh   = get_test_fm(grid_vh,true,p_ref);
    auto vh_in   = set_input_params("remap_vertical_horizontal",io_comm,t0.to_string(),p_ref);
    AtmosphereInput test_input(vh_in,fm_vh);
    test_input.read_variables();

    // Check the "test" metadata, which should match the field name
    // Note: the FieldAtPressureLevel diag should get the attribute from its input field,
    //       so the valuf for "Y_int"_at_XPa should be "Y_int"
    std::string att_val;
    const auto& filename = vh_in.get<std::string>("Filename");
    for (auto& fname : fnames) {
      att_val = scorpio::get_attribute<std::string>(filename,fname,"test");
      REQUIRE (att_val==fname);
    }
    std::string f_at_lev_name = "Y_int_at_" + std::to_string(p_ref) + "Pa";
    att_val = scorpio::get_attribute<std::string>(filename,f_at_lev_name,"test");
    REQUIRE (att_val=="Y_int");
    test_input.finalize();

    // Test vertically + horizontally remapped output.
    // This test is a combination of the vertical test and horizontal test above.
    // There should be maksing in the vertical in all locations where the target pressure
    // is lower/higher than the min/max of the surface pressure, just like in the vertical test.  This should
    // also translate to more masking in the horizontal reamapping.  So we must check for potential
    // masking for all variables rather than just the Y_int_at_XPa variable for the horizontal interpolation.
    //
    // NOTE: For scorpio_output.cpp the mask value for vertical remapping is std::numeric_limits<Real>::max()/10.0
    const auto& Yf_f_vh = fm_vh->get_field("Y_flat");
    const auto& Ys_f_vh = fm_vh->get_field("Y_int_at_"+std::to_string(p_ref)+"Pa");
    const auto& Ym_f_vh = fm_vh->get_field("Y_mid");
    const auto& Yi_f_vh = fm_vh->get_field("Y_int");
    const auto& Vm_f_vh = fm_vh->get_field("V_mid");
    const auto& Vi_f_vh = fm_vh->get_field("V_int");

    const auto& Yf_v_vh = Yf_f_vh.get_view<Real*,Host>();
    const auto& Ys_v_vh = Ys_f_vh.get_view<Real*,Host>();
    const auto& Ym_v_vh = Ym_f_vh.get_view<Real**,Host>();
    const auto& Yi_v_vh = Yi_f_vh.get_view<Real**,Host>();
    const auto& Vm_v_vh = Vm_f_vh.get_view<Real***,Host>();
    const auto& Vi_v_vh = Vi_f_vh.get_view<Real***,Host>();

    for (int ii=0; ii<ncols_tgt_l; ii++) {
      const int col1 = 2*ii;
      const int col2 = 2*ii+1;
      REQUIRE(approx(Yf_v_vh(ii),Yf_v(col1)*wgt + Yf_v(col2)*(1.0-wgt)));
      for (int jj=0; jj<nlevs_tgt; jj++) {
        auto p_jj = p_tgt[jj];
        const Real mid_mask_1 = (p_jj<=pm_v(col1,nlevs_src-1) && p_jj>=pm_v(col1,0));
        const Real mid_mask_2 = (p_jj<=pm_v(col2,nlevs_src-1) && p_jj>=pm_v(col2,0));
        const Real int_mask_1 = (p_jj<=pi_v(col1,nlevs_src)   && p_jj>=pi_v(col1,0));
        const Real int_mask_2 = (p_jj<=pi_v(col2,nlevs_src)   && p_jj>=pi_v(col2,0));
        Real test_mid;
        Real test_int;
        if (mid_mask_1 + mid_mask_2 > 0.0) {
          test_mid = (mid_mask_1*calculate_output(p_jj,col1,0)*wgt + mid_mask_2*calculate_output(p_jj,col2,0)*(1-wgt))/(mid_mask_1*wgt + mid_mask_2*(1-wgt));
        } else {
          // This point is completely masked out, assign masked value
          test_mid = mask_val;
        }
        if (int_mask_1 + int_mask_2 > 0.0) {
          test_int = (int_mask_1*calculate_output(p_jj,col1,0)*wgt + int_mask_2*calculate_output(p_jj,col2,0)*(1-wgt))/(int_mask_1*wgt + int_mask_2*(1-wgt));
        } else {
          // This point is completely masked out, assign masked value
          test_int = mask_val;
        }
        REQUIRE(approx(Ym_v_vh(ii,jj), test_mid));
        REQUIRE(approx(Yi_v_vh(ii,jj), test_int));
        for (int cc=0; cc<2; cc++) {
          if (mid_mask_1 + mid_mask_2 > 0.0) {
            test_mid = (mid_mask_1*calculate_output(p_jj,col1,cc+1)*wgt + mid_mask_2*calculate_output(p_jj,col2,cc+1)*(1-wgt))/(mid_mask_1*wgt + mid_mask_2*(1-wgt));
          } else {
            // This point is completely masked out, assign masked value
            test_mid = mask_val;
          }
          if (int_mask_1 + int_mask_2 > 0.0) {
            test_int = (int_mask_1*calculate_output(p_jj,col1,cc+1)*wgt + int_mask_2*calculate_output(p_jj,col2,cc+1)*(1-wgt))/(int_mask_1*wgt + int_mask_2*(1-wgt));
          } else {
            // This point is completely masked out, assign masked value
            test_int = mask_val;
          }
          REQUIRE(approx(Vm_v_vh(ii,cc,jj), test_mid));
          REQUIRE(approx(Vi_v_vh(ii,cc,jj), test_int));
        }
      }
      // For the pressured sliced variable we expect it to match the solution from horizontal mapping only so we use the same syntax.
      Real Ys_exp = 0.0;
      Real Ys_wgt = 0.0;
      bool found  = false;
      if (p_ref<=pi_v(col1,nlevs_src) && p_ref>=pi_v(col1,0)) {
        found = true;
        Ys_exp += calculate_output(p_ref,col1,0)*wgt;
        Ys_wgt += wgt;
      }
      if (p_ref<=pi_v(col2,nlevs_src) && p_ref>=pi_v(col2,0)) {
        found = true;
        Ys_exp += calculate_output(p_ref,col2,0)*(1.0-wgt);
        Ys_wgt += (1.0 - wgt);
      }
      if (found) {
        Ys_exp /= Ys_wgt;
      } else {
        Ys_exp = mask_val;
      }
      REQUIRE(approx(Ys_v_vh(ii), Ys_exp));
    }
    print ("    -> vertical + horizontal remap ... done\n",io_comm);
  }
  // ------------------------------------------------------------------------------------------------------
  // All Done
  print (" -> Test Remapped Output ... done\n",io_comm);
  scorpio::finalize_subsystem();

}
/*==========================================================================================================*/
Real set_pressure(const Real p_top, const Real p_bot, const int nlevs, const int level)
{
  const Real dp = (p_bot - p_top) / (nlevs-1);
  return (p_top + dp*level);
}
/*==========================================================================================================*/
Real calculate_output(const Real pressure, const int col, const int cmp)
{
  const int p_slp = 100.0;
  const int x_slp = 2.0;
  const int c_slp = 7.0;
  return p_slp*pressure + x_slp*col + c_slp*cmp;
}
/*==========================================================================================================*/
std::shared_ptr<GridsManager> get_test_gm(const ekat::Comm& io_comm, const Int num_gcols, const Int num_levs)
{
  auto gm = create_mesh_free_grids_manager(io_comm,0,0,num_levs,num_gcols);
  gm->build_grids();
  return gm;
}
/*==========================================================================================================*/
std::shared_ptr<FieldManager> get_test_fm(std::shared_ptr<const AbstractGrid> grid, const bool midonly, const int p_ref)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;
  using FL = FieldLayout;
  using FR = FieldRequest;

  // Create a fm
  auto fm = std::make_shared<FieldManager>(grid);

  const int num_lcols = grid->get_num_local_dofs();
  const int num_levs = grid->get_num_vertical_levels();

  const FieldTag itag = midonly ? LEV : ILEV;
  const int      isz  = midonly ? num_levs : num_levs + 1;

  // Create some fields for this fm
  std::vector<FieldTag> tag_h  = {COL};
  std::vector<FieldTag> tag_2d_m = {COL,LEV};
  std::vector<FieldTag> tag_2d_i = {COL,itag};
  std::vector<FieldTag> tag_3d_m = {COL,CMP,LEV};
  std::vector<FieldTag> tag_3d_i = {COL,CMP,itag};

  std::vector<Int>     dims_h  = {num_lcols};
  std::vector<Int>     dims_2d_m = {num_lcols,num_levs};
  std::vector<Int>     dims_2d_i = {num_lcols,isz};
  std::vector<Int>     dims_3d_m = {num_lcols,2,num_levs};
  std::vector<Int>     dims_3d_i = {num_lcols,2,isz};

  const std::string& gn = grid->name();

  FieldIdentifier fid_ps("p_surf",FL{tag_h,dims_h},Pa,gn);
  FieldIdentifier fid_pm("p_mid", FL{tag_2d_m,dims_2d_m},Pa,gn);
  FieldIdentifier fid_pi("p_int", FL{tag_2d_i,dims_2d_i},Pa/m,gn);
  FieldIdentifier fid_Yf("Y_flat",FL{tag_h,dims_h},m,gn);
  FieldIdentifier fid_Ym("Y_mid", FL{tag_2d_m,dims_2d_m},m,gn);
  FieldIdentifier fid_Yi("Y_int", FL{tag_2d_i,dims_2d_i},m,gn);
  FieldIdentifier fid_Vm("V_mid", FL{tag_3d_m,dims_3d_m},m,gn);
  FieldIdentifier fid_Vi("V_int", FL{tag_3d_i,dims_3d_i},m,gn);

  // Register fields with fm
  // Make sure packsize isn't bigger than the packsize for this machine, but not so big that we end up with only 1 pack.
  fm->registration_begins();
  fm->register_field(FR{fid_ps,"output"});
  fm->register_field(FR{fid_pm,"output",Pack::n});
  fm->register_field(FR{fid_pi,"output",Pack::n});
  fm->register_field(FR{fid_Yf,"output"});
  fm->register_field(FR{fid_Ym,"output",Pack::n});
  fm->register_field(FR{fid_Yi,"output",Pack::n});
  fm->register_field(FR{fid_Vm,"output",Pack::n});
  fm->register_field(FR{fid_Vi,"output",Pack::n});
  if (p_ref>=0) {
    FieldIdentifier fid_di("Y_int_at_"+std::to_string(p_ref)+"Pa", FL{tag_h,dims_h},m,gn);
    fm->register_field(FR{fid_di,"output"});
  }
  fm->registration_ends();

  // Initialize these fields
  auto f_ps = fm->get_field(fid_ps);
  auto f_pm = fm->get_field(fid_pm);
  auto f_pi = fm->get_field(fid_pi);
  auto f_Yf = fm->get_field(fid_Yf);
  auto f_Ym = fm->get_field(fid_Ym);
  auto f_Yi = fm->get_field(fid_Yi);
  auto f_Vm = fm->get_field(fid_Vm);
  auto f_Vi = fm->get_field(fid_Vi);

  // Set some string to be written to file as attribute to the variables
  for (const std::string fname : {"Y_flat","Y_mid","Y_int","V_mid","V_int"}) {
    auto& f = fm->get_field(fname);
    auto& str_atts = f.get_header().get_extra_data<stratts_t>("io: string attributes");
    str_atts["test"] = fname;
  }
 // Update timestamp
  util::TimeStamp time ({2000,1,1},{0,0,0});
  fm->init_fields_time_stamp(time);

  // Sync back to device
  f_ps.sync_to_dev();
  f_pm.sync_to_dev();
  f_pi.sync_to_dev();
  f_Yf.sync_to_dev();
  f_Ym.sync_to_dev();
  f_Yi.sync_to_dev();
  f_Vm.sync_to_dev();
  f_Vi.sync_to_dev();
  if (p_ref>=0) {
    auto f_di = fm->get_field("Y_int_at_"+std::to_string(p_ref)+"Pa");
    f_di.sync_to_dev();
  }

  return fm;
}
/*==========================================================================================================*/
ekat::ParameterList set_output_params(const std::string& name, const std::string& remap_filename, const int p_ref, const bool vert_remap, const bool horiz_remap)
{
  using vos_type = std::vector<std::string>;
  ekat::ParameterList params;

  params.set<std::string>("filename_prefix",name);
  params.set<std::string>("Averaging Type","Instant");
  params.set<int>("Max Snapshots Per File",1);
  params.set<std::string>("Floating Point Precision","real");
  auto& oc = params.sublist("output_control");
  oc.set<int>("Frequency",1);
  oc.set<std::string>("frequency_units","nsteps");

  vos_type fields_out = {"Y_flat", "Y_mid", "Y_int", "V_mid", "V_int"};
  if (p_ref>=0) {
    fields_out.push_back("Y_int_at_"+std::to_string(p_ref)+"Pa");
  }
  if (!vert_remap && !horiz_remap) {
    fields_out.push_back("p_surf");
    fields_out.push_back("p_mid");
    fields_out.push_back("p_int");
  }
  params.set<vos_type>("Field Names",fields_out);

  if (vert_remap) {
    params.set<std::string>("vertical_remap_file",remap_filename); // TODO, make this work for general np=?
  }
  if (horiz_remap) {
    params.set<std::string>("horiz_remap_file",remap_filename); // TODO, make this work for general np=?
  }

  return params;
}
/*==========================================================================================================*/
ekat::ParameterList set_input_params(const std::string& name, ekat::Comm& comm, const std::string& tstamp, const int p_ref)
{
  using vos_type = std::vector<std::string>;
  ekat::ParameterList in_params("Input Parameters");
  std::string filename = name + ".INSTANT.nsteps_x1.np" + std::to_string(comm.size()) + "." + tstamp + ".nc";
  in_params.set<std::string>("Filename",filename);
  vos_type fields_in =  {"Y_flat", "Y_mid", "Y_int", "V_mid", "V_int"};
  if (p_ref>=0) {
    fields_in.push_back("Y_int_at_"+std::to_string(p_ref)+"Pa");
  }

  in_params.set<vos_type>("Field Names", fields_in);
  in_params.set<std::string>("Floating Point Precision","real");
  return in_params;
}
/*==========================================================================================================*/

} //namespace
