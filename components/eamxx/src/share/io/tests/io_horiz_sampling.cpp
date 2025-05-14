#include <catch2/catch.hpp>
#include <memory>

#include "diagnostics/register_diagnostics.hpp"

#include "share/io/eamxx_output_manager.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/io/eamxx_scorpio_interface.hpp"
#include "share/field/field_utils.hpp"
#include "share/util/eamxx_setup_random_test.hpp"

#include "share/grid/point_grid.hpp"

namespace scream {

Field create_f (const std::string& name,
                const FieldLayout layout,
                const std::string& grid_name)
{
  const auto nondim = ekat::units::Units::nondimensional();
  FieldIdentifier fid(name,layout,nondim,grid_name);
  Field f(fid);
  f.allocate_view();
  return f;
}

ekat::ParameterList output_params(const std::string& map_file)
{
  using strvec_t = std::vector<std::string>;

  ekat::ParameterList params;
  params.set<std::string>("filename_prefix","horiz_sampling");
  params.set<std::string>("averaging_type","instant");
  params.set<std::string>("floating_point_precision","real");
  auto& oc = params.sublist("output_control");
  oc.set<int>("frequency",1);
  oc.set<std::string>("frequency_units","nsteps");
  params.set<strvec_t>("field_names",{"s2d","s3d"});
  params.set<std::string>("horiz_remap_file",map_file);

  return params;
}

void print (const std::string& msg, const ekat::Comm& comm) {
  if (comm.am_i_root()) {
    printf("%s",msg.c_str());
  }
}

TEST_CASE("io_remap_test","io_remap_test")
{
  using gid_type = AbstractGrid::gid_type;

  // Init scorpio
  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);

  util::TimeStamp t0 ({2000,1,1},{0,0,0});

  // Random number generation
  using RPDF  = std::uniform_real_distribution<Real>;
  auto engine = setup_random_test(&comm);
  RPDF pdf(0, 1);

  // Create src grid
  const std::string& gname = "point_grid";
  const int ngcols_src = 10*comm.size();
  const int nlevs = 4;
  const auto src_grid = create_point_grid (gname,ngcols_src,nlevs,comm);
  const int  nlcols_src = src_grid->get_num_local_dofs();
  const auto gids_src_h = src_grid->get_dofs_gids().get_view<const gid_type*,Host>();

  // Create remap file. The mapping strategy is simply to take every other column,
  // but making sure to NOT pick the 1st one, so that the min col GID in map file is 2
  print (" -> Create remap file ... \n",comm);
  const int ngcols_tgt = ngcols_src / 2;
  const int nlcols_tgt = nlcols_src / 2;
  std::vector<int> col(nlcols_tgt), row(nlcols_tgt);
  std::vector<Real> S(nlcols_tgt,1.0);
  for (int i=0; i<nlcols_tgt; ++i) {
    row[i] = 1 + i + nlcols_tgt*comm.rank();
    col[i] = 1 + gids_src_h[2*i] + 1;
  }

  // Write remap data to file
  const std::string remap_filename = "horiz_sampling_remap_np"+std::to_string(comm.size())+".nc";
  scorpio::register_file(remap_filename, scorpio::FileMode::Write);

  scorpio::define_dim(remap_filename,"n_a",ngcols_src);
  scorpio::define_dim(remap_filename,"n_b",ngcols_tgt);
  scorpio::define_dim(remap_filename,"n_s",ngcols_tgt);

  scorpio::define_var(remap_filename,"col",   {"n_s"},"int");
  scorpio::define_var(remap_filename,"row",   {"n_s"},"int");
  scorpio::define_var(remap_filename,"S",     {"n_s"},"real");

  // Linear decomposition
  scorpio::set_dim_decomp(remap_filename,"n_s",comm.rank()*nlcols_tgt,nlcols_tgt);
  scorpio::enddef(remap_filename);

  scorpio::write_var(remap_filename,"row",   row.data());
  scorpio::write_var(remap_filename,"col",   col.data());
  scorpio::write_var(remap_filename,"S",     S.data());

  scorpio::release_file(remap_filename);
  print (" -> Create remap file ... done\n",comm);

  // Create random source data
  print (" -> Create source data ... \n",comm);

  // Create the fields and randomize
  auto s2d_src = create_f("s2d",src_grid->get_2d_scalar_layout(),gname);
  auto s3d_src = create_f("s3d",src_grid->get_3d_scalar_layout(true),gname);
  randomize(s2d_src,engine,pdf);
  randomize(s3d_src,engine,pdf);

  // Stuff fields in a FieldManager, since that's what OuputManager wants
  auto fm = std::make_shared<FieldManager> (src_grid,RepoState::Closed);
  fm->add_field(s2d_src);
  fm->add_field(s3d_src);
  fm->init_fields_time_stamp(t0);
  print (" -> Create source data ... done\n",comm);

  print (" -> Write output ... \n",comm);
  double dt = 1.5;
  OutputManager om;
  auto params = output_params(remap_filename);
  om.initialize (comm, params, t0, false);
  om.setup(fm,{gname});

  om.init_timestep(t0,dt);
  om.run(t0+dt);
  om.finalize();
  print (" -> Write output ... done\n",comm);

  print (" -> Check output ... \n",comm);

  // Read output file
  std::string filename = "horiz_sampling.INSTANT.nsteps_x1.np" + std::to_string(comm.size()) + "." + t0.to_string() + ".nc";
  auto tgt_grid = create_point_grid(gname + "_tgt",ngcols_tgt,nlevs,comm);
  auto s2d_tgt = create_f("s2d",tgt_grid->get_2d_scalar_layout(),gname+"_tgt");
  auto s3d_tgt = create_f("s3d",tgt_grid->get_3d_scalar_layout(true),gname+"_tgt");

  AtmosphereInput reader(filename,tgt_grid,{s2d_tgt,s3d_tgt});
  reader.read_variables();
  reader.finalize(); // manually finalize, or scorpio cleanup will complain about a file still open

  // Check values
  auto s2d_src_h = s2d_src.get_view<const Real* ,Host>();
  auto s3d_src_h = s3d_src.get_view<const Real**,Host>();
  auto s2d_tgt_h = s2d_tgt.get_view<const Real* ,Host>();
  auto s3d_tgt_h = s3d_tgt.get_view<const Real**,Host>();
  for (int i=0; i<nlcols_tgt; ++i) {
    REQUIRE (s2d_tgt_h(i)==s2d_src_h(2*i+1));
    for (int k=0; k<nlevs; ++k) {
      REQUIRE (s3d_tgt_h(i,k)==s3d_src_h(2*i+1,k));
    }
  }
  print (" -> Check output ... done\n",comm);
  
  // Cleanup scorpio
  scorpio::finalize_subsystem();
}

} //namespace scream
