#include <catch2/catch.hpp>

#include "share/grid/remap/horizontal_remap_utility.hpp"

#include "share/io/scorpio_input.hpp"
#include "share/io/scream_output_manager.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/grid/point_grid.hpp"

#include "share/field/field_manager.hpp"

#include "share/util/scream_setup_random_test.hpp"

#include "ekat/util/ekat_units.hpp"
#include "ekat/util/ekat_test_utils.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"
namespace {

using namespace scream;

using gid_type = AbstractGrid::gid_type;
using Pack = ekat::Pack<Real,SCREAM_PACK_SIZE>;

using KT = KokkosTypes<DefaultDevice>;

template <typename S>
using view_1d = typename KT::template view_1d<S>;

template <typename S>
using view_2d = typename KT::template view_2d<S>;

template <typename S>
using view_3d = typename KT::template view_3d<S>;

template <typename S>
using view_1d_host = typename KT::template view_1d<S>::HostMirror;

std::shared_ptr<GridsManager> get_test_gm(const ekat::Comm& comm, const Int num_gcols, const Int num_levs);

void run(std::mt19937_64& engine, const ekat::Comm& comm, const gid_type src_min_dof)
{
/* Testing for the HorizontalMap structure and using it to apply a remap to data.
 * This test has three main parts:
 *   1. Construction of the remapping data, construction of the source data to be
 *      remapped, and the establishment of a baseline of remapped data.
 *   2. The construction of a HorizontalMap where the remap data is gathered directly from
 *      the views derived in part 1.  Apply the remap and compare against baseline.
 *   3. After writing the remapping data and source data to file, we construct a
 *      HorizontalMap reampper with data gathered from the file.  We then apply the remap
 *      and compare against the baseline.
 *
 * By breaking the test into these parts we are able to make sure that each test
 * generate random data, i.e. we are not getting lucky with static numbers.  We are
 * also able to test the two main implementations of the horizontal remap structure;
 * the generation of weights using a map constructed online, and the use of remapping
 * weights constructed offline and save to file.
 *
 * INPUT ARGS
 *  - engine: This is just the random seed provided by the EAMxx testing setup
 *  - src_min_dof: This is the degree-of-freedom index for the first column.  This
 *                 simulates that some remapping algorithms may treat the first column
 *                 with index=0 and others with index=1.  Making this an input allows
 *                 the test to check both cases.
 */

  using IPDF = std::uniform_int_distribution<Int>;
  using RPDF = std::uniform_real_distribution<Real>;
  constexpr Real tol  = std::numeric_limits<Real>::epsilon()*10;

//---------------------- Test configuration
  int num_src_cols = 20;  // Number of columns from source data
  int num_tgt_cols = 10;  // Number of columns on target grid
  int num_levels   = 2*SCREAM_PACK_SIZE+1;   // Dummy variable needed for initialization of a grid.

  // Construct a target grid for remapped data
  auto gm_tgt           = get_test_gm(comm, num_tgt_cols, num_levels);
  auto grid_tgt         = gm_tgt->get_grid("Point Grid");     // Retrieve the actual grid from the grids manager.
  auto num_loc_tgt_cols = grid_tgt->get_num_local_dofs();     // Number of columns (dofs) on this rank.
  auto tgt_min_dof      = grid_tgt->get_global_min_dof_gid(); // The starting index of the EAMxx grid.
  auto dofs_gids        = grid_tgt->get_dofs_gids().get_view<const gid_type*>();          // A view of the dofs on this rank with their global id.
  auto dofs_gids_h      = grid_tgt->get_dofs_gids().get_view<const gid_type*,Host>();

//----------------------
// Step 1: Generate a random remapping (source col, target col, weight) and random source data.
//         Store these in local 1D views and create a remap file.
// Note, randomization is done entirely on ROOT so that we only have one version of the remap
// data.
//
// OUT - y_baseline, num_src_to_tgt_cols, n_s, map_src_cols, map_tgt_cols, map_wgts

  // When we create random source data we again only do this on root so that we
  // only have one set of data.
  std::vector<Real> vec_src_data(num_src_cols);
  if (comm.am_i_root()) {
    RPDF pdf_src_data(-1,1);               // Randomizer for generating source data to be mapped from.
    for (int ii=0; ii<num_src_cols; ii++) {
      vec_src_data[ii] = pdf_src_data(engine);
    }
  } // if am_i_root
  comm.barrier();
  comm.broadcast(vec_src_data.data(),vec_src_data.size(),comm.root_rank());
  comm.barrier();
  view_1d<Real> source_remap_data("",num_src_cols);
  auto source_remap_data_h = Kokkos::create_mirror_view(source_remap_data);
  for (int ii=0; ii<num_src_cols; ii++) {
    source_remap_data_h(ii) = vec_src_data[ii];
  }
  Kokkos::deep_copy(source_remap_data, source_remap_data_h);

  // Because we are building the remap  on the fly, use vectors to set up the remap data
  // We can also construct a baseline set of remapped data for testing.
  Int                     n_s_local=0;      // Total number of remapping src->tgt pairs.
  std::vector<Int>        vec_src_col, vec_tgt_col;
  std::vector<Real>       vec_wgt;
  std::map<gid_type,Real> y_baseline;
  std::map<gid_type,Int>  num_src_to_tgt_cols;  // A record of how many source columns map to a specific target.  Will be used later for testing segments.
  IPDF pdf_map_src_num(2,num_src_cols); // Randomizer for number of source columns mapping to a target column
  IPDF pdf_map_src(0,num_src_cols-1);   // Randomizer for which source columns map to a target column
  RPDF pdf_map_wgt(0,1);                // Randomizer for mapping weights
  // Cycle through all tgt_cols and set up a remap
  for (int tcol=0; tcol<num_loc_tgt_cols; tcol++) {
    Int nmap = pdf_map_src_num(engine); // Number of source columns that will map to this target column
    n_s_local += nmap;
    num_src_to_tgt_cols.emplace(dofs_gids_h(tcol),nmap);
    // We need to treat random weights specially, to ensure sum(wgt) = 1;
    std::vector<Real> temp_wgt;
    std::vector<Int>  temp_src;  
    Real wgt_sum = 0.0;
    Real y_sum   = 0.0;
    for (int ii=0; ii<nmap; ii++) {
      Int src_idx = pdf_map_src(engine);
      // In order to mimic a typical remap file which may not be 0 based we offset the raw
      // column dof indices by tgt_min_dof and src_min_dof.
      vec_src_col.push_back(src_idx+src_min_dof);
      vec_tgt_col.push_back(dofs_gids_h(tcol)-tgt_min_dof+src_min_dof);
      Real wgt = pdf_map_wgt(engine);
      wgt_sum += wgt;
      temp_wgt.push_back(wgt);
      y_sum += source_remap_data_h(src_idx) * wgt;
    }
    y_sum /= wgt_sum;  // account for the normalization of the weights
    y_baseline.emplace(dofs_gids_h(tcol),y_sum);
    // Normalize weights
    for (int ii=0; ii<nmap; ii++) {
      temp_wgt[ii] /= wgt_sum;
    }
    // Append normalized weights to vector of overall weights
    vec_wgt.insert(vec_wgt.end(), temp_wgt.begin(), temp_wgt.end());
  }
  // Gather data from all ranks
  Int n_s;
  Int* n_s_all_ranks = new Int[comm.size()]; 
  comm.all_gather(&n_s_local,n_s_all_ranks,1);
  comm.all_reduce(&n_s_local,&n_s,1,MPI_SUM);
  Int n_s_offset = 0; // The offset for this rank, used for writing the remap to output
  for (int ii=0; ii<comm.rank(); ii++) {
    n_s_offset += n_s_all_ranks[ii];
  }
  // Now copy these vectors to views
  view_1d_host<Int>  map_src_cols_h(vec_src_col.data(),vec_wgt.size());
  view_1d_host<Int>  map_tgt_cols_h(vec_tgt_col.data(),vec_wgt.size());
  view_1d_host<Real> map_wgts_h(vec_wgt.data(),vec_wgt.size());
  auto map_src_cols = Kokkos::create_mirror_view(map_src_cols_h);
  auto map_tgt_cols = Kokkos::create_mirror_view(map_tgt_cols_h);
  auto map_wgts     = Kokkos::create_mirror_view(map_wgts_h    );
  Kokkos::deep_copy(map_src_cols, map_src_cols_h);
  Kokkos::deep_copy(map_tgt_cols, map_tgt_cols_h);
  Kokkos::deep_copy(map_wgts    , map_wgts_h    );

//----------------------
  // Step 2: Use HorizontalMap with remap values and source data, and compare against baseline.
  std::map<gid_type,Real> y_from_views;
  HorizontalMap remap_from_views(comm,"From Views", dofs_gids, tgt_min_dof);
  for (int iseg=0; iseg<num_loc_tgt_cols; iseg++) {
    gid_type seg_dof = dofs_gids_h(iseg);
    Int seg_len = num_src_to_tgt_cols.at(seg_dof+tgt_min_dof); // Need to offset dof index because this is a look up.
    // The following is done on HOST given the way the random data has been organized
    view_1d<gid_type> source_dofs("",seg_len);
    view_1d<Real>     weights("",seg_len);    
    auto source_dofs_h = Kokkos::create_mirror_view(source_dofs);
    auto weights_h     = Kokkos::create_mirror_view(weights);
    int idx = 0;
    for (int ii=0; ii<map_src_cols_h.extent_int(0); ii++) {
      if (map_tgt_cols_h(ii)-src_min_dof==seg_dof) {
        source_dofs_h(idx) = map_src_cols_h(ii)-src_min_dof;  // Make 0-based dof indexing
        weights_h(idx)     = map_wgts_h(ii);
        idx ++;
      }
    }
    Kokkos::deep_copy(source_dofs,source_dofs_h);
    Kokkos::deep_copy(weights,weights_h);
    // Create a HorizontalMapSegment for this mapping
    HorizontalMapSegment seg(seg_dof,seg_len, source_dofs, weights);
    remap_from_views.add_remap_segment(seg);
  }
  remap_from_views.set_unique_source_dofs();
  remap_from_views.check();
//  remap_from_views.print();  // Un-comment this line when debugging a failure.  Should give a comprehensive printing of the map
  // We need to extract just that source data that is related to this ranks need.
  const auto& unique_dofs_from_views = remap_from_views.get_unique_source_dofs();
  Int num_unique_dofs_from_views = unique_dofs_from_views.extent(0);
  view_1d<Real> x_data_from_views("",num_unique_dofs_from_views);
  Kokkos::parallel_for("", num_unique_dofs_from_views, KOKKOS_LAMBDA (const int& ii) {
    const int idx = unique_dofs_from_views(ii);
    x_data_from_views(ii) = source_remap_data(idx);
  });
  // Apply the remap
  view_1d<Real> y_data_from_views("",num_loc_tgt_cols);
  remap_from_views.apply_remap(x_data_from_views,y_data_from_views);
  // The remapped values should match the y_baseline
  auto y_data_from_views_h = Kokkos::create_mirror_view(y_data_from_views);
  Kokkos::deep_copy(y_data_from_views_h,y_data_from_views);
  for (int ii=0; ii<num_loc_tgt_cols; ii++) {
    gid_type dof = dofs_gids_h(ii);  // Offset to 0-based
    auto y_base = y_baseline[dof];
    REQUIRE(std::abs(y_data_from_views_h(ii)-y_base)<tol);
  }
//----------------------
  // Step 3: Write source data and remap data to a file.
  std::string filename = "horizontal_remap_test_data_" + std::to_string(comm.size()) + "_ranks_" + std::to_string(src_min_dof) + "_case.nc";
  // Initialize the pio_subsystem for this test:
  MPI_Fint fcomm = MPI_Comm_c2f(comm.mpi_comm()); 
  scorpio::eam_init_pio_subsystem(fcomm);   // Gather the initial PIO subsystem data creater by component coupler
  // Register the output file
  scorpio::register_file(filename,scorpio::Write);
  //   - Dimensions
  scorpio::register_dimension(filename,"n_s","n_s",n_s,true);
  scorpio::register_dimension(filename,"n_a","n_a",num_src_cols,false);
  scorpio::register_dimension(filename,"n_b","n_b",num_tgt_cols,false);
  scorpio::register_dimension(filename,"ncol","ncol",num_src_cols,true);
  //   - Variables
  std::string remap_decomp_tag_r = "n_s_real";
  std::string data_decomp_tag_r  = "data_real";
  std::string remap_decomp_tag_i = "n_s_int";
  std::vector<std::string> vec_of_remap_dims = {"n_s"};
  std::vector<std::string> vec_of_data_dims  = {"ncol"};
  std::string int_type = "int";
  scorpio::register_variable(filename,"col","col","unitless",vec_of_remap_dims,"int","int",remap_decomp_tag_i);
  scorpio::register_variable(filename,"row","row","unitless",vec_of_remap_dims,"int","int",remap_decomp_tag_i);
  scorpio::register_variable(filename,"S","S",    "unitless",vec_of_remap_dims,"real","real",remap_decomp_tag_r);
  scorpio::register_variable(filename,"src_data","src_data", "m",vec_of_data_dims,"real","real",data_decomp_tag_r);
  //   - DOFs
  std::vector<int64_t> var_dof(n_s_local);
  std::iota(var_dof.begin(),var_dof.end(),n_s_offset);
  scorpio::set_dof(filename,"col",var_dof.size(),var_dof.data());
  scorpio::set_dof(filename,"row",var_dof.size(),var_dof.data());
  scorpio::set_dof(filename,"S",var_dof.size(),var_dof.data());
  var_dof.resize(unique_dofs_from_views.size());
  auto unique_dofs_from_views_h = Kokkos::create_mirror_view(unique_dofs_from_views);
  Kokkos::deep_copy(unique_dofs_from_views_h,unique_dofs_from_views);
  for (size_t ii=0; ii<var_dof.size(); ii++) {
    var_dof[ii] = unique_dofs_from_views_h(ii);
  }
  scorpio::set_dof(filename,"src_data",var_dof.size(),var_dof.data());
  scorpio::eam_pio_enddef(filename);
  // Write data
  auto x_data_from_views_h = Kokkos::create_mirror_view(x_data_from_views);
  Kokkos::deep_copy(x_data_from_views_h,x_data_from_views);
  scorpio::grid_write_data_array(filename,"col",map_src_cols_h.data(),vec_wgt.size());
  scorpio::grid_write_data_array(filename,"row",map_tgt_cols_h.data(),vec_wgt.size());
  scorpio::grid_write_data_array(filename,"S",map_wgts_h.data(),vec_wgt.size());
  scorpio::grid_write_data_array(filename,"src_data",x_data_from_views_h.data(),num_unique_dofs_from_views);
  // All done writing 
  scorpio::eam_pio_closefile(filename);

//----------------------
  // Step 4: Load remap information and source data from the file created in step 3
  //         and use this data to test remapping from file.  Compare remapped data
  //         against the baseline. 
  HorizontalMap remap_from_file(comm,"From File", dofs_gids, tgt_min_dof);
  remap_from_file.set_remap_segments_from_file(filename);
  remap_from_file.set_unique_source_dofs();
  remap_from_file.check();
  auto unique_dofs_from_file = remap_from_file.get_unique_source_dofs();
  auto unique_dofs_from_file_h = Kokkos::create_mirror_view(unique_dofs_from_file);
  Kokkos::deep_copy(unique_dofs_from_file_h,unique_dofs_from_file);
  // Read source data at unique points
  scorpio::register_file(filename,scorpio::Read);
  scorpio::get_variable(filename,"src_data","src_data",vec_of_data_dims,"real",data_decomp_tag_r);
  var_dof.resize(unique_dofs_from_file.size());
  for (size_t ii=0; ii<var_dof.size(); ii++) {
    var_dof[ii] = unique_dofs_from_file_h(ii);
  }
  scorpio::set_dof(filename,"src_data",var_dof.size(),var_dof.data());
  scorpio::set_decomp(filename);
  view_1d<Real> x_data_from_file("",var_dof.size());
  auto x_data_from_file_h = Kokkos::create_mirror_view(x_data_from_file);
  scorpio::grid_read_data_array(filename,"src_data",0,x_data_from_file.data(),x_data_from_file.size()); 
  scorpio::eam_pio_closefile(filename);
  Kokkos::deep_copy(x_data_from_file_h,x_data_from_file);
  // Apply remap using the data
  view_1d<Real> y_data_from_file("",num_loc_tgt_cols);
  remap_from_file.apply_remap(x_data_from_file,y_data_from_file);
  // The remapped values should match the y_baseline
  auto y_data_from_file_h = Kokkos::create_mirror_view(y_data_from_file);
  Kokkos::deep_copy(y_data_from_file_h,y_data_from_file);
  for (int ii=0; ii<num_loc_tgt_cols; ii++) {
    gid_type dof = dofs_gids_h(ii);
    auto y_base = y_baseline[dof];
    REQUIRE(std::abs(y_data_from_file_h(ii)-y_base)<tol);
  } 

  // All Done with testing
  scorpio::eam_pio_finalize();

  // Test to check that a 2D and 3D views work
  Int num_bands = 4;
  view_2d<Real> x_2d_data("",unique_dofs_from_file.size(),num_levels);
  view_2d<Real> y_2d_data("",num_loc_tgt_cols,num_levels);
  view_3d<Real> x_3d_data("",unique_dofs_from_file.size(),num_bands,num_levels);
  view_3d<Real> y_3d_data("",num_loc_tgt_cols,num_bands,num_levels);
  auto x_2d_data_h = Kokkos::create_mirror_view(x_2d_data);
  auto x_3d_data_h = Kokkos::create_mirror_view(x_3d_data);
  for (size_t ii=0;ii<unique_dofs_from_file.size();ii++) {
    for (int kk=0; kk<num_levels; kk++) {
      x_2d_data_h(ii,kk) = x_data_from_file_h(ii)*(kk+1);
      for (int nn=0; nn<num_bands; nn++) {
        x_3d_data_h(ii,nn,kk) = x_data_from_file_h(ii)*(kk+1)*(nn+1);
      }
    }
  }
  Kokkos::deep_copy(x_2d_data,x_2d_data_h);
  remap_from_file.apply_remap(x_2d_data,y_2d_data);
  Kokkos::deep_copy(x_3d_data,x_3d_data_h);
  remap_from_file.apply_remap(x_3d_data,y_3d_data);
  auto y_2d_data_h = Kokkos::create_mirror_view(y_2d_data);
  auto y_3d_data_h = Kokkos::create_mirror_view(y_3d_data);
  Kokkos::deep_copy(y_2d_data_h,y_2d_data);
  Kokkos::deep_copy(y_3d_data_h,y_3d_data);
  for (int ii=0; ii<num_loc_tgt_cols; ii++) {
    gid_type dof = dofs_gids_h(ii);
    auto y_base = y_baseline[dof];
    for (int kk=0; kk<num_levels; kk++) {
      REQUIRE(std::abs(y_2d_data_h(ii,kk)-(kk+1)*y_base)<tol*10);
      for (int nn=0; nn<num_bands; nn++) {
        REQUIRE(std::abs(y_3d_data_h(ii,nn,kk)-(nn+1)*(kk+1)*y_base)<tol*100);
      }
    }
  } 
  
  
} // end function run

//===============================================================================
TEST_CASE("horizontal_remap_test", "[horizontal_remap_test]"){

  using scream::Real;

  ekat::Comm comm (MPI_COMM_WORLD);
  auto engine = scream::setup_random_test();

  if (comm.am_i_root()) {
    printf(" -> Testing horizontal remapping for minimum source dof = 0...");
  }
  run(engine, comm, 0);
  if (comm.am_i_root()) {
    printf("ok!\n");
  }

  if (comm.am_i_root()) {
    printf(" -> Testing horizontal remapping for minimum source dof = 1...");
  }
  run(engine, comm, 1);
  if (comm.am_i_root()) {
    printf("ok!\n");
  }

} // TEST_CASE
//===============================================================================

TEST_CASE("horizontal_remap_units", "") {

  using namespace scream;
  using IPDF = std::uniform_int_distribution<int>;
  using RPDF = std::uniform_real_distribution<Real>;
  auto engine = scream::setup_random_test();
  ekat::Comm comm (MPI_COMM_WORLD);

  // Test the remap_segment structure
  {
    IPDF pdf_seg_len(2,100);
    RPDF pdf_seg_wgt(0,1);
    Int len = pdf_seg_len(engine);
    HorizontalMapSegment test_segment(0,len);
    auto seg_weights = test_segment.get_weights();
    ekat::genRandArray(seg_weights,engine,pdf_seg_wgt);
    Real wgt_normalize = 0;
    auto weights_h = Kokkos::create_mirror_view(seg_weights);
    Kokkos::deep_copy(weights_h,seg_weights);
    for (int ii=0; ii<len; ii++) {
      wgt_normalize += weights_h(ii);
    }
    for (int ii=0; ii<len; ii++) {
      weights_h(ii) /= wgt_normalize;
    };
    Kokkos::deep_copy(seg_weights,weights_h);
    REQUIRE(test_segment.check());
    // Make it fail by adjusting one of the weights so that they no longer add to 1.0
    weights_h(0) += 0.1;
    Kokkos::deep_copy(seg_weights,weights_h);
    REQUIRE(!test_segment.check());
    weights_h(0) -= 0.2;
    Kokkos::deep_copy(seg_weights,weights_h);
    REQUIRE(!test_segment.check());
  }


  // Create a map with two remap segments for some simple map structure testing.
  {
    IPDF pdf_seg_len(2,100);
    RPDF pdf_seg_wgt(0,1);
    // Create a HorizontalMap to test
    HorizontalMap test_map(comm,"Test Map");
    Int num_dofs = 2;
    view_1d<gid_type> dof_gids("",num_dofs);
    Kokkos::parallel_for("", num_dofs, KOKKOS_LAMBDA (const int& ii) {
      dof_gids(ii) = ii;
    });
    test_map.set_dof_gids(dof_gids,0);
    // Create a remap segment
    Int len = pdf_seg_len(engine);
    HorizontalMapSegment test_seg(1,len);
    auto test_weights = test_seg.get_weights();
    ekat::genRandArray(test_weights,engine,pdf_seg_wgt);
    Real wgt_normalize = 0;
    auto weights_h = Kokkos::create_mirror_view(test_weights);
    Kokkos::deep_copy(weights_h,test_weights);
    for (int ii=0; ii<len; ii++) {
      wgt_normalize += weights_h(ii);
    }
    for (int ii=0; ii<len; ii++) {
      weights_h(ii) /= wgt_normalize;
    };
    Kokkos::deep_copy(test_weights,weights_h);
    // Add segment to map
    test_map.add_remap_segment(test_seg);
    test_map.check();
    REQUIRE(test_map.get_num_of_segments()==1);
    auto check_seg = test_map.get_map_segments()[0];
    REQUIRE(check_seg.get_length()==len);

    // Set a wrong weight, to ensure the check method can spot it
    auto check_wgts = check_seg.get_weights();
    weights_h = Kokkos::create_mirror_view(check_wgts);
    Kokkos::deep_copy(weights_h,check_wgts);
    Real wgt_to_add = weights_h(0)/2;
    weights_h(0) = 1-weights_h(0);
    Kokkos::deep_copy(check_wgts,weights_h);
    REQUIRE_THROWS(test_map.check());
    // Restore w/2, so that adding a complementary weight segment should
    // restore the valid map
    weights_h(0) = wgt_to_add;
    Kokkos::deep_copy(check_wgts,weights_h);

    // Create a new segment that has the same DOF as the first segment
    HorizontalMapSegment new_seg(1,1);
    auto new_wgts = new_seg.get_weights();
    weights_h = Kokkos::create_mirror_view(new_wgts);
    weights_h(0)  = wgt_to_add;
    Kokkos::deep_copy(new_wgts,weights_h);
    test_map.add_remap_segment(new_seg);
    // Because the new segment is the same DOF, add_segment should have
    // appended it to the first segment.
    REQUIRE(test_map.get_num_of_segments()==1);
    test_map.check();
    // Create a completely new segment for a different DOF and check
    HorizontalMapSegment newer_seg(2,1);
    auto newer_wgts = newer_seg.get_weights();
    Kokkos::deep_copy(newer_wgts,1.0);
    test_map.add_remap_segment(newer_seg);
    REQUIRE(test_map.get_num_of_segments()==2);
    test_map.check();
  }

  // Test retrieving unique columns from a full map 
  {
    // Create a HorizontalMap to test
    HorizontalMap test_map(comm,"Test Map 2");
    // Create a set of segments each with overlapping source dofs
    IPDF pdf_seg_len(4,10),
         pdf_num_seg(5,15);
    Int num_of_segs = pdf_num_seg(engine);
    gid_type seg_start = 0;
    gid_type dof_end;
    Int seg_len;
    for (int iseg=0; iseg<num_of_segs; iseg++) {
      seg_len = pdf_seg_len(engine);
      dof_end = seg_start + seg_len;
      HorizontalMapSegment seg(iseg,seg_len);
      auto source_dofs = seg.get_source_dofs();
      Kokkos::parallel_for("", seg_len, KOKKOS_LAMBDA (const int& ii) {
        source_dofs(ii) = seg_start + ii;
      });
      test_map.add_remap_segment(seg);
      seg_start += seg_len-2;
    }
    Int num_dofs = dof_end;
    view_1d<gid_type> dof_gids("",num_dofs);
    Kokkos::parallel_for("", num_dofs, KOKKOS_LAMBDA (const int& ii) {
      dof_gids(ii) = ii;
    });
    test_map.set_dof_gids(dof_gids,0);
    // All of the segments combined should overlap by 2 on either side.  The set of unique
    // columns for all segments should be 0,1,...,seg_start+seg_len-1.  So we check that.
    test_map.set_unique_source_dofs();
    const auto& unique_dofs = test_map.get_unique_source_dofs();
    REQUIRE(unique_dofs.extent_int(0) == dof_end);
    Kokkos::parallel_for("", dof_end, KOKKOS_LAMBDA (const int& ii) {
      EKAT_KERNEL_REQUIRE(unique_dofs(ii) == ii);
    });
  }

} // TEST_CASE horizontal remap
/*===================================================================================================*/
std::shared_ptr<GridsManager> get_test_gm(const ekat::Comm& comm, const Int num_gcols, const Int num_levs)
{
/* Simple routine to construct and return a grids manager given number of columns and levels */
  ekat::ParameterList gm_params;
  gm_params.set("number_of_global_columns",num_gcols);
  gm_params.set("number_of_vertical_levels",num_levs);
  auto gm = create_mesh_free_grids_manager(comm,gm_params);
  gm->build_grids();
  return gm;
} // end get_test_gm
/*==========================================================================================================*/

} // end namespace
