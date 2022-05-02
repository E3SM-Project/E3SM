#ifndef SPA_FUNCTIONS_IMPL_HPP
#define SPA_FUNCTIONS_IMPL_HPP

#include "share/scream_types.hpp"
#include "share/io/scream_scorpio_interface.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/grid/point_grid.hpp"
#include "physics/share/physics_constants.hpp"

#include "ekat/kokkos/ekat_subview_utils.hpp"
#include "ekat/util/ekat_lin_interp.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"

#include <numeric>

/*-----------------------------------------------------------------
 * The main SPA routines used to convert SPA data into a format that
 * is usable by the rest of the atmosphere processes.
 *
 * SPA or Simple Prescribed Aerosols provides a way to prescribe
 * aerosols for an atmospheric simulation using pre-computed data.
 * 
 * The data is typically provided at a frequency of monthly, and
 * does not necessarily have to be on the same horizontal or vertical
 * domain as the atmospheric simulation.
 *
 * In order to accomodate coarse temporal resolution and a potentially
 * different spatial resolution it is necessary to perform a series
 * of interpolations, which make up the main body of the SPA routines.
 *
 * The interpolations can be broken into three categories.
 * 1. Horizontal Interpolation: TODO - not done yet.
 * The SPA data set does not have to be provided on the same grid as
 * the atmospheric simulation.  Whenever SPA data is loaded, it is
 * interpolated horizontally onto the simulation grid to provide
 * forcing at every location.  This can be done by,
 *   a) Preloaded remapping wieghts which are applied at every
 *      horizontal column.
 *   b) Online calculation of remapping weights given a set of lat/lon
 *      for the source data and comparing it with the lat/lon of each
 *      column in the simulation.  TODO: THIS HAS NOT BEEN IMPLEMENTED YET
 *
 * 2. Temporal Interpolation:
 * As noted above, the SPA data is provided at some fixed frequency.  Typically
 * as monthly data.  As a result, the data must be interpolated to the current
 * time of the simulation at each time step.  Temporal interpolation follows
 * a basic linear interpolation and is performed for all SPA data at all columns
 * and levels.
 * Note: There is also a temporal interpolation of the surface pressure for the SPA
 * data, which is used in the vertical reconstruction of the pressure profile.
 *
 * 3. Vertical Interpolation:
 * Given that the SPA data has been generated elsewhere it is very likely that
 * the vertical pressure profiles of the data won't match the simulation pressure
 * profiles.  The vertical SPA data structure must be remapped onto the simulation
 * pressure profile.
 * This is done using the EKAT linear interpolation code, see /externals/ekat/util/ekat_lin_interp.hpp
 * The SPA pressure profiles are calculated using the surface pressure which was
 * temporally interpolated in the last step and the set of hybrid coordinates (hyam and hybm)
 * that are used in EAM to construct the physics pressure profiles.
 * The SPA data is then projected onto the simulation pressure profile (pmid)
 * using EKAT linear interpolation. 
-----------------------------------------------------------------*/

namespace scream {
namespace spa {

/*-----------------------------------------------------------------*/
// The main SPA routine which handles projecting SPA data onto the
// horizontal columns and vertical pressure profiles of the atmospheric
// state.
// Inputs:
//   time_state: A structure defined in spa_functions.hpp which handles
//     the current temporal state of the simulation.
//   p_tgt: the vertical pressure profile for the atmospheric simulation state
//   data_beg: A structure defined in spa_functions.hpp which handles the full
//     set of SPA data for the beginning of the month.
//   data_end: Similar to data_beg, but for SPA data for the end of the month.
//   data_out: A structure defined in spa_functions.hpp which handles the full
//     set of SPA data projected onto the pressure profile of the current atmosphere
//     state.  This is the data that will be passed to other processes.
//   ncols_tgt, nlevs_tgt: The number of columns and levels in the simulation grid.
//     (not to be confused with the number of columns and levels used for the SPA data, 
//      which can be different.)
//   nswbands, nlwbands: The number of shortwave (sw) and longwave (lw) aerosol bands 
//     for the data that will be passed to radiation.
template <typename S, typename D>
void SPAFunctions<S,D>
::spa_main(
  const SPATimeState& time_state,
  const view_2d<const Spack>& p_tgt,
  const view_2d<      Spack>& p_src,
  const SPAInput&   data_beg,
  const SPAInput&   data_end,
  const SPAInput&   data_tmp,
  const SPAOutput&  data_out)
{
  // Beg/End/Tmp month must have all sizes matching
  EKAT_REQUIRE_MSG (
      data_end.data.nswbands==data_beg.data.nswbands &&
      data_end.data.nswbands==data_tmp.data.nswbands &&
      data_end.data.nlwbands==data_beg.data.nlwbands &&
      data_end.data.nlwbands==data_tmp.data.nlwbands,
      "Error! SPAInput data structs must have the same number of SW/LW bands.\n");
  EKAT_REQUIRE_MSG (
      data_end.data.ncols==data_beg.data.ncols &&
      data_end.data.ncols==data_tmp.data.ncols &&
      data_end.data.nlevs==data_beg.data.nlevs &&
      data_end.data.nlevs==data_tmp.data.nlevs,
      "Error! SPAInput data structs must have the same number of columns/levels.\n");

  // Output must have same number of bands
  EKAT_REQUIRE_MSG (
      data_end.data.nswbands==data_out.nswbands &&
      data_end.data.nlwbands==data_out.nlwbands,
      "Error! SPAInput and SPAOutput data structs must have the same number of SW/LW bands.\n");

  // Horiz interpolation can be expensive, and does not depend on the particular time of
  // the month, so it can be done ONCE per month, *outside* spa_main (when updating
  // the beg/end states, reading them from file).
  EKAT_REQUIRE_MSG (
      data_end.data.ncols==data_out.ncols,
      "Error! Horizontal interpolation is performed *before* calling spa_main,\n"
      "       SPAInput and SPAOutput data structs must have the same number columns.\n");

  // Step 1. Perform time interpolation
  perform_time_interpolation(time_state,data_beg,data_end,data_tmp);

  // Step 2. Compute source pressure levels
  compute_source_pressure_levels(data_tmp.PS, p_src, data_beg.hyam, data_beg.hybm);

  // Step 3. Perform vertical interpolation
  perform_vertical_interpolation(p_src, p_tgt, data_tmp.data, data_out);
}

/*-----------------------------------------------------------------*/
template <typename S, typename D>
void SPAFunctions<S,D>
::perform_time_interpolation(
  const SPATimeState& time_state,
  const SPAInput&  data_beg,
  const SPAInput&  data_end,
  const SPAInput&  data_out)
{
  // NOTE: we *assume* data_beg and data_end have the *same* hybrid v coords.
  //       IF this ever ceases to be the case, you can interp those too.

  using ExeSpace = typename KT::ExeSpace;
  using ESU = ekat::ExeSpaceUtils<ExeSpace>;

  // Gather time stamp info
  auto& t_now = time_state.t_now;
  auto& t_beg = time_state.t_beg_month;
  auto& delta_t = time_state.days_this_month;

  // Makes no sense to have different number of bands
  EKAT_REQUIRE(data_end.data.nswbands==data_beg.data.nswbands);
  EKAT_REQUIRE(data_end.data.nlwbands==data_beg.data.nlwbands);

  // At this stage, begin/end must have the same dimensions
  EKAT_REQUIRE(data_end.data.ncols==data_beg.data.ncols);
  EKAT_REQUIRE(data_end.data.nlevs==data_beg.data.nlevs);

  // We can ||ize over columns as well as over variables and bands
  const int num_vars = 1+data_beg.data.nswbands*3+data_beg.data.nlwbands;
  const int outer_iters = data_beg.data.ncols*num_vars;
  const int num_vert_packs = ekat::PackInfo<Spack::n>::num_packs(data_beg.data.nlevs);
  const auto policy = ESU::get_default_team_policy(outer_iters, num_vert_packs);

  auto delta_t_fraction = (t_now-t_beg) / delta_t;

  Kokkos::parallel_for("spa_time_interp_loop", policy,
    KOKKOS_LAMBDA(const MemberType& team) {

    // The policy is over ncols*num_vars, so retrieve icol/ivar
    const int icol = team.league_rank() / num_vars;
    const int ivar = team.league_rank() % num_vars;

    // Compute ps out only once
    if (ivar==0) {
      // PS is a 2d var, so we need to make one team member handle it.
      Kokkos::single(Kokkos::PerTeam(team),[&]{
          data_out.PS(icol) = linear_interp(data_beg.PS(icol),data_end.PS(icol),delta_t_fraction);
      });
    }

    // Get column of beg/end/out variable
    auto var_beg = get_var_column (data_beg.data,icol,ivar);
    auto var_end = get_var_column (data_end.data,icol,ivar);
    auto var_out = get_var_column (data_out.data,icol,ivar);

    Kokkos::parallel_for (Kokkos::TeamThreadRange(team,num_vert_packs),
                          [&] (const int& k) {
      var_out(k) = linear_interp(var_beg(k),var_end(k),delta_t_fraction);
    });
  });
  Kokkos::fence();
}

template<typename S, typename D>
void SPAFunctions<S,D>::
compute_source_pressure_levels(
  const view_1d<const Real>& ps_src,
  const view_2d<      Spack>& p_src,
  const view_1d<const Spack>& hyam,
  const view_1d<const Spack>& hybm)
{
  using ExeSpace = typename KT::ExeSpace;
  using ESU = ekat::ExeSpaceUtils<ExeSpace>;
  using C = scream::physics::Constants<Real>;

  constexpr auto P0 = C::P0;

  const int ncols = ps_src.extent(0);
  const int num_vert_packs = p_src.extent(1);
  const auto policy = ESU::get_default_team_policy(ncols, num_vert_packs);

  Kokkos::parallel_for("spa_compute_p_src_loop", policy,
    KOKKOS_LAMBDA (const MemberType& team) {
    const int icol = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,num_vert_packs),
                         [&](const int k) {
      p_src(icol,k) = ps_src(icol) * hybm(k)  + P0 * hyam(k);
    });
  });
}

template<typename S, typename D>
void SPAFunctions<S,D>::
perform_vertical_interpolation(
  const view_2d<const Spack>& p_src,
  const view_2d<const Spack>& p_tgt,
  const SPAData& input,
  const SPAData& output)
{
  using ExeSpace = typename KT::ExeSpace;
  using ESU = ekat::ExeSpaceUtils<ExeSpace>;
  using LIV = ekat::LinInterp<Real,Spack::n>;

  // Makes no sense to have different number of bands
  EKAT_REQUIRE(input.nswbands==output.nswbands);
  EKAT_REQUIRE(input.nlwbands==output.nlwbands);

  // At this stage, begin/end must have the same horiz dimensions
  EKAT_REQUIRE(input.ncols==output.ncols);

  const int ncols     = input.ncols;
  const int nlevs_src = input.nlevs;
  const int nlevs_tgt = output.nlevs;

  LIV vert_interp(ncols,nlevs_src,nlevs_tgt);

  // We can ||ize over columns as well as over variables and bands
  const int num_vars = 1+input.nswbands*3+input.nlwbands;
  const int num_vert_packs = ekat::PackInfo<Spack::n>::num_packs(nlevs_tgt);
  const auto policy_setup = ESU::get_default_team_policy(ncols, num_vert_packs);

  // Setup the linear interpolation object
  Kokkos::parallel_for("spa_vert_interp_setup_loop", policy_setup,
    KOKKOS_LAMBDA(typename LIV::MemberType const& team) {

    const int icol = team.league_rank();
  
    // Setup
    vert_interp.setup(team, ekat::subview(p_src,icol),
                            ekat::subview(p_tgt,icol));
  });
  Kokkos::fence();

  // Now use the interpolation object in || over all variables.
  const int outer_iters = ncols*num_vars;
  const auto policy_interp = ESU::get_default_team_policy(outer_iters, num_vert_packs);
  Kokkos::parallel_for("spa_vert_interp_loop", policy_interp,
    KOKKOS_LAMBDA(typename LIV::MemberType const& team) {

    const int icol = team.league_rank() / num_vars;
    const int ivar = team.league_rank() % num_vars;

    const auto x1 = ekat::subview(p_src,icol);
    const auto x2 = ekat::subview(p_tgt,icol);

    const auto y1 = get_var_column(input, icol,ivar);
    const auto y2 = get_var_column(output,icol,ivar);

    vert_interp.lin_interp(team, x1, x2, y1, y2, icol);
  });
  Kokkos::fence();
}

/*-----------------------------------------------------------------*/
// Function to read the weights for conducting horizontal remapping
// from a file.
template <typename S, typename D>
void SPAFunctions<S,D>
::set_remap_weights_one_to_one(
    const Int                ncols_scream,
    gid_type                 min_dof,
    const view_1d<gid_type>& dofs_gids,
          SPAHorizInterp&    spa_horiz_interp
  )
{
  // There may be cases where the SPA data is defined on the same grid as the simulation
  // and thus no remapping is required.  This simple routine establishes a 1-1 horizontal
  // mapping.
  auto num_local_cols = dofs_gids.size();
  spa_horiz_interp.length            = num_local_cols; 
  spa_horiz_interp.source_grid_ncols = ncols_scream;
  spa_horiz_interp.weights           = view_1d<Real>("",spa_horiz_interp.length);
  spa_horiz_interp.source_grid_loc   = view_1d<gid_type> ("",spa_horiz_interp.length);
  spa_horiz_interp.target_grid_loc   = view_1d<gid_type> ("",spa_horiz_interp.length);
  Kokkos::deep_copy(spa_horiz_interp.weights,1.0);
  Kokkos::parallel_for("", num_local_cols, KOKKOS_LAMBDA(const int& ii) {
    spa_horiz_interp.target_grid_loc(ii) = ii;
    // Note we are interested in the vector index, not the actual global-id 
    // Here we want the index in a the source data vector corresponding to this column
    // which needs to be offset by the minimum degree of freedom in the whole grid.
    // That way the first global-id will map to the 0th entry in the source grid data.
    spa_horiz_interp.source_grid_loc(ii) = dofs_gids(ii) - min_dof;
  });
  // Determine the set of unique columns in this remapping
  spa_horiz_interp.set_unique_cols();
} // END set_remap_weights_one_to_one
/*-----------------------------------------------------------------*/
// Function to read the weights for conducting horizontal remapping
// from a file.
template <typename S, typename D>
void SPAFunctions<S,D>
::get_remap_weights_from_file(
    const std::string&       remap_file_name,
    const Int                ncols_scream,
    gid_type                 min_dof,
    const view_1d<gid_type>& dofs_gids,
          SPAHorizInterp&    spa_horiz_interp
  )
{
  // Note, the remap file doesn't follow a conventional grid setup so
  // here we manually go through all of the input steps rather than
  // use the scorpio_input class.

  // Open input file: 
  scorpio::register_file(remap_file_name,scorpio::Read);

  // Gather the size of the remap data from file.
  // NOTE, we are currently assuming the remap file was generated by NCO and thus follows
  // NCO conventions.
  // As such the dimensions have rather general names.  Here is an explanation:
  //   n_s is the total number of nonzero remap weights.
  //   n_a is the total number of columns in the SOURCE grid (i.e. in an ne2 -> ne4 remap this would have 218 columns)
  //   n_b is the total number of columns in the TARGET grid (i.e. in an ne2 -> ne4 remap this would have 866 columns)
  // Conceptually, if we wrote the remapping as y = W*x then W would be a matrix with n_a columns, n_b rows and n_s total nonzero entries.
  // The data is stored in sparse format, so we have the following variables:
  //   S   is a vector of length n_s that represents all of the nonzero remapping wieghts
  //   col is a vector of length n_s that contains the indices of the SOURCE grid column associated with the appropriate weight in S
  //   row is a vector of length n_s that contains the indices of the TARGET grid column associated with the appropriate weight in S
  // Thus considering the matrix W described above.  If W[i,j] = w, and this was the n'th nonzero weight in W then
  //   S[n]   = w
  //   col[N] = j
  //   row[N] = i
  // TODO: provide infrastructure for using different horizontal remap files.
  spa_horiz_interp.length = scorpio::get_dimlen_c2f(remap_file_name.c_str(),"n_s");
  // And the number of columns that should be in the data source file
  spa_horiz_interp.source_grid_ncols = scorpio::get_dimlen_c2f(remap_file_name.c_str(),"n_a");
  // Check that the target grid size matches the remap file
  Int target_ncols = scorpio::get_dimlen_c2f(remap_file_name.c_str(),"n_b");
  EKAT_REQUIRE_MSG(target_ncols==ncols_scream,"ERROR: SPA get_remap_weights_from_file, remap target domain does not match simulation domain size");

  // Construct local arrays to read data into
  view_1d<Real> S_global("weights",spa_horiz_interp.length);
  view_1d<Int>  row_global("row",spa_horiz_interp.length); 
  view_1d<Int>  col_global("col",spa_horiz_interp.length); 
  auto S_global_h   = Kokkos::create_mirror_view(S_global);
  auto col_global_h = Kokkos::create_mirror_view(col_global); // Note, in remap files col -> row (src -> tgt)
  auto row_global_h = Kokkos::create_mirror_view(row_global);

  // Setup the scorpio structures needed for input
  // Register variables for input
  std::vector<std::string> vec_of_dims = {"n_s"};
  std::string r_decomp = "Real-n_s";
  std::string i_decomp = "Int-n_s";
  scorpio::get_variable(remap_file_name, "S", "S", vec_of_dims.size(), vec_of_dims, PIO_REAL, r_decomp);
  scorpio::get_variable(remap_file_name, "row", "row", vec_of_dims.size(), vec_of_dims, PIO_INT, i_decomp);
  scorpio::get_variable(remap_file_name, "col", "col", vec_of_dims.size(), vec_of_dims, PIO_INT, i_decomp);
  // Set the dof's to read in variables, since we will have all mpi ranks read in the full set of data the dof's are the whole array
  std::vector<int> var_dof(spa_horiz_interp.length);
  std::iota(var_dof.begin(),var_dof.end(),0);
  scorpio::set_dof(remap_file_name,"S",var_dof.size(),var_dof.data());
  scorpio::set_dof(remap_file_name,"row",var_dof.size(),var_dof.data());
  scorpio::set_dof(remap_file_name,"col",var_dof.size(),var_dof.data());
  scorpio::set_decomp(remap_file_name);
  
  // Now read all of the input
  scorpio::grid_read_data_array(remap_file_name,"S",0,S_global_h.data()); 
  scorpio::grid_read_data_array(remap_file_name,"row",0,row_global_h.data()); 
  scorpio::grid_read_data_array(remap_file_name,"col",0,col_global_h.data()); 

  // Finished, close the file
  scorpio::eam_pio_closefile(remap_file_name);

  // Retain only the information needed on this rank. 
  auto dofs_gids_h = Kokkos::create_mirror_view(dofs_gids);
  Kokkos::deep_copy(dofs_gids_h,dofs_gids);
  std::vector<int> local_idx;
  std::vector<int> global_idx;
  for (int idx=0;idx<spa_horiz_interp.length;idx++) {
    int dof = row_global_h(idx) - (1-min_dof); // Note, the dof ids may start with 0 or 1.  In the data they certainly start with 1.  This maps 1 -> true min dof.
    for (int id=0;id<dofs_gids_h.extent_int(0);id++) {
      if (dof == dofs_gids_h(id)) {
        global_idx.push_back(idx);
        local_idx.push_back(id);
        break;
      }
    }
  }
  // Now that we have the full list of indexs in the global remap data that correspond to local columns we can construct
  // the spa_horiz_weights data.   Note: This is an important step when running with multiple MPI ranks.
  spa_horiz_interp.length          = local_idx.size();
  spa_horiz_interp.weights         = view_1d<Real>("",local_idx.size());
  spa_horiz_interp.source_grid_loc = view_1d<Int>("",local_idx.size());
  spa_horiz_interp.target_grid_loc = view_1d<Int>("",local_idx.size());
  auto weights_h         = Kokkos::create_mirror_view(spa_horiz_interp.weights);
  auto source_grid_loc_h = Kokkos::create_mirror_view(spa_horiz_interp.source_grid_loc);
  auto target_grid_loc_h = Kokkos::create_mirror_view(spa_horiz_interp.target_grid_loc);
  for (size_t idx=0;idx<local_idx.size();idx++) {
      int ii = global_idx[idx];
      weights_h(idx)         = S_global_h(ii);
      target_grid_loc_h(idx) = local_idx[idx];
      // Note that the remap column location starts with 1.
      // Here we want the index in a the source data vector corresponding to this column
      // which needs to start with 0 since cpp starts with 0.
      source_grid_loc_h(idx) = col_global_h(ii) - 1;
  }
  Kokkos::deep_copy(spa_horiz_interp.weights        , weights_h        );
  Kokkos::deep_copy(spa_horiz_interp.source_grid_loc, source_grid_loc_h);
  Kokkos::deep_copy(spa_horiz_interp.target_grid_loc, target_grid_loc_h);
  // Determine the set of unique columns in this remapping
  spa_horiz_interp.set_unique_cols();
}  // END get_remap_weights_from_file
/*-----------------------------------------------------------------*/
/* Note: In this routine the SPA source data is padded in the vertical
 * to facilitate the proper behavior at the boundaries when doing the
 * vertical interpolation in spa_main.
 *
 * The vertical interpolation routine will map the source data (y1)
 * defined at the pressure levels (x1) onto the target data vector (y2)
 * defined on the target pressure levels (x2) following the equation:
 *
 * y2 = y1(k) + ( y1(k+1)-y1(k) ) * ( x2-x1(k) )/( x1(k+1) - x1(k) )
 *
 * where k and k+1 are the pressure levels that bound the target pressure
 * level x2.
 *
 * We pad the left-hand-side (lhs) of the data with 0.0 and the right-hand-side (rhs)
 * of the data with a copy of the last data entry.  Similarly, we pad the
 * lhs of the hybrid coordinate system to ensure that the lhs of the
 * source pressure levels is 0.0 and the rhs is big enough to be larger than
 * the highest pressure in the target pressure levels.
 * 
 * The padding sets the lhs in a vertical column of data to 0.0.
 * This has the effect of ramping the top-of-model target data down to 0.0
 * when the top target pressure level is smaller than the top source pressure,
 * as follow:
 *      y2 = y1(1) * x2/x1(1)
 * Since x2 <= x1(1) this will slowly taper down to 0.0.
 *
 * The rhs value is set to match the actual rhs value of the data, such that
 *      y1[N+1] = y1(N)
 * where N is the number of points in the source data.
 *
 * Thus when the bottom target pressure level is larger than the bottom source
 * pressure level the interpolation is:
 *     y2 = y1[N]
 *
 * Thus, each column is padded as follows:
 *   src_data = [ 0.0, actual_spa_data, actual_spa_data[-1] ]
 *
 * Note, the hyam and hybm vectors are also padded so that the source data
 * pressure profile that is constructed in spa_main is
 * a) the current length, and
 * b) ensures that whatever the target pressure profile is, it's within the
 *    the range of the source data.                                        
 */
template<typename S, typename D>
void SPAFunctions<S,D>
::update_spa_data_from_file(
    const std::string&          spa_data_file_name,
    const Int                   time_index,
    const Int                   nswbands,
    const Int                   nlwbands,
          SPAHorizInterp&       spa_horiz_interp,
          SPAInput&             spa_data)
{
  // Ensure all ranks are operating independently when reading the file, so there's a copy on all ranks
  auto comm = spa_horiz_interp.m_comm;

  // We have enough info to start opening the file
  std::vector<std::string> fnames = {"hyam","hybm","PS","CCN3","AER_G_SW","AER_SSA_SW","AER_TAU_SW","AER_TAU_LW"};
  ekat::ParameterList spa_data_in_params;
  spa_data_in_params.set("Field Names",fnames);
  spa_data_in_params.set("Filename",spa_data_file_name);
  spa_data_in_params.set("Skip_Grid_Checks",true);  // We need to skip grid checks because multiple ranks may want the same column of source data.
  AtmosphereInput spa_data_input(comm,spa_data_in_params);

  // Note that the SPA data follows a conventional GLL grid format, albeit at a different resolution than
  // the simulation.  For simplicity we can use the scorpio_input object class but we must construct a
  // local grid to match the size of the SPA data file.

  // To construct the grid we need to determine the number of columns and levels in the data file.
  Int ncol = scorpio::get_dimlen_c2f(spa_data_file_name.c_str(),"ncol");
  const int source_data_nlevs = scorpio::get_dimlen_c2f(spa_data_file_name.c_str(),"lev");
  Int num_local_cols = spa_horiz_interp.num_unique_cols;
  // Check that padding matches source size:
  EKAT_REQUIRE(source_data_nlevs+2 == spa_data.data.nlevs);
  // while we have the file open, check that the dimensions map the simulation and the horizontal interpolation structure
  EKAT_REQUIRE_MSG(nswbands==scorpio::get_dimlen_c2f(spa_data_file_name.c_str(),"swband"),"ERROR update_spa_data_from_file: Number of SW bands in simulation doesn't match the SPA data file");
  EKAT_REQUIRE_MSG(nlwbands==scorpio::get_dimlen_c2f(spa_data_file_name.c_str(),"lwband"),"ERROR update_spa_data_from_file: Number of LW bands in simulation doesn't match the SPA data file");
  EKAT_REQUIRE_MSG(ncol==spa_horiz_interp.source_grid_ncols,"ERROR update_spa_data_from_file: Number of columns in remap data (" 
                        + std::to_string(ncol) + " doesn't match the SPA data file (" + std::to_string(spa_horiz_interp.source_grid_ncols) + ").");

  // Construct local arrays to read data into
  // Note, all of the views being created here are meant to hold the local "coarse" resolution
  // data that will need to be horizontally interpolated to the simulation grid using the remap
  // data.  For example, 
  //   We will first define the data surface pressure PS_v and read that from file.
  //   then we will use the horizontal interpolation structure, spa_horiz_interp, to
  //   interpolate PS_v onto the simulation grid: PS_v -> spa_data.PS
  //   and so on for the other variables.
  typename view_1d<Real>::HostMirror hyam_v_h("hyam",source_data_nlevs);
  typename view_1d<Real>::HostMirror hybm_v_h("hybm",source_data_nlevs);
  typename view_1d<Real>::HostMirror PS_v_h("PS",num_local_cols);
  typename view_2d<Real>::HostMirror CCN3_v_h("CCN3",num_local_cols,source_data_nlevs);
  typename view_3d<Real>::HostMirror AER_G_SW_v_h("AER_G_SW",num_local_cols,nswbands,source_data_nlevs);
  typename view_3d<Real>::HostMirror AER_SSA_SW_v_h("AER_SSA_SW",num_local_cols,nswbands,source_data_nlevs);
  typename view_3d<Real>::HostMirror AER_TAU_SW_v_h("AER_TAU_SW",num_local_cols,nswbands,source_data_nlevs);
  typename view_3d<Real>::HostMirror AER_TAU_LW_v_h("AER_TAU_LW",num_local_cols,nlwbands,source_data_nlevs);

  // Construct the grid needed for input:
  auto grid = std::make_shared<PointGrid>("grid",num_local_cols,source_data_nlevs,comm);
  PointGrid::dofs_list_type dof_gids("",num_local_cols);
  auto dof_gids_h = Kokkos::create_mirror_view(dof_gids);
  for (const auto& nn : spa_horiz_interp.source_local_col_map) {
    dof_gids_h(nn.second) = nn.first;
  }
  Kokkos::deep_copy(dof_gids,dof_gids_h);
  grid->set_dofs(dof_gids);

  using namespace ShortFieldTagsNames;
  FieldLayout scalar1d_layout { {LEV}, {source_data_nlevs} };
  FieldLayout scalar2d_layout_mid { {COL}, {num_local_cols} };
  FieldLayout scalar3d_layout_mid { {COL,LEV}, {num_local_cols, source_data_nlevs} };
  FieldLayout scalar3d_swband_layout { {COL,SWBND, LEV}, {num_local_cols, nswbands, source_data_nlevs} }; 
  FieldLayout scalar3d_lwband_layout { {COL,LWBND, LEV}, {num_local_cols, nlwbands, source_data_nlevs} };
  std::map<std::string,view_1d_host<Real>> host_views;
  std::map<std::string,FieldLayout>  layouts;
  // Define each input variable we need
  host_views["hyam"] = view_1d_host<Real>(hyam_v_h.data(),hyam_v_h.size());
  layouts.emplace("hyam", scalar1d_layout);
  host_views["hybm"] = view_1d_host<Real>(hybm_v_h.data(),hybm_v_h.size());
  layouts.emplace("hybm", scalar1d_layout);
  //
  host_views["PS"] = view_1d_host<Real>(PS_v_h.data(),PS_v_h.size());
  layouts.emplace("PS", scalar2d_layout_mid);
  //
  host_views["CCN3"] = view_1d_host<Real>(CCN3_v_h.data(),CCN3_v_h.size());
  layouts.emplace("CCN3",scalar3d_layout_mid);
  //
  host_views["AER_G_SW"] = view_1d_host<Real>(AER_G_SW_v_h.data(),AER_G_SW_v_h.size());
  layouts.emplace("AER_G_SW",scalar3d_swband_layout);
  //
  host_views["AER_SSA_SW"] = view_1d_host<Real>(AER_SSA_SW_v_h.data(),AER_SSA_SW_v_h.size());
  layouts.emplace("AER_SSA_SW",scalar3d_swband_layout);
  //
  host_views["AER_TAU_SW"] = view_1d_host<Real>(AER_TAU_SW_v_h.data(),AER_TAU_SW_v_h.size());
  layouts.emplace("AER_TAU_SW",scalar3d_swband_layout);
  //
  host_views["AER_TAU_LW"] = view_1d_host<Real>(AER_TAU_LW_v_h.data(),AER_TAU_LW_v_h.size());
  layouts.emplace("AER_TAU_LW",scalar3d_lwband_layout);
  //
  
  // Now that we have all the variables defined we can use the scorpio_input class to grab the data.
  spa_data_input.init(grid,host_views,layouts);
  spa_data_input.read_variables(time_index);
  spa_data_input.finalize();
 
  // Now that we have the data we can map the data onto the target data.
  auto hyam_h       = Kokkos::create_mirror_view(spa_data.hyam);
  auto hybm_h       = Kokkos::create_mirror_view(spa_data.hybm);
  auto ps_h         = Kokkos::create_mirror_view(spa_data.PS);
  auto ccn3_h       = Kokkos::create_mirror_view(spa_data.data.CCN3);
  auto aer_g_sw_h   = Kokkos::create_mirror_view(spa_data.data.AER_G_SW);
  auto aer_ssa_sw_h = Kokkos::create_mirror_view(spa_data.data.AER_SSA_SW);
  auto aer_tau_sw_h = Kokkos::create_mirror_view(spa_data.data.AER_TAU_SW);
  auto aer_tau_lw_h = Kokkos::create_mirror_view(spa_data.data.AER_TAU_LW);
  Kokkos::deep_copy(hyam_h,0.0);
  Kokkos::deep_copy(hybm_h,0.0);
  Kokkos::deep_copy(ps_h,0.0);
  Kokkos::deep_copy(ccn3_h,0.0);
  Kokkos::deep_copy(aer_g_sw_h,0.0);
  Kokkos::deep_copy(aer_ssa_sw_h,0.0);
  Kokkos::deep_copy(aer_tau_sw_h,0.0);
  Kokkos::deep_copy(aer_tau_lw_h,0.0);

  auto weights_h         = Kokkos::create_mirror_view(spa_horiz_interp.weights);
  auto source_grid_loc_h = Kokkos::create_mirror_view(spa_horiz_interp.source_grid_loc);
  auto target_grid_loc_h = Kokkos::create_mirror_view(spa_horiz_interp.target_grid_loc);
  Kokkos::deep_copy(weights_h,         spa_horiz_interp.weights);
  Kokkos::deep_copy(source_grid_loc_h, spa_horiz_interp.source_grid_loc);
  Kokkos::deep_copy(target_grid_loc_h, spa_horiz_interp.target_grid_loc);
  for (int idx=0;idx<spa_horiz_interp.length;idx++) {
    auto src_wgt = weights_h(idx);
    int  src_col = spa_horiz_interp.source_local_col_map[source_grid_loc_h(idx)];
    int  tgt_col = target_grid_loc_h(idx);
    // PS is defined only over columns
    ps_h(tgt_col) += PS_v_h(src_col)*src_wgt;
    // CCN3 and all AER variables have levels
    for (int kk=0; kk<source_data_nlevs; kk++) {
      // Note, all variables we map to are packed, while all the data we just loaded as
      // input are in real N-D views.  So we need to set the pack and index of the actual
      // data ahead by one value.
      // Note, we want to pad the actual source data such that
      //   Y[0]   = 0.0, note this is handled by the deep copy above
      //   Y[k+1] = y[k], k = 0,source_data_nlevs (y is the data from file)
      //   Y[N+2] = y[N-1], N = source_data_nlevs
      int pack = (kk+1) / Spack::n; 
      int kidx = (kk+1) % Spack::n;
      ccn3_h(tgt_col,pack)[kidx] += CCN3_v_h(src_col,kk)*src_wgt;
      for (int n=0; n<nswbands; n++) {
        aer_g_sw_h(tgt_col,n,pack)[kidx]   += AER_G_SW_v_h(src_col,n,kk)*src_wgt;
        aer_ssa_sw_h(tgt_col,n,pack)[kidx] += AER_SSA_SW_v_h(src_col,n,kk)*src_wgt;
        aer_tau_sw_h(tgt_col,n,pack)[kidx] += AER_TAU_SW_v_h(src_col,n,kk)*src_wgt;
      }
      for (int n=0; n<nlwbands; n++) {
        aer_tau_lw_h(tgt_col,n,pack)[kidx] += AER_TAU_LW_v_h(src_col,n,kk)*src_wgt;
      }
    }
    int kk = source_data_nlevs-1;
    int pack = (kk+2) / Spack::n; 
    int kidx = (kk+2) % Spack::n;
    ccn3_h(tgt_col,pack)[kidx] += CCN3_v_h(src_col,kk)*src_wgt;
    for (int n=0; n<nswbands; n++) {
      aer_g_sw_h(tgt_col,n,pack)[kidx]   += AER_G_SW_v_h(src_col,n,kk)*src_wgt;
      aer_ssa_sw_h(tgt_col,n,pack)[kidx] += AER_SSA_SW_v_h(src_col,n,kk)*src_wgt;
      aer_tau_sw_h(tgt_col,n,pack)[kidx] += AER_TAU_SW_v_h(src_col,n,kk)*src_wgt;
    }
    for (int n=0; n<nlwbands; n++) {
      aer_tau_lw_h(tgt_col,n,pack)[kidx] += AER_TAU_LW_v_h(src_col,n,kk)*src_wgt;
    }
  }
  // We also need to pad the hyam and hybm views with
  //   hya/b[0] = 0.0, note this is handled by deep copy above
  //   hya/b[N+2] = BIG number so always bigger than likely pmid for target
  for (int kk=0; kk<source_data_nlevs; kk++) {
    int pack = (kk+1) / Spack::n; 
    int kidx = (kk+1) % Spack::n;
    hyam_h(pack)[kidx] = hyam_v_h(kk);
    hybm_h(pack)[kidx] = hybm_v_h(kk);
  }
  const int pack = (source_data_nlevs+1) / Spack::n;
  const int kidx = (source_data_nlevs+1) % Spack::n;
  hyam_h(pack)[kidx] = 1e5; 
  hybm_h(pack)[kidx] = 0.0;
  
  Kokkos::deep_copy(spa_data.hyam,hyam_h);
  Kokkos::deep_copy(spa_data.hybm,hybm_h);
  Kokkos::deep_copy(spa_data.PS,ps_h);
  Kokkos::deep_copy(spa_data.data.CCN3,ccn3_h);
  Kokkos::deep_copy(spa_data.data.AER_G_SW,aer_g_sw_h);
  Kokkos::deep_copy(spa_data.data.AER_SSA_SW,aer_ssa_sw_h);
  Kokkos::deep_copy(spa_data.data.AER_TAU_SW,aer_tau_sw_h);
  Kokkos::deep_copy(spa_data.data.AER_TAU_LW,aer_tau_lw_h);

} // END update_spa_data_from_file

/*-----------------------------------------------------------------*/
template<typename S, typename D>
void SPAFunctions<S,D>
::update_spa_timestate(
  const std::string&     spa_data_file_name,
  const Int              nswbands,
  const Int              nlwbands,
  const util::TimeStamp& ts,
        SPAHorizInterp&  spa_horiz_interp,
        SPATimeState&    time_state, 
        SPAInput&        spa_beg,
        SPAInput&        spa_end)
{

  // Now we check if we have to update the data that changes monthly
  // NOTE:  This means that SPA assumes monthly data to update.  Not
  //        any other frequency.
  const auto month = ts.get_month();
  if (month != time_state.current_month or !time_state.inited) {
    // Update the SPA time state information
    time_state.current_month = month;
    time_state.t_beg_month = util::TimeStamp({0,month,1}, {0,0,0}).frac_of_year_in_days();
    time_state.days_this_month = util::days_in_month(ts.get_year(),month);
    // Update the SPA forcing data for this month and next month
    // Start by copying next months data to this months data structure.  
    // NOTE: If the timestep is bigger than monthly this could cause the wrong values
    //       to be assigned.  A timestep greater than a month is very unlikely so we
    //       will proceed.
    update_spa_data_from_file(spa_data_file_name,time_state.current_month,nswbands,nlwbands,spa_horiz_interp,spa_beg);
    Int next_month = time_state.current_month==12 ? 1 : time_state.current_month+1;
    update_spa_data_from_file(spa_data_file_name,next_month,nswbands,nlwbands,spa_horiz_interp,spa_end);
    // If time state was not initialized it is now:
    time_state.inited = true;
  }

} // END updata_spa_timestate

template<typename S,typename D>
KOKKOS_INLINE_FUNCTION
auto SPAFunctions<S,D>::
get_var_column (const SPAData& data, const int icol, const int ivar)
  -> view_1d<Spack> 
{
  if (ivar==0) {
    return ekat::subview(data.CCN3,icol);
  } else {
    // NOTE: if we ever get 2+ long-wave vars, you will have to do 
    //       something like
    //    if (jvar<num_sw_vars) {..} else { // compute lw var index }
    int jvar   = (ivar-1) / data.nswbands;
    int swband = (ivar-1) % data.nswbands;
    int lwband = (ivar-1) - 3*data.nswbands;
    switch (jvar) {
      case 0: return ekat::subview(data.AER_G_SW,icol,swband);
      case 1: return ekat::subview(data.AER_SSA_SW,icol,swband); 
      case 2: return ekat::subview(data.AER_TAU_SW,icol,swband); 
      default: return ekat::subview(data.AER_TAU_LW,icol,lwband); 

    }
  }
}

template<typename S,typename D>
template<typename ScalarX,typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarX SPAFunctions<S,D>::
linear_interp(const ScalarX& x0, const ScalarX& x1, const ScalarT& t)
{
  return (1 - t)*x0 + t*x1;
}
/*-----------------------------------------------------------------*/

} // namespace spa
} // namespace scream

#endif // SPA_FUNCTIONS_IMPL_HPP
