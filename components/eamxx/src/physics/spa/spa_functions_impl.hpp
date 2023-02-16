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

#include "share/util/scream_timing.hpp"
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

  EKAT_REQUIRE_MSG (delta_t_fraction>=0 && delta_t_fraction<=1,
      "Error! Convex interpolation with coefficient out of [0,1].\n"
      "  t_now  : " + std::to_string(t_now) + "\n"
      "  t_beg  : " + std::to_string(t_beg) + "\n"
      "  delta_t: " + std::to_string(delta_t) + "\n");

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

    Kokkos::parallel_for (Kokkos::TeamVectorRange(team,num_vert_packs),
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
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team,num_vert_packs),
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
// Function to set the remap and weights for a one-to-one mapping.
// This is used when the SPA data and the simulation grid are the
// same and no remapping is needed.
// Note: This function should be called only once during SPA::init
template <typename S, typename D>
void SPAFunctions<S,D>
::set_remap_weights_one_to_one(
    gid_type                       min_dof,
    const view_1d<const gid_type>& dofs_gids,
          SPAHorizInterp&          spa_horiz_interp
  )
{
  // There may be cases where the SPA data is defined on the same grid as the simulation
  // and thus no remapping is required.  This simple routine establishes a 1-1 horizontal
  // mapping
  const int num_local_cols = dofs_gids.size();
  auto& spa_horiz_map = spa_horiz_interp.horiz_map;
  spa_horiz_map = HorizontalMap(spa_horiz_interp.m_comm,"SPA 1-1 Remap",dofs_gids,min_dof);
  view_1d<gid_type> src_dofs("",1);
  view_1d<Real>     src_wgts("",1);
  auto dofs_gids_h = Kokkos::create_mirror_view(dofs_gids);
  Kokkos::deep_copy(dofs_gids_h,dofs_gids);
  for (int ii=0;ii<num_local_cols;ii++) {
    gid_type seg_dof = dofs_gids_h(ii)-min_dof;
    Kokkos::deep_copy(src_dofs,seg_dof);
    Kokkos::deep_copy(src_wgts,   1.0);
    HorizontalMapSegment seg(seg_dof,1,src_dofs,src_wgts);
    spa_horiz_map.add_remap_segment(seg);
  }
  spa_horiz_map.set_unique_source_dofs();

} // END set_remap_weights_one_to_one
/*-----------------------------------------------------------------*/
// Function to gather the remap and weights from a specific file.
// This is used when the SPA data is not on the same grid as the
// simulation, and so some remapping is needed.  In this case the
// remapping is provided by file.
// Note: This function should only be called once during SPA::init
template <typename S, typename D>
void SPAFunctions<S,D>
::get_remap_weights_from_file(
    const std::string&             remap_file_name,
    const gid_type                 min_dof,
    const view_1d<const gid_type>& dofs_gids,
          SPAHorizInterp&          spa_horiz_interp
  )
{
  start_timer("EAMxx::SPA::get_remap_weights_from_file");
  auto& spa_horiz_map = spa_horiz_interp.horiz_map;
  spa_horiz_map = HorizontalMap(spa_horiz_interp.m_comm,"SPA File Remap", dofs_gids, min_dof);
  spa_horiz_map.set_remap_segments_from_file(remap_file_name);
  spa_horiz_map.set_unique_source_dofs();
  stop_timer("EAMxx::SPA::get_remap_weights_from_file");

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
    const int                   time_index, // zero-based
    const int                   nswbands,
    const int                   nlwbands,
          SPAHorizInterp&       spa_horiz_interp,
          SPAInput&             spa_data)
{
  start_timer("EAMxx::SPA::update_spa_data_from_file");
  // Ensure all ranks are operating independently when reading the file, so there's a copy on all ranks
  auto comm = spa_horiz_interp.m_comm;

  // Use HorizontalMap to define the set of source column data we need to load
  auto& spa_horiz_map = spa_horiz_interp.horiz_map;
  auto unique_src_dofs = spa_horiz_map.get_unique_source_dofs();
  const int num_local_cols = spa_horiz_map.get_num_unique_dofs();
  scorpio::register_file(spa_data_file_name,scorpio::Read);
  const int source_data_nlevs = scorpio::get_dimlen_c2f(spa_data_file_name.c_str(),"lev");
  scorpio::eam_pio_closefile(spa_data_file_name);

  // Construct local arrays to read data into
  // Note, all of the views being created here are meant to hold the source resolution
  // data that will need to be horizontally interpolated to the simulation grid using the remap
  // data.  For example, 
  //   We will first define the data surface pressure PS_v and read that from file.
  //   then we will use the horizontal interpolation structure, spa_horiz_interp, to
  //   interpolate PS_v onto the simulation grid: PS_v -> spa_data.PS
  //   and so on for the other variables.
  start_timer("EAMxx::SPA::update_spa_data_from_file::read_data");
  std::vector<std::string> fnames = {"hyam","hybm","PS","CCN3","AER_G_SW","AER_SSA_SW","AER_TAU_SW","AER_TAU_LW"};
  ekat::ParameterList spa_data_in_params;
  spa_data_in_params.set("Field Names",fnames);
  spa_data_in_params.set("Filename",spa_data_file_name);
  spa_data_in_params.set("Skip_Grid_Checks",true);  // We need to skip grid checks because multiple ranks may want the same column of source data.
  AtmosphereInput spa_data_input(comm,spa_data_in_params);
  
  // Construct the grid needed for input:
  auto grid = std::make_shared<PointGrid>("grid",num_local_cols,source_data_nlevs,comm);
  Kokkos::deep_copy(grid->get_dofs_gids().template get_view<gid_type*>(),unique_src_dofs);
  grid->get_dofs_gids().sync_to_host();

  // Check that padding matches source size:
  EKAT_REQUIRE(source_data_nlevs+2 == spa_data.data.nlevs);
  EKAT_REQUIRE_MSG(nswbands==scorpio::get_dimlen_c2f(spa_data_file_name.c_str(),"swband"),"ERROR update_spa_data_from_file: Number of SW bands in simulation doesn't match the SPA data file");
  EKAT_REQUIRE_MSG(nlwbands==scorpio::get_dimlen_c2f(spa_data_file_name.c_str(),"lwband"),"ERROR update_spa_data_from_file: Number of LW bands in simulation doesn't match the SPA data file");

  // Constuct views to read source data in from file
  typename view_1d<Real>::HostMirror hyam_v_h("hyam",source_data_nlevs);
  typename view_1d<Real>::HostMirror hybm_v_h("hybm",source_data_nlevs);
  view_1d<Real> PS_v("PS",num_local_cols);
  view_2d<Real> CCN3_v("CCN3",num_local_cols,source_data_nlevs);
  view_3d<Real> AER_G_SW_v("AER_G_SW",num_local_cols,nswbands,source_data_nlevs);
  view_3d<Real> AER_SSA_SW_v("AER_SSA_SW",num_local_cols,nswbands,source_data_nlevs);
  view_3d<Real> AER_TAU_SW_v("AER_TAU_SW",num_local_cols,nswbands,source_data_nlevs);
  view_3d<Real> AER_TAU_LW_v("AER_TAU_LW",num_local_cols,nlwbands,source_data_nlevs);

  auto PS_v_h         = Kokkos::create_mirror_view(PS_v);
  auto CCN3_v_h       = Kokkos::create_mirror_view(CCN3_v);      
  auto AER_G_SW_v_h   = Kokkos::create_mirror_view(AER_G_SW_v);
  auto AER_SSA_SW_v_h = Kokkos::create_mirror_view(AER_SSA_SW_v);
  auto AER_TAU_SW_v_h = Kokkos::create_mirror_view(AER_TAU_SW_v);
  auto AER_TAU_LW_v_h = Kokkos::create_mirror_view(AER_TAU_LW_v);

  // Set up input structure to read data from file.
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
  stop_timer("EAMxx::SPA::update_spa_data_from_file::read_data");
  start_timer("EAMxx::SPA::update_spa_data_from_file::apply_remap");
  // Copy data from host back to the device views.
  Kokkos::deep_copy(PS_v,PS_v_h);
  Kokkos::deep_copy(CCN3_v      , CCN3_v_h);      
  Kokkos::deep_copy(AER_G_SW_v  , AER_G_SW_v_h);
  Kokkos::deep_copy(AER_SSA_SW_v, AER_SSA_SW_v_h);
  Kokkos::deep_copy(AER_TAU_SW_v, AER_TAU_SW_v_h);
  Kokkos::deep_copy(AER_TAU_LW_v, AER_TAU_LW_v_h);

  // Apply the remap to this data
  spa_horiz_map.apply_remap(PS_v,spa_data.PS); // Note PS is not padded, so remap can be applied right away
  // For padded data we need create temporary arrays to store the direct remapped data, then we can add
  // padding.
  int tgt_ncol = spa_data.data.ncols;
  int tgt_nlev = spa_data.data.nlevs-2;  // Note, the spa data already accounts for padding in the nlevs, so we subtract 2
  view_2d<Real> CCN3_unpad("",tgt_ncol,tgt_nlev);
  view_3d<Real> AER_G_SW_unpad("",tgt_ncol,nswbands,tgt_nlev);
  view_3d<Real> AER_SSA_SW_unpad("",tgt_ncol,nswbands,tgt_nlev);
  view_3d<Real> AER_TAU_SW_unpad("",tgt_ncol,nswbands,tgt_nlev);
  view_3d<Real> AER_TAU_LW_unpad("",tgt_ncol,nlwbands,tgt_nlev);
  // Apply remap to "unpadded" data
  spa_horiz_map.apply_remap(CCN3_v,CCN3_unpad);
  spa_horiz_map.apply_remap(AER_G_SW_v, AER_G_SW_unpad);
  spa_horiz_map.apply_remap(AER_SSA_SW_v, AER_SSA_SW_unpad);
  spa_horiz_map.apply_remap(AER_TAU_SW_v, AER_TAU_SW_unpad);
  spa_horiz_map.apply_remap(AER_TAU_LW_v, AER_TAU_LW_unpad);
  stop_timer("EAMxx::SPA::update_spa_data_from_file::apply_remap");
  start_timer("EAMxx::SPA::update_spa_data_from_file::copy_and_pad");
  // Copy unpadded data to SPA data structure, add padding.
  // Note, all variables we map to are packed, while all the data we just loaded as
  // input are in real N-D views.  So we need to set the pack and index of the actual
  // data ahead by one value.
  // Note, we want to pad the actual source data such that
  //   Y[0]   = 0.0, note this is handled by the deep copy above
  //   Y[k+1] = y[k], k = 0,source_data_nlevs (y is the data from file)
  //   Y[N+2] = y[N-1], N = source_data_nlevs
  Kokkos::deep_copy(spa_data.data.CCN3,0.0);
  Kokkos::deep_copy(spa_data.data.AER_G_SW,0.0);
  Kokkos::deep_copy(spa_data.data.AER_SSA_SW,0.0);
  Kokkos::deep_copy(spa_data.data.AER_TAU_SW,0.0);
  Kokkos::deep_copy(spa_data.data.AER_TAU_LW,0.0);
  // 2D vars - CCN3
  Kokkos::parallel_for("", tgt_ncol*tgt_nlev, KOKKOS_LAMBDA (const int& idx) {
    int icol  = idx / tgt_nlev;
    int klev1 = idx % tgt_nlev;
    int kpack = (klev1+1) / Spack::n;
    int klev2 = (klev1+1) % Spack::n;
    spa_data.data.CCN3(icol,kpack)[klev2] = CCN3_unpad(icol,klev1);
    if (klev1 == tgt_nlev-1) {
      int kpack = (tgt_nlev+1) / Spack::n;
      int klev2 = (tgt_nlev+1) % Spack::n;
      spa_data.data.CCN3(icol,kpack)[klev2] = CCN3_unpad(icol,klev1);
    }
  });
  Kokkos::fence();
  // 3D vars - AER_G_SW, AER_SSA_SW, AER_TAU_SW, AER_TAU_LW
  Kokkos::parallel_for("", tgt_ncol*tgt_nlev*nswbands, KOKKOS_LAMBDA (const int& idx) {
    int icol  = idx / (tgt_nlev*nswbands);
    int nband = (idx-icol*tgt_nlev*nswbands) / tgt_nlev;
    int klev1 = idx % tgt_nlev;
    int kpack = (klev1+1) / Spack::n;
    int klev2 = (klev1+1) % Spack::n;
    spa_data.data.AER_G_SW(icol,nband,kpack)[klev2] = AER_G_SW_unpad(icol,nband,klev1);
    spa_data.data.AER_SSA_SW(icol,nband,kpack)[klev2] = AER_SSA_SW_unpad(icol,nband,klev1);
    spa_data.data.AER_TAU_SW(icol,nband,kpack)[klev2] = AER_TAU_SW_unpad(icol,nband,klev1);
    if (klev1 == tgt_nlev-1) {
      int kpack = (tgt_nlev+1) / Spack::n;
      int klev2 = (tgt_nlev+1) % Spack::n;
      spa_data.data.AER_G_SW(icol,nband,kpack)[klev2] = AER_G_SW_unpad(icol,nband,klev1);
      spa_data.data.AER_SSA_SW(icol,nband,kpack)[klev2] = AER_SSA_SW_unpad(icol,nband,klev1);
      spa_data.data.AER_TAU_SW(icol,nband,kpack)[klev2] = AER_TAU_SW_unpad(icol,nband,klev1);
    }
  });
  Kokkos::parallel_for("", tgt_ncol*tgt_nlev*nlwbands, KOKKOS_LAMBDA (const int& idx) {
    int icol  = idx / (tgt_nlev*nlwbands);
    int nband = (idx-icol*tgt_nlev*nlwbands) / tgt_nlev;
    int klev1 = idx % tgt_nlev;
    int kpack = (klev1+1) / Spack::n;
    int klev2 = (klev1+1) % Spack::n;
    spa_data.data.AER_TAU_LW(icol,nband,kpack)[klev2] = AER_TAU_LW_unpad(icol,nband,klev1);
    if (klev1 == tgt_nlev-1) {
      int kpack = (tgt_nlev+1) / Spack::n;
      int klev2 = (tgt_nlev+1) % Spack::n;
      spa_data.data.AER_TAU_LW(icol,nband,kpack)[klev2] = AER_TAU_LW_unpad(icol,nband,klev1);
    }
  });
  stop_timer("EAMxx::SPA::update_spa_data_from_file::copy_and_pad");
  // The hybrid coordinates are just vertical data, so not remapped.  We still need to make
  // a padded version of the data.
  //   hya/b[0] = 0.0, note this is handled by deep copy above
  //   hya/b[N+2] = BIG number so always bigger than likely pmid for target
  auto hyam_h       = Kokkos::create_mirror_view(spa_data.hyam);
  auto hybm_h       = Kokkos::create_mirror_view(spa_data.hybm);
  Kokkos::deep_copy(hyam_h,0.0);
  Kokkos::deep_copy(hybm_h,0.0);
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
  stop_timer("EAMxx::SPA::update_spa_data_from_file");

} // END update_spa_data_from_file

/*-----------------------------------------------------------------*/
template<typename S, typename D>
void SPAFunctions<S,D>
::update_spa_timestate(
  const std::string&     spa_data_file_name,
  const int              nswbands,
  const int              nlwbands,
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
    time_state.t_beg_month = util::TimeStamp({ts.get_year(),month,1}, {0,0,0}).frac_of_year_in_days();
    time_state.days_this_month = util::days_in_month(ts.get_year(),month);
    // Update the SPA forcing data for this month and next month
    // Start by copying next months data to this months data structure.  
    // NOTE: If the timestep is bigger than monthly this could cause the wrong values
    //       to be assigned.  A timestep greater than a month is very unlikely so we
    //       will proceed.
    // NOTE: we use zero-based time indexing here.
    update_spa_data_from_file(spa_data_file_name,time_state.current_month-1,nswbands,nlwbands,spa_horiz_interp,spa_beg);
    int next_month = time_state.current_month==12 ? 1 : time_state.current_month+1;
    update_spa_data_from_file(spa_data_file_name,next_month-1,nswbands,nlwbands,spa_horiz_interp,spa_end);
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
