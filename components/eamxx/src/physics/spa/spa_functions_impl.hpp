#ifndef SPA_FUNCTIONS_IMPL_HPP
#define SPA_FUNCTIONS_IMPL_HPP

#include "physics/share/physics_constants.hpp"
#include "share/grid/remap/coarsening_remapper.hpp"
#include "share/grid/remap/refining_remapper_p2p.hpp"
#include "share/grid/remap/identity_remapper.hpp"
#include "share/io/scream_scorpio_interface.hpp"
#include "share/util/scream_timing.hpp"
#include "share/scream_types.hpp"

#include <ekat/kokkos/ekat_subview_utils.hpp>
#include <ekat/kokkos/ekat_kokkos_utils.hpp>
#include <ekat/util/ekat_lin_interp.hpp>
#include <ekat/ekat_pack_utils.hpp>
#include <ekat/ekat_pack_kokkos.hpp>

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

template <typename S, typename D>
std::shared_ptr<AbstractRemapper>
SPAFunctions<S,D>::
create_horiz_remapper (
    const std::shared_ptr<const AbstractGrid>& model_grid,
    const std::string& spa_data_file,
    const std::string& map_file,
    const bool use_iop)
{
  using namespace ShortFieldTagsNames;

  scorpio::register_file(spa_data_file,scorpio::Read,0);
  const int nlevs_data = scorpio::get_dimlen(spa_data_file,"lev");
  const int ncols_data = scorpio::get_dimlen(spa_data_file,"ncol");
  const int nswbands   = scorpio::get_dimlen(spa_data_file,"swband");
  const int nlwbands   = scorpio::get_dimlen(spa_data_file,"lwband");
  scorpio::eam_pio_closefile(spa_data_file);

  // We could use model_grid directly if using same num levels,
  // but since shallow clones are cheap, we may as well do it (less lines of code)
  auto horiz_interp_tgt_grid = model_grid->clone("spa_horiz_interp_tgt_grid",true);
  horiz_interp_tgt_grid->reset_num_vertical_lev(nlevs_data);

  const int ncols_model = model_grid->get_num_global_dofs();
  std::shared_ptr<AbstractRemapper> remapper;
  if (ncols_data==ncols_model
      or
      use_iop /*IOP class defines it's own remapper for file data*/) {
    remapper = std::make_shared<IdentityRemapper>(horiz_interp_tgt_grid,IdentityRemapper::SrcAliasTgt);
  } else {
    EKAT_REQUIRE_MSG (ncols_data<=ncols_model,
      "Error! We do not allow to coarsen spa data to fit the model. We only allow\n"
      "       spa data to be at the same or coarser resolution as the model.\n");
    // We must have a valid map file
    EKAT_REQUIRE_MSG (map_file!="",
        "ERROR: Spa data is on a different grid than the model one,\n"
        "       but spa_remap_file is missing from SPA parameter list.");

    remapper = std::make_shared<RefiningRemapperP2P>(horiz_interp_tgt_grid,map_file);
  }

  remapper->registration_begins();

  const auto tgt_grid = remapper->get_tgt_grid();

  const auto layout_2d   = tgt_grid->get_2d_scalar_layout();
  const auto layout_ccn3 = tgt_grid->get_3d_scalar_layout(true);
  const auto layout_sw   = tgt_grid->get_3d_vector_layout(true,SWBND,nswbands);
  const auto layout_lw   = tgt_grid->get_3d_vector_layout(true,LWBND,nlwbands);
  const auto nondim = ekat::units::Units::nondimensional();

  Field ps          (FieldIdentifier("PS",        layout_2d,  nondim,tgt_grid->name()));
  Field ccn3        (FieldIdentifier("CCN3",      layout_ccn3,nondim,tgt_grid->name()));
  Field aero_g_sw   (FieldIdentifier("AER_G_SW",  layout_sw,  nondim,tgt_grid->name()));
  Field aero_ssa_sw (FieldIdentifier("AER_SSA_SW",layout_sw,  nondim,tgt_grid->name()));
  Field aero_tau_sw (FieldIdentifier("AER_TAU_SW",layout_sw,  nondim,tgt_grid->name()));
  Field aero_tau_lw (FieldIdentifier("AER_TAU_LW",layout_lw,  nondim,tgt_grid->name()));
  ps.allocate_view();
  ccn3.allocate_view();
  aero_g_sw.allocate_view();
  aero_ssa_sw.allocate_view();
  aero_tau_sw.allocate_view();
  aero_tau_lw.allocate_view();

  remapper->register_field_from_tgt (ps);
  remapper->register_field_from_tgt (ccn3);
  remapper->register_field_from_tgt (aero_g_sw);
  remapper->register_field_from_tgt (aero_ssa_sw);
  remapper->register_field_from_tgt (aero_tau_sw);
  remapper->register_field_from_tgt (aero_tau_lw);

  remapper->registration_ends();

  return remapper;
}

template <typename S, typename D>
std::shared_ptr<AtmosphereInput>
SPAFunctions<S,D>::
create_spa_data_reader (
    const std::shared_ptr<AbstractRemapper>& horiz_remapper,
    const std::string& spa_data_file)
{
  std::vector<Field> io_fields;
  for (int i=0; i<horiz_remapper->get_num_fields(); ++i) {
    io_fields.push_back(horiz_remapper->get_src_field(i));
  }
  const auto io_grid = horiz_remapper->get_src_grid();
  return std::make_shared<AtmosphereInput>(spa_data_file,io_grid,io_fields,true);
}

template <typename S, typename D>
std::shared_ptr<typename SPAFunctions<S,D>::IOPReader>
SPAFunctions<S,D>::
create_spa_data_reader (
    iop_ptr_type& iop,
    const std::shared_ptr<AbstractRemapper>& horiz_remapper,
    const std::string& spa_data_file)
{
  std::vector<Field> io_fields;
  for (int i=0; i<horiz_remapper->get_num_fields(); ++i) {
    io_fields.push_back(horiz_remapper->get_src_field(i));
  }
  const auto io_grid = horiz_remapper->get_src_grid();
  return std::make_shared<IOPReader>(iop, spa_data_file, io_fields, io_grid);
}

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
    std::shared_ptr<AtmosphereInput>& scorpio_reader,
    std::shared_ptr<IOPReader>&       iop_reader,
    const util::TimeStamp&            ts,
    const int                         time_index, // zero-based
    AbstractRemapper&                 spa_horiz_interp,
    SPAInput&                         spa_input)
{
  using namespace ShortFieldTagsNames;
  using ESU = ekat::ExeSpaceUtils<typename DefaultDevice::execution_space>;
  using Member = typename KokkosTypes<DefaultDevice>::MemberType;

  start_timer("EAMxx::SPA::update_spa_data_from_file");

  // 1. Read from file
  start_timer("EAMxx::SPA::update_spa_data_from_file::read_data");
  if (iop_reader) {
    iop_reader->read_variables(time_index, ts);
  } else {
    scorpio_reader->read_variables(time_index);
  }
  stop_timer("EAMxx::SPA::update_spa_data_from_file::read_data");

  // 2. Run the horiz remapper (it is a do-nothing op if spa data is on same grid as model)
  start_timer("EAMxx::SPA::update_spa_data_from_file::horiz_remap");
  spa_horiz_interp.remap(/*forward = */ true);
  stop_timer("EAMxx::SPA::update_spa_data_from_file::horiz_remap");

  // 3. Copy from the tgt field of the remapper into the spa_data, padding data if necessary
  start_timer("EAMxx::SPA::update_spa_data_from_file::copy_and_pad");
  // Recall, the fields are registered in the order: ps, ccn3, g_sw, ssa_sw, tau_sw, tau_lw
  auto ps          = spa_horiz_interp.get_tgt_field (0).get_view<const Real*>();
  auto ccn3        = spa_horiz_interp.get_tgt_field (1).get_view<const Real**>();
  auto aero_g_sw   = spa_horiz_interp.get_tgt_field (2).get_view<const Real***>();
  auto aero_ssa_sw = spa_horiz_interp.get_tgt_field (3).get_view<const Real***>();
  auto aero_tau_sw = spa_horiz_interp.get_tgt_field (4).get_view<const Real***>();
  auto aero_tau_lw = spa_horiz_interp.get_tgt_field (5).get_view<const Real***>();

  const auto& sw_layout = spa_horiz_interp.get_tgt_field (2).get_header().get_identifier().get_layout();
  const auto& lw_layout = spa_horiz_interp.get_tgt_field (5).get_header().get_identifier().get_layout();

  const int ncols    = sw_layout.dim(COL);
  const int nlevs    = sw_layout.dim(LEV);
  const int nswbands = sw_layout.dim(SWBND);
  const int nlwbands = lw_layout.dim(LWBND);

  Kokkos::deep_copy(spa_input.PS,ps);
  Kokkos::fence();
  auto spa_data_ccn3        = ekat::scalarize(spa_input.data.CCN3);
  auto spa_data_aero_g_sw   = ekat::scalarize(spa_input.data.AER_G_SW);
  auto spa_data_aero_ssa_sw = ekat::scalarize(spa_input.data.AER_SSA_SW);
  auto spa_data_aero_tau_sw = ekat::scalarize(spa_input.data.AER_TAU_SW);
  auto spa_data_aero_tau_lw = ekat::scalarize(spa_input.data.AER_TAU_LW);

  auto copy_and_pad = KOKKOS_LAMBDA (const Member& team) {
    int icol = team.league_rank();

    auto copy_col = [&](const int k) {
      spa_data_ccn3(icol,k+1) = ccn3(icol,k);
      for (int isw=0; isw<nswbands; ++isw) {
        spa_data_aero_g_sw(icol,isw,k+1)   = aero_g_sw(icol,isw,k);
        spa_data_aero_ssa_sw(icol,isw,k+1) = aero_ssa_sw(icol,isw,k);
        spa_data_aero_tau_sw(icol,isw,k+1) = aero_tau_sw(icol,isw,k);
      }
      for (int ilw=0; ilw<nlwbands; ++ilw) {
        spa_data_aero_tau_lw(icol,ilw,k+1) = aero_tau_lw(icol,ilw,k);
      }
    };
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team,nlevs),copy_col);

    // Set the first/last entries of the spa data, so that linear interp
    // can extrapolate if the p_tgt is outside the p_src bounds
    Kokkos::single(Kokkos::PerTeam(team),[&]{
      spa_data_ccn3(icol,0) = 0;
      for (int isw=0; isw<nswbands; ++isw) {
        spa_data_aero_g_sw(icol,isw,0)   = 0;
        spa_data_aero_ssa_sw(icol,isw,0) = 0;
        spa_data_aero_tau_sw(icol,isw,0) = 0;
      }
      for (int ilw=0; ilw<nlwbands; ++ilw) {
        spa_data_aero_tau_lw(icol,ilw,0) = 0;
      }
      spa_data_ccn3(icol,nlevs+1) = ccn3(icol,nlevs-1);
      for (int isw=0; isw<nswbands; ++isw) {
        spa_data_aero_g_sw(icol,isw,nlevs+1)   = aero_g_sw(icol,isw,nlevs-1);
        spa_data_aero_ssa_sw(icol,isw,nlevs+1) = aero_ssa_sw(icol,isw,nlevs-1);
        spa_data_aero_tau_sw(icol,isw,nlevs+1) = aero_tau_sw(icol,isw,nlevs-1);
      }
      for (int ilw=0; ilw<nlwbands; ++ilw) {
        spa_data_aero_tau_lw(icol,ilw,nlevs+1) = aero_tau_lw(icol,ilw,nlevs-1);
      }
    });
  };
  auto policy = ESU::get_default_team_policy(ncols,nlevs);
  Kokkos::parallel_for("", policy, copy_and_pad);
  Kokkos::fence();
  stop_timer("EAMxx::SPA::update_spa_data_from_file::copy_and_pad");

  stop_timer("EAMxx::SPA::update_spa_data_from_file");
} // END update_spa_data_from_file

/*-----------------------------------------------------------------*/
template<typename S, typename D>
void SPAFunctions<S,D>
::update_spa_timestate(
    std::shared_ptr<AtmosphereInput>& scorpio_reader,
    std::shared_ptr<IOPReader>&       iop_reader,
    const util::TimeStamp&            ts,
    AbstractRemapper&                 spa_horiz_interp,
    SPATimeState&                     time_state,
    SPAInput&                         spa_beg,
    SPAInput&                         spa_end)
{
  // Now we check if we have to update the data that changes monthly
  // NOTE:  This means that SPA assumes monthly data to update.  Not
  //        any other frequency.
  const auto month = ts.get_month() - 1; // Make it 0-based
  if (month != time_state.current_month) {
    // Update the SPA time state information
    time_state.current_month = month;
    time_state.t_beg_month = util::TimeStamp({ts.get_year(),month+1,1}, {0,0,0}).frac_of_year_in_days();
    time_state.days_this_month = util::days_in_month(ts.get_year(),month+1);

    // Copy spa_end'data into spa_beg'data, and read in the new spa_end
    std::swap(spa_beg,spa_end);

    // Update the SPA forcing data for this month and next month
    // Start by copying next months data to this months data structure.
    // NOTE: If the timestep is bigger than monthly this could cause the wrong values
    //       to be assigned.  A timestep greater than a month is very unlikely so we
    //       will proceed.
    int next_month = (time_state.current_month + 1) % 12;
    update_spa_data_from_file(scorpio_reader,iop_reader,ts,next_month,spa_horiz_interp,spa_end);
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
