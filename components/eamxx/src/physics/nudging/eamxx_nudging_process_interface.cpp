#include "eamxx_nudging_process_interface.hpp"

#include "share/util/eamxx_universal_constants.hpp"
#include "share/grid/remap/refining_remapper_p2p.hpp"
#include "share/util/eamxx_utils.hpp"
#include "share/io/eamxx_scorpio_interface.hpp"

#include <ekat/util/ekat_lin_interp.hpp>
#include <ekat/util/ekat_math_utils.hpp>
#include <ekat/kokkos/ekat_kokkos_utils.hpp>

namespace scream
{

// =========================================================================================
Nudging::Nudging (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  m_datafiles  = filename_glob(m_params.get<std::vector<std::string>>("nudging_filenames_patterns"));
  m_timescale = m_params.get<int>("nudging_timescale",0);

  m_fields_nudge = m_params.get<std::vector<std::string>>("nudging_fields");
  m_use_weights   = m_params.get<bool>("use_nudging_weights",false);
  m_skip_vert_interpolation   = m_params.get<bool>("skip_vert_interpolation",false);
  // If we are doing horizontal refine-remapping, we need to get the mapfile from user
  m_refine_remap_file = m_params.get<std::string>(
      "nudging_refine_remap_mapfile", "no-file-given");
  m_refine_remap_vert_cutoff = m_params.get<Real>(
      "nudging_refine_remap_vert_cutoff", 0.0);
  auto src_pres_type = m_params.get<std::string>("source_pressure_type","TIME_DEPENDENT_3D_PROFILE");
  if (src_pres_type=="TIME_DEPENDENT_3D_PROFILE") {
    m_src_pres_type = TIME_DEPENDENT_3D_PROFILE;
  } else if (src_pres_type=="STATIC_1D_VERTICAL_PROFILE") {
    m_src_pres_type = STATIC_1D_VERTICAL_PROFILE;
    // Check for a designated source pressure file, default to first nudging data source if not given.
    m_static_vertical_pressure_file = m_params.get<std::string>("source_pressure_file",m_datafiles[0]);
    EKAT_REQUIRE_MSG(m_skip_vert_interpolation == false,
                   "Error! It makes no sense to not interpolate if src press is uniform and constant ");
  } else {
    EKAT_ERROR_MSG("ERROR! Nudging::parameter_list - unsupported source_pressure_type provided.  Current options are [TIME_DEPENDENT_3D_PROFILE,STATIC_1D_VERTICAL_PROFILE].  Please check");
  }
  int first_file_levs = scorpio::get_dimlen(m_datafiles[0],"lev");
  for (const auto& file : m_datafiles) {
      int current_file_levs = scorpio::get_dimlen(file,"lev");
      EKAT_REQUIRE_MSG(current_file_levs == first_file_levs,
                       "Error! Inconsistent 'lev' dimension found in nudging data files.");
  }
  // use nudging weights
  if (m_use_weights)
    m_weights_file = m_params.get<std::string>("nudging_weights_file");

  // TODO: Add some warning messages here.
  // 1. if m_timescale is <= 0 we will do direct replacement.
  // 2. if m_fields_nudge is empty or =NONE then we will skip nudging altogether.
}

// =========================================================================================
void Nudging::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;

  m_grid = grids_manager->get_grid("Physics");
  const auto& grid_name = m_grid->name();
  m_num_cols = m_grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = m_grid->get_num_vertical_levels();  // Number of levels per column

  FieldLayout scalar3d_layout_mid = m_grid->get_3d_scalar_layout(true);
  FieldLayout horiz_wind_layout = m_grid->get_3d_vector_layout(true,2);

  constexpr int ps = 1;
  add_field<Required>("p_mid", scalar3d_layout_mid, Pa, grid_name, ps);

  /* ----------------------- WARNING --------------------------------*/
  /* The following is a HACK to get things moving, we don't want to
   * add all fields as "updated" long-term.  A separate stream of work
   * is adapting the infrastructure to allow for a generic "add_field" call
   * to be used here which we can then setup using the m_fields_nudge variable
   * For now we check if a field is intended to be nudged via the m_fields_nudge
   * vector, if it is we register it.  For now we are limited to just T_mid, qv,
   * U and V
   */
  if (ekat::contains(m_fields_nudge,"T_mid")) {
    add_field<Updated>("T_mid", scalar3d_layout_mid, K, grid_name, ps);
  }
  if (ekat::contains(m_fields_nudge,"qv")) {
    add_tracer<Updated>("qv", m_grid, kg/kg, ps);
  }
  if (ekat::contains(m_fields_nudge,"U") or ekat::contains(m_fields_nudge,"V")) {
    add_field<Updated>("horiz_winds",   horiz_wind_layout,   m/s,     grid_name, ps);
  }

  /* ----------------------- WARNING --------------------------------*/

  //Now need to read in the file
  if (m_src_pres_type == TIME_DEPENDENT_3D_PROFILE) {
    m_num_src_levs = scorpio::get_dimlen(m_datafiles[0],"lev");
  } else {
    m_num_src_levs = scorpio::get_dimlen(m_static_vertical_pressure_file,"lev");
  }
  if (m_skip_vert_interpolation) {
    EKAT_REQUIRE_MSG(m_num_src_levs == m_num_levs,
                   "Error! skip_vert_interpolation requires the vertical level to be "
                   << " the same as model vertical level ");
  }

  /* Check for consistency between nudging files, map file, and remapper */

  // Number of columns globally
  auto num_cols_global = m_grid->get_num_global_dofs();

  // Get the information from the first nudging data file
  int num_cols_src = scorpio::get_dimlen(m_datafiles[0],"ncol");

  if (num_cols_src != num_cols_global) {
    // If differing cols, check if remap file is provided
    EKAT_REQUIRE_MSG(m_refine_remap_file != "no-file-given",
                     "Error! Nudging::set_grids - the number of columns in the nudging data file "
                     << std::to_string(num_cols_src) << " does not match the number of columns in the "
                     << "model grid " << std::to_string(num_cols_global) << ".  Please check the "
                     << "nudging data file and/or the model grid.");
    // If remap file is provided, check if it is consistent with the nudging data file
    // First get the data from the mapfile
    int num_cols_remap_a = scorpio::get_dimlen(m_refine_remap_file,"n_a");
    int num_cols_remap_b = scorpio::get_dimlen(m_refine_remap_file,"n_b");
    // Then, check if n_a (source) and n_b (target) are consistent
    EKAT_REQUIRE_MSG(num_cols_remap_a == num_cols_src,
                     "Error! Nudging::set_grids - the number of columns in the nudging data file "
                     << std::to_string(num_cols_src) << " does not match the number of columns in the "
                     << "mapfile " << std::to_string(num_cols_remap_a) << ".  Please check the "
                     << "nudging data file and/or the mapfile.");
    EKAT_REQUIRE_MSG(num_cols_remap_b == num_cols_global,
                     "Error! Nudging::set_grids - the number of columns in the model grid "
                     << std::to_string(num_cols_global) << " does not match the number of columns in the "
                     << "mapfile " << std::to_string(num_cols_remap_b) << ".  Please check the "
                     << "model grid and/or the mapfile.");
    EKAT_REQUIRE_MSG(m_use_weights == false,
                     "Error! Nudging::set_grids - it seems that the user intends to use both nuding "
                     << "from coarse data as well as weighted nudging simultaneously. This is not supported. "
                     << "If the user wants to use both at their own risk, the user should edit the source code "
                     << "by deleting this error message.");
    // If we get here, we are good to go!
    m_refine_remap = true;
  } else {
    // If the number of columns is the same, we don't need to do any remapping,
    // but print a warning if the user provided a mapfile
    if (m_refine_remap_file != "no-file-given") {
      m_atm_logger->warn(
           "[Nudging::set_grids] Warning! Map file provided, but it is not needed.\n"
           "  - num cols in nudging data file: " + std::to_string(num_cols_src) + "\n"
           "  - num cols in model grid       : " + std::to_string(num_cols_global) + "\n"
           " Please, make sure the nudging data file and/or model grid are correct.\n"
           " The map file is only needed if the above two numbers differ.");
    }
    // If the user gives us the vertical cutoff, warn them
    if (m_refine_remap_vert_cutoff > 0.0) {
      m_atm_logger->warn(
          "[Nudging::set_grids] Warning! Non-zero vertical cutoff provided, but it is not needed\n"
          " - vertical cutoff: " + std::to_string(m_refine_remap_vert_cutoff) + "\n"
          " Please, check your settings. This parameter is only needed if we are remapping.");
    }
    // Set m_refine_remap to false
    m_refine_remap = false;
  }
}
// =========================================================================================
void Nudging::apply_tendency(Field& state, const Field& nudge, const Real dt) const
{
  // Calculate the weight to apply the tendency
  const Real dtend = dt / m_timescale;

  using cview_2d = decltype(state.get_view<const Real**>());

  auto state_view = state.get_view<Real**>();
  auto nudge_view = nudge.get_view<Real**>();
  cview_2d w_view, pmid_view;

  if (m_use_weights) {
    auto weights = get_helper_field("nudging_weights");
    w_view = weights.get_view<const Real**>();
  }
  if (m_refine_remap_vert_cutoff>0) {
    pmid_view = get_field_in("p_mid").get_view<const Real**>();
  }

  auto use_weights = m_use_weights;
  auto cutoff = m_refine_remap_vert_cutoff;
  auto policy = Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {m_num_cols, m_num_levs});
  auto update = KOKKOS_LAMBDA(const int& i, const int& j) {
    if (cutoff>0 and pmid_view(i,j)>=cutoff) {
      return;
    }

    auto tend = nudge_view(i,j) - state_view(i,j);
    if (use_weights) {
      tend *= w_view(i,j);
    }
    state_view(i,j) += dtend * tend;
  };
  Kokkos::parallel_for(policy,update);
}
// =============================================================================================================
void Nudging::initialize_impl (const RunType /* run_type */)
{
  using namespace ShortFieldTagsNames;

  // The first thing we do is time interpolation.
  // The second thing we do is horiz interpolation. The reason for doing horizontal
  // before vertical is that we do not have a coarse p_mid (which would be needed
  // as tgt pressure during vert remap). To get it, we'd have to "remap back" p_mid,
  // but that seems overly complicated. For horiz remap we do not need anything
  // on the target grid.

  // The "intermediate" grid, is the grid after horiz remap, and before vert remap
  auto grid_tmp = m_grid->clone("after_horiz_before_vert",true);
  grid_tmp->reset_num_vertical_lev(m_num_src_levs);

  if (m_refine_remap) {
    // P2P remapper
    m_horiz_remapper = std::make_shared<RefiningRemapperP2P>(grid_tmp, m_refine_remap_file);
  } else {
    // We set up an IdentityRemapper, specifying that tgt is an alias
    // of src, so that the remap method will do nothing
    auto r = std::make_shared<IdentityRemapper>(grid_tmp);
    r->set_aliasing(IdentityRemapper::TgtAliasSrc);
    m_horiz_remapper = r;
  }

  // Now that we have the remapper, we can grab the grid where the input data lives
  auto grid_ext = m_horiz_remapper->get_src_grid();

  // Initialize the time interpolator and horiz remapper
  m_time_interp = util::TimeInterpolation(grid_ext, m_datafiles);
  m_time_interp.set_logger(m_atm_logger,"[EAMxx::Nudging] Reading nudging data");

  // NOTE: we are ASSUMING all fields are 3d and scalar!
  const auto layout_ext = grid_ext->get_3d_scalar_layout(true);
  const auto layout_tmp = grid_tmp->get_3d_scalar_layout(true);
  const auto layout_atm = m_grid->get_3d_scalar_layout(true);
  m_horiz_remapper->registration_begins();
  for (auto name : m_fields_nudge) {
    std::string name_ext = name + "_ext";
    std::string name_tmp = name + "_tmp";

    // First copy of the field: what's read from file, and time-interpolated.
    auto field_ext = create_helper_field(name_ext, layout_ext, grid_ext->name());

    // Second copy of the field: after horiz interp (alias "ext" if no remap)
    Field field_tmp;
    if (m_refine_remap) {
      field_tmp = create_helper_field(name_tmp, layout_tmp, grid_tmp->name());
    } else {
      field_tmp = field_ext.alias(name_tmp);
      m_helper_fields[name_tmp] = field_tmp;
    }

    // Add the field to the time interpolator
    m_time_interp.add_field(field_ext.alias(name), true);

    // Register the fields with the remapper
    m_horiz_remapper->register_field(field_ext, field_tmp);

    if (m_timescale>0) {
      // Third copy of the field: after vert interpolation.
      // We cannot store directly in get_field_out(name),
      // since we need to back out tendencies
      create_helper_field(name, layout_atm, m_grid->name());
    } else {
      // We do not need to back out any tendency; the input data is used
      // to directly replace the atm state
      m_helper_fields[name] = get_field_out_wrap(name);
    }
  }

  // A helper field, where we copy each field after horiz remap, padding it
  // at top/bot, to allow vert lin interp to extrapolate outside the bounds of p_mid
  FieldLayout layout_padded ({COL,LEV},{m_num_cols,m_num_src_levs+2});
  create_helper_field("padded_field",layout_padded,"");

  if (m_src_pres_type == TIME_DEPENDENT_3D_PROFILE && !m_skip_vert_interpolation) {
    // If the pressure profile is 3d and time-dep, we need to interpolate (in time/horiz)
    auto pmid_ext = create_helper_field("p_mid_ext", layout_ext, grid_ext->name());
    m_time_interp.add_field(pmid_ext.alias("p_mid"),true);
    Field pmid_tmp;
    if (m_refine_remap) {
      pmid_tmp = create_helper_field("p_mid_tmp", layout_tmp, grid_tmp->name());
    } else {
      pmid_tmp = pmid_ext.alias("p_mid_tmp");
      m_helper_fields["p_mid_tmp"] = pmid_tmp;
    }
    m_horiz_remapper->register_field(pmid_ext,pmid_tmp);
    create_helper_field("padded_p_mid_tmp",layout_padded,"");
  } else if (m_src_pres_type == STATIC_1D_VERTICAL_PROFILE) {
    // For static 1D profile, we can read p_mid now
    auto pmid_ext = create_helper_field("p_mid_ext", grid_ext->get_vertical_layout(true), grid_ext->name());
    AtmosphereInput src_input(m_static_vertical_pressure_file,grid_ext,{pmid_ext.alias("p_levs")},true);
    src_input.read_variables(-1);

    // For static 1d profile, p_mid_tmp is an alias of p_mid_ext
    m_helper_fields["p_mid_tmp"] = pmid_ext.alias("p_mid_tmp");

    // The padded p_mid is also 1d
    FieldLayout pmid1d_padded_layout({COL},{m_num_src_levs+2});
    create_helper_field("padded_p_mid_tmp",pmid1d_padded_layout,"");
  }

  // Close the registration
  m_time_interp.initialize_data_from_files();
  m_horiz_remapper->registration_ends();

  // load nudging weights from file
  // NOTE: the regional nudging use the same grid as the run, no need to
  // do the interpolation.
  if (m_use_weights) {
    auto nudging_weights = create_helper_field("nudging_weights", layout_atm, m_grid->name());
    AtmosphereInput src_weights_input(m_weights_file, m_grid, {nudging_weights},true);
    src_weights_input.read_variables();
  }
}

// =========================================================================================
void Nudging::run_impl (const double dt)
{
  using KT            = KokkosTypes<DefaultDevice>;
  using RangePolicy   = typename KT::RangePolicy;
  using MemberType    = typename KT::MemberType;
  using ESU           = ekat::ExeSpaceUtils<typename KT::ExeSpace>;
  using PackT         = ekat::Pack<Real,1>;
  using view_1d       = KT::view_1d<PackT>;
  using view_2d       = KT::view_2d<PackT>;

  // Perform time interpolation
  m_time_interp.perform_time_interpolation(end_of_step_ts());

  // If the input data contains "masked" values (sometimes also called "filled" values),
  // the horiz remapping would smear them around. To prevent that, we need to "cure"
  // these values. Masked values can only happen at top/bot of the model (with top
  // being not common), and they must be a contiguous set of entries. So to cure them,
  // we simply set all bot/top masked entries equal to the first non-masked value
  // from the bot/top respectively. This corresponds to a constant extrapolation.
  // NOTE: we need to do a tol check, since time interpolation may not return fillValue,
  //       even if both f(t_beg)/f(t_end) are equal to fillValue (due to rounding).
  // NOTE: if f(t_beg)==fillValue!=f(t_end), or viceversa, the time-interpolated value can
  //       substantially differ from fillValue. Here, we assume it didn't happen.
  auto correct_masked_values = [&](const Field f) {
    const auto fl = f.get_header().get_identifier().get_layout();
    const auto v  = f.get_view<Real**>();

    Real var_fill_value = constants::DefaultFillValue<Real>().value;
    // Query the helper field for the fill value, if not present use default
    if (f.get_header().has_extra_data("mask_value")) {
      var_fill_value = f.get_header().get_extra_data<Real>("mask_value");
    }

    const int ncols = fl.dim(0);
    const int nlevs = fl.dim(1);
    const auto thresh = std::abs(var_fill_value)*0.0001;
    auto lambda = KOKKOS_LAMBDA(const int icol) {
      int first_good = nlevs;
      int last_good = -1;
      for (int k=0; k<nlevs; ++k) {
        if (std::abs(v(icol,k)-var_fill_value)>thresh) {
          // This entry is substantially different from var_fill_value, so it's good
          first_good = ekat::impl::min(first_good,k);
          last_good  = ekat::impl::max(last_good,k);
        }
      }
      EKAT_KERNEL_REQUIRE_MSG (first_good<nlevs and last_good>=0,
          "[Nudging] Error! Could not locate a non-masked entry in a column.\n");

      // Fix near TOM
      for (int k=0; k<first_good; ++k) {
        v(icol,k) = v(icol,first_good);
      }
      // Fix near surf
      for (int k=last_good+1; k<nlevs; ++k) {
        v(icol,k) = v(icol,last_good);
      }
    };

    Kokkos::parallel_for(RangePolicy(0,ncols),lambda);
  };

  // Correct before horiz remap
  for (const auto& name: m_fields_nudge) {
    const auto f  = get_helper_field(name+"_ext");
    correct_masked_values(f);
  }

  // Perform horizontal remap (if needed)
  m_horiz_remapper->remap_fwd();

  // bypass copy_and_pad and vert_interp for skip_vert_interpolation:
  if (m_skip_vert_interpolation) {
    for (const auto& name : m_fields_nudge) {
      auto tmp_state_field = get_helper_field(name+"_tmp");

      if (m_timescale > 0) {
        auto atm_state_field = get_field_out_wrap(name);
        apply_tendency(atm_state_field,tmp_state_field,dt);
      }
    }
    return;
  }

  // Copy remapper tgt fields into padded views, to allow extrapolation at top/bot,
  // then call remapping routines

  const int ncols = m_num_cols;
  const int nlevs_src = m_num_src_levs;
  auto copy_and_pad = [&](const Field from, const Field to, const bool is_pmid) {
    auto from_view = from.get_view<const Real**>();
    auto to_view = to.get_view<Real**>();
    auto fl = from.get_header().get_identifier().get_layout();

    auto copy_3d = KOKKOS_LAMBDA (const MemberType& team) {
      int icol = team.league_rank();

      auto copy_col = [&](const int k) {
        to_view(icol,k+1) = from_view(icol,k);
      };
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team,nlevs_src),copy_col);

      // Set the first/last entries of data, so that linear interp
      // can extrapolate if the p_tgt is outside the p_src bounds
      Kokkos::single(Kokkos::PerTeam(team),[&]{
        to_view(icol,0) = 0; // Does this make sense for *every field*?
        if (is_pmid) {
          // For pmid, we put a very large value, so that any p_mid_tgt
          // that is larger than input p_mid bnds will end up in the
          // last interval.
          to_view(icol,nlevs_src+1) = 1e7;
        } else {
          // For data, we set last entry equal to second-to-last.
          // This will cause constant extrapolation outside of
          // the input p_mid bounds
          to_view(icol,nlevs_src+1) = from_view(icol,nlevs_src-1);
        }
      });
    };

    auto policy = ESU::get_default_team_policy(ncols,nlevs_src);
    Kokkos::parallel_for("", policy, copy_3d);
  };

  // First, copy/pad p_mid, and extract the right copy (1d vs 3d)
  if (m_src_pres_type==TIME_DEPENDENT_3D_PROFILE) {
    copy_and_pad (get_helper_field("p_mid_tmp"),get_helper_field("padded_p_mid_tmp"),true);
  } else {
    // pmid is a 1d view. Just pad by hand
    auto from = get_helper_field("p_mid_tmp");
    auto to   = get_helper_field("padded_p_mid_tmp");
    auto from_v = from.get_view<const Real*>();
    auto to_v   = to.get_view<Real*>();
    auto lambda = KOKKOS_LAMBDA(const int& lev) {
      to_v(lev+1) = from_v(lev);
      if (lev==0) {
        to_v(0) = 0;
      } else if (lev==nlevs_src-1) {
        to_v(lev+2) = to_v(lev+1);
      }
    };
    Kokkos::parallel_for(RangePolicy(0,nlevs_src),lambda);
  }
  const auto& p_mid_v = get_field_in("p_mid").get_view<const PackT**>();
  view_2d p_mid_tmp_3d;
  view_1d p_mid_tmp_1d;
  bool src_pmid_3d;
  if (m_src_pres_type == TIME_DEPENDENT_3D_PROFILE) {
    p_mid_tmp_3d = get_helper_field("padded_p_mid_tmp").get_view<PackT**>();
    src_pmid_3d = true;
  } else {
    p_mid_tmp_1d = get_helper_field("padded_p_mid_tmp").get_view<PackT*>();
    src_pmid_3d = false;
  }

  // Setup the linear interpolation object
  using LI = ekat::LinInterp<Real,1>;
  const int nlevs_tgt = m_num_levs;
  LI vert_interp(ncols,nlevs_src+2,nlevs_tgt);
  const auto policy_vinterp = ESU::get_default_team_policy(ncols, nlevs_tgt);
  auto p_tgt = get_field_in("p_mid").get_view<const PackT**>();
  Kokkos::parallel_for("nudging_vert_interp_setup_loop", policy_vinterp,
    KOKKOS_LAMBDA(const MemberType& team) {

    const int icol = team.league_rank();

    // Setup
    if (src_pmid_3d) {
      vert_interp.setup(team, ekat::subview(p_mid_tmp_3d,icol),
                              ekat::subview(p_tgt,icol));
    } else {
      vert_interp.setup(team, p_mid_tmp_1d,
                              ekat::subview(p_tgt,icol));
    }
  });
  Kokkos::fence();

  // Then loop over fields, and do copy_and_pad + vremap
  auto padded_field = get_helper_field("padded_field");
  for (const auto& name : m_fields_nudge) {
    copy_and_pad(get_helper_field(name+"_tmp"),padded_field,false);

    auto field_after_vinterp = get_helper_field(name);
    auto view_in  = padded_field.get_view<const PackT**>();
    auto view_out = field_after_vinterp.get_view<PackT**>();
    // Now use the interpolation object in || over all variables.
    auto vinterp = KOKKOS_LAMBDA(const MemberType& team) {
      const int icol = team.league_rank();

      view_1d x1;
      if (src_pmid_3d) {
        x1 = ekat::subview(p_mid_tmp_3d,icol);
      } else {
        x1 = p_mid_tmp_1d;
      }
      auto x2 = ekat::subview(p_tgt,icol);

      auto y1 = ekat::subview(view_in, icol);
      auto y2 = ekat::subview(view_out,icol);

      vert_interp.lin_interp(team, x1, x2, y1, y2, icol);
    };
    Kokkos::parallel_for("nudging_vert_interp_loop", policy_vinterp, vinterp);
    Kokkos::fence();

    // If timescale==0, the call get_helper_field(name) returns the same
    // fields as get_field_out_wrap(name) they are alias, so nothing to do.
    // If timescale>0, then we need to back out a tendency.
    if (m_timescale > 0) {
      auto atm_state_field = get_field_out_wrap(name);
      apply_tendency(atm_state_field,field_after_vinterp,dt);
    }
  }
}

// =========================================================================================
void Nudging::finalize_impl()
{
  m_time_interp.finalize();
}
// =========================================================================================
Field Nudging::create_helper_field (const std::string& name,
                                             const FieldLayout& layout,
                                             const std::string& grid_name,
                                             const int ps)
{
  using namespace ekat::units;

  // For helper fields we don't bother w/ units, so we set them to non-dimensional
  FieldIdentifier id(name,layout,Units::nondimensional(),grid_name);

  // Create the field. Init with NaN's, so we spot instances of uninited memory usage
  Field f(id);
  f.get_header().get_alloc_properties().request_allocation(ps);
  f.allocate_view();
  f.deep_copy(ekat::ScalarTraits<Real>::invalid());

  m_helper_fields[name] = f;
  return m_helper_fields[name];
}

// =========================================================================================
Field Nudging::get_field_out_wrap(const std::string& field_name) {
  if (field_name == "U" or field_name == "V") {
    auto hw = get_field_out("horiz_winds");
    if (field_name == "U") {
      return hw.get_component(0);
    } else {
      return hw.get_component(1);
    }
  } else {
    return get_field_out(field_name);
  }
}

} // namespace scream
