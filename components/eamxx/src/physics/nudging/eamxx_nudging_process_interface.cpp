#include "eamxx_nudging_process_interface.hpp"

#include "share/field/field_reader.hpp"
#include "share/util/eamxx_utils.hpp"
#include "share/util/eamxx_universal_constants.hpp"
#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"

#include <ekat_math_utils.hpp>

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
  // use nudging weights
  if (m_use_weights)
    m_weights_file = m_params.get<std::string>("nudging_weights_file");

  // TODO: Add some warning messages here.
  // 1. if m_timescale is <= 0 we will do direct replacement.
  // 2. if m_fields_nudge is empty or =NONE then we will skip nudging altogether.
}

// =========================================================================================
void Nudging::create_requests()
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  m_grid = m_grids_manager->get_grid("physics");
  const auto& grid_name = m_grid->name();
  m_num_cols = m_grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = m_grid->get_num_vertical_levels();  // Number of levels per column

  FieldLayout scalar3d_layout_mid = m_grid->get_3d_scalar_layout(LEV);
  FieldLayout horiz_wind_layout = m_grid->get_3d_vector_layout(LEV,2);

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

  auto num_cols_global = m_grid->get_num_global_dofs();
  int num_cols_src = scorpio::get_dimlen(m_datafiles[0],"ncol");

  if (num_cols_src != num_cols_global) {
    EKAT_REQUIRE_MSG(m_refine_remap_file != "no-file-given",
                     "Error! Nudging::create_requests - the number of columns in the nudging data file "
                     << std::to_string(num_cols_src) << " does not match the number of columns in the "
                     << "model grid " << std::to_string(num_cols_global) << ".  Please check the "
                     << "nudging data file and/or the model grid.");
    EKAT_REQUIRE_MSG(m_use_weights == false,
                     "Error! Nudging::create_requests - it seems that the user intends to use both nudging "
                     << "from coarse data as well as weighted nudging simultaneously. This is not supported. "
                     << "If the user wants to use both at their own risk, the user should edit the source code "
                     << "by deleting this error message.");
    m_refine_remap = true;
  } else {
    if (m_refine_remap_file != "no-file-given") {
      m_atm_logger->warn(
           "[Nudging::create_requests] Warning! Map file provided, but it is not needed.\n"
           "  - num cols in nudging data file: " + std::to_string(num_cols_src) + "\n"
           "  - num cols in model grid       : " + std::to_string(num_cols_global) + "\n"
           " Please, make sure the nudging data file and/or model grid are correct.\n"
           " The map file is only needed if the above two numbers differ.");
    }
    if (m_refine_remap_vert_cutoff > 0.0) {
      m_atm_logger->warn(
          "[Nudging::create_requests] Warning! Non-zero vertical cutoff provided, but it is not needed\n"
          " - vertical cutoff: " + std::to_string(m_refine_remap_vert_cutoff) + "\n"
          " Please, check your settings. This parameter is only needed if we are remapping.");
    }
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
  const auto layout_atm = m_grid->get_3d_scalar_layout(LEV);

  std::vector<Field> nudge_fields;
  for (const auto& name : m_fields_nudge) {
    nudge_fields.push_back(create_helper_field(name,layout_atm,m_grid->name()));
  }

  m_data_interpolation = std::make_shared<DataInterpolation>(m_grid,nudge_fields);
  m_data_interpolation->set_name("Nudging");
  m_data_interpolation->set_logger(m_atm_logger);
  m_data_interpolation->setup_time_database(m_datafiles,util::TimeLine::Linear,DataInterpolation::Linear);
  m_data_interpolation->create_horiz_remappers(m_refine_remap ? m_refine_remap_file : "");

  DataInterpolation::VertRemapData vremap_data;
  if (m_skip_vert_interpolation) {
    vremap_data.vr_type = DataInterpolation::None;
  } else if (m_src_pres_type == TIME_DEPENDENT_3D_PROFILE) {
    vremap_data.vr_type = DataInterpolation::Dynamic3D;
    vremap_data.pname = "p_mid";
    vremap_data.pmid  = get_field_in("p_mid");
  } else {
    vremap_data.vr_type = DataInterpolation::Static1D;
    vremap_data.pname = "p_levs";
    vremap_data.pmid  = get_field_in("p_mid");
    if (m_static_vertical_pressure_file != m_datafiles[0]) {
      vremap_data.pfile = m_static_vertical_pressure_file;
    }
  }
  m_data_interpolation->create_vert_remapper(vremap_data);
  m_data_interpolation->init_data_interval(start_of_step_ts());

  if (m_use_weights) {
    auto nudging_weights = create_helper_field("nudging_weights", layout_atm, m_grid->name());
    read_fields(m_weights_file,{nudging_weights},m_grid->get_partitioned_dim_gids(),m_comm);
  }
}

// =========================================================================================
void Nudging::run_impl (const double dt)
{
  m_data_interpolation->run(end_of_step_ts());

  for (const auto& name : m_fields_nudge) {
    auto nudge_field = get_helper_field(name);
    correct_fill_values(nudge_field);

    auto state_field = get_field_out_wrap(name);

    if (m_timescale <= 0) {
      state_field.deep_copy(nudge_field);
    } else {
      apply_tendency(state_field,nudge_field,dt);
    }
  }
}

// =========================================================================================
void Nudging::finalize_impl()
{}

// =========================================================================================
void Nudging::correct_fill_values (const Field& f) const
{
  using namespace ShortFieldTagsNames;
  const auto& fl = f.get_header().get_identifier().get_layout();

  const bool has_vert_dim = fl.has_tag(LEV) or fl.has_tag(ILEV);
  if (not fl.has_tag(COL) or not has_vert_dim) {
    return;
  }

  using KT = KokkosTypes<DefaultDevice>;
  using RangePolicy = typename KT::RangePolicy;

  constexpr Real fill_value = constants::fill_value<Real>;
  const auto thresh = std::abs(fill_value)*0.0001;

  const int ncols = fl.dim(COL);
  const int nlevs = fl.dim(fl.rank()-1);

  if (fl.rank()==2) {
    const auto v = f.get_view<Real**>();
    Kokkos::parallel_for(RangePolicy(0,ncols), KOKKOS_LAMBDA(const int icol) {
      int first_good = nlevs;
      int last_good  = -1;
      for (int k=0; k<nlevs; ++k) {
        if (Kokkos::abs(v(icol,k)-fill_value)>thresh) {
          first_good = ekat::impl::min(first_good,k);
          last_good  = ekat::impl::max(last_good,k);
        }
      }
      if (first_good==nlevs or last_good<0) return; // all fill — skip
      for (int k=0; k<first_good; ++k)       v(icol,k) = v(icol,first_good);
      for (int k=last_good+1; k<nlevs; ++k)  v(icol,k) = v(icol,last_good);
    });
  } else if (fl.rank()==3) {
    const auto v = f.get_view<Real***>();
    const int ncomps = fl.get_vector_dim();
    Kokkos::parallel_for(RangePolicy(0,ncols*ncomps), KOKKOS_LAMBDA(const int idx) {
      const int icol = idx / ncomps;
      const int icmp = idx % ncomps;
      int first_good = nlevs;
      int last_good  = -1;
      for (int k=0; k<nlevs; ++k) {
        if (Kokkos::abs(v(icol,icmp,k)-fill_value)>thresh) {
          first_good = ekat::impl::min(first_good,k);
          last_good  = ekat::impl::max(last_good,k);
        }
      }
      if (first_good==nlevs or last_good<0) return; // all fill — skip
      for (int k=0; k<first_good; ++k)       v(icol,icmp,k) = v(icol,icmp,first_good);
      for (int k=last_good+1; k<nlevs; ++k)  v(icol,icmp,k) = v(icol,icmp,last_good);
    });
  }
}
// =========================================================================================
Field Nudging::create_helper_field (const std::string& name,
                                             const FieldLayout& layout,
                                             const std::string& grid_name,
                                             const int ps)
{
  // For helper fields we don't bother w/ units, so we set them to non-dimensional
  FieldIdentifier id(name,layout,ekat::units::none,grid_name);

  // Create the field. Init with NaN's, so we spot instances of uninited memory usage
  Field f(id);
  f.get_header().get_alloc_properties().request_allocation(ps);
  f.allocate_view();
  f.deep_copy(ekat::invalid<Real>());

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
