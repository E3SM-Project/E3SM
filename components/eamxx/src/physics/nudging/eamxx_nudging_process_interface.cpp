#include "eamxx_nudging_process_interface.hpp"
#include "share/util/scream_universal_constants.hpp"

namespace scream
{

// =========================================================================================
Nudging::Nudging (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  m_datafiles  = m_params.get<std::vector<std::string>>("nudging_filename");
  m_timescale = m_params.get<int>("nudging_timescale",0);
  m_fields_nudge = m_params.get<std::vector<std::string>>("nudging_fields");
  m_use_weights   = m_params.get<bool>("use_nudging_weights",false);
  auto src_pres_type = m_params.get<std::string>("source_pressure_type","TIME_DEPENDENT_3D_PROFILE");
  if (src_pres_type=="TIME_DEPENDENT_3D_PROFILE") {
    m_src_pres_type = TIME_DEPENDENT_3D_PROFILE;
  } else if (src_pres_type=="STATIC_1D_VERTICAL_PROFILE") {
    m_src_pres_type = STATIC_1D_VERTICAL_PROFILE;
    // Check for a designated source pressure file, default to first nudging data source if not given.
    m_static_vertical_pressure_file = m_params.get<std::string>("source_pressure_file",m_datafiles[0]);
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
void Nudging::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  m_grid = grids_manager->get_grid("Physics");
  const auto& grid_name = m_grid->name();
  m_num_cols = m_grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = m_grid->get_num_vertical_levels();  // Number of levels per column

  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols, m_num_levs} };
  FieldLayout horiz_wind_layout { {COL,CMP,LEV}, {m_num_cols,2,m_num_levs} };

  constexpr int ps = 1;
  auto Q = kg/kg;
  Q.set_string("kg/kg");
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
    add_field<Updated>("qv",    scalar3d_layout_mid, Q, grid_name, "tracers", ps);
  }
  if (ekat::contains(m_fields_nudge,"U") or ekat::contains(m_fields_nudge,"V")) {
    add_field<Updated>("horiz_winds",   horiz_wind_layout,   m/s,     grid_name, ps);
  }

  /* ----------------------- WARNING --------------------------------*/

  //Now need to read in the file
  if (m_src_pres_type == TIME_DEPENDENT_3D_PROFILE) {
    scorpio::register_file(m_datafiles[0],scorpio::Read);
    m_num_src_levs = scorpio::get_dimlen(m_datafiles[0],"lev");
    scorpio::eam_pio_closefile(m_datafiles[0]);
  } else {
    scorpio::register_file(m_static_vertical_pressure_file,scorpio::Read);
    m_num_src_levs = scorpio::get_dimlen(m_static_vertical_pressure_file,"lev");
    scorpio::eam_pio_closefile(m_static_vertical_pressure_file);
  }
}
// =========================================================================================
void Nudging::apply_tendency(Field& base, const Field& next, const Real dt)
{
  // Calculate the weight to apply the tendency
  const Real dtend = dt/Real(m_timescale);
  EKAT_REQUIRE_MSG(dtend>=0,"Error! Nudging::apply_tendency - timescale tendency of " << std::to_string(dt) 
		  << " / " << std::to_string(m_timescale) << " = " << std::to_string(dtend) 
		  << " is invalid.  Please check the timescale and/or dt");
  // Now apply the tendency.
  Field tend = base.clone();
  // Use update internal to set tendency, will be (1.0*next - 1.0*base), note tend=base at this point.
  tend.update(next,Real(1.0),Real(-1.0));
  base.update(tend,dtend,Real(1.0));
}
// =========================================================================================
void Nudging::apply_weighted_tendency(Field& base, const Field& next, const Field& weights, const Real dt)
{
  // Calculate the weight to apply the tendency
  const Real dtend = dt/Real(m_timescale);
  EKAT_REQUIRE_MSG(dtend>=0,"Error! Nudging::apply_tendency - timescale tendency of " << std::to_string(dt)
                  << " / " << std::to_string(m_timescale) << " = " << std::to_string(dtend)
                  << " is invalid.  Please check the timescale and/or dt");
  // Now apply the tendency.
  Field tend = base.clone();

  // Use update internal to set tendency, will be (weights*next - weights*base), note tend=base at this point.
  auto base_view = base.get_view<const Real**>();
  auto tend_view = tend.get_view<      Real**>();
  auto next_view = next.get_view<      Real**>();
  auto w_view    = weights.get_view<   Real**>();

  const int num_cols       = base_view.extent(0);
  const int num_vert_packs = base_view.extent(1);
  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {num_cols, num_vert_packs}), KOKKOS_LAMBDA(int i, int j) {
    tend_view(i,j) = next_view(i,j)*w_view(i,j) - base_view(i,j)*w_view(i,j);
  });
  base.update(tend, dtend, Real(1.0));
}
// =============================================================================================================
void Nudging::initialize_impl (const RunType /* run_type */)
{
  using namespace ShortFieldTagsNames;

  // Initialize the time interpolator
  auto grid_ext = m_grid->clone(m_grid->name(), false);
  grid_ext->reset_num_vertical_lev(m_num_src_levs);
  FieldLayout scalar2d_layout_mid { {LEV}, {m_num_src_levs} };
  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols, m_num_src_levs} };
  m_time_interp = util::TimeInterpolation(grid_ext, m_datafiles);

  constexpr int ps = SCREAM_PACK_SIZE;
  const auto& grid_name = m_grid->name();
  if (m_src_pres_type == TIME_DEPENDENT_3D_PROFILE) {
    create_helper_field("p_mid_ext", scalar3d_layout_mid, grid_name, ps);
    auto pmid_ext = get_helper_field("p_mid_ext");
    m_time_interp.add_field(pmid_ext.alias("p_mid"),true);
  } else if (m_src_pres_type == STATIC_1D_VERTICAL_PROFILE) {
    // Load p_levs from source data file
    ekat::ParameterList in_params;
    in_params.set("Filename",m_static_vertical_pressure_file);
    in_params.set("Skip_Grid_Checks",true);  // We need to skip grid checks because multiple ranks may want the same column of source data.
    std::map<std::string,view_1d_host<Real>> host_views;
    std::map<std::string,FieldLayout>  layouts;
    create_helper_field("p_mid_ext", scalar2d_layout_mid, grid_name, ps);
    auto pmid_ext = get_helper_field("p_mid_ext");
    auto pmid_ext_v = pmid_ext.get_view<Real*,Host>();
    in_params.set<std::vector<std::string>>("Field Names",{"p_levs"});
    host_views["p_levs"] = pmid_ext_v;
    layouts.emplace("p_levs",scalar2d_layout_mid);
    AtmosphereInput src_input(in_params,grid_ext,host_views,layouts);
    src_input.read_variables(-1);
    src_input.finalize();
    pmid_ext.sync_to_dev();
  }
  for (auto name : m_fields_nudge) {
    std::string name_ext = name + "_ext";
    // Helper fields that will temporarily store the target state, which can then
    // be used to back out a nudging tendency
    auto field  = get_field_out_wrap(name);
    auto layout = field.get_header().get_identifier().get_layout();
    create_helper_field(name,     layout,              grid_name, ps);
    create_helper_field(name_ext, scalar3d_layout_mid, grid_name, ps);
    auto field_ext = get_helper_field(name_ext);
    m_time_interp.add_field(field_ext.alias(name),true);
  }
  m_time_interp.initialize_data_from_files();

  // load nudging weights from file
  // NOTE: the regional nudging use the same grid as the run, no need to
  // do the interpolation.
  if (m_use_weights)
  {
    FieldLayout scalar3d_layout_grid { {COL,LEV}, {m_num_cols, m_num_levs} };	  
    create_helper_field("nudging_weights", scalar3d_layout_grid, grid_name, ps);
    std::vector<Field> fields;
    auto nudging_weights = get_helper_field("nudging_weights");
    fields.push_back(nudging_weights);
    AtmosphereInput src_weights_input(m_weights_file, grid_ext, fields);
    src_weights_input.read_variables();
    src_weights_input.finalize();
    nudging_weights.sync_to_dev();
  }
}

// =========================================================================================
void Nudging::run_impl (const double dt)
{
  using namespace scream::vinterp;

  // Have to add dt because first time iteration is at 0 seconds where you will
  // not have any data from the field. The timestamp is only iterated at the
  // end of the full step in scream.
  auto ts = timestamp()+dt;

  // Perform time interpolation
  m_time_interp.perform_time_interpolation(ts);

  // Process data and nudge the atmosphere state
  const auto& p_mid_v = get_field_in("p_mid").get_view<const mPack**>();
  view_Nd<mPack,2> p_mid_ext_p;
  view_Nd<mPack,1> p_mid_ext_1d;
  if (m_src_pres_type == TIME_DEPENDENT_3D_PROFILE) {
    p_mid_ext_p = get_helper_field("p_mid_ext").get_view<mPack**>();
  } else if (m_src_pres_type == STATIC_1D_VERTICAL_PROFILE) {
    p_mid_ext_1d   = get_helper_field("p_mid_ext").get_view<mPack*>();
  }

  for (auto name : m_fields_nudge) {
    auto atm_state_field = get_field_out_wrap(name);
    auto int_state_field = get_helper_field(name);
    auto ext_state_field = get_helper_field(name+"_ext");
    auto ext_state_view  = ext_state_field.get_view<mPack**>();
    auto atm_state_view  = atm_state_field.get_view<mPack**>();  // TODO: Right now assume whatever field is defined on COLxLEV
    auto int_state_view  = int_state_field.get_view<mPack**>();
    auto int_mask_view = m_buffer.int_mask_view;
    // Masked values in the source data can lead to strange behavior in the vertical interpolation.
    // We pre-process the data and map any masked values (sometimes called "filled" values) to the
    // nearest un-masked value.
    // Here we are updating the ext_state_view, which is the time interpolated values taken from the nudging
    // data.
    Real var_fill_value = constants::DefaultFillValue<Real>().value;
    // Query the helper field for the fill value, if not present use default
    if (ext_state_field.get_header().has_extra_data("mask_value")) {
      var_fill_value = ext_state_field.get_header().get_extra_data<float>("mask_value");
    }
    const int num_cols           = ext_state_view.extent(0);
    const int num_vert_packs     = ext_state_view.extent(1);
    const int num_src_levs       = m_num_src_levs;
    const auto policy = ESU::get_default_team_policy(num_cols, num_vert_packs);
    Kokkos::parallel_for("correct_for_masked_values", policy,
       	       KOKKOS_LAMBDA(MemberType const& team) {
      const int icol = team.league_rank();
      auto ext_state_view_1d = ekat::subview(ext_state_view,icol);
      Real fill_value;
      int  fill_idx = -1;
      // Scan top to surf and backfill all values near TOM that are masked.
      for (int kk=0; kk<num_src_levs; ++kk) {
        const auto ipack = kk / mPack::n;
	const auto iidx  = kk % mPack::n;
        // Check if this index is masked
	if (ext_state_view_1d(ipack)[iidx]!=var_fill_value) {
	  fill_value = ext_state_view_1d(ipack)[iidx];
	  fill_idx = kk;
	  for (int jj=0; jj<kk; ++jj) {
            const auto jpack = jj / mPack::n;
	    const auto jidx  = jj % mPack::n;
	    ext_state_view_1d(jpack)[jidx] = fill_value;
	  }
	  break;
	}
      }
      // Now fill the rest, the fill_idx should be non-negative.  If it isn't that means
      // we have a column that is fully masked 
      for (int kk=fill_idx+1; kk<num_src_levs; ++kk) {
        const auto ipack = kk / mPack::n;
	const auto iidx  = kk % mPack::n;
        // Check if this index is masked
	if (!(ext_state_view_1d(ipack)[iidx]==var_fill_value)) {
	  fill_value = ext_state_view_1d(ipack)[iidx];
	} else {
	  ext_state_view_1d(ipack)[iidx] = fill_value;
	}
      }
    });

    // Vertical Interpolation onto atmosphere state pressure levels
    if (m_src_pres_type == TIME_DEPENDENT_3D_PROFILE) {
      perform_vertical_interpolation<Real,1,2>(p_mid_ext_p,
                                               p_mid_v,
                                               ext_state_view,
                                               int_state_view,
                                               int_mask_view,
                                               m_num_src_levs,
                                               m_num_levs);
    } else if (m_src_pres_type == STATIC_1D_VERTICAL_PROFILE) {
      perform_vertical_interpolation<Real,1,2>(p_mid_ext_1d,
                                               p_mid_v,
                                               ext_state_view,
                                               int_state_view,
                                               int_mask_view,
                                               m_num_src_levs,
                                               m_num_levs);
    }

    // Check that none of the nudging targets are masked, if they are, set value to
    // nearest unmasked value above.
    // NOTE: We use an algorithm whichs scans from TOM to the surface.
    //       If TOM is masked we keep scanning until we hit an unmasked value,
    //       we then set all masked values above to the unmasked value.
    //       We continue scanning towards the surface until we hit an unmasked value, we
    //       then assign that masked value the most recent unmasked value, until we hit the
    //       surface.
    // Here we change the int_state_view which represents the vertically interpolated fields onto
    // the simulation grid.
    const int num_levs = m_num_levs;
    Kokkos::parallel_for("correct_for_masked_values", policy,
       	       KOKKOS_LAMBDA(MemberType const& team) {
      const int icol = team.league_rank();
      auto int_mask_view_1d  = ekat::subview(int_mask_view,icol);
      auto int_state_view_1d = ekat::subview(int_state_view,icol);
      Real fill_value;
      int  fill_idx = -1;
      // Scan top to surf and backfill all values near TOM that are masked.
      for (int kk=0; kk<num_levs; ++kk) {
        const auto ipack = kk / mPack::n;
	const auto iidx  = kk % mPack::n;
        // Check if this index is masked
	if (!int_mask_view_1d(ipack)[iidx]) {
	  fill_value = int_state_view_1d(ipack)[iidx];
	  fill_idx = kk;
	  for (int jj=0; jj<fill_idx; ++jj) {
            const auto jpack = jj / mPack::n;
	    const auto jidx  = jj % mPack::n;
	    int_state_view_1d(jpack)[jidx] = fill_value;
	  }
	  break;
	}
      }
      // Now fill the rest, the fill_idx should be non-negative.  If it isn't that means
      // we have a column that is fully masked
      for (int kk=fill_idx+1; kk<num_levs; ++kk) {
        const auto ipack = kk / mPack::n;
	const auto iidx  = kk % mPack::n;
        // Check if this index is masked
	if (!int_mask_view_1d(ipack)[iidx]) {
	  fill_value = int_state_view_1d(ipack)[iidx];
	} else {
	  int_state_view_1d(ipack)[iidx] = fill_value;
	}
      }
    });

    // Apply the nudging tendencies to the ATM state
    if (m_timescale <= 0) {
      // We do direct replacement
      Kokkos::deep_copy(atm_state_view,int_state_view);
    } else {
      // Back out a tendency and apply it.
      if (m_use_weights) {
        // get nudging weights field
        // NOTES: do we really need the vertical interpolation for nudging weights? Since we are going to 
        //        use the same grids as the case by providing the nudging weights file.
        //        I would not apply the vertical interpolation here, but it depends...
        //
        auto nudging_weights_field = get_helper_field("nudging_weights");
        // appply the nudging tendencies to the ATM states
        apply_weighted_tendency(atm_state_field, int_state_field, nudging_weights_field, dt);
      } else {
	 apply_tendency(atm_state_field, int_state_field, dt);
      }      
    }
  }
}

// =========================================================================================
void Nudging::finalize_impl()
{
  m_time_interp.finalize();
}
// =========================================================================================
void Nudging::create_helper_field (const std::string& name,
                                             const FieldLayout& layout,
                                             const std::string& grid_name,
                                             const int ps)
{
  using namespace ekat::units;
  // For helper fields we don't bother w/ units, so we set them to non-dimensional
  FieldIdentifier id(name,layout,Units::nondimensional(),grid_name);

  // Create the field. Init with NaN's, so we spot instances of uninited memory usage
  Field f(id);
  if (ps>=0) {
    f.get_header().get_alloc_properties().request_allocation(ps);
  } else {
    f.get_header().get_alloc_properties().request_allocation();
  }
  f.allocate_view();
  f.deep_copy(ekat::ScalarTraits<Real>::invalid());

  m_helper_fields[name] = f;
}

// =========================================================================================
size_t Nudging::requested_buffer_size_in_bytes() const {
  return m_buffer.num_2d_midpoint_mask_views*m_num_cols*m_num_levs*sizeof(mMask);
}

// =========================================================================================
void Nudging::init_buffers(const ATMBufferManager& buffer_manager) {
  EKAT_REQUIRE_MSG(buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(),
                   "Error, Nudging::init_buffers! Buffers size not sufficient.\n");
  mMask* mem = reinterpret_cast<mMask*>(buffer_manager.get_memory());

  m_buffer.int_mask_view = decltype(m_buffer.int_mask_view)(mem,m_num_cols,m_num_levs);
  mem += m_buffer.int_mask_view.size();

  size_t used_mem = (reinterpret_cast<Real*>(mem) - buffer_manager.get_memory())*sizeof(Real);
  EKAT_REQUIRE_MSG(used_mem==requested_buffer_size_in_bytes(), "Error: Nudging::init_buffers! Used memory != requested memory.");
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
