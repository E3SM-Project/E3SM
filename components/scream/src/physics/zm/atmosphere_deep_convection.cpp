#include "physics/zm/scream_zm_interface.hpp"
#include "physics/zm/atmosphere_deep_convection.hpp"
#include "physics/zm/zm_inputs_initializer.cpp"

#include "ekat/ekat_assert.hpp"

namespace scream

{



ZMDeepConvection::ZMDeepConvection (const ekat::Comm& comm,const ekat::ParameterList& params )
  : m_zm_comm (comm)
  , m_zm_params(params)
{
  m_initializer = create_field_initializer<ZMInputsInitializer>();
  
}

void ZMDeepConvection::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace std;
  using namespace ekat;
  using namespace units;

  auto Q = kg/kg;
  auto nondim = m/m;
  Q.set_string("kg/kg");

  constexpr int NVL = 72;  /* TODO THIS NEEDS TO BE CHANGED TO A CONFIGURABLE */
  constexpr int QSZ =  35;  /* TODO THIS NEEDS TO BE CHANGED TO A CONFIGURABLE */
  auto grid = grids_manager->get_grid("Physics");
  const int num_dofs = grid->get_num_local_dofs();
  const int nc = num_dofs;

  using namespace ShortFieldTagsNames;

  FieldLayout scalar3d_layout_mid { {COL,VL}, {nc,NVL} }; // Note that C++ and Fortran read array dimensions in reverse
  FieldLayout scalar3d_layout_int { {COL,VL}, {nc,NVL+1} }; // Note that C++ and Fortran read array dimensions in reverse
  FieldLayout vector3d_layout_mid{ {COL,CMP,VL}, {nc,QSZ,NVL} };
  FieldLayout tracers_layout { {COL,VAR,VL}, {nc,QSZ,NVL} };
  FieldLayout scalar2d_layout{ {COL}, {nc} };
 
  std::vector<FieldLayout> layout_opts = {scalar3d_layout_mid, scalar3d_layout_int,
					vector3d_layout_mid, tracers_layout, scalar2d_layout};
 
  set_grid_opts();
  
  for ( auto i = opt_map.begin(); i != opt_map.end(); ++i) {
    m_required_fields.emplace((i->second).name, layout_opts[((i->second).field_idx)], Q, grid->name());
    if ( (i->second).isOut == true ) {
      m_computed_fields.emplace((i->second).name, layout_opts[((i->second).field_idx)], Q, grid->name());
    }
  }
  
}

// =========================================================================================
void ZMDeepConvection::initialize_impl (const util::TimeStamp& t0)
{
  zm_init_f90 (*m_raw_ptrs_in["limcnv_in"], m_raw_ptrs_in["no_deep_pbl_in"]);

  m_current_ts = t0;

  using strvec = std::vector<std::string>;
  const strvec& allowed_to_init = m_zm_params.get<strvec>("Initializable Inputs",strvec(0));
  const bool can_init_all = m_zm_params.get<bool>("Can Initialize All Inputs", false);
  const bool init_all_or_none = m_zm_params.get<bool>("Must Init All Inputs Or None", true);
  const strvec& initable = can_init_all ? zm_inputs : allowed_to_init;
  
  if (initable.size()>0) {
    bool all_inited = true, all_uninited = true;
    for (const auto& name : initable) {
      const auto& f = m_zm_fields_in.at(name);
      const auto& track = f.get_header().get_tracking();
      if (track.get_init_type()==InitType::None) {
        // Nobody claimed to init this field. P3InputsInitializer will take care of it
        m_initializer->add_me_as_initializer(f);
        all_uninited &= true;
        all_inited &= false;
      } else {
        all_uninited &= false;
        all_inited &= true;
      }
    }
  }
}
// =========================================================================================
void ZMDeepConvection::run_impl (const Real dt)
{
  std::vector<const Real*> in;
  std::vector<Real*> out;

  // Copy inputs to host. Copy also outputs, cause we might "update" them, rather than overwrite them.
  for (auto& it : m_zm_fields_in) {
    Kokkos::deep_copy(m_zm_host_views_in.at(it.first),it.second.get_view());
  }
  for (auto& it : m_zm_fields_out) {
    Kokkos::deep_copy(m_zm_host_views_out.at(it.first),it.second.get_view());
  }

  Real** temp = &m_raw_ptrs_out["fracis"];
  Real*** fracis = &temp;
 
  zm_main_f90(*m_raw_ptrs_out["lchnk"], *m_raw_ptrs_out["ncol"], m_raw_ptrs_out["t"], 
  	      m_raw_ptrs_out["qh"], m_raw_ptrs_out["prec"], m_raw_ptrs_out["jctop"], 
              m_raw_ptrs_out["jcbot"], m_raw_ptrs_out["pblh"], m_raw_ptrs_out["zm"], 
              m_raw_ptrs_out["geos"], m_raw_ptrs_out["zi"], m_raw_ptrs_out["qtnd"], 
              m_raw_ptrs_out["heat"], m_raw_ptrs_out["pap"], m_raw_ptrs_out["paph"], 
              m_raw_ptrs_out["dpp"], *m_raw_ptrs_out["delt"], m_raw_ptrs_out["mcon"], 
	      m_raw_ptrs_out["cme"], m_raw_ptrs_out["cape"], m_raw_ptrs_out["tpert"], 
              m_raw_ptrs_out["dlf"], m_raw_ptrs_out["plfx"], m_raw_ptrs_out["zdu"], 
	      m_raw_ptrs_out["rprd"], m_raw_ptrs_out["mu"], m_raw_ptrs_out["md"], 
              m_raw_ptrs_out["du"], m_raw_ptrs_out["eu"], m_raw_ptrs_out["ed"], 
	      m_raw_ptrs_out["dp"], m_raw_ptrs_out["dsubcld"], m_raw_ptrs_out["jt"], 
	      m_raw_ptrs_out["maxg"], m_raw_ptrs_out["ideep"], *m_raw_ptrs_out["lengath"], 
              m_raw_ptrs_out["ql"], m_raw_ptrs_out["rliq"], m_raw_ptrs_out["landfrac"], 
              m_raw_ptrs_out["hu_nm1"], m_raw_ptrs_out["cnv_nm1"], m_raw_ptrs_out["tm1"], 
              m_raw_ptrs_out["qm1"], &m_raw_ptrs_out["t_star"], &m_raw_ptrs_out["q_star"], 
              m_raw_ptrs_out["dcape"], m_raw_ptrs_out["q"], &m_raw_ptrs_out["tend_s"], 
              &m_raw_ptrs_out["tend_q"], &m_raw_ptrs_out["cld"], m_raw_ptrs_out["snow"], 
              m_raw_ptrs_out["ntprprd"], m_raw_ptrs_out["ntsnprd"], 
              &m_raw_ptrs_out["flxprec"], &m_raw_ptrs_out["flxsnow"], 
              *m_raw_ptrs_out["ztodt"], m_raw_ptrs_out["pguall"], m_raw_ptrs_out["pgdall"], 
              m_raw_ptrs_out["icwu"], *m_raw_ptrs_out["ncnst"], fracis); 
  m_current_ts += dt;
  for (int i = 0; i < zm_inputs.size(); i++){
  	m_zm_fields_out.at(zm_inputs[i]).get_header().get_tracking().update_time_stamp(m_current_ts);
  }
}
// =========================================================================================
void ZMDeepConvection::finalize_impl()
{
  zm_finalize_f90 ();
}
// =========================================================================================
void ZMDeepConvection::register_fields (FieldRepository<Real, device_type>& field_repo) const {
   for (auto& fid : m_required_fields) {
     field_repo.register_field(fid);
   }
   for (auto& fid : m_computed_fields) {
     field_repo.register_field(fid);
   }
 }

void ZMDeepConvection::set_required_field_impl (const Field<const Real, device_type>& f) {
  // @Meredith: Diff between add_me_as_a_customer and get_tracking().add_customer? 
  
  // Store a copy of the field. We need this in order to do some tracking checks
  // at the beginning of the run call. Other than that, there would be really
  // no need to store a scream field here; we could simply set the view ptr
  // in the Homme's view, and be done with it.
  const auto& name = f.get_header().get_identifier().name();
  m_zm_fields_in.emplace(name,f);
  m_zm_host_views_in[name] = Kokkos::create_mirror_view(f.get_view());
  m_raw_ptrs_in[name] = m_zm_host_views_in[name].data();
  // Add myself as customer to the field
  add_me_as_customer(f);

}

void ZMDeepConvection::set_computed_field_impl (const Field<      Real, device_type>& f) {
  // Store a copy of the field. We need this in order to do some tracking updates
  // at the end of the run call. Other than that, there would be really
  // no need to store a scream field here; we could simply set the view ptr
  // in the Homme's view, and be done with it.
  const auto& name = f.get_header().get_identifier().name();
  m_zm_fields_out.emplace(name,f);
  m_zm_host_views_out[name] = Kokkos::create_mirror_view(f.get_view());
  m_raw_ptrs_out[name] = m_zm_host_views_out[name].data();

  // Add myself as provider for the field
  add_me_as_provider(f);
  }




} // namespace scream
