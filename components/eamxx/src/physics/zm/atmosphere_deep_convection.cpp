#include "physics/zm/scream_zm_interface.hpp"
#include "physics/zm/atmosphere_deep_convection.hpp"

#include "ekat/ekat_assert.hpp"

namespace {
// A helper struct and fcn;
struct GridOpts {
  GridOpts(const std::string& n, bool out, const ekat::units::Units& u, int idx)
    : name(n), isOut(out), unit(u), field_idx(idx) {}
  std::string name;
  bool isOut;
  const ekat::units::Units& unit;
  int field_idx;
};

void set_grid_opts(std::map<std::string, GridOpts>& opt_map);
}

namespace scream
{

ZMDeepConvection::ZMDeepConvection (const ekat::Comm& comm,const ekat::ParameterList& params )
  : AtmosphereProcess(comm, params)
{
  // Nothing to do here
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

  FieldLayout scalar3d_layout_mid { {COL,LEV}, {nc,NVL} }; // Note that C++ and Fortran read array dimensions in reverse
  FieldLayout scalar3d_layout_int { {COL,ILEV}, {nc,NVL+1} }; // Note that C++ and Fortran read array dimensions in reverse
  FieldLayout vector3d_layout_mid{ {COL,CMP,LEV}, {nc,QSZ,NVL} };
  FieldLayout tracers_layout { {COL,CMP,LEV}, {nc,QSZ,NVL} };
  FieldLayout scalar2d_layout{ {COL}, {nc} };

  std::vector<FieldLayout> layout_opts = {scalar3d_layout_mid, scalar3d_layout_int,
					vector3d_layout_mid, tracers_layout, scalar2d_layout};

  std::map<std::string, GridOpts> opt_map;
  set_grid_opts(opt_map);

  for ( auto i = opt_map.begin(); i != opt_map.end(); ++i) {
    add_required_field((i->second).name, layout_opts[((i->second).field_idx)], Q, grid->name());
    if ( (i->second).isOut == true ) {
      add_computed_field((i->second).name, layout_opts[((i->second).field_idx)], Q, grid->name());
    }
  }

}

// =========================================================================================
void ZMDeepConvection::initialize_impl (const RunType /* run_type */)
{
  zm_init_f90 (*m_raw_ptrs_in["limcnv_in"], m_raw_ptrs_in["no_deep_pbl_in"]);
}
// =========================================================================================
void ZMDeepConvection::run_impl (const double dt)
{
  std::vector<const Real*> in;
  std::vector<Real*> out;

  // Copy inputs to host. Copy also outputs, cause we might "update" them, rather than overwrite them.
  for (auto& it : m_zm_fields_in) {
    it.second.sync_to_host();
  }
  for (auto& it : m_zm_fields_out) {
    it.second.sync_to_host();
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
              m_raw_ptrs_out["dcape"], m_raw_ptrs_out["qv"], &m_raw_ptrs_out["tend_s"],
              &m_raw_ptrs_out["tend_q"], &m_raw_ptrs_out["cld"], m_raw_ptrs_out["snow"],
              m_raw_ptrs_out["ntprprd"], m_raw_ptrs_out["ntsnprd"],
              &m_raw_ptrs_out["flxprec"], &m_raw_ptrs_out["flxsnow"],
              *m_raw_ptrs_out["ztodt"], m_raw_ptrs_out["pguall"], m_raw_ptrs_out["pgdall"],
              m_raw_ptrs_out["icwu"], *m_raw_ptrs_out["ncnst"], fracis);
  auto ts = timestamp();
  ts += dt;
  for (auto& it : m_zm_fields_out) {
    it.second.get_header().get_tracking().update_time_stamp(ts);
  }
}
// =========================================================================================
void ZMDeepConvection::finalize_impl()
{
  zm_finalize_f90 ();
}
// =========================================================================================
void ZMDeepConvection::
register_fields (const std::map<std::string,std::shared_ptr<FieldManager<Real>>>& field_mgrs) const {
  auto& field_mgr = *field_mgrs.at("Physics");
   for (const auto& fid : get_required_fields()) {
     field_mgr.register_field(fid);
   }
   for (const auto& fid : get_computed_fields()) {
     field_mgr.register_field(fid);
   }
 }

void ZMDeepConvection::set_required_field_impl (const Field<const Real>& f) {
  // @Meredith: Diff between add_me_as_a_customer and get_tracking().add_customer?

  // Store a copy of the field. We need this in order to do some tracking checks
  // at the beginning of the run call. Other than that, there would be really
  // no need to store a scream field here; we could simply set the view ptr
  // in the Homme's view, and be done with it.
  const auto& name = f.get_header().get_identifier().name();
  m_zm_fields_in.emplace(name,f);
  m_zm_host_views_in[name] = f.get_view<Host>();
  m_raw_ptrs_in[name] = m_zm_host_views_in[name].data();

  // Add myself as customer to the field
  add_me_as_customer(f);

}

void ZMDeepConvection::set_computed_field_impl (const Field<      Real>& f) {
  // Store a copy of the field. We need this in order to do some tracking updates
  // at the end of the run call. Other than that, there would be really
  // no need to store a scream field here; we could simply set the view ptr
  // in the Homme's view, and be done with it.
  const auto& name = f.get_header().get_identifier().name();
  m_zm_fields_out.emplace(name,f);
  m_zm_host_views_out[name] = f.get_view<Host>();
  m_raw_ptrs_out[name] = m_zm_host_views_out[name].data();

  // Add myself as provider for the field
  add_me_as_provider(f);
}

} // namespace scream

namespace {
//Initializes struct GridOpts fields
void set_grid_opts_helper(std::map<std::string, GridOpts>& opt_map,
                          const std::string& n, bool out, int field_idx,
                          const ekat::units::Units& unit = ekat::units::Units::nondimensional())
{
  opt_map.emplace(std::piecewise_construct,
                  std::forward_as_tuple(n),
                  std::forward_as_tuple(n,out,unit,field_idx));
}

void set_grid_opts(std::map<std::string, GridOpts>& opt_map){
  //Layout options are set as an int
  //to be passed into the GridOpts struct
  const int SCALAR_3D_MID = 0;
  const int SCALAR_3D_INT = 1;
  const int VECTOR_3D_MID = 2;
  const int LINEAR = 4;

  const auto Pa = ekat::units::Pa;

  // Sets the value of the grid opt struct categories
  set_grid_opts_helper(opt_map, "limcnv_in",       true, 0);
  set_grid_opts_helper(opt_map, "no_deep_pbl_in",  true, VECTOR_3D_MID);
  set_grid_opts_helper(opt_map, "lchnk",           true, SCALAR_3D_MID, Pa); //temperature(K)
  set_grid_opts_helper(opt_map, "ncol",            true, SCALAR_3D_MID, Pa); //temperature(K)
  set_grid_opts_helper(opt_map, "t",               true, SCALAR_3D_MID, Pa); //temperature(K)
  set_grid_opts_helper(opt_map, "qh",              true, SCALAR_3D_MID);
  set_grid_opts_helper(opt_map, "prec",            true, LINEAR);
  set_grid_opts_helper(opt_map, "jctop",           true, LINEAR);
  set_grid_opts_helper(opt_map, "jcbot",           true, LINEAR);
  set_grid_opts_helper(opt_map, "pblh",            true, LINEAR);
  set_grid_opts_helper(opt_map, "zm",              true, SCALAR_3D_MID);
  set_grid_opts_helper(opt_map, "geos",            true, LINEAR);
  set_grid_opts_helper(opt_map, "zi",              true, SCALAR_3D_INT);
  set_grid_opts_helper(opt_map, "qtnd",            true, SCALAR_3D_MID);
  set_grid_opts_helper(opt_map, "heat",            true, SCALAR_3D_MID);
  set_grid_opts_helper(opt_map, "pap",             true, SCALAR_3D_MID);
  set_grid_opts_helper(opt_map, "paph",            true, SCALAR_3D_INT);
  set_grid_opts_helper(opt_map, "dpp",             true, SCALAR_3D_MID);
  set_grid_opts_helper(opt_map, "delt",            true, SCALAR_3D_MID);
  set_grid_opts_helper(opt_map, "mcon",            true, SCALAR_3D_MID);
  set_grid_opts_helper(opt_map, "cme",             true, SCALAR_3D_MID);
  set_grid_opts_helper(opt_map, "cape",            true, LINEAR);
  set_grid_opts_helper(opt_map, "tpert",           true, LINEAR);
  set_grid_opts_helper(opt_map, "dlf",             true, SCALAR_3D_MID);
  set_grid_opts_helper(opt_map, "pflx",            true, SCALAR_3D_MID);
  set_grid_opts_helper(opt_map, "zdu",             true, SCALAR_3D_MID);
  set_grid_opts_helper(opt_map, "rprd",            true, SCALAR_3D_MID);
  set_grid_opts_helper(opt_map, "mu",              true, SCALAR_3D_MID);
  set_grid_opts_helper(opt_map, "md",              true, SCALAR_3D_MID);
  set_grid_opts_helper(opt_map, "du",              true, SCALAR_3D_MID);
  set_grid_opts_helper(opt_map, "eu",              true, SCALAR_3D_MID);
  set_grid_opts_helper(opt_map, "ed",              true, SCALAR_3D_MID);
  set_grid_opts_helper(opt_map, "dp",              true, SCALAR_3D_MID);
  set_grid_opts_helper(opt_map, "dsubcld",         true, LINEAR);
  set_grid_opts_helper(opt_map, "jt",              true, LINEAR);
  set_grid_opts_helper(opt_map, "maxg",            true, LINEAR);
  set_grid_opts_helper(opt_map, "ideep",           true, LINEAR);
  set_grid_opts_helper(opt_map, "lengath",         true, SCALAR_3D_MID);
  set_grid_opts_helper(opt_map, "ql",              true, SCALAR_3D_MID);
  set_grid_opts_helper(opt_map, "rliq",            true, SCALAR_3D_MID);
  set_grid_opts_helper(opt_map, "landfrac",        true, LINEAR);
  set_grid_opts_helper(opt_map, "hu_nm1",          true, SCALAR_3D_MID);
  set_grid_opts_helper(opt_map, "cnv_nm1",         true, SCALAR_3D_MID);
  set_grid_opts_helper(opt_map, "tm1",             true, SCALAR_3D_MID);
  set_grid_opts_helper(opt_map, "qm1",             true, SCALAR_3D_MID);
  set_grid_opts_helper(opt_map, "t_star",          true, SCALAR_3D_MID);
  set_grid_opts_helper(opt_map, "q_star",          true, SCALAR_3D_MID);
  set_grid_opts_helper(opt_map, "dcape",           true, LINEAR);
  set_grid_opts_helper(opt_map, "qv",              true, SCALAR_3D_MID);
  set_grid_opts_helper(opt_map, "tend_s",          true, LINEAR);
  set_grid_opts_helper(opt_map, "tend_q",          true, LINEAR);
  set_grid_opts_helper(opt_map, "cld",             true, LINEAR);
  set_grid_opts_helper(opt_map, "snow",            true, LINEAR);
  set_grid_opts_helper(opt_map, "ntprprd",         true, SCALAR_3D_MID);
  set_grid_opts_helper(opt_map, "ntsnprd",         true, SCALAR_3D_MID);
  set_grid_opts_helper(opt_map, "flxprec",         true, SCALAR_3D_MID);
  set_grid_opts_helper(opt_map, "flxsnow",         true, SCALAR_3D_MID);
  set_grid_opts_helper(opt_map, "ztodt",           true, SCALAR_3D_MID);
  set_grid_opts_helper(opt_map, "pguall",          true, VECTOR_3D_MID);
  set_grid_opts_helper(opt_map, "pgdall",          true, VECTOR_3D_MID);
  set_grid_opts_helper(opt_map, "icwu",            true, VECTOR_3D_MID);
  set_grid_opts_helper(opt_map, "ncnst",           true, VECTOR_3D_MID);
  set_grid_opts_helper(opt_map, "fracis",          true, VECTOR_3D_MID);
}
} // anonymous namespace

