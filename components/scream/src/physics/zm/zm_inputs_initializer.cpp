#include "physics/zm/zm_inputs_initializer.hpp"
#include "physics/zm/scream_zm_interface.hpp"

#include <array>
#include <string>
#include <typeinfo>
#include <string>

namespace scream
{

using namespace std;
std::vector<string> zm_inputs = {"limcnv_in", "no_deep_pbl_in","lchnk", "ncol", "t", "qh", "prec", 
			"jctop", "jcbot", "pblh", "zm", "geos", "zi", "qtnd", "heat", 
			"pap", "paph", "dpp", "delt", "mcon", "cme", "cape", "tpert",
			"dlf", "pflx", "zdu", "rprd", "mu", "md", "du", "eu", "ed", 
			"dp", "dsubcld", "jt", "maxg", "ideep", "lengath", "ql", 
			"rliq", "landfrac", "hu_nm1", "cnv_nm1", "tm1", "qm1", "t_star", 
			"q_star", "dcape", "q", "tend_s", "tend_q", "cld", "snow", 
			"ntprprd", "ztodt", "ntsnprd", "flxprec", "flxsnow", "pguall", 
			"pgdall", "icwu", "ncnst", "fracis"};


const int INPUT_SIZE = 63;


//Layout options are set as an int 
//to be passed into the GridOpts struct
const int SCALAR_3D_MID = 0;
const int SCALAR_3D_INT = 1;
const int VECTOR_3D_MID = 2;
const int TRACERS = 3;
const int LINEAR = 4;

using namespace std;
using namespace scream;
using namespace ekat;
using namespace units;

struct GridOpts{
  string name;
  bool isOut;
  const Units* unit;
  int field_idx;
};

unordered_map<string, GridOpts> opt_map;

//Initializes struct GridOpts fields
void set_grid_opts_helper(GridOpts O, string n, bool out, const Units* unit, int field_idx)
{
  O.name = n;
  O.isOut = out;
  O.field_idx = field_idx;
  O.unit = unit;
  opt_map.insert({O.name, O});
}

void set_grid_opts(){

  // Declares the field structs 
  GridOpts limcnv_in;
  GridOpts no_deep_pbl_in;
  GridOpts lchnk;
  GridOpts ncol;
  GridOpts t;
  GridOpts qh;
  GridOpts prec;
  GridOpts jctop;
  GridOpts jcbot;
  GridOpts pblh;
  GridOpts zm;
  GridOpts geos;
  GridOpts zi;
  GridOpts qtnd;
  GridOpts heat;
  GridOpts pap;
  GridOpts paph;
  GridOpts dpp;
  GridOpts delt;
  GridOpts mcon;
  GridOpts cme;
  GridOpts cape;
  GridOpts tpert;
  GridOpts dlf;
  GridOpts pflx;
  GridOpts zdu;
  GridOpts rprd;
  GridOpts mu;
  GridOpts md;
  GridOpts du;
  GridOpts eu;
  GridOpts ed;
  GridOpts dp;
  GridOpts dsubcld;
  GridOpts jt;
  GridOpts maxg;
  GridOpts ideep;
  GridOpts lengath;
  GridOpts ql;
  GridOpts rliq;
  GridOpts landfrac;
  GridOpts hu_nm1;
  GridOpts cnv_nm1;
  GridOpts tm1;
  GridOpts qm1;
  GridOpts t_star;
  GridOpts q_star;
  GridOpts dcape;
  GridOpts q;
  GridOpts tend_s;
  GridOpts tend_q;
  GridOpts cld;
  GridOpts snow;
  GridOpts ntprprd;
  GridOpts ntsnprd;
  GridOpts flxprec;
  GridOpts flxsnow;
  GridOpts ztodt;
  GridOpts pguall;
  GridOpts pgdall;
  GridOpts icwu;
  GridOpts ncnst;
  GridOpts fracis;
  
  // Sets the value of the grid opt struct categories
  set_grid_opts_helper(limcnv_in, "limcnv_in", true, NULL, NULL);
  set_grid_opts_helper(no_deep_pbl_in, "no_deep_pbl_in", true, NULL, VECTOR_3D_MID);
  set_grid_opts_helper(lchnk, "lchnk", true, &Pa, SCALAR_3D_MID); //temperature(K)
  set_grid_opts_helper(ncol, "ncol", true, &Pa, SCALAR_3D_MID); //temperature(K)
  set_grid_opts_helper(t, "t", true, &Pa, SCALAR_3D_MID); //temperature(K)
  set_grid_opts_helper(qh, "qh", true, NULL, SCALAR_3D_MID);
  set_grid_opts_helper(prec, "prec", true, NULL, LINEAR);
  set_grid_opts_helper(jctop, "jctop", true, NULL, LINEAR);
  set_grid_opts_helper(jcbot, "jcbot", true, NULL, LINEAR);
  set_grid_opts_helper(pblh, "pblh", true, NULL, LINEAR);
  set_grid_opts_helper(zm, "zm", true, NULL, SCALAR_3D_MID);
  set_grid_opts_helper(geos, "geos", true, NULL, LINEAR);
  set_grid_opts_helper(zi, "zi", true, NULL, SCALAR_3D_INT);
  set_grid_opts_helper(qtnd, "qtnd", true, NULL, SCALAR_3D_MID);
  set_grid_opts_helper(heat, "heat", true, NULL, SCALAR_3D_MID);
  set_grid_opts_helper(pap, "pap", true, NULL, SCALAR_3D_MID);
  set_grid_opts_helper(paph, "paph", true, NULL, SCALAR_3D_INT);
  set_grid_opts_helper(dpp, "dpp", true, NULL, SCALAR_3D_MID);
  set_grid_opts_helper(delt, "delt", true, NULL, SCALAR_3D_MID);
  set_grid_opts_helper(mcon, "mcon", true, NULL, SCALAR_3D_MID);
  set_grid_opts_helper(cme, "cme", true, NULL, SCALAR_3D_MID);
  set_grid_opts_helper(cape, "cape", true, NULL, LINEAR);
  set_grid_opts_helper(tpert, "tpert", true, NULL, LINEAR);
  set_grid_opts_helper(dlf, "dlf", true, NULL, SCALAR_3D_MID);
  set_grid_opts_helper(pflx, "pflx", true, NULL, SCALAR_3D_MID);
  set_grid_opts_helper(zdu, "zdu", true, NULL, SCALAR_3D_MID);
  set_grid_opts_helper(rprd, "rprd", true, NULL, SCALAR_3D_MID);
  set_grid_opts_helper(mu, "mu", true, NULL, SCALAR_3D_MID);
  set_grid_opts_helper(md, "md", true, NULL, SCALAR_3D_MID);
  set_grid_opts_helper(du, "du", true, NULL, SCALAR_3D_MID);
  set_grid_opts_helper(eu, "eu", true, NULL, SCALAR_3D_MID);
  set_grid_opts_helper(ed, "ed", true, NULL, SCALAR_3D_MID);
  set_grid_opts_helper(dp, "dp", true, NULL, SCALAR_3D_MID);
  set_grid_opts_helper(dsubcld, "dsubcld", true, NULL, LINEAR);
  set_grid_opts_helper(jt, "jt", true, NULL, LINEAR);
  set_grid_opts_helper(maxg, "maxg", true, NULL, LINEAR);
  set_grid_opts_helper(ideep, "ideep", true, NULL, LINEAR);
  set_grid_opts_helper(lengath, "lengath", true, NULL, SCALAR_3D_MID);
  set_grid_opts_helper(ql, "ql", true, NULL, SCALAR_3D_MID);
  set_grid_opts_helper(rliq, "rliq", true, NULL, SCALAR_3D_MID);
  set_grid_opts_helper(landfrac, "landfrac", true, NULL, LINEAR);
  set_grid_opts_helper(hu_nm1, "hu_nm1", true, NULL, SCALAR_3D_MID);
  set_grid_opts_helper(cnv_nm1, "cnv_nm1", true, NULL, SCALAR_3D_MID);
  set_grid_opts_helper(tm1, "tm1", true, NULL, SCALAR_3D_MID);
  set_grid_opts_helper(qm1, "qm1", true, NULL, SCALAR_3D_MID);
  set_grid_opts_helper(t_star, "t_star", true, NULL, SCALAR_3D_MID);
  set_grid_opts_helper(q_star, "q_star", true, NULL, SCALAR_3D_MID);
  set_grid_opts_helper(dcape, "dcape", true, NULL, LINEAR);
  set_grid_opts_helper(q, "q", true, NULL, SCALAR_3D_MID);
  set_grid_opts_helper(tend_s, "tend_s", true, NULL, LINEAR);
  set_grid_opts_helper(tend_q, "tend_q", true, NULL, LINEAR);
  set_grid_opts_helper(cld, "cld", true, NULL, LINEAR);
  set_grid_opts_helper(snow, "snow", true, NULL, LINEAR);
  set_grid_opts_helper(ntprprd, "ntprprd", true, NULL, SCALAR_3D_MID);
  set_grid_opts_helper(ntsnprd, "ntsnprd", true, NULL, SCALAR_3D_MID);
  set_grid_opts_helper(flxprec, "flxprec", true, NULL, SCALAR_3D_MID);
  set_grid_opts_helper(flxsnow, "flxsnow", true, NULL, SCALAR_3D_MID);
  set_grid_opts_helper(ztodt, "ztodt", true, NULL, SCALAR_3D_MID);
  set_grid_opts_helper(pguall, "pguall", true, NULL, VECTOR_3D_MID);
  set_grid_opts_helper(pgdall, "pgdall", true, NULL, VECTOR_3D_MID);
  set_grid_opts_helper(icwu, "icwu", true, NULL, VECTOR_3D_MID);
  set_grid_opts_helper(ncnst, "ncnst", true, NULL, VECTOR_3D_MID);
  set_grid_opts_helper(fracis, "fracis", true, NULL, VECTOR_3D_MID);
} 

void ZMInputsInitializer::add_field (const field_type &f)
{
  const auto& id = f.get_header().get_identifier();
  
  m_fields.emplace(id.name(),f);
  m_fields_id.insert(id);
}

void ZMInputsInitializer::
add_field (const field_type &f, const field_type& f_ref,
           const remapper_ptr_type& remapper)
{
  if (m_remapper) {
    // Sanity check
    EKAT_REQUIRE_MSG (m_remapper->get_src_grid()->name()==remapper->get_src_grid()->name(),
      "Error! A remapper was already set in P3InputsInitializer, but its src grid differs from"
      "       the grid of the input remapper of this call.\n");
  } else {
    m_remapper = remapper;
    m_remapper->registration_begins();
  }

  const auto& id = f.get_header().get_identifier();
  const auto& id_ref = f_ref.get_header().get_identifier();

  // To the AD, we only expose the fact that we init f_ref...
  m_fields_id.insert(id_ref);

  // ...but P3 only knows how to init f...
  m_fields.emplace(id.name(),f);

  // ...hence, we remap to f_ref.
  m_remapper->register_field(f, f_ref);
}

void ZMInputsInitializer :: initialize_fields(){

  //check if zm inputs have been registered as fields
  int count = 0;
  for (size_t j = 0; j < zm_inputs.size(); j++){
    count += m_fields.count(zm_inputs[j]);
    if (count==0){
      return;
    }
  }
  

    EKAT_REQUIRE_MSG (count==INPUT_SIZE,
    "Error! ZMInputsInitializer is expected to init _______________.\n"
    "Only " + std::to_string(count) + " of those have been found.\n"
    "Please, check the atmosphere processes you are using,"
    "and make sure they agree on who's initializing each field.\n");



  for (size_t i = 0; i < zm_inputs.size(); i++){
    //Get and store device view using input name
    Kokkos::View<Real*, Kokkos::LayoutRight, DefaultDevice> d_v = m_fields.at(zm_inputs[i]).get_view();
    //Create and store host mirrors using device views
    Kokkos::View<scream::Real*, Kokkos::LayoutRight, HostDevice> h_m = Kokkos::create_mirror_view(d_v);
    //Create and store host mirrors raw pointers
    // Real* r_p = h_m.data();
    //Deep copy back to device
    Kokkos::deep_copy(d_v, h_m);
  }
}


} // namespace scream


