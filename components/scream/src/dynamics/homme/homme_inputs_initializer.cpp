#include "dynamics/homme/homme_inputs_initializer.hpp"

#include "dynamics/homme//interface/scream_homme_interface.hpp"
#include "dynamics/homme/homme_dimensions.hpp"

#include "Hommexx_Session.hpp"
#include "ElementsState.hpp"
#include "Tracers.hpp"
#include "ElementsForcing.hpp"
#include "Types.hpp"

namespace scream {

void HommeInputsInitializer::initialize_fields () {

  // We cheat a little bit here. It is way easier for Homme to init
  // fields during the whole prim initialization sequence.
  // However, such sequence happens during calls to HommeDynamics'
  // ctor and initialize method. The AD attempts to init the
  // inputs *after* all atm process have been inited.
  // So here's what we do: if the user asks to init inputs,
  // we simply go ahead and to it, by calling this method inside
  // HommeDynamics' initialize call. When later the AD calls this
  // initializer, we simply skip everything, since we already inited
  // all the fields.
  // Note: if some other atm proc is trying to init any of the fields
  //       inited by this class, we the AD should throw an error anyways.

  if (!m_fields_inited) {
    // Sanity check. This class should init *all* homme's inputs:
    // FM, FVTheta, FPHI, FQ, v, dp, w, phi, vtheta, q, qdp
    // auto count = m_names.count("v")
    //            + m_names.count("w")
    //            + m_names.count("dp")
    //            + m_names.count("vtheta")
    //            + m_names.count("phi")
    //            + m_names.count("qdp")
    //            + m_names.count("q")
    //            + m_names.count("fm")
    //            + m_names.count("fq")
    //            + m_names.count("fphi")
    //            + m_names.count("ft")
    //            + m_names.count("fvtheta");
    // EKAT_REQUIRE_MSG (count==12, "Error! Homme should init *all* input, not just a few.\n");

    // Compute initial conditions in f90, then copy all data to the
    // C++ elements views. Notice that this does *not( reallocate
    // views; it simply writes over them. Well, it actuall does
    // reallocate some views in Diagnostics, but those views are private
    // to the Diagnostics class (unlike ElementsStates, which share the
    // same pointer with the Scream fields!).
    prim_set_test_initial_conditions_f90 ();
  }
}

void HommeInputsInitializer::add_field (const field_type& f) {
  // We don't really need to store the field, since we can access all views
  // from Homme's data structures. But we need to store the fields
  // identifiers, since we need to expose them to the AD.
  const auto& id = f.get_header().get_identifier();
  m_fids.insert(id);
  m_names.insert(id.name());
}

} // namespace scream
