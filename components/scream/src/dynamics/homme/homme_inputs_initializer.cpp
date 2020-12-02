#include "dynamics/homme/homme_inputs_initializer.hpp"

#include "dynamics/homme//interface/scream_homme_interface.hpp"

namespace scream {

void HommeInputsInitializer::initialize_fields () {

  // We cheat a little bit here. It is way easier for Homme to init
  // fields during the whole prim initialization sequence.
  // However, such sequence happens during calls to HommeDynamics'
  // ctor and initialize method. The AD attempts to init the
  // inputs *after* all atm process have been inited.
  // So here's what we do: if the user asks to init inputs,
  // we simply go ahead and do it, by calling this method inside
  // HommeDynamics' initialize call. When later the AD calls this
  // initializer, we simply skip everything, since we already inited
  // all the fields.
  // Note: if some other atm proc is trying to init any of the fields
  //       inited by this class, the AD should throw an error anyways.

  if (!m_fields_inited) {
    // TODO: add checks to ensure we're in charge of initialization for *all*
    //       the dynamics inputs (otherwise there may be some inconsistencies).

    // Compute initial conditions in f90, then copy all data to the
    // C++ elements views. Notice that this does *not( reallocate
    // views; it simply writes over them. Well, it actuall does
    // reallocate some views in Diagnostics, but those views are private
    // to the Diagnostics class (unlike ElementsStates, which share the
    // same pointer with the Scream fields!).
    prim_set_test_initial_conditions_f90 ();

    if (m_remapper) {
      // Copy inited fields to the ref grid, then dispose of the remapper.
      m_remapper->registration_ends();
      m_remapper->remap(true);
      m_remapper = nullptr;
    }
  }
}

void HommeInputsInitializer::add_field (const field_type& f) {
  // We don't really need to store the field, since we can access all views
  // from Homme's data structures. But we need to store the fields
  // identifiers, since we need to expose them to the AD.
  const auto& id = f.get_header().get_identifier();
  m_fids.insert(id);
}

void HommeInputsInitializer::
add_field (const field_type& f, const field_type& f_ref,
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

  // We don't really need to store the field, since we can access all views
  // from Homme's data structures. But we need to store the fields
  // identifiers, since we need to expose them to the AD.
  // Also, the remapper needs to store the fields
  const auto& id = f.get_header().get_identifier();
  m_fids.insert(id);

  const auto& id_ref = f_ref.get_header().get_identifier();

  // To the AD, we only expose the fact that we init f_ref...
  m_fids.insert(id_ref);

  // ...but homme only knows how to init f,
  // hence, we remap to f_ref.
  m_remapper->register_field(f, f_ref);
}

} // namespace scream
