#ifndef SCREAM_ABSTRACT_REMAPPER_HPP
#define SCREAM_ABSTRACT_REMAPPER_HPP

#include "share/field/field.hpp"
#include "share/grid/abstract_grid.hpp"

#include "ekat/util/ekat_factory.hpp"
#include "ekat/util/ekat_string_utils.hpp"
#include "ekat/std_meta/ekat_std_utils.hpp"
#include "ekat/ekat_parameter_list.hpp"

namespace scream
{

// An abstract interface for a remapper

// A remapper is basically a functor, that, given two fields,
// copies the first into the second, or viceversa. The copy must
// account for different layouts and/or different mpi distributions.
// This concept can be extended to remaps that involve interpolation,
// but as of now (07/2019) it is not the intent and projected use
// of this class in the scream framework
class AbstractRemapper
{
public:
  using grid_ptr_type   = std::shared_ptr<const AbstractGrid>;

  AbstractRemapper () = default;
  AbstractRemapper (const grid_ptr_type& src_grid,
                    const grid_ptr_type& tgt_grid);

  virtual ~AbstractRemapper () = default;

  // ----- Registration/setup methods ---- //

  // Call this before you begin registering fields with this remapper.
  void registration_begins ();

  // This method registers a source field to be remapped to a target field.
  void register_field (const Field& src, const Field& tgt);

  // Like the above, but create tgt/src internally
  virtual void register_field_from_src (const Field& src);
  virtual void register_field_from_tgt (const Field& tgt);

  // Call this to indicate that field registration is complete.
  void registration_ends ();


  //  ------- Getter methods ------- //
  RepoState get_state () const { return m_state; }

  int get_num_fields () const { return m_num_fields; }

  grid_ptr_type get_src_grid () const { return m_src_grid; }
  grid_ptr_type get_tgt_grid () const { return m_tgt_grid; }

  const Field& get_src_field (const int i) const;
  const Field& get_tgt_field (const int i) const;

  bool fwd_allowed () const { return m_fwd_allowed; }
  bool bwd_allowed () const { return m_bwd_allowed; }

  // ------ Runtime methods ------- //

  void remap_fwd ();
  void remap_bwd ();

  FieldIdentifier create_src_fid (const FieldIdentifier& tgt_fid) const;
  FieldIdentifier create_tgt_fid (const FieldIdentifier& src_fid) const;

  // These are public, so that InverseRemapper can call them on the stored remapper
  // We provide reasonable default implementations, but derived classes can specialize
  virtual FieldLayout create_src_layout (const FieldLayout& tgt_layout) const;
  virtual FieldLayout create_tgt_layout (const FieldLayout& src_layout) const;
  virtual bool compatible_layouts (const FieldLayout& src, const FieldLayout& tgt) const;
  virtual bool is_valid_src_layout (const FieldLayout& layout) const {
    return m_src_grid->is_valid_layout(layout);
  }
  virtual bool is_valid_tgt_layout (const FieldLayout& layout) const {
    return m_tgt_grid->is_valid_layout(layout);
  }

protected:

  virtual FieldLayout create_layout (const FieldLayout& from_layout,
                                     const grid_ptr_type& to_grid) const;

  void set_grids (const grid_ptr_type& src_grid,
                  const grid_ptr_type& tgt_grid);

  // Override this method to insert logic executed when the field registration
  // process ends.
  virtual void registration_ends_impl () { /* Default to do-nothing */ }

  // Override this method to implement the forward remapping process using
  // the protected data members of this class.
  virtual void remap_fwd_impl () {
    EKAT_ERROR_MSG ("Error! Missing override of remap_fwd_impl for this remapper.\n");
  }

  // Override this method to implement the backward/inverse remapping process
  // using the protected data members of this class.
  virtual void remap_bwd_impl () {
    EKAT_ERROR_MSG ("Error! Missing override of remap_bwd_impl for this remapper.\n");
  }

  // By default, assume both are allowed.
  bool m_fwd_allowed = true;
  bool m_bwd_allowed = true;

  // The state of the remapper
  RepoState     m_state = RepoState::Clean;

  // The grids associated with the src and tgt fields
  grid_ptr_type m_src_grid;
  grid_ptr_type m_tgt_grid;

  // The number of fields to remap
  int           m_num_fields = 0;


  std::vector<Field> m_src_fields;
  std::vector<Field> m_tgt_fields;

private:
  // If source/target is read-only, we set the corresponding var to true.
  // This allow better error msg in case fwd/bwd remap is not allowed.
  bool m_has_read_only_src_fields = false;
  bool m_has_read_only_tgt_fields = false;
};

// A short name for an AbstractRemapper factory
using RemapperFactory =
    ekat::Factory<AbstractRemapper,
                  ekat::CaseInsensitiveString,
                  std::shared_ptr<AbstractRemapper>,
                  const ekat::ParameterList&>;

} // namespace scream

#endif // SCREAM_ABSTRACT_REMAPPER_HPP
