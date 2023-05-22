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
  using field_type      = Field;
  using identifier_type = typename field_type::identifier_type;
  using layout_type     = typename identifier_type::layout_type;
  using grid_type       = AbstractGrid;
  using grid_ptr_type   = std::shared_ptr<const grid_type>;
  using ci_string       = ekat::CaseInsensitiveString;

  AbstractRemapper () = default;
  AbstractRemapper (const grid_ptr_type& src_grid,
                    const grid_ptr_type& tgt_grid);

  virtual ~AbstractRemapper () = default;

  // Call this before you begin registering fields with this remapper.
  void registration_begins ();

  // This method registers a source field to be remapped to a target field.
  void register_field (const field_type& src, const field_type& tgt);

  // This method registers a source field to be remapped to a target field
  // using fields associated with the given field identifiers.
  void register_field (const identifier_type& src, const identifier_type& tgt);

  // Like the above, but figure out tgt using create_tgt_fid
  void register_field_from_src (const identifier_type& src) {
    register_field(src,create_tgt_fid(src));
  }

  // Like the above, but figure out src using create_src_fid
  void register_field_from_tgt (const identifier_type& tgt) {
    register_field(create_src_fid(tgt),tgt);
  }

  // Call this to indicate that field registration is complete.
  void registration_ends ();

  // The user is allowed to only provide identifiers in the registration phase.
  // In that case, fields have to be bound after registration is complete, and
  // before any call to remap.
  void bind_field (const field_type& src, const field_type& tgt);

  RepoState get_state () const { return m_state; }

  // The actual remap routine.
  void remap (const bool forward);

  // Getter methods
  grid_ptr_type get_src_grid () const { return m_src_grid; }
  grid_ptr_type get_tgt_grid () const { return m_tgt_grid; }

  // Returns the source field identifier for a specified field index (assigned
  // during field registration).
  const identifier_type& get_src_field_id (const int ifield) const {
    EKAT_REQUIRE_MSG(ifield>=0 && ifield<m_num_registered_fields,
                       "Error! Field index out of bounds.\n");
    return do_get_src_field_id(ifield);
  }

  // Returns the target field identifier for a specified field index (assigned
  // during field registration).
  const identifier_type& get_tgt_field_id (const int ifield) const {
    EKAT_REQUIRE_MSG(ifield>=0 && ifield<m_num_registered_fields,
                       "Error! Field index out of bounds.\n");
    return do_get_tgt_field_id(ifield);
  }

  // Returns the source field for the given field index.
  const field_type& get_src_field (const int ifield) const {
    EKAT_REQUIRE_MSG(m_state==RepoState::Closed,
                       "Error! Cannot call 'get_src_field' until registration has ended.\n");
    EKAT_REQUIRE_MSG(ifield>=0 && ifield<m_num_registered_fields,
                       "Error! Field index out of bounds.\n");
    return do_get_src_field(ifield);
  }

  // Returns the target field for the given field index.
  const field_type& get_tgt_field (const int ifield) const {
    EKAT_REQUIRE_MSG(m_state==RepoState::Closed,
                       "Error! Cannot call 'get_tgt_field' until registration has ended.\n");
    EKAT_REQUIRE_MSG(ifield>=0 && ifield<m_num_registered_fields,
                       "Error! Field index out of bounds.\n");
    return do_get_tgt_field(ifield);
  }

  virtual FieldLayout create_src_layout (const FieldLayout& tgt_layout) const = 0;
  virtual FieldLayout create_tgt_layout (const FieldLayout& src_layout) const = 0;

  FieldIdentifier create_src_fid (const FieldIdentifier& tgt_fid) const {
    EKAT_REQUIRE_MSG (tgt_fid.get_grid_name()==m_tgt_grid->name(),
        "Error! Input FieldIdentifier has the wrong grid name:\n"
        "   - input tgt fid grid name: " + tgt_fid.get_grid_name() + "\n"
        "   - remapper tgt grid name:  " + m_tgt_grid->name() + "\n");

    const auto& name = tgt_fid.name();
    const auto& layout = create_src_layout(tgt_fid.get_layout());
    const auto& units = tgt_fid.get_units();

    return FieldIdentifier(name,layout,units,m_src_grid->name());
  }

  FieldIdentifier create_tgt_fid (const FieldIdentifier& src_fid) const {
    EKAT_REQUIRE_MSG (src_fid.get_grid_name()==m_src_grid->name(),
        "Error! Input FieldIdentifier has the wrong grid name:\n"
        "   - input src fid grid name: " + src_fid.get_grid_name() + "\n"
        "   - remapper src grid name:  " + m_src_grid->name() + "\n");

    const auto& name = src_fid.name();
    const auto& layout = create_tgt_layout(src_fid.get_layout());
    const auto& units = src_fid.get_units();

    return FieldIdentifier(name,layout,units,m_tgt_grid->name());
  }

  bool has_src_field (const identifier_type& fid) const {
    for (int i=0; i<m_num_registered_fields; ++i) {
      if (get_src_field_id(i) == fid) {
        return true;
      }
    }
    return false;
  }

  bool has_tgt_field (const identifier_type& fid) const {
    for (int i=0; i<m_num_registered_fields; ++i) {
      if (get_tgt_field_id(i) == fid) {
        return true;
      }
    }
    return false;
  }

  int get_num_fields () const {
    EKAT_REQUIRE_MSG(m_state!=RepoState::Open,
      "Error! Cannot call 'get_num_fields' during the registration phase.\n"
      "       This number is set at 'registration_ends' time.\n"
      " Note: you can call 'num_registered_fields' and 'num_bound_fields' though.\n");
    return m_num_fields;
  }
  int get_num_registered_fields () const {
    return m_num_registered_fields;
  }
  int get_num_bound_fields () const {
    return m_num_bound_fields;
  }

  virtual bool compatible_layouts (const layout_type& src,
                                   const layout_type& tgt) const {
    // By default, the only compatible layouts are identical
    return src==tgt;
  }

protected:

  void set_grids (const grid_ptr_type& src_grid,
                  const grid_ptr_type& tgt_grid);

  virtual const identifier_type& do_get_src_field_id (const int ifield) const = 0;
  virtual const identifier_type& do_get_tgt_field_id (const int ifield) const = 0;
  virtual const field_type& do_get_src_field (const int ifield) const = 0;
  virtual const field_type& do_get_tgt_field (const int ifield) const = 0;

  // Override this method to insert logic executed at the beginning of the field
  // registration process.
  virtual void do_registration_begins () = 0;

  // Override this method to insert logic executed when a pair of source/target
  // field identifiers are registered.
  virtual void do_register_field (const identifier_type& src, const identifier_type& tgt) = 0;

  // Override this method to insert logic executed when a pair of source/target
  // fields are bound to an index within the remapper.
  virtual void do_bind_field (const int ifield, const field_type& src, const field_type& tgt) = 0;

  // Override this method to insert logic executed when the field registration
  // process ends.
  virtual void do_registration_ends () = 0;

  // Override this method to implement the forward remapping process using
  // the protected data members of this class.
  virtual void do_remap_fwd () = 0;

  // Override this method to implement the backward/inverse remapping process
  // using the protected data members of this class.
  virtual void do_remap_bwd () = 0;

  // This helper function searches for a pair of source and target field
  // identifiers in the remapper's set of bound fields. It returns -1 if fields
  // were not registered at all (via the field identifier). If they were
  // registered, it returns the index, regardless of whether the fields have
  // already been bound.
  int find_field (const identifier_type& src,
                  const identifier_type& tgt) {
    int ifield = -1;
    for (int i=0; i<m_num_registered_fields; ++i) {
      if (src==get_src_field_id(i) &&
          tgt==get_tgt_field_id(i)) {
        ifield = i;
        break;
      }
    }
    return ifield;
  }

  // These flags will be updated once all fields are bound.
  // By default, assume both are allowed. During bind_field, if source/target
  // field is read-only, then bwd/fwd is NOT allowed
  bool m_fwd_allowed = true;
  bool m_bwd_allowed = true;

  // The state of the remapper
  RepoState     m_state = RepoState::Clean;

  // The grids associated with the src and tgt fields
  grid_ptr_type m_src_grid;
  grid_ptr_type m_tgt_grid;

  // The number of fields to remap, and the number of fields currently registered.
  // The latter is guaranteed to be equal to the former only when registration is
  // not undergoing. During registration, m_num_fields=0<=m_num_registered_fields.
  int           m_num_fields = 0;
  int           m_num_registered_fields = 0;

  // This vector maps the indices of registered fields to booleans that indicate
  // whether these fields have been bound to the remapper. This vector is
  // necessary because the binding of fields is separate from their
  // regіstration: recall that one may register a field using its identifier,
  // and bind the actual field later.
  // NOTE: vector<bool> is a strange beast, and doesn't necessarily behave as
  // expected. Use caution when manipulating this member, and don't rely on the
  // usual assumptions about how the boolean elements are stored. In particular,
  // avoid range for loops and iterators. Instead, use op[] when you need to
  // access its entries.
  std::vector<bool>   m_fields_are_bound;
  int                 m_num_bound_fields = 0;
};

// A short name for an AbstractRemapper factory
using RemapperFactory =
    ekat::Factory<AbstractRemapper,
                  ekat::CaseInsensitiveString,
                  std::shared_ptr<AbstractRemapper>,
                  const ekat::ParameterList&>;

} // namespace scream

#endif // SCREAM_ABSTRACT_REMAPPER_HPP
