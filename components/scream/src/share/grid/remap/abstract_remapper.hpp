#ifndef SCREAM_ABSTRACT_REMAPPER_HPP
#define SCREAM_ABSTRACT_REMAPPER_HPP

#include "share/field/field.hpp"
#include "share/field/field_utils.hpp"
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

template<typename RealType>
class AbstractRemapper
{
public:
  using real_type       = RealType;
  using field_type      = Field<real_type>;
  using identifier_type = typename field_type::identifier_type;
  using layout_type     = typename identifier_type::layout_type;
  using grid_type       = AbstractGrid;
  using grid_ptr_type   = std::shared_ptr<const grid_type>;
  using ci_string       = ekat::CaseInsensitiveString;

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

  // This method unregisters source and target fields associated with the given
  // identifiers, indicating that they are no longer to be remapped.
  void unregister_field (const identifier_type& src, const identifier_type& tgt);

  // Call this to indicate that field registration is complete.
  void registration_ends ();

  // The user is allowed to only provide identifiers in the registration phase.
  // In that case, fields have to be bound after registration is complete, and
  // before any call to remap.
  void bind_field (const field_type& src, const field_type& tgt);

  RepoState get_state () const { return m_state; }

  // The actual remap routine.
  void remap (const bool forward) const {
    EKAT_REQUIRE_MSG(m_state!=RepoState::Open,
                       "Error! Cannot perform remapping at this time.\n"
                       "       Did you forget to call 'registration_ends'?\n");

    EKAT_REQUIRE_MSG(m_num_bound_fields==m_num_fields,
                       "Error! Not all fields have been set in the remapper.\n"
                       "       In particular, field " +
                       std::to_string(std::distance(m_fields_are_bound.begin(),std::find(m_fields_are_bound.begin(),m_fields_are_bound.end(),false))) +
                       " has not been bound.\n");

    if (m_state!=RepoState::Clean) {
      if (forward) {
        do_remap_fwd ();
      } else {
        do_remap_bwd ();
      }
    }
  }

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

  // Override this method to insert logic executed when a field is unregistered.
  virtual void do_unregister_field (const int ifield) = 0;

  // Override this method to insert logic executed when the field registration
  // process ends.
  virtual void do_registration_ends () = 0;

  // Override this method to implement the forward remapping process using
  // the protected data members of this class.
  virtual void do_remap_fwd () const = 0;

  // Override this method to implement the backward/inverse remapping process
  // using the protected data members of this class.
  virtual void do_remap_bwd () const = 0;

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

  // The state of the remapper
  RepoState     m_state;

  // The grids associated with the src and tgt fields
  grid_ptr_type m_src_grid;
  grid_ptr_type m_tgt_grid;

  // The number of fields to remap, and the number of fields currently registered.
  // The latter is guaranteed to be equal to the former only when registration is
  // not undergoing. During registration, m_num_fields=0<=m_num_registered_fields.
  int           m_num_fields;
  int           m_num_registered_fields;

  // This vector maps the indices of registered fields to booleans that indicate
  // whether these fields have been bound to the remapper. This vector is
  // necessary because the binding of fields is separate from their
  // regіstration: recall that one may register a field using its identifier,
  // and bind the actual field later.
  // NOTE: vector<bool> is a strange beast, and doesn't necessarily behave as
  // expected. Use caution when manipulating this member, and don't rely on the
  // usual assumptions about how the boolean elements are stored.
  std::vector<bool>   m_fields_are_bound;
  int                 m_num_bound_fields;
};

template<typename RealType>
AbstractRemapper<RealType>::
AbstractRemapper (const grid_ptr_type& src_grid,
                  const grid_ptr_type& tgt_grid)
 : m_state                 (RepoState::Clean)
 , m_src_grid              (src_grid)
 , m_tgt_grid              (tgt_grid)
 , m_num_fields            (0)
 , m_num_registered_fields (0)
 , m_fields_are_bound      (0)
 , m_num_bound_fields      (0)
{
  EKAT_REQUIRE_MSG(static_cast<bool>(src_grid), "Error! Invalid source grid pointer.\n");
  EKAT_REQUIRE_MSG(static_cast<bool>(tgt_grid), "Error! Invalid target grid pointer.\n");
}

template<typename RealType>
void AbstractRemapper<RealType>::
registration_begins () {
  EKAT_REQUIRE_MSG(m_state==RepoState::Clean,
                       "Error! Cannot start registration on a non-clean repo.\n"
                       "       Did you call 'registration_begins' already?\n");

  do_registration_begins();

  m_state = RepoState::Open;
}

template<typename RealType>
void AbstractRemapper<RealType>::
register_field (const identifier_type& src, const identifier_type& tgt) {
  EKAT_REQUIRE_MSG(m_state!=RepoState::Clean,
                       "Error! Cannot register fields in the remapper at this time.\n"
                       "       Did you forget to call 'registration_begins' ?");
  EKAT_REQUIRE_MSG(m_state!=RepoState::Closed,
                       "Error! Cannot register fields in the remapper at this time.\n"
                       "       Did you accidentally call 'registration_ends' already?");

  EKAT_REQUIRE_MSG(src.get_grid_name()==m_src_grid->name(),
                       "Error! Source field stores the wrong grid.\n");
  EKAT_REQUIRE_MSG(tgt.get_grid_name()==m_tgt_grid->name(),
                       "Error! Target field stores the wrong grid.\n");

  EKAT_REQUIRE_MSG(compatible_layouts(src.get_layout(),tgt.get_layout()),
                     "Error! Source and target layouts are not compatible.\n");

  do_register_field (src,tgt);

  m_fields_are_bound.push_back(false);
  ++m_num_registered_fields;
}

template<typename RealType>
void AbstractRemapper<RealType>::
register_field (const field_type& src, const field_type& tgt) {
  register_field(src.get_header().get_identifier(),
                 tgt.get_header().get_identifier());
  bind_field(src,tgt);
}

template<typename RealType>
void AbstractRemapper<RealType>::
unregister_field (const identifier_type& src, const identifier_type& tgt) {
  EKAT_REQUIRE_MSG(m_state==RepoState::Open,
                     "Error! You can only un-register fields during the registration phase.\n");

  const int ifield = find_field(src,tgt);
  EKAT_REQUIRE_MSG(ifield>=0,
                     "Error! The src/tgt pair of fields \n"
                     "         " + src.get_id_string() + "\n"
                     "         " + tgt.get_id_string() + "\n"
                     "       was not registered.\n");

  const bool was_bound = m_fields_are_bound[ifield];
  do_unregister_field(ifield);

  --m_num_registered_fields;
  if (was_bound) {
    --m_num_bound_fields;
  }
  m_fields_are_bound.erase(m_fields_are_bound.begin()+ifield);
}

template<typename RealType>
void AbstractRemapper<RealType>::
bind_field (const field_type& src, const field_type& tgt) {
  EKAT_REQUIRE_MSG(m_state!=RepoState::Clean,
                     "Error! Cannot bind fields in the remapper at this time.\n"
                     "       Did you forget to call 'registration_begins' ?");

  const auto& src_fid = src.get_header().get_identifier();
  const auto& tgt_fid = tgt.get_header().get_identifier();

  // Try to locate the pair of fields
  const int ifield = find_field(src_fid, tgt_fid);
  EKAT_REQUIRE_MSG(ifield>=0,
                     "Error! The src/tgt field pair\n"
                     "         " + src_fid.get_id_string() + "\n"
                     "         " + tgt_fid.get_id_string() + "\n"
                     "       was not registered. Please, register fields before binding them.\n");

  EKAT_REQUIRE_MSG(src.is_allocated(), "Error! Source field is not yet allocated.\n");
  EKAT_REQUIRE_MSG(tgt.is_allocated(), "Error! Source field is not yet allocated.\n");

  EKAT_REQUIRE_MSG(!m_fields_are_bound[ifield],
                     "Error! Field already bound.\n");

  do_bind_field(ifield,src,tgt);

  m_fields_are_bound[ifield] = true;
  ++m_num_bound_fields;
}

template<typename RealType>
void AbstractRemapper<RealType>::
registration_ends () {
  EKAT_REQUIRE_MSG(m_state!=RepoState::Closed,
                       "Error! Cannot call registration_ends at this time.\n"
                       "       Did you accidentally call 'registration_ends' already?");

  m_num_fields = m_num_registered_fields;

  do_registration_ends();

  m_state = RepoState::Closed;
}

// A short name for an AbstractRemapper factory
template<typename RealType>
using RemapperFactory =
    ekat::Factory<AbstractRemapper<RealType>,
                  ekat::CaseInsensitiveString,
                  std::shared_ptr<AbstractRemapper<RealType> >,
                  const ekat::ParameterList&>;

} // namespace scream

#endif // SCREAM_ABSTRACT_REMAPPER_HPP
