#ifndef SCREAM_ABSTRACT_REMAPPER_HPP
#define SCREAM_ABSTRACT_REMAPPER_HPP

#include "share/field/field.hpp"
#include "share/field/field_utils.hpp"
#include "share/grid/abstract_grid.hpp"
#include "ekat/util/scream_factory.hpp"
#include "ekat/util/string_utils.hpp"
#include "ekat/util/scream_std_utils.hpp"
#include "ekat/scream_parameter_list.hpp"

namespace scream
{

// An abstract interface for a remapper

// A remapper is basically a functor, that, given two fields
// copies the first into the second, or viceversa. The copy must
// account for different layouts and/or different mpi distributions.
// This concept can be extended to remaps that involve interpolation,
// but as of now (07/2019) it is not the intent and projected use
// of this class in the scream framework

template<typename ScalarType, typename DeviceType>
class AbstractRemapper
{
public:
  using scalar_type     = ScalarType;
  using device_type     = DeviceType;
  using field_type      = Field<scalar_type,device_type>;
  using identifier_type = typename field_type::identifier_type;
  using layout_type     = typename identifier_type::layout_type;
  using grid_type       = AbstractGrid;
  using grid_ptr_type   = std::shared_ptr<const grid_type>;

  AbstractRemapper (const grid_ptr_type& src_grid,
                    const grid_ptr_type& tgt_grid);

  virtual ~AbstractRemapper () = default;

  void registration_begins ();
  void register_field (const field_type& src, const field_type& tgt);
  void register_field (const identifier_type& src, const identifier_type& tgt);
  void unregister_field (const identifier_type& src, const identifier_type& tgt);
  void registration_ends ();

  // The user is allowed to only provide identifiers in the registration phase.
  // In that case, fields have to be bound after registration is complete, and
  // before any call to remap.
  void bind_field (const field_type& src, const field_type& tgt);

  RepoState get_state () const { return m_state; }

  // The actual remap routine.
  void remap (const bool forward) const {
    scream_require_msg(m_state!=RepoState::Open,
                       "Error! Cannot perform remapping at this time.\n"
                       "       Did you forget to call 'registration_ends'?\n");

    scream_require_msg(m_num_bound_fields==m_num_fields,
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

  const identifier_type& get_src_field_id (const int ifield) const {
    scream_require_msg(ifield>=0 && ifield<m_num_registered_fields,
                       "Error! Field index out of bounds.\n");
    return do_get_src_field_id(ifield);
  }

  const identifier_type& get_tgt_field_id (const int ifield) const {
    scream_require_msg(ifield>=0 && ifield<m_num_registered_fields,
                       "Error! Field index out of bounds.\n");
    return do_get_tgt_field_id(ifield);
  }

  const field_type& get_src_field (const int ifield) const {
    scream_require_msg(m_state==RepoState::Closed,
                       "Error! Cannot call 'get_src_field' until registration has ended.\n");
    scream_require_msg(ifield>=0 && ifield<m_num_registered_fields,
                       "Error! Field index out of bounds.\n");
    return do_get_src_field(ifield);
  }

  const field_type& get_tgt_field (const int ifield) const {
    scream_require_msg(m_state==RepoState::Closed,
                       "Error! Cannot call 'get_tgt_field' until registration has ended.\n");
    scream_require_msg(ifield>=0 && ifield<m_num_registered_fields,
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
    scream_require_msg(m_state!=RepoState::Open,
                       "Error! Cannot call 'get_num_fields' durin the registration phase.\n");
    return m_num_fields;
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

  virtual void do_registration_begins () = 0;
  virtual void do_register_field (const identifier_type& src, const identifier_type& tgt) = 0;
  virtual void do_bind_field (const int ifield, const field_type& src, const field_type& tgt) = 0;
  virtual void do_unregister_field (const int ifield) = 0;
  virtual void do_registration_ends () = 0;

  virtual void do_remap_fwd () const = 0;
  virtual void do_remap_bwd () const = 0;

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

  // Whether fields have been provided. Recall that the user may register
  // fields passing only an identifier, and actually bind fields only
  // at a later moment.
  std::vector<bool>   m_fields_are_bound;
  int                 m_num_bound_fields;
};

template<typename ScalarType, typename DeviceType>
AbstractRemapper<ScalarType,DeviceType>::
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
  scream_require_msg(static_cast<bool>(src_grid), "Error! Invalid source grid pointer.\n");
  scream_require_msg(static_cast<bool>(tgt_grid), "Error! Invalid target grid pointer.\n");
}

template<typename ScalarType, typename DeviceType>
void AbstractRemapper<ScalarType,DeviceType>::
registration_begins () {
  scream_require_msg(m_state==RepoState::Clean,
                       "Error! Cannot start registration on a non-clean repo.\n"
                       "       Did you call 'registration_begins' already?\n");

  do_registration_begins();

  m_state = RepoState::Open;
}

template<typename ScalarType, typename DeviceType>
void AbstractRemapper<ScalarType,DeviceType>::
register_field (const identifier_type& src, const identifier_type& tgt) {
  scream_require_msg(m_state!=RepoState::Clean,
                       "Error! Cannot register fields in the remapper at this time.\n"
                       "       Did you forget to call 'registration_begins' ?");
  scream_require_msg(m_state!=RepoState::Closed,
                       "Error! Cannot register fields in the remapper at this time.\n"
                       "       Did you accidentally call 'registration_ends' already?");

  scream_require_msg(src.get_grid_name()==m_src_grid->name(),
                       "Error! Source field stores the wrong grid.\n");
  scream_require_msg(tgt.get_grid_name()==m_tgt_grid->name(),
                       "Error! Target field stores the wrong grid.\n");

  scream_require_msg(compatible_layouts(src.get_layout(),tgt.get_layout()),
                     "Error! Source and target layouts are not compatible.\n");

  do_register_field (src,tgt);

  m_fields_are_bound.push_back(false);
  ++m_num_registered_fields;
}

template<typename ScalarType, typename DeviceType>
void AbstractRemapper<ScalarType,DeviceType>::
register_field (const field_type& src, const field_type& tgt) {
  register_field(src.get_header().get_identifier(),
                 tgt.get_header().get_identifier());
  bind_field(src,tgt);
}

template<typename ScalarType, typename DeviceType>
void AbstractRemapper<ScalarType,DeviceType>::
unregister_field (const identifier_type& src, const identifier_type& tgt) {
  scream_require_msg(m_state==RepoState::Open,
                     "Error! You can only un-register fields during the registration phase.\n");

  const int ifield = find_field(src,tgt);
  scream_require_msg(ifield>=0,
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

template<typename ScalarType, typename DeviceType>
void AbstractRemapper<ScalarType,DeviceType>::
bind_field (const field_type& src, const field_type& tgt) {
  scream_require_msg(m_state!=RepoState::Clean,
                     "Error! Cannot bind fields in the remapper at this time.\n"
                     "       Did you forget to call 'registration_begins' ?");

  const auto& src_fid = src.get_header().get_identifier();
  const auto& tgt_fid = tgt.get_header().get_identifier();

  // Try to locate the pair of fields
  const int ifield = find_field(src_fid, tgt_fid);
  scream_require_msg(ifield>=0,
                     "Error! The src/tgt field pair\n"
                     "         " + src_fid.get_id_string() + "\n"
                     "         " + tgt_fid.get_id_string() + "\n"
                     "       was not registered. Please, register fields before binding them.\n");

  scream_require_msg(src.is_allocated(), "Error! Source field is not yet allocated.\n");
  scream_require_msg(tgt.is_allocated(), "Error! Source field is not yet allocated.\n");

  scream_require_msg(!m_fields_are_bound[ifield],
                     "Error! Field already bound.\n");

  do_bind_field(ifield,src,tgt);

  m_fields_are_bound[ifield] = true;
  ++m_num_bound_fields;
}

template<typename ScalarType, typename DeviceType>
void AbstractRemapper<ScalarType,DeviceType>::
registration_ends () {
  scream_require_msg(m_state!=RepoState::Closed,
                       "Error! Cannot call registration_ends at this time.\n"
                       "       Did you accidentally call 'registration_ends' already?");

  m_num_fields = m_num_registered_fields;

  do_registration_ends();

  m_state = RepoState::Closed;
}

// A short name for an AbstractRemapper factory
template<typename ScalarType, typename Device>
using RemapperFactory =
    util::Factory<AbstractRemapper<ScalarType,Device>,
                  util::CaseInsensitiveString,
                  std::shared_ptr<AbstractRemapper<ScalarType,Device>>,
                  const ParameterList&>;

} // namespace scream

#endif // SCREAM_ABSTRACT_REMAPPER_HPP
