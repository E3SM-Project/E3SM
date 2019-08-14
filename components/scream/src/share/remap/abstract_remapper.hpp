#ifndef SCREAM_ABSTRACT_REMAPPER_HPP
#define SCREAM_ABSTRACT_REMAPPER_HPP

#include "share/field/field.hpp"
#include "share/util/factory.hpp"
#include "share/util/string_utils.hpp"
#include "share/parameter_list.hpp"
#include "share/remap/remap_utils.hpp"
#include "share/grid/abstract_grid.hpp"

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

  AbstractRemapper (const grid_ptr_type src_grid,
                    const grid_ptr_type tgt_grid);

  virtual ~AbstractRemapper () = default;

  // If the number of fields is known, you can call this to pre-allocate
  // some memory. This can just be a 'good guess'. The actual number of
  // registered fields is computed during registration_complete.
  // TODO: should this be removed alltogether, and just go with dynamic
  //       sizing of arrays?
  void set_num_fields (const int num_fields);

  void registration_begins ();
  void register_field (const field_type& src, const field_type& tgt);
  void register_field (const identifier_type& src, const identifier_type& tgt);
  void registration_complete ();

  // The user is allowed to only provide identifiers in the registration phase.
  // In that case, fields have to be bound after registration is complete, and
  // before any call to remap.
  void bind_field (const int ifield, const field_type& src, const field_type& tgt);

  RepoState get_state () const { return m_state; }

  // The actual remap routine.
  void remap (const bool forward) const {
    scream_require_msg(m_state!=RepoState::Open,
                       "Error! Cannot perform remapping at this time.\n"
                       "       Did you forget to call 'registration_complete'?\n");

    scream_require_msg(m_all_fields_are_bound,
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

  const std::set<identifier_type>& get_src_fields_ids () const {
    return m_src_fields_ids;
  }
  const std::set<identifier_type>& get_tgt_fields_ids () const {
    return m_tgt_fields_ids;
  }

  bool field_is_bound (const int ifield) const { return m_fields_are_bound[ifield]; }

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

  int get_num_fields () const {
    scream_require_msg(m_state==RepoState::Closed,
                       "Error! Cannot call 'get_num_fields' until registration has ended.\n"
                       "       Before then, this number may be incorrect.\n");
    return m_num_fields;
  }

protected:
  virtual const identifier_type& do_get_src_field_id (const int ifield) const = 0;
  virtual const identifier_type& do_get_tgt_field_id (const int ifield) const = 0;
  virtual const field_type& do_get_src_field (const int ifield) const = 0;
  virtual const field_type& do_get_tgt_field (const int ifield) const = 0;

  virtual void do_registration_begins () = 0;
  virtual void do_register_field (const identifier_type& src, const identifier_type& tgt) = 0;
  virtual void do_bind_field (const int ifield, const field_type& src, const field_type& tgt) = 0;
  virtual void do_registration_complete () = 0;

  virtual void do_remap_fwd () const = 0;
  virtual void do_remap_bwd () const = 0;

  // The state of the remapper
  RepoState     m_state;

  // The grids associated with the src and tgt fields
  grid_ptr_type m_src_grid;
  grid_ptr_type m_tgt_grid;

  // The number of fields to remap
  int           m_num_fields;

  // The number of fields currently registered. This is guaranteed to be equal
  // to m_num_fields only when registration is not undergoing. During registration,
  // we have either 0<=m_num_registered_fields<=m_num_fields (if `set_num_fields`
  // was called), or m_num_fields=0<=m_num_registered_fields.
  int           m_num_registered_fields;

  // Whether fields have been provided. Recall that the user may register
  // fields passing only an identifier, and actually bind fields only
  // at a later moment.
  std::vector<bool>   m_fields_are_bound;
  bool                m_all_fields_are_bound;

  std::set<identifier_type>  m_src_fields_ids;
  std::set<identifier_type>  m_tgt_fields_ids;
};

template<typename ScalarType, typename DeviceType>
AbstractRemapper<ScalarType,DeviceType>::
AbstractRemapper (const grid_ptr_type src_grid,
                  const grid_ptr_type tgt_grid)
 : m_state                 (RepoState::Clean)
 , m_src_grid              (src_grid)
 , m_tgt_grid              (tgt_grid)
 , m_num_fields            (0)
 , m_num_registered_fields (0)
 , m_fields_are_bound      (0)
 , m_all_fields_are_bound  (true)
{
  // Nothing to do here
}

template<typename ScalarType, typename DeviceType>
void AbstractRemapper<ScalarType,DeviceType>::
set_num_fields (const int num_fields) {
  scream_require_msg(m_state==RepoState::Clean,
                     "Error! Cannot re-set the number of fields after registration has beginned.\n");
  scream_require_msg(num_fields>=0,
                     "Error! Invalid number of fields specified.\n");

  m_num_fields = num_fields;
  m_fields_are_bound.resize(m_num_fields,false);
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
                       "       Did you accidentally call 'registration_complete' already?");

  scream_require_msg(src.get_grid_name()==m_src_grid->name(),
                       "Error! Source field stores the wrong grid.\n");
  scream_require_msg(tgt.get_grid_name()==m_tgt_grid->name(),
                       "Error! Target field stores the wrong grid.\n");

  const auto src_lt = get_layout_type(src.get_layout());
  const auto tgt_lt = get_layout_type(tgt.get_layout());

  scream_require_msg(src_lt==tgt_lt, "Error! Source and target layouts do not match.\n");

  m_fields_are_bound.push_back(false);

  do_register_field (src,tgt);

  ++m_num_registered_fields;
  m_all_fields_are_bound = false;
}

template<typename ScalarType, typename DeviceType>
void AbstractRemapper<ScalarType,DeviceType>::
register_field (const field_type& src, const field_type& tgt) {
  register_field(src.get_header().get_identifier(),
                 tgt.get_header().get_identifier());
  bind_field(m_num_registered_fields-1,src,tgt);
}

template<typename ScalarType, typename DeviceType>
void AbstractRemapper<ScalarType,DeviceType>::
bind_field (const int ifield, const field_type& src, const field_type& tgt) {
  scream_require_msg(m_state!=RepoState::Clean,
                       "Error! Cannot bind fields in the remapper at this time.\n"
                       "       Did you forget to call 'registration_begins' ?");
  scream_require_msg(ifield>=0 && ifield<m_num_registered_fields,
                     "Error! Field index (" + std::to_string(ifield) + ") out of bounds.\n");
  scream_require_msg(!m_fields_are_bound[ifield],
                     "Error! Fields already bound for index " + std::to_string(ifield) + ".\n");

  scream_require_msg(src.get_header().get_identifier()==get_src_field_id(ifield),
                       "Error! Source field has the wrong identifier.\n");
  scream_require_msg(tgt.get_header().get_identifier()==get_tgt_field_id(ifield),
                       "Error! Target field has the wrong identifier.\n");

  scream_require_msg(src.is_allocated(), "Error! Source field is not yet allocated.\n");
  scream_require_msg(tgt.is_allocated(), "Error! Target field is not yet allocated.\n");

  do_bind_field(ifield,src,tgt);

  m_fields_are_bound[ifield] = true;

  // Assume all good
  m_all_fields_are_bound = true;
  for (auto bound : m_fields_are_bound) {
    m_all_fields_are_bound &= bound;
  }
}

template<typename ScalarType, typename DeviceType>
void AbstractRemapper<ScalarType,DeviceType>::
registration_complete () {
  scream_require_msg(m_state!=RepoState::Closed,
                       "Error! Cannot call registration_complete at this time.\n"
                       "       Did you accidentally call 'registration_complete' already?");

  m_num_fields = m_num_registered_fields;

  do_registration_complete();

  for (int ifield=0; ifield<m_num_fields; ++ifield) {
    auto it_bool = m_src_fields_ids.insert(get_src_field_id(ifield));
    scream_require_msg(it_bool.second, "Error! The remapper has the same field registered twice.\n");
    it_bool = m_tgt_fields_ids.insert(get_tgt_field_id(ifield));
    scream_require_msg(it_bool.second, "Error! The remapper has the same field registered twice.\n");
  }

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
