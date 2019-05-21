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

// A remapper is basically a functor, that, given two fieldsm
// copies the first into the second, or viceversa. The copy must
// account for different layouts and/or different mpi distributions.

template<typename ScalarType, typename DeviceType>
class AbstractRemapper
{
public:
  using scalar_type       = ScalarType;
  using device_type       = DeviceType;
  using field_type        = Field<scalar_type,device_type>;
  using identifier_type   = typename field_type::identifier_type;
  using layout_type       = typename identifier_type::layout_type;

  AbstractRemapper (const std::shared_ptr<AbstractGrid> src_grid,
                    const std::shared_ptr<AbstractGrid> tgt_grid);

  virtual ~AbstractRemapper () = default;

  void set_num_fields (const int num_fields);
  void register_field (const field_type& src, const field_type& tgt);
  void registration_complete ();

  // The actual remap routine.
  void remap (const bool forward) const {
    error::runtime_check(m_state!=RepoState::Open, "Error! Cannot perform remapping at this time.\n"
                                                   "       Did you forget to call 'registration_complete'?\n");
    if (m_state!=RepoState::Clean) {
      if (forward) {
        do_remap_fwd ();
      } else {
        do_remap_bwd ();
      }
    }
  }

  const layout_type& get_src_layout (const int ifield) const {
    error::runtime_check(ifield>=0 && ifield<get_num_fields(), "Error! Field index out of bounds.\n");
    error::runtime_check(m_state==RepoState::Closed, "Error! Cannot query fields layout at this time.\n"
                                                   "         Did you forget to call 'registration_complete'?\n");
    return do_get_src_layout(ifield);
  }
  const layout_type& get_tgt_layout (const int ifield) const {
    error::runtime_check(ifield>=0 && ifield<get_num_fields(), "Error! Field index out of bounds.\n");
    error::runtime_check(m_state==RepoState::Closed, "Error! Cannot query fields layout at this time.\n"
                                                   "         Did you forget to call 'registration_complete'?\n");
    return do_get_tgt_layout(ifield);
  }

  int get_num_fields () const { return m_num_fields; }

protected:
  virtual const layout_type& do_get_src_layout (const int ifield) const = 0;
  virtual const layout_type& do_get_tgt_layout (const int ifield) const = 0;

  virtual void do_registration_start () = 0;
  virtual void do_register_field (const field_type& src, const field_type& tgt) = 0;
  virtual void do_registration_complete () = 0;

  virtual void do_remap_fwd () const = 0;
  virtual void do_remap_bwd () const = 0;

  // The state of the remapper
  RepoState     m_state;

  // The number of fields to remap
  int           m_num_fields;

  // The grids associated with the src and tgt fields
  std::shared_ptr<AbstractGrid> m_src_grid;
  std::shared_ptr<AbstractGrid> m_tgt_grid;
};

template<typename ScalarType, typename DeviceType>
AbstractRemapper<ScalarType,DeviceType>::
AbstractRemapper (const std::shared_ptr<AbstractGrid> src_grid,
                  const std::shared_ptr<AbstractGrid> tgt_grid)
 : m_state      (RepoState::Clean)
 , m_num_fields (0)
 , m_src_grid   (src_grid)
 , m_tgt_grid   (tgt_grid)
{
  // Nothing to do here
}

template<typename ScalarType, typename DeviceType>
void AbstractRemapper<ScalarType,DeviceType>::
set_num_fields (const int num_fields) {
  error::runtime_check(m_state==RepoState::Clean,
                       "Error! Cannot re-set the number of fields after it has been set.\n");

  m_num_fields = num_fields;

  do_registration_start();

  m_state = RepoState::Open;
}

template<typename ScalarType, typename DeviceType>
void AbstractRemapper<ScalarType,DeviceType>::
register_field (const field_type& src, const field_type& tgt) {
  error::runtime_check(m_state!=RepoState::Clean,
                       "Error! Cannot register fields in the remapper at this time.\n"
                       "       Did you forget to call 'set_num_fields' ?");
  error::runtime_check(m_state!=RepoState::Closed,
                       "Error! Cannot register fields in the remapper at this time.\n"
                       "       Did you accidentally call 'registration_complete' already?");

  error::runtime_check(src.get_header().get_identifier().get_grid_name()==m_src_grid->name(),
                       "Error! Source field stores the wrong grid.\n");
  error::runtime_check(src.get_header().get_identifier().get_grid_name()==m_src_grid->name(),
                       "Error! Target field stores the wrong grid.\n");

  do_register_field(src,tgt);
}

template<typename ScalarType, typename DeviceType>
void AbstractRemapper<ScalarType,DeviceType>::
registration_complete () {
  error::runtime_check(m_state!=RepoState::Closed,
                       "Error! Cannot call registration_complete at this time.\n"
                       "       Did you accidentally call 'registration_complete' already?");

  do_registration_complete();

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
