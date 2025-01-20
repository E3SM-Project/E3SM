#ifndef EAMXX_IDENTITY_REMAPPER_HPP
#define EAMXX_IDENTITY_REMAPPER_HPP

#include "share/grid/remap/abstract_remapper.hpp"

namespace scream
{

/*
 *  A remapper representing an 'identity'.
 *
 *  This remapper ensures that after remap operations, src and tgt fields
 *  contain the same data. For memory/performance reasons, the user can
 *  request that src and tgt fields alias each other, in which case remap
 *  calls are in fact no-ops. But the fields must be registered via one of
 *  register_field_from_tgt/src methods.
 */

class IdentityRemapper : public AbstractRemapper
{
public:

  enum Aliasing {
    SrcAliasTgt,
    TgtAliasSrc,
    NoAliasing
  };

  IdentityRemapper (const grid_ptr_type grid,
                    const Aliasing aliasing = NoAliasing);

  ~IdentityRemapper () = default;

  void set_aliasing (const Aliasing aliasing);

  void register_field_from_src (const Field& src) override;
  void register_field_from_tgt (const Field& tgt) override;

protected:

  void registration_ends_impl () override;
  void remap_fwd_impl () override;
  void remap_bwd_impl () override;

  Aliasing m_aliasing;
};

} // namespace scream

#endif // EAMXX_IDENTITY_REMAPPER_HPP
