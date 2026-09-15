#include "dynamics/homme/eamxx_homme_context.hpp"

namespace scream
{

EamxxHommeContext& EamxxHommeContext::singleton () {
  static EamxxHommeContext c;
  return c;
}

} // namespace scream
