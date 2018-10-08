#include "field_repository.hpp"
#include "field_repository_impl.hpp"

#include "scream_types.hpp"

namespace scream {

template class FieldRepository<ExecMemSpace>;

#ifdef CUDA_BUILD
// HostMemSpace differs from ExecMemSpace *only* for cuda builds.
template class FieldRepository<HostMemSpace>;
#endif

} // namespace scream
