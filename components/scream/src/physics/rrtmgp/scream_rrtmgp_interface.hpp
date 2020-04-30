#ifndef SCREAM_RRTMGP_INTERFACE_HPP
#define SCREAM_RRTMGP_INTERFACE_HPP

#include "share/scream_assert.hpp"
#include "share/util/scream_utils.hpp"

namespace scream {
    extern "C" {
        void rrtmgp_init_f90();
        void rrtmgp_main_f90();
        void rrtmgp_finalize_f90();
    }  // extern "C"
}  // namespace scream

#endif  // SCREAM_RRTMGP_INTERFACE_HPP
