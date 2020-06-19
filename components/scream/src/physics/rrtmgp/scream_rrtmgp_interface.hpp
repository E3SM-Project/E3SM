#ifndef SCREAM_RRTMGP_INTERFACE_HPP
#define SCREAM_RRTMGP_INTERFACE_HPP

#include "ekat/scream_assert.hpp"
#include "ekat/util/scream_utils.hpp"

namespace scream {
    namespace rrtmgp {
        /* 
         * These may or may not be needed to interface with the fortran code.
         * If we jump directly to the C++ version of RRTMGP, or use the C++
         * API for RRTMGP, we should be able to do away with these.
         */
        extern "C" {
            void rrtmgp_init_f90();
            void rrtmgp_main_f90();
            void rrtmgp_finalize_f90();
        }  // extern "C"
        /*
         * Assuming we can jump directly to using a C++ API...
         */
        extern void rrtmgp_init();
        extern void rrtmgp_main();
        extern void rrtmgp_finalize();
    } // namespace rrtmgp
}  // namespace scream

#endif  // SCREAM_RRTMGP_INTERFACE_HPP
