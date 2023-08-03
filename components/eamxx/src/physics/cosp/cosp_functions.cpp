#include "cosp_functions.hpp"
namespace scream {
    namespace CospFunc {
        void initialize(int ncol, int nsubcol, int nlay) {
            cosp_c2f_init(ncol, nsubcol, nlay);
        };
        void finalize() {
            cosp_c2f_final();
        };
    }
}
