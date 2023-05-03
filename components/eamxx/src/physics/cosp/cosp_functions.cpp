#include "cosp_functions.hpp"
namespace scream {
    namespace CospFunc {
        void initialize(int ncol, int nsubcol, int nlay) {
            std::cout << "brhdebug: call cosp_c2f_init()" << std::endl;
            cosp_c2f_init(ncol, nsubcol, nlay);
        };
        void finalize() {
            cosp_c2f_final();
        };
    }
}
