#include "cosp_functions.hpp"
namespace CospFunc {
    void initialize(int ncol, int nsubcol, int nlay) {
        cosp_c2f_init(ncol, nsubcol, nlay);
    };
    void main(int ncol, int nsubcol, int nlay) {
        cosp_c2f_main();
    };
    void finalize() {
        cosp_c2f_final();
    };
}
