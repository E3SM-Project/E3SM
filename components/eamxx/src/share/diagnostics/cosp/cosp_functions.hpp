#ifndef SCREAM_COSP_FUNCTIONS_HPP
#define SCREAM_COSP_FUNCTIONS_HPP
#include "share/core/eamxx_types.hpp"

// Fortran COSP interface always uses double precision (c_double)
extern "C" void cosp_c2f_init(int ncol, int nsubcol, int nlay);
extern "C" void cosp_c2f_final();
extern "C" void cosp_c2f_run(const int ncol, const int nsubcol, const int nlay, const int ntau, const int nctp, const int ncth,
    const double emsfc_lw, const double* sunlit, const double* skt,
    const double* T_mid, const double* p_mid, const double* p_int, const double* z_mid, const double* qv, const double* qc, const double* qi,
    const double* cldfrac, const double* reff_qc, const double* reff_qi, const double* dtau067, const double* dtau105,
    double* isccp_cldtot, double* isccp_ctptau, double* modis_ctptau, double* misr_cthtau);

namespace scream {

    namespace CospFunc {
        using lview_host_1d = typename ekat::KokkosTypes<HostDevice>::template lview<double*  >;
        using lview_host_2d = typename ekat::KokkosTypes<HostDevice>::template lview<double** >;
        using lview_host_3d = typename ekat::KokkosTypes<HostDevice>::template lview<double***>;
        template <typename S>
        using view_1d = typename ekat::KokkosTypes<HostDevice>::template view_1d<S>;
        template <typename S>
        using view_2d = typename ekat::KokkosTypes<HostDevice>::template view_2d<S>;
        template <typename S>
        using view_3d = typename ekat::KokkosTypes<HostDevice>::template view_3d<S>;

        inline void initialize(int ncol, int nsubcol, int nlay) {
            cosp_c2f_init(ncol, nsubcol, nlay);
        };
        inline void finalize() {
            cosp_c2f_final();
        };
        inline void main(
                const int ncol, const int nsubcol, const int nlay, const int ntau, const int nctp, const int ncth, const Real emsfc_lw,
                const view_1d<const Real>& sunlit , const view_1d<const Real>& skt,
                const view_2d<const Real>& T_mid  , const view_2d<const Real>& p_mid  ,
                const view_2d<const Real>& p_int,  const view_2d<const Real>& z_mid,
                const view_2d<const Real>& qv     , const view_2d<const Real>& qc,
                const view_2d<const Real>& qi, const view_2d<const Real>& cldfrac,
                const view_2d<const Real>& reff_qc, const view_2d<const Real>& reff_qi,
                const view_2d<const Real>& dtau067, const view_2d<const Real>& dtau105,
                const view_1d<Real>& isccp_cldtot , const view_3d<Real>& isccp_ctptau,
                const view_3d<Real>& modis_ctptau, const view_3d<Real>& misr_cthtau) {

            // Fortran COSP always uses double precision. Create double-precision
            // layoutLeft host temporaries for the Fortran interface.
            lview_host_1d sunlit_dbl("sunlit_dbl", ncol), skt_dbl("skt_dbl", ncol);
            lview_host_1d isccp_cldtot_dbl("isccp_cldtot_dbl", ncol);
            lview_host_2d
                  T_mid_dbl("T_mid_dbl", ncol, nlay), p_mid_dbl("p_mid_dbl", ncol, nlay), p_int_dbl("p_int_dbl", ncol, nlay+1),
                  z_mid_dbl("z_mid_dbl", ncol, nlay), qv_dbl("qv_dbl", ncol, nlay), qc_dbl("qc_dbl", ncol, nlay), qi_dbl("qi_dbl", ncol, nlay),
                  cldfrac_dbl("cldfrac_dbl", ncol, nlay),
                  reff_qc_dbl("reff_qc_dbl", ncol, nlay), reff_qi_dbl("reff_qi_dbl", ncol, nlay),
                  dtau067_dbl("dtau067_dbl", ncol, nlay), dtau105_dbl("dtau105_dbl", ncol, nlay);
            lview_host_3d isccp_ctptau_dbl("isccp_ctptau_dbl", ncol, ntau, nctp);
            lview_host_3d modis_ctptau_dbl("modis_ctptau_dbl", ncol, ntau, nctp);
            lview_host_3d misr_cthtau_dbl("misr_cthtau_dbl", ncol, ntau, ncth);

            // Copy inputs to double-precision layoutLeft host views
            for (int i = 0; i < ncol; i++) {
                sunlit_dbl(i) = static_cast<double>(sunlit(i));
                skt_dbl(i)    = static_cast<double>(skt(i));
                for (int j = 0; j < nlay; j++) {
                    T_mid_dbl(i,j)   = static_cast<double>(T_mid(i,j));
                    p_mid_dbl(i,j)   = static_cast<double>(p_mid(i,j));
                    z_mid_dbl(i,j)   = static_cast<double>(z_mid(i,j));
                    qv_dbl(i,j)      = static_cast<double>(qv(i,j));
                    qc_dbl(i,j)      = static_cast<double>(qc(i,j));
                    qi_dbl(i,j)      = static_cast<double>(qi(i,j));
                    cldfrac_dbl(i,j) = static_cast<double>(cldfrac(i,j));
                    reff_qc_dbl(i,j) = static_cast<double>(reff_qc(i,j));
                    reff_qi_dbl(i,j) = static_cast<double>(reff_qi(i,j));
                    dtau067_dbl(i,j) = static_cast<double>(dtau067(i,j));
                    dtau105_dbl(i,j) = static_cast<double>(dtau105(i,j));
                }
            }
            for (int i = 0; i < ncol; i++) {
                for (int j = 0; j < nlay+1; j++) {
                    p_int_dbl(i,j) = static_cast<double>(p_int(i,j));
                }
            }

            // Call COSP Fortran wrapper (always double precision)
            const double emsfc_lw_dbl = static_cast<double>(emsfc_lw);
            cosp_c2f_run(ncol, nsubcol, nlay, ntau, nctp, ncth,
                    emsfc_lw_dbl, sunlit_dbl.data(), skt_dbl.data(), T_mid_dbl.data(), p_mid_dbl.data(), p_int_dbl.data(),
                    z_mid_dbl.data(), qv_dbl.data(), qc_dbl.data(), qi_dbl.data(),
                    cldfrac_dbl.data(), reff_qc_dbl.data(), reff_qi_dbl.data(), dtau067_dbl.data(), dtau105_dbl.data(),
                    isccp_cldtot_dbl.data(), isccp_ctptau_dbl.data(), modis_ctptau_dbl.data(), misr_cthtau_dbl.data());

            // Copy double-precision outputs back to Real views
            for (int i = 0; i < ncol; i++) {
                isccp_cldtot(i) = static_cast<Real>(isccp_cldtot_dbl(i));
                for (int j = 0; j < ntau; j++) {
                    for (int k = 0; k < nctp; k++) {
                        isccp_ctptau(i,j,k) = static_cast<Real>(isccp_ctptau_dbl(i,j,k));
                        modis_ctptau(i,j,k) = static_cast<Real>(modis_ctptau_dbl(i,j,k));
                    }
                    for (int k = 0; k < ncth; k++) {
                        misr_cthtau(i,j,k) = static_cast<Real>(misr_cthtau_dbl(i,j,k));
                    }
                }
            }
        }
    }
}
#endif  /* SCREAM_COSP_FUNCTIONS_HPP */
