#ifndef SCREAM_COSP_FUNCTIONS_HPP
#define SCREAM_COSP_FUNCTIONS_HPP
#include "share/core/eamxx_types.hpp"

#include <vector>

using scream::Real;

/*
 * C++ side of the COSP diagnostic's fortran bridge.
 *
 * This is a sibling of physics/cosp/cosp_functions.hpp, not a replacement: the
 * Cosp *process* keeps its own copy, and the two coexist. The fortran behind
 * this one lives in cosp_diag_c2f.F90, whose entry points are renamed
 * (cosp_c2f_*) precisely so that both can be linked into one executable.
 * Both drive the same external COSP fortran library.
 *
 * The one behavioral difference: COSP's fortran interface is unconditionally
 * double precision (it declares every argument as real(kind=c_double)), while
 * eamxx may be built single precision. Everything crossing the boundary is
 * therefore staged through double buffers and cast on the way in and out, so
 * that the diagnostic works at either precision.
 *
 * The casting is close to free. The 2d/3d fields already had to be copied into
 * LayoutLeft buffers to be readable from fortran, so the cast rides along with
 * a copy that was happening anyway; the 1d ones are a handful of columns. In a
 * double precision build every cast is an identity.
 */
using cosp_real = double;

extern "C" void cosp_c2f_init(int ncol, int nsubcol, int nlay);
extern "C" void cosp_c2f_final();
extern "C" void cosp_c2f_get_bins(cosp_real* tau_centers, cosp_real* tau_edges,
                                       cosp_real* prs_centers, cosp_real* prs_edges,
                                       cosp_real* cth_centers, cosp_real* cth_edges);
extern "C" void cosp_c2f_run(const int ncol, const int nsubcol, const int nlay, const int ntau, const int nctp, const int ncth,
    const cosp_real emsfc_lw, const cosp_real* sunlit, const cosp_real* skt,
    const cosp_real* T_mid, const cosp_real* p_mid, const cosp_real* p_int, const cosp_real* z_mid, const cosp_real* qv, const cosp_real* qc, const cosp_real* qi,
    const cosp_real* cldfrac, const cosp_real* reff_qc, const cosp_real* reff_qi, const cosp_real* dtau067, const cosp_real* dtau105,
    cosp_real* isccp_cldtot, cosp_real* isccp_ctptau, cosp_real* modis_ctptau, cosp_real* misr_cthtau);

namespace scream {

    namespace CospFunc {
        using lview_host_1d = typename ekat::KokkosTypes<HostDevice>::template lview<cosp_real*  >;
        using lview_host_2d = typename ekat::KokkosTypes<HostDevice>::template lview<cosp_real** >;
        using lview_host_3d = typename ekat::KokkosTypes<HostDevice>::template lview<cosp_real***>;
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

        // Bin centers and edges of the tau/pressure/height dimensions, straight
        // from mod_cosp_config. The fortran edges arrays are (2,nbins) in
        // column-major order, hence the 2*n sizes: the flat layout is
        // [lower_0, upper_0, lower_1, upper_1, ...].
        inline void get_bins(const int ntau, const int nprs, const int ncth,
                             Real* tau_centers, Real* tau_edges,
                             Real* prs_centers, Real* prs_edges,
                             Real* cth_centers, Real* cth_edges) {
            std::vector<cosp_real> tau_c(ntau), tau_e(2*ntau),
                                   prs_c(nprs), prs_e(2*nprs),
                                   cth_c(ncth), cth_e(2*ncth);
            cosp_c2f_get_bins(tau_c.data(), tau_e.data(),
                                   prs_c.data(), prs_e.data(),
                                   cth_c.data(), cth_e.data());

            auto to_real = [](const std::vector<cosp_real>& src, Real* dst) {
                for (size_t i = 0; i < src.size(); ++i) {
                    dst[i] = static_cast<Real>(src[i]);
                }
            };
            to_real(tau_c, tau_centers); to_real(tau_e, tau_edges);
            to_real(prs_c, prs_centers); to_real(prs_e, prs_edges);
            to_real(cth_c, cth_centers); to_real(cth_e, cth_edges);
        };

        inline void main(
                const Int ncol, const Int nsubcol, const Int nlay, const Int ntau, const Int nctp, const Int ncth, const Real emsfc_lw,
                const view_1d<const Real>& sunlit , const view_1d<const Real>& skt,
                const view_2d<const Real>& T_mid  , const view_2d<const Real>& p_mid  ,
                const view_2d<const Real>& p_int,  const view_2d<const Real>& z_mid,
                const view_2d<const Real>& qv     , const view_2d<const Real>& qc,
                const view_2d<const Real>& qi, const view_2d<const Real>& cldfrac,
                const view_2d<const Real>& reff_qc, const view_2d<const Real>& reff_qi,
                const view_2d<const Real>& dtau067, const view_2d<const Real>& dtau105,
                const view_1d<Real>& isccp_cldtot , const view_3d<Real>& isccp_ctptau,
                const view_3d<Real>& modis_ctptau, const view_3d<Real>& misr_cthtau) {

            // Make host copies, permuting the layout and widening to double as
            // needed. Note sunlit/skt/isccp_cldtot are staged too: they are
            // contiguous already, but they still have to change precision.
            lview_host_1d
                  sunlit_h("sunlit_h", ncol), skt_h("skt_h", ncol),
                  isccp_cldtot_h("isccp_cldtot_h", ncol);
            lview_host_2d
                  T_mid_h("T_mid_h", ncol, nlay), p_mid_h("p_mid_h", ncol, nlay), p_int_h("p_int_h", ncol, nlay+1),
                  z_mid_h("z_mid_h", ncol, nlay), qv_h("qv_h", ncol, nlay), qc_h("qc_h", ncol, nlay), qi_h("qi_h", ncol, nlay),
                  cldfrac_h("cldfrac_h", ncol, nlay),
                  reff_qc_h("reff_qc_h", ncol, nlay), reff_qi_h("reff_qi_h", ncol, nlay),
                  dtau067_h("dtau_067_h", ncol, nlay), dtau105_h("dtau105_h", ncol, nlay);
            lview_host_3d isccp_ctptau_h("isccp_ctptau_h", ncol, ntau, nctp);
            lview_host_3d modis_ctptau_h("modis_ctptau_h", ncol, ntau, nctp);
            lview_host_3d misr_cthtau_h("misr_cthtau_h", ncol, ntau, ncth);

            // Copy to layoutLeft host views
            for (int i = 0; i < ncol; i++) {
                sunlit_h(i) = sunlit(i);
                skt_h(i)    = skt(i);
                for (int j = 0; j < nlay; j++) {
                    T_mid_h(i,j) = T_mid(i,j);
                    p_mid_h(i,j) = p_mid(i,j);
                    z_mid_h(i,j) = z_mid(i,j);
                    qv_h(i,j) = qv(i,j);
                    qc_h(i,j) = qc(i,j);
                    qi_h(i,j) = qi(i,j);
                    cldfrac_h(i,j) = cldfrac(i,j);
                    reff_qc_h(i,j) = reff_qc(i,j);
                    reff_qi_h(i,j) = reff_qi(i,j);
                    dtau067_h(i,j) = dtau067(i,j);
                    dtau105_h(i,j) = dtau105(i,j);
                }
            }
            for (int i = 0; i < ncol; i++) {
                for (int j = 0; j < nlay+1; j++) {
                    p_int_h(i,j) = p_int(i,j);
                }
            }

            // Call COSP wrapper
            cosp_c2f_run(ncol, nsubcol, nlay, ntau, nctp, ncth,
                    emsfc_lw, sunlit_h.data(), skt_h.data(), T_mid_h.data(), p_mid_h.data(), p_int_h.data(),
                    z_mid_h.data(), qv_h.data(), qc_h.data(), qi_h.data(),
                    cldfrac_h.data(), reff_qc_h.data(), reff_qi_h.data(), dtau067_h.data(), dtau105_h.data(),
                    isccp_cldtot_h.data(), isccp_ctptau_h.data(), modis_ctptau_h.data(), misr_cthtau_h.data());

            // Copy outputs back to layoutRight views
            for (int i = 0; i < ncol; i++) {
                isccp_cldtot(i) = isccp_cldtot_h(i);
                for (int j = 0; j < ntau; j++) {
                    for (int k = 0; k < nctp; k++) {
                        isccp_ctptau(i,j,k) = isccp_ctptau_h(i,j,k);
                        modis_ctptau(i,j,k) = modis_ctptau_h(i,j,k);
                    }
                    for (int k = 0; k < ncth; k++) {
                        misr_cthtau(i,j,k) = misr_cthtau_h(i,j,k);
                    }
                }
            }
        }
    }
}
#endif  /* SCREAM_COSP_FUNCTIONS_HPP */
