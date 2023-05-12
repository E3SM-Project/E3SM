#ifndef SCREAM_COSP_FUNCTIONS_HPP
#define SCREAM_COSP_FUNCTIONS_HPP
#include "share/scream_types.hpp"
#include "ekat/ekat_pack_kokkos.hpp"
#include "share/util/scream_deep_copy.hpp"
using scream::Real;
extern "C" void cosp_c2f_init(int ncol, int nsubcol, int nlay);
extern "C" void cosp_c2f_final();
extern "C" void cosp_c2f_run(int ncol, int nsubcol, int nlay, Real emsfc_lw, Real* sunlit, Real* skt,
    Real* T_mid, Real* p_mid, Real* p_int, Real* qv,
    Real* cldfrac, Real* reff_qc, Real* reff_qi, Real* dtau067, Real* dtau105,
    Real* isccp_cldtot);

namespace scream {

    namespace CospFunc {
        // views for single- and multi-column data
        //using view_1d_int   = typename KT::template view_1d<int>;
        //using view_1d       = typename KT::template view_1d<Real>;
        //using view_1d_const = typename KT::template view_1d<const Real>;
        //using view_2d       = typename KT::template view_2d<Real>;
        //using view_2d_const = typename KT::template view_2d<const Real>;

        using KT  = ekat::KokkosTypes<DefaultDevice>;
        using lview_host_1d = typename ekat::KokkosTypes<HostDevice>::template lview<Real* >;
        using lview_host_2d = typename ekat::KokkosTypes<HostDevice>::template lview<Real**>;
        template <typename S>
        using view_1d = typename KT::template view_1d<S>;
        template <typename S>
        using view_2d = typename KT::template view_2d<S>;

        template <typename S>
        using SmallPack = ekat::Pack<S,SCREAM_SMALL_PACK_SIZE>;
        using Spack = SmallPack<Real>;
        using Pack  = ekat::Pack<Real,Spack::n>;

        void initialize(int ncol, int nsubcol, int nlay) {
            cosp_c2f_init(ncol, nsubcol, nlay);
        };
        void finalize() {
            cosp_c2f_final();
        };
        void main(
                Int ncol, Int nsubcol, Int nlay, Real emsfc_lw,
                view_1d<const Real>& sunlit , view_1d<const Real>& skt,
                view_2d<const Pack>& T_mid  , view_2d<const Pack>& p_mid  , view_2d<const Pack>& p_int,
                view_2d<const Pack>& qv     , view_2d<const Pack>& cldfrac,
                view_2d<const Pack>& reff_qc, view_2d<const Pack>& reff_qi,
                view_2d<const Pack>& dtau067, view_2d<const Pack>& dtau105,
                view_1d<Real>& isccp_cldtot) {

            // Make host copies and permute data as needed
            lview_host_2d
                  T_mid_h("T_mid_h", ncol, nlay), p_mid_h("p_mid_h", ncol, nlay), p_int_h("p_int_h", ncol, nlay+1),
                  qv_h("qv_h", ncol, nlay), cldfrac_h("cldfrac_h", ncol, nlay),
                  reff_qc_h("reff_qc_h", ncol, nlay), reff_qi_h("reff_qi_h", ncol, nlay),
                  dtau067_h("dtau_067_h", ncol, nlay), dtau105_h("dtau105_h", ncol, nlay);
            auto sunlit_h = Kokkos::create_mirror_view(sunlit);
            auto skt_h = Kokkos::create_mirror_view(skt);
            auto isccp_cldtot_h = create_mirror_view(isccp_cldtot);

            {
                Kokkos::deep_copy(sunlit_h, sunlit);
                Kokkos::deep_copy(skt_h, skt);
            }

            {
                std::vector<view_2d<const Pack>> device_views = {T_mid, p_mid, qv, cldfrac, reff_qc, reff_qi, dtau067, dtau105};
                ekat::device_to_host({T_mid_h.data(), p_mid_h.data(), qv_h.data(), cldfrac_h.data(), reff_qc_h.data(), reff_qi_h.data(), dtau067_h.data(), dtau105_h.data()}, ncol, nlay, device_views, true);
            }

            {
                std::vector<view_2d<const Pack>> device_views = {p_int};
                ekat::device_to_host({p_int_h.data()}, Int(ncol), Int(nlay+1), device_views, true);
            }

            // Subsample here?

            // Call COSP wrapper
            cosp_c2f_run(ncol, nsubcol, nlay, emsfc_lw,
                    sunlit_h.data(), skt_h.data(), T_mid_h.data(), p_mid_h.data(), p_int_h.data(),
                    qv_h.data(),
                    cldfrac_h.data(), reff_qc_h.data(), reff_qi_h.data(), dtau067_h.data(), dtau105_h.data(),
                    isccp_cldtot_h.data());

            // Copy outputs back to device
            {
                Kokkos::deep_copy(isccp_cldtot, isccp_cldtot_h);
            }
        }
    }
}
#endif  /* SCREAM_COSP_FUNCTIONS_HPP */
