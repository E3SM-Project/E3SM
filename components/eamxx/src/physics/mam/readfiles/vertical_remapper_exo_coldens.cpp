#include "vertical_remapper_exo_coldens.hpp"
#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"
#include <mam4xx/mam4.hpp>

namespace scream
{
VerticalRemapperExoColdensMAM4::
VerticalRemapperExoColdensMAM4 (const grid_ptr_type& src_grid,
                  const grid_ptr_type& tgt_grid)
:VerticalRemapper(src_grid,
                  tgt_grid,
                  true,
                  false) {
}

void VerticalRemapperExoColdensMAM4::remap_fwd_impl ()
{
  for (int i=0; i<m_num_fields; ++i) {
    const auto& f_src    = m_src_fields[i];
          auto& f_tgt    = m_tgt_fields[i];
    apply_vertical_interpolation(f_src, f_tgt);

  }
}

void VerticalRemapperExoColdensMAM4::set_delta_pressure(const std::string& file_name, const Field& pint )
{
    using view_1d_host    = typename KT::view_1d<Real>::HostMirror;

    const auto pint_host = pint.get_view<const Real **, Host>();
    const Real ptop_ref=pint_host(0,0);
    constexpr Real hPa2Pa = 100;
    scorpio::register_file(file_name,scorpio::FileMode::Read);
    const int n_exo_levs = scorpio::get_dimlen(file_name,"lev");
    view_1d_host lev_src("lev", n_exo_levs);
    scorpio::read_var(file_name,"lev",lev_src.data());
    scorpio::release_file(file_name);

    m_exo_ki=-1;
    m_exo_delp=0.0;
    if (ptop_ref <= lev_src(0)*hPa2Pa) { // levs(1) in Fortran is levs[0] in C++
        m_exo_ki = 0;
        m_exo_delp = 0.0;
    } else {
    for (int ki = 1; ki < n_exo_levs; ++ki) { // ki starts from 2 in Fortran, so we start from 1 in C++
      if (ptop_ref <= lev_src(ki) * hPa2Pa) {
        m_exo_delp = std::log(ptop_ref / lev_src(ki - 1)/ hPa2Pa) / std::log(lev_src(ki) / lev_src(ki - 1));
        m_exo_ki=ki;
        break; // exit the loop
        }
      }
    }

    EKAT_REQUIRE_MSG (m_exo_ki >= int(0),
      "[VerticalRemapperExoColdensMAM4] Error! m_exo_ki is negative.\n"
      "  m_exo_ki : " + std::to_string(m_exo_ki) + "\n"
    );
  }
/* Invokes MAM4XX routines for vertical interpolation.*/
void VerticalRemapperExoColdensMAM4::
apply_vertical_interpolation(const Field& f_src, const Field& f_tgt) const
{
  const auto datain = f_src.get_view<Real **>();
  const auto dataout =  f_tgt.get_view<Real **>();

  const int ncols = m_src_grid->get_num_local_dofs();

  using TPF = ekat::TeamPolicyFactory<KT::ExeSpace>;

  const auto policy =
      TPF::get_default_team_policy(ncols, 1);
  using Team = Kokkos::TeamPolicy<KT::ExeSpace>::member_type;

  const Real delp=m_exo_delp;
  const int ki = m_exo_ki;
  const int kl = ki == 0 ? 0: m_exo_ki-1;
  Kokkos::parallel_for(
      "vert_interpolation_exo_coldens", policy,
      KOKKOS_LAMBDA(const Team &team) {
        const int icol = team.league_rank();
        dataout(icol, 0) = datain(icol, kl)
                 + delp * (datain(icol,ki)
                 - datain(icol,kl));
  });
}

} // namespace scream
