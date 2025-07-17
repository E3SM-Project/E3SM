#include "vertical_remapper_mam4.hpp"
#include "share/io/eamxx_scorpio_interface.hpp"
#include <mam4xx/mam4.hpp>

namespace scream
{
VerticalRemapperMAM4::
VerticalRemapperMAM4 (const grid_ptr_type& src_grid,
                  const grid_ptr_type& tgt_grid,
                  const VertRemapType& vremp_type)
:VerticalRemapper(src_grid,
                  tgt_grid,
                  true,
                  false) {
  m_vremap_type=vremp_type;
}

void VerticalRemapperMAM4::remap_fwd_impl ()
{
  for (int i=0; i<m_num_fields; ++i) {
    const auto& f_src    = m_src_fields[i];
          auto& f_tgt    = m_tgt_fields[i];
    apply_vertical_interpolation(f_src, f_tgt, m_src_pmid, m_tgt_pmid);

  }
}
/* The DataInterpolation class only uses this method if the Custom type is employed.
   We only use this function for elevated emissions. */
void VerticalRemapperMAM4::
set_target_pressure (const Field& p)
{
  m_tgt_pmid=p;
}
/* It reads altitude from an NC file and sets m_src_pmid as altitude_int_src.
* DataInterpolation assumes that pressure is the variable for interpolation.
* Here, we use m_src_pmid to pass altitude, but note that we are using
* the MAM4XX interpolation for elevated emissions.*/
void VerticalRemapperMAM4::
set_source_pressure (const std::string& file_name )
{
  if (m_vremap_type == MAM4_ELEVATED_EMISSIONS) {
    auto layout = m_src_grid->get_vertical_layout(false);
    auto mbar = ekat::units::Units(100*ekat::units::Pa,"mbar");
    Field altitude_int_src(FieldIdentifier("altitude_int_field",layout,mbar,m_src_grid->name()));
    altitude_int_src.allocate_view();
    scorpio::register_file(file_name,scorpio::FileMode::Read);
    scorpio::read_var(file_name,"altitude_int",altitude_int_src.get_view<Real*,Host>().data());
    altitude_int_src.sync_to_dev();
    scorpio::release_file(file_name);
    m_src_pmid=altitude_int_src;
  }
}
/* Invokes MAM4XX routines for vertical interpolation.*/
void VerticalRemapperMAM4::
apply_vertical_interpolation(const Field& f_src, const Field& f_tgt,
                             const Field& p_src, const Field& p_tgt) const
{
  const auto p_tgt_c = p_tgt.get_view<const Real **>();
  const auto datain = f_src.get_view<Real **>();
  const auto dataout =  f_tgt.get_view<Real **>();

  const auto& f_tgt_l = f_tgt.get_header().get_identifier().get_layout();
  const int ncols = m_src_grid->get_num_local_dofs();
  const int nlevs_tgt = f_tgt_l.dims().back();

  const auto& f_src_l = f_src.get_header().get_identifier().get_layout();
  const int levsiz = f_src_l.dims().back();

  const auto policy =
      ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(ncols, nlevs_tgt);
  using Team = Kokkos::TeamPolicy<KT::ExeSpace>::member_type;
  if (m_vremap_type== MAM4_PSRef) {

    const int unit_factor_pin=1;
    const auto p_src_c = p_src.get_view<const Real **>();

    Kokkos::parallel_for(
      "vert_interp_psref", policy,
      KOKKOS_LAMBDA(const Team &team) {
        const int icol = team.league_rank();
        const auto pin_at_icol     = ekat::subview(p_src_c, icol);
        const auto pmid_at_icol    = ekat::subview(p_tgt_c, icol);
        const auto datain_at_icol  = ekat::subview(datain, icol);
        const auto dataout_at_icol = ekat::subview(dataout, icol);
        mam4::vertical_interpolation::vert_interp(
            team, levsiz, nlevs_tgt, pin_at_icol, pmid_at_icol, datain_at_icol,
            dataout_at_icol, unit_factor_pin);
          });
  } else if (m_vremap_type== MAM4_ZONAL)
  {
    // unit conversion from mbar->pascals
    const int unit_factor_pin=100;
    const auto p_src_c = p_src.get_view<const Real *>();

    Kokkos::parallel_for(
      "vert_interp_mam4_zonal", policy,
      KOKKOS_LAMBDA(const Team &team) {
        const int icol = team.league_rank();
        const auto pmid_at_icol    = ekat::subview(p_tgt_c, icol);
        const auto datain_at_icol  = ekat::subview(datain, icol);
        const auto dataout_at_icol = ekat::subview(dataout, icol);
        mam4::vertical_interpolation::vert_interp(
            team, levsiz, nlevs_tgt, p_src_c, pmid_at_icol, datain_at_icol,
            dataout_at_icol, unit_factor_pin);
          });
  } else if (m_vremap_type == MAM4_ELEVATED_EMISSIONS)
  {
      const auto src_x = p_src.get_view<const Real *>();
      const int pverp=nlevs_tgt+1;
      //FIXME: get this values from grid.
      constexpr int nlev = mam4::nlev;
      constexpr Real m2km    = 1e-3;
      Kokkos::parallel_for(
      "vert_interpolation_elevated_emissions", policy,
      KOKKOS_LAMBDA(const Team &team) {
        const int icol = team.league_rank();
        const auto datain_at_icol  = ekat::subview(datain, icol);
        const auto dataout_at_icol = ekat::subview(dataout, icol);

        // FIXME: Try to avoid copy of trg_x by modifying rebin
        // trg_x
        Real trg_x[nlev + 1];
        // I am trying to do this:
        // model_z(1:pverp) = m2km * state(c)%zi(i,pverp:1:-1)
        for(int i = 0; i < pverp; ++i) {
          trg_x[pverp - i - 1] = m2km * p_tgt_c(icol, i);
        }
        team.team_barrier();
        mam4::vertical_interpolation::rebin(team, levsiz, nlevs_tgt,
                                            src_x, trg_x, datain_at_icol, dataout_at_icol);
      });
  }

}

} // namespace scream
