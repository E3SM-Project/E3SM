#include "vertical_remapper_elevated_emissions_mam4.hpp"
#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"
#include <ekat_team_policy_utils.hpp>
#include <mam4xx/mam4.hpp>

namespace scream
{

VerticalRemapperElevatedEmissionsMAM4::VerticalRemapperElevatedEmissionsMAM4(
  const grid_ptr_type& src_grid,
  const grid_ptr_type& tgt_grid)
{
  set_name("Vertical " + tgt_grid->name());

  m_bwd_allowed = false;

  EKAT_REQUIRE_MSG(
    src_grid->get_2d_scalar_layout().congruent(tgt_grid->get_2d_scalar_layout()),
    "Error! Source and target grid can only differ for their number of levels.\n");

  set_grids(src_grid, tgt_grid);
}

void VerticalRemapperElevatedEmissionsMAM4::remap_fwd_impl ()
{
  for (int i=0; i<m_num_fields; ++i) {
    const auto& f_src = m_src_fields[i];
          auto& f_tgt = m_tgt_fields[i];
    apply_vertical_interpolation(f_src, f_tgt);
  }
}

/* Reads altitude_int from an NC file and stores it as m_alt_int_src.
 * altitude_int is a 1D field (same for all columns) in km. */
void VerticalRemapperElevatedEmissionsMAM4::
set_source_interface_height (const std::string& file_name)
{
  using namespace ShortFieldTagsNames;
  auto layout = m_src_grid->get_vertical_layout(ILEV);
  auto mbar = ekat::units::Units(100*ekat::units::Pa,"mbar");
  Field altitude_int_src(FieldIdentifier("altitude_int_field",layout,mbar,m_src_grid->name()));
  altitude_int_src.allocate_view();
  scorpio::register_file(file_name,scorpio::FileMode::Read);
  scorpio::read_var(file_name,"altitude_int",altitude_int_src.get_view<Real*,Host>().data());
  altitude_int_src.sync_to_dev();
  scorpio::release_file(file_name);
  m_alt_int_src = altitude_int_src;
}

/* Stores the model interface height field (z_mam4_int) as m_z_iface_tgt. */
void VerticalRemapperElevatedEmissionsMAM4::
set_target_interface_height (const Field& z_iface)
{
  m_z_iface_tgt = z_iface;
}

/* Invokes the MAM4xx rebin routine for vertical interpolation from altitude
 * levels (km) to model interface heights. */
void VerticalRemapperElevatedEmissionsMAM4::
apply_vertical_interpolation(const Field& f_src, const Field& f_tgt) const
{
  const auto p_tgt_c = m_z_iface_tgt.get_view<const Real **>();
  const auto src_x   = m_alt_int_src.get_view<const Real *>();
  const auto datain  = f_src.get_view<Real **>();
  const auto dataout = f_tgt.get_view<Real **>();

  const auto& f_tgt_l = f_tgt.get_header().get_identifier().get_layout();
  const int ncols     = m_src_grid->get_num_local_dofs();
  const int nlevs_tgt = f_tgt_l.dims().back();

  const auto& f_src_l = f_src.get_header().get_identifier().get_layout();
  const int levsiz    = f_src_l.dims().back();

  const int pverp = nlevs_tgt + 1;
  //FIXME: get this value from grid.
  constexpr int nlev  = mam4::nlev;
  constexpr Real m2km = 1e-3;

  using TPF  = ekat::TeamPolicyFactory<KT::ExeSpace>;
  using Team = Kokkos::TeamPolicy<KT::ExeSpace>::member_type;
  const auto policy = TPF::get_default_team_policy(ncols, nlevs_tgt);

  Kokkos::parallel_for(
    "vert_interpolation_elevated_emissions", policy,
    KOKKOS_LAMBDA(const Team &team) {
      const int icol = team.league_rank();
      const auto datain_at_icol  = ekat::subview(datain, icol);
      const auto dataout_at_icol = ekat::subview(dataout, icol);

      // FIXME: Try to avoid copy of trg_x by modifying rebin
      // Reverse z_mam4_int from top-to-bottom to bottom-to-top (as required
      // by rebin) and convert from m to km.
      // This mirrors: model_z(1:pverp) = m2km * state(c)%zi(i,pverp:1:-1)
      Real trg_x[nlev + 1];
      for(int i = 0; i < pverp; ++i) {
        trg_x[pverp - i - 1] = m2km * p_tgt_c(icol, i);
      }
      team.team_barrier();
      mam4::vertical_interpolation::rebin(team, levsiz, nlevs_tgt,
                                          src_x, trg_x, datain_at_icol, dataout_at_icol);
    });
}

} // namespace scream
