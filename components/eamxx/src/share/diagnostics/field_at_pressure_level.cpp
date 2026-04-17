#include "field_at_pressure_level.hpp"
#include "share/util/eamxx_universal_constants.hpp"

#include <ekat_std_utils.hpp>
#include <ekat_upper_bound.hpp>
#include <ekat_units.hpp>

namespace scream
{

// =========================================================================================
FieldAtPressureLevel::
FieldAtPressureLevel (const ekat::Comm& comm, const ekat::ParameterList& params)
 : AtmosphereDiagnostic(comm,params)
{
  m_field_name = m_params.get<std::string>("field_name");

  const auto units = m_params.get<std::string>("pressure_units");
  EKAT_REQUIRE_MSG (units=="mb" or units=="hPa" or units=="Pa",
      "Error! Invalid units for FieldAtPressureLevel.\n"
      " - input units: " + units + "\n"
      " - valid units: 'mb', 'hPa', 'Pa'\n");

  // Figure out the pressure value, and convert to Pa if needed
  auto p_value = m_params.get<std::string>("pressure_value");

  if (units=="mb" || units=="hPa") {
    m_pressure_level = std::stod(p_value)*100;
  } else {
    m_pressure_level = std::stod(p_value);
  }

  m_diag_name = m_field_name + "_at_" + p_value + units;

  // Determine interpolation type (default: log-linear)
  const auto vit = m_params.get<std::string>("vert_interpolation_type","log-linear");
  EKAT_REQUIRE_MSG (vit=="linear" or vit=="log-linear",
      "Error! Invalid value for 'vert_interpolation_type' in FieldAtPressureLevel.\n"
      " - input value: " + vit + "\n"
      " - valid values: linear, log-linear\n");
  m_log_linear = (vit=="log-linear");
  m_log_pressure_level = std::log(m_pressure_level);
}

void FieldAtPressureLevel::
create_requests ()
{
  const auto& gname = m_params.get<std::string>("grid_name");
  add_field<Required>(m_field_name,gname);

  // We don't know yet which one we need
  add_field<Required>("p_mid",gname);
  add_field<Required>("p_int",gname);
}

void FieldAtPressureLevel::
initialize_impl (const RunType /*run_type*/)
{
  const auto& f = get_field_in(m_field_name);
  const auto& fid = f.get_header().get_identifier();

  // Sanity checks
  using namespace ShortFieldTagsNames;
  const auto& layout = fid.get_layout();
  EKAT_REQUIRE_MSG (layout.rank()>=2 && layout.rank()<=3,
      "Error! Field rank not supported by FieldAtPressureLevel.\n"
      " - field name: " + fid.name() + "\n"
      " - field layout: " + layout.to_string() + "\n"
      "NOTE: if you requested something like 'field_horiz_avg_at_Y',\n"
      "      you can avoid this error by requesting 'fieldX_at_Y_horiz_avg' instead.\n");
  const auto tag = layout.tags().back();
  EKAT_REQUIRE_MSG (tag==LEV || tag==ILEV,
      "Error! FieldAtPressureLevel diagnostic expects a layout ending with 'LEV'/'ILEV' tag.\n"
      " - field name  : " + fid.name() + "\n"
      " - field layout: " + layout.to_string() + "\n");

  // All good, create the diag output
  auto d_fid = fid.clone(m_diag_name).reset_layout(layout.clone().strip_dim(tag));
  m_diagnostic_output = Field(d_fid);
  m_diagnostic_output.allocate_view();

  m_pressure_name = tag==LEV ? "p_mid" : "p_int";

  // Add a field representing the mask as extra data to the diagnostic field.
  m_diagnostic_output.create_valid_mask();
  m_diagnostic_output.get_header().set_may_be_filled(true);

  using stratts_t = std::map<std::string,std::string>;

  // Propagate any io string attribute from input field to diag field
  const auto& src = get_fields_in().front();
  const auto& src_atts = src.get_header().get_extra_data<stratts_t>("io: string attributes");
        auto& dst_atts = m_diagnostic_output.get_header().get_extra_data<stratts_t>("io: string attributes");
  for (const auto& [name, val] : src_atts) {
    dst_atts[name] = val;
  }
}

// =========================================================================================
void FieldAtPressureLevel::compute_diagnostic_impl()
{
  using KT = KokkosTypes<DefaultDevice>;
  using MemberType = typename KT::MemberType;

  //This is 2D source pressure
  const Field& p_src = get_field_in(m_pressure_name);
  const auto p_src_v = p_src.get_view<const Real**>();
  const Field& f = get_field_in(m_field_name);

  // The setup for interpolation varies depending on the rank of the input field:
  const int rank = f.rank();

  const auto& pl = p_src.get_header().get_identifier().get_layout();
  const int ncols = pl.dim(0);
  const int nlevs = pl.dim(1);

  auto p_tgt = m_pressure_level;
  auto log_p_tgt = m_log_pressure_level;
  auto log_linear = m_log_linear;
  constexpr auto fval = constants::fill_value<Real>;
  if (rank==2) {
    auto policy = KT::RangePolicy(0,ncols);
    auto diag = m_diagnostic_output.get_view<Real*>();
    auto mask = m_diagnostic_output.get_valid_mask().get_view<int*>();
    auto f_v  = f.get_view<const Real**>();
    Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const int icol) {
      auto x1 = ekat::subview(p_src_v,icol);
      auto y1 = ekat::subview(f_v,icol);
      auto beg = x1.data();
      auto end = beg + nlevs;
      auto last = beg + (nlevs-1);
      if (p_tgt<*beg or p_tgt>*last) {
        diag(icol) = fval;
        mask(icol) = 0;
      } else {
        auto ub = ekat::upper_bound(beg,end,p_tgt);     
        auto k1 = ub - beg;
        if (k1==0) {
          // Corner case: p_tgt==y1(0)
          diag(icol) = y1(0);
        } else if (k1==nlevs) {
          // Corner case: p_tgt==y1(nlevs-1)
          diag(icol) = y1(nlevs-1);
        } else {
          // General case: interpolate between k1-1 and k1
          if (log_linear) {
            auto lx0 = Kokkos::log(x1(k1-1));
            auto lx1 = Kokkos::log(x1(k1));
            diag(icol) = y1(k1-1) + (y1(k1)-y1(k1-1))/(lx1 - lx0) * (log_p_tgt - lx0);
          } else {
            diag(icol) = y1(k1-1) + (y1(k1)-y1(k1-1))/(x1(k1) - x1(k1-1)) * (p_tgt-x1(k1-1));
          }
        }
        mask(icol) = 1;
      }
    });
  } else if (rank==3) {
    const int ndims = f.get_header().get_identifier().get_layout().get_vector_dim();
    auto policy = KT::TeamPolicy(ncols,ndims);
    auto diag = m_diagnostic_output.get_view<Real**>();
    auto mask = m_diagnostic_output.get_valid_mask().get_view<int**>();
    auto f_v  = f.get_view<const Real***>();
    Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const MemberType& team) {
      int icol = team.league_rank();
      auto x1 = ekat::subview(p_src_v,icol);
      auto beg = x1.data();
      auto end = beg + nlevs;
      auto last = beg + (nlevs-1);
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team,ndims),[&](const int idim) {
        if (p_tgt<*beg or p_tgt>*last) {
          diag(icol,idim) = fval; // TODO: don't bother setting an arbitrary value
          mask(icol,idim) = 0;
        } else {
          auto y1 = ekat::subview(f_v,icol,idim);
          auto ub = ekat::upper_bound(beg,end,p_tgt);     
          auto k1 = ub - beg;
          if (k1==0) {
            // Corner case: p_tgt==y1(0)
            diag(icol,idim) = y1(0);
          } else if (k1==nlevs) {
            // Corner case: p_tgt==y1(nlevs-1)
            diag(icol,idim) = y1(nlevs-1);
          } else {
            // General case: interpolate between k1-1 and k1
            if (log_linear) {
              auto lx0 = Kokkos::log(x1(k1-1));
              auto lx1 = Kokkos::log(x1(k1));
              diag(icol,idim) = y1(k1-1) + (y1(k1)-y1(k1-1))/(lx1 - lx0) * (log_p_tgt - lx0);
            } else {
              diag(icol,idim) = y1(k1-1) + (y1(k1)-y1(k1-1))/(x1(k1) - x1(k1-1)) * (p_tgt-x1(k1-1));
            }
          }
          mask(icol,idim) = 1;
        }
      });
    });
  } else {
    EKAT_ERROR_MSG("Error! field at pressure level only supports fields ranks 2 and 3 \n");
  }

}

} //namespace scream
