#ifndef SCREAM_SHOC_INPUTS_INITIALIZER_HPP
#define SCREAM_SHOC_INPUTS_INITIALIZER_HPP

#include "share/field/field_initializer.hpp"
#include "physics/share/physics_functions.hpp"
#include "physics/shoc/shoc_constants.hpp"
#include "physics/shoc/shoc_functions.hpp"

namespace scream {

class SHOCInputsInitializer : public FieldInitializer
{
  using SHF            = shoc::Functions<Real, DefaultDevice>;
  using Spack          = typename SHF::Spack;
  using view_1d_scalar = typename SHF::view_1d<Real>;
  using view_2d        = typename SHF::view_2d<Spack>;

public:

  virtual ~SHOCInputsInitializer () = default;

  // The name of the initializer
  std::string name () const { return "SHOCInputsInitializer"; }

  // Initialize fields
  void initialize_fields ();

  const std::set<FieldIdentifier>& get_inited_fields () const {
    return m_fields_id;
  }

  // Linearly interpolates data at a specific elevation using elevation and
  // data arrays.
  KOKKOS_FUNCTION
  Real interpolate_data(const int N,
                        const view_1d_scalar& ref_elevations,
                        const view_1d_scalar& ref_data,
                        Real z)
  {
    // Find lower bound index s.t. ref_elevations(indx) >= z
    int index = 0;
    for (int i=0; i<N; ++i) {
      if (ref_elevations(i) >= z) break;
      ++index;
    }

    if (index == 0)
      return ref_data(0);
    else if (index < N) {
      const Real a = (z - ref_elevations(index-1)) /
                     (ref_elevations(index) - ref_elevations(index-1));
      return (1.0 - a) * ref_data(index-1) + a * ref_data(index);
    }
    else {
      // Don't extrapolate off the end of the table.
      return ref_data(N-1);
    }
  }

  // Calculates hydrostatic pressure for a specific column given elevation
  // data.
  KOKKOS_FUNCTION
  void compute_column_pressure(Int i, Int nlevs,
                               const Real rair, const Real cpair, const Real ggr,
                               const Real p0,
                               const view_1d_scalar z_ref, const view_1d_scalar theta_ref,
                               const view_2d& z, const view_2d& pres) {
    const Real surface_pressure = 1015e2;

    const Real c0 = rair / cpair;
    const Real c1 = -ggr * pow(p0, c0) / rair;
    const Real p_s = surface_pressure;

    // Move up the column, computing the pressures at each elevation.
    for (Int k=nlevs-1; k>=0; --k ) {
      const int view_k = k/Spack::n;
      const int view_kp1 = (k+1)/Spack::n;
      const int pack_k = k%Spack::n;
      const int pack_kp1 = (k+1)%Spack::n;

      Real z0 = (k  == nlevs-1) ? 0.0 : z(i, view_kp1)[pack_kp1];
      Real z1 = z(i, view_k)[pack_k];
      Real th0 = interpolate_data(5, z_ref, theta_ref, z0);
      Real th1 = interpolate_data(5, z_ref, theta_ref, z1);
      Real p0 = (k  == nlevs-1) ? p_s : pres(i, view_kp1)[pack_kp1];
      if (std::abs(th0 - th1) < 1e-14 * th0) {
        pres(i, view_k)[pack_k] = pow(pow(p0, c0) + c0*c1*(z1 - z0)/th0, 1.0/c0);
      }
      else {
        Real ra = (z1 - z0)/(th1 - th0);
        pres(i, view_k)[pack_k] = pow(pow(p0, c0) + c0*c1*ra*log(th1/th0), 1.0/c0);
      }
    }
  }

protected:

  void add_field (const field_type& f);
  void add_field (const field_type& f, const field_type& f_ref,
                  const remapper_ptr_type& remapper);

  std::map<std::string,const field_type>  m_fields;

  std::set<FieldIdentifier> m_fields_id;

  std::shared_ptr<AbstractRemapper<Real>> m_remapper;

};

} // namespace scream

#endif // SCREAM_SHOC_INPUTS_INITIALIZER_HPP
