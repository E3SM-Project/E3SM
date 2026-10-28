#ifndef SCREAM_PHYSICS_SHARE_HPP
#define SCREAM_PHYSICS_SHARE_HPP

#include "share/core/eamxx_types.hpp"
#include "share/util/eamxx_time_stamp.hpp"
#include <ekat_parameter_list.hpp>

namespace scream {

namespace physics {

// Host geometry for one evaluation timestamp and orbit context; recreate after either changes.
// Constructor angles are radians, calendar_day is one-based, and equation_of_time is seconds.
class SolarGeometry {
public:
	SolarGeometry(double calendar_day, double declination, double eccentricity_factor,
								double equation_of_time, bool use_equation_of_time)
		: m_calendar_day(calendar_day), m_declination(declination),
			m_eccentricity_factor(eccentricity_factor), m_equation_of_time(equation_of_time),
			m_use_equation_of_time(use_equation_of_time) {}

	double eccentricity_factor() const { return m_eccentricity_factor; }
  // Zero width is instantaneous; positive width averages forward over [t,t+dt].
  // Callers retain their own windows (radiation cadence, MAM timestep, etc.).
	double cosine_zenith(double latitude_radians, double longitude_radians,
											 double averaging_seconds) const;

private:
	double m_calendar_day, m_declination, m_eccentricity_factor, m_equation_of_time;
	bool m_use_equation_of_time;
};

// Installed driver geometry owns calendar/orbit/EOT policy, even with EOT off.
// Standalone fallback uses the timestamp calendar and the supplied orbital parameters.
SolarGeometry solar_geometry(const util::TimeStamp& timestamp, const ekat::ParameterList& params,
														 int orbital_year_setting, double eccentricity,
														 double obliquity, double perihelion);

}

extern "C" {

Real scream_pow(Real base, Real exp);
Real scream_sqrt(Real base);
Real scream_cbrt(Real base);
Real scream_gamma(Real input);
Real scream_log(Real input);
Real scream_log10(Real input);
Real scream_exp(Real input);
Real scream_expm1(Real input);
Real scream_tanh(Real input);
Real scream_cos(Real input);
Real scream_erf(Real input);

}

} // namespace scream

#endif // SCREAM_PHYSICS_SHARE_HPP
