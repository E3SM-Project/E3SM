#ifndef PYUNITS_HPP
#define PYUNITS_HPP

#include <ekat/util/ekat_units.hpp>

#include <pybind11/pybind11.h>

namespace scream {

struct PyUnits {
  ekat::units::Units units;

  PyUnits ()
   : units(ekat::RationalConstant(0))
  {
    // Nothing to do here
  }
  PyUnits (const ekat::units::Units& src)
   : units(src)
  {
    // Nothing to do here
  }

  std::string str() {
    return units.to_string();
  }

  PyUnits (const PyUnits&) = default;
};

namespace short_units {

PyUnits one(ekat::units::Units::nondimensional());
PyUnits m(ekat::units::m);
PyUnits s(ekat::units::s);
PyUnits kg(ekat::units::kg);
PyUnits A(ekat::units::A);
PyUnits mol(ekat::units::mol);
PyUnits cd(ekat::units::cd);
PyUnits K(ekat::units::K);

PyUnits N (ekat::units::N);
PyUnits J (ekat::units::J);
PyUnits W (ekat::units::W);
PyUnits Pa (ekat::units::Pa);
}

} // namespace scream

#endif // PYUNITS_HPP
