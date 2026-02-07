/**
 * @file eatm_atm_emulator.hpp
 * @brief Concrete atmosphere emulator for MCT coupling.
 *
 * Inherits from emulator::Emulator.  Owns the grid decomposition
 * and synthetic lat-lon data that the Fortran MCT layer queries.
 */

#ifndef EATM_ATM_EMULATOR_HPP
#define EATM_ATM_EMULATOR_HPP

#include "emulator.hpp"
#include <vector>

namespace eatm {

class AtmEmulator : public emulator::Emulator {
public:
  /**
   * @brief Construct a new AtmEmulator.
   *
   * @param fcomm  Fortran MPI communicator handle (integer)
   * @param compid Component id assigned by the driver
   * @param nxg    Number of longitude points
   * @param nyg    Number of latitude points
   */
  AtmEmulator(int fcomm, int compid, int nxg, int nyg);

  /**
   * @brief Compute a 1-D contiguous decomposition and fill
   *        synthetic grid arrays (lon, lat, area, mask, frac).
   *
   * Must be called after construction and before any MCT
   * queries from the Fortran side.
   */
  void compute_decomposition();

  // ---- accessors used by the C binding layer ----
  int  lsize()  const { return m_lsize; }
  int  nxg()    const { return m_nxg; }
  int  nyg()    const { return m_nyg; }

  const std::vector<int>&    gindex() const { return m_gindex; }
  const std::vector<double>& lons()   const { return m_lons; }
  const std::vector<double>& lats()   const { return m_lats; }
  const std::vector<double>& areas()  const { return m_areas; }
  const std::vector<double>& masks()  const { return m_masks; }
  const std::vector<double>& fracs()  const { return m_fracs; }

protected:
  // Emulator virtual overrides
  void init_impl()       override;
  void run_impl(int dt)  override;
  void final_impl()      override;

private:
  int m_fcomm;    ///< Fortran communicator handle
  int m_comp_id;  ///< driver-assigned component id
  int m_nxg;      ///< global longitude count
  int m_nyg;      ///< global latitude count
  int m_lsize;    ///< local grid size on this PE

  std::vector<int>    m_gindex; ///< 1-based global indices
  std::vector<double> m_lons;   ///< longitude  (degrees)
  std::vector<double> m_lats;   ///< latitude   (degrees)
  std::vector<double> m_areas;  ///< cell area  (steradians)
  std::vector<double> m_masks;  ///< land/ocean mask (1.0)
  std::vector<double> m_fracs;  ///< grid fraction   (1.0)
};

} // namespace eatm

#endif // EATM_ATM_EMULATOR_HPP
