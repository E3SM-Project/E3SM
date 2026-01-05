

# File atm\_field\_manager.cpp

[**File List**](files.md) **>** [**components**](dir_409f97388efe006bc3438b95e9edef48.md) **>** [**emulator\_comps**](dir_cd6ef227c082afa5b90fe3621cc9f093.md) **>** [**eatm**](dir_54689134e1a693092e83f56806593839.md) **>** [**src**](dir_1c3b735e18de9b9534f50214e18facf2.md) **>** [**impl**](dir_6975f7b28201ba7a9e865ff30c48a340.md) **>** [**atm\_field\_manager.cpp**](atm__field__manager_8cpp.md)

[Go to the documentation of this file](atm__field__manager_8cpp.md)


```C++


#include "atm_field_manager.hpp"

namespace emulator {
namespace impl {

void AtmFieldManager::init_field_map() {
  if (!m_field_map.empty())
    return;

  // Register hardcoded fields
  m_field_map["shf"] = &shf;
  m_field_map["cflx"] = &cflx;
  m_field_map["lhf"] = &lhf;
  m_field_map["wsx"] = &wsx;
  m_field_map["wsy"] = &wsy;
  m_field_map["lwup"] = &lwup;
  m_field_map["asdir"] = &asdir;
  m_field_map["aldir"] = &aldir;
  m_field_map["asdif"] = &asdif;
  m_field_map["aldif"] = &aldif;
  m_field_map["ts"] = &ts;
  m_field_map["surface_temperature"] = &ts; // Alias
  m_field_map["sst"] = &sst;
  m_field_map["snowhland"] = &snowhland;
  m_field_map["snowhice"] = &snowhice;
  m_field_map["tref"] = &tref;
  m_field_map["qref"] = &qref;
  m_field_map["Q2m"] = &qref; // Alias
  m_field_map["u10"] = &u10;
  m_field_map["UGRD10m"] = &u10; // Alias
  m_field_map["u10withgusts"] = &u10withgusts;
  m_field_map["icefrac"] = &icefrac;
  m_field_map["sea_ice_fraction"] = &icefrac; // Alias
  m_field_map["ocnfrac"] = &ocnfrac;
  m_field_map["ocean_fraction"] = &ocnfrac; // Alias
  m_field_map["lndfrac"] = &lndfrac;
  m_field_map["land_fraction"] = &lndfrac; // Alias

  m_field_map["zbot"] = &zbot;
  m_field_map["ubot"] = &ubot;
  m_field_map["vbot"] = &vbot;
  m_field_map["tbot"] = &tbot;
  m_field_map["ptem"] = &ptem;
  m_field_map["shum"] = &shum;
  m_field_map["dens"] = &dens;
  m_field_map["pbot"] = &pbot;
  m_field_map["PRESsfc"] = &pbot; // Alias
  m_field_map["pslv"] = &pslv;
  m_field_map["lwdn"] = &lwdn;
  m_field_map["DLWRFsfc"] = &lwdn; // Alias
  m_field_map["rainc"] = &rainc;
  m_field_map["rainl"] = &rainl;
  m_field_map["snowc"] = &snowc;
  m_field_map["snowl"] = &snowl;
  m_field_map["swndr"] = &swndr;
  m_field_map["swvdr"] = &swvdr;
  m_field_map["swndf"] = &swndf;
  m_field_map["swvdf"] = &swvdf;
  m_field_map["swnet"] = &swnet;
  m_field_map["DSWRFtoa"] =
      &swnet; // Use swnet as placeholder for solar? Or dynamic?
}

std::vector<double> *AtmFieldManager::get_field_ptr(const std::string &name) {
  init_field_map();

  // Check hardcoded map first
  auto it = m_field_map.find(name);
  if (it != m_field_map.end()) {
    return it->second;
  }

  // Check dynamic fields
  auto dyn_it = dynamic_fields.find(name);
  if (dyn_it != dynamic_fields.end()) {
    return &(dyn_it->second);
  }

  return nullptr;
}

void AtmFieldManager::register_dynamic_field(const std::string &name) {
  init_field_map();
  if (m_field_map.find(name) != m_field_map.end())
    return; // Already hardcoded

  if (dynamic_fields.find(name) == dynamic_fields.end()) {
    dynamic_fields[name] = std::vector<double>();
    if (m_allocated && m_ncols > 0) {
      dynamic_fields[name].resize(m_ncols);
    }
  }
}

void AtmFieldManager::allocate(int ncols) {
  if (m_allocated) {
    deallocate();
  }

  m_ncols = ncols;

  // AI model buffers - sized dynamically later if needed, but safe fast alloc
  // here net_inputs/outputs sizing moved to emulator_atm logic or just cleared
  net_inputs.clear();
  net_outputs.clear();

  // Import fields (x2a) - allocated only, no defaults
  shf.resize(ncols);
  cflx.resize(ncols);
  lhf.resize(ncols);
  wsx.resize(ncols);
  wsy.resize(ncols);
  lwup.resize(ncols);
  asdir.resize(ncols);
  aldir.resize(ncols);
  asdif.resize(ncols);
  aldif.resize(ncols);
  ts.resize(ncols);
  sst.resize(ncols);
  snowhland.resize(ncols);
  snowhice.resize(ncols);
  tref.resize(ncols);
  qref.resize(ncols);
  u10.resize(ncols);
  u10withgusts.resize(ncols);
  icefrac.resize(ncols);
  ocnfrac.resize(ncols);
  lndfrac.resize(ncols);

  // Export fields (a2x) - allocated only, filled from IC file
  zbot.resize(ncols);
  ubot.resize(ncols);
  vbot.resize(ncols);
  tbot.resize(ncols);
  ptem.resize(ncols);
  shum.resize(ncols);
  dens.resize(ncols);
  pbot.resize(ncols);
  pslv.resize(ncols);
  lwdn.resize(ncols);
  rainc.resize(ncols);
  rainl.resize(ncols);
  snowc.resize(ncols);
  snowl.resize(ncols);
  swndr.resize(ncols);
  swvdr.resize(ncols);
  swndf.resize(ncols);
  swvdf.resize(ncols);
  swnet.resize(ncols);

  // Allocate dynamic fields
  for (auto &kv : dynamic_fields) {
    kv.second.resize(ncols);
  }

  m_allocated = true;
}

void AtmFieldManager::deallocate() {
  net_inputs.clear();
  net_outputs.clear();

  shf.clear();
  cflx.clear();
  lhf.clear();
  wsx.clear();
  wsy.clear();
  lwup.clear();
  asdir.clear();
  aldir.clear();
  asdif.clear();
  aldif.clear();
  ts.clear();
  sst.clear();
  snowhland.clear();
  snowhice.clear();
  tref.clear();
  qref.clear();
  u10.clear();
  u10withgusts.clear();
  icefrac.clear();
  ocnfrac.clear();
  lndfrac.clear();

  zbot.clear();
  ubot.clear();
  vbot.clear();
  tbot.clear();
  ptem.clear();
  shum.clear();
  dens.clear();
  pbot.clear();
  pslv.clear();
  lwdn.clear();
  rainc.clear();
  rainl.clear();
  snowc.clear();
  snowl.clear();
  swndr.clear();
  swvdr.clear();
  swndf.clear();
  swvdf.clear();
  swnet.clear();

  for (auto &kv : dynamic_fields) {
    kv.second.clear();
  }

  m_allocated = false;
  m_ncols = 0;
}

void AtmFieldManager::set_defaults(int ncols) {
  for (int i = 0; i < ncols; ++i) {
    zbot[i] = 10.0;     // 10m reference height [m]
    ubot[i] = 5.0;      // Zonal wind [m/s]
    vbot[i] = 0.0;      // Meridional wind [m/s]
    tbot[i] = 288.0;    // ~15°C [K]
    ptem[i] = 288.0;    // Potential temp [K]
    shum[i] = 0.008;    // Specific humidity [kg/kg]
    dens[i] = 1.225;    // Air density [kg/m3]
    pbot[i] = 101325.0; // Surface pressure [Pa]
    pslv[i] = 101325.0; // Sea level pressure [Pa]
    lwdn[i] = 300.0;    // Downward LW [W/m2]
    rainc[i] = 0.0;     // Convective rain [kg/m2/s]
    rainl[i] = 0.0;     // Large-scale rain [kg/m2/s]
    snowc[i] = 0.0;     // Convective snow [kg/m2/s]
    snowl[i] = 0.0;     // Large-scale snow [kg/m2/s]
    swnet[i] = 200.0;   // Net SW [W/m2]

    // Partition swnet into components
    swvdr[i] = swnet[i] * 0.28; // Visible direct
    swndr[i] = swnet[i] * 0.31; // NIR direct
    swvdf[i] = swnet[i] * 0.24; // Visible diffuse
    swndf[i] = swnet[i] * 0.17; // NIR diffuse
  }
}

} // namespace impl
} // namespace emulator
```


