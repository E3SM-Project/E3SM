

# File atm\_coupling.cpp

[**File List**](files.md) **>** [**components**](dir_409f97388efe006bc3438b95e9edef48.md) **>** [**emulator\_comps**](dir_cd6ef227c082afa5b90fe3621cc9f093.md) **>** [**eatm**](dir_54689134e1a693092e83f56806593839.md) **>** [**src**](dir_1c3b735e18de9b9534f50214e18facf2.md) **>** [**impl**](dir_6975f7b28201ba7a9e865ff30c48a340.md) **>** [**atm\_coupling.cpp**](atm__coupling_8cpp.md)

[Go to the documentation of this file](atm__coupling_8cpp.md)


```C++


#include "atm_coupling.hpp"

namespace emulator {
namespace impl {

void AtmCouplingIndices::initialize(CouplingFieldsBase &fields) {
  // Export indices (a2x)
  Sa_z = fields.get_export_index("Sa_z");
  Sa_u = fields.get_export_index("Sa_u");
  Sa_v = fields.get_export_index("Sa_v");
  Sa_tbot = fields.get_export_index("Sa_tbot");
  Sa_ptem = fields.get_export_index("Sa_ptem");
  Sa_shum = fields.get_export_index("Sa_shum");
  Sa_dens = fields.get_export_index("Sa_dens");
  Sa_pbot = fields.get_export_index("Sa_pbot");
  Sa_pslv = fields.get_export_index("Sa_pslv");
  Faxa_lwdn = fields.get_export_index("Faxa_lwdn");
  Faxa_rainc = fields.get_export_index("Faxa_rainc");
  Faxa_rainl = fields.get_export_index("Faxa_rainl");
  Faxa_snowc = fields.get_export_index("Faxa_snowc");
  Faxa_snowl = fields.get_export_index("Faxa_snowl");
  Faxa_swndr = fields.get_export_index("Faxa_swndr");
  Faxa_swvdr = fields.get_export_index("Faxa_swvdr");
  Faxa_swndf = fields.get_export_index("Faxa_swndf");
  Faxa_swvdf = fields.get_export_index("Faxa_swvdf");
  Faxa_swnet = fields.get_export_index("Faxa_swnet");

  // Import indices (x2a)
  Sx_t = fields.get_import_index("Sx_t");
  So_t = fields.get_import_index("So_t");
  Faxx_sen = fields.get_import_index("Faxx_sen");
  Faxx_lat = fields.get_import_index("Faxx_lat");
  Faxx_taux = fields.get_import_index("Faxx_taux");
  Faxx_tauy = fields.get_import_index("Faxx_tauy");
  Faxx_lwup = fields.get_import_index("Faxx_lwup");
  Faxx_evap = fields.get_import_index("Faxx_evap");
  Sx_avsdr = fields.get_import_index("Sx_avsdr");
  Sx_anidr = fields.get_import_index("Sx_anidr");
  Sx_avsdf = fields.get_import_index("Sx_avsdf");
  Sx_anidf = fields.get_import_index("Sx_anidf");
  Sl_snowh = fields.get_import_index("Sl_snowh");
  Si_snowh = fields.get_import_index("Si_snowh");
  Sx_tref = fields.get_import_index("Sx_tref");
  Sx_qref = fields.get_import_index("Sx_qref");
  Sx_u10 = fields.get_import_index("Sx_u10");
  Sf_ifrac = fields.get_import_index("Sf_ifrac");
  Sf_ofrac = fields.get_import_index("Sf_ofrac");
  Sf_lfrac = fields.get_import_index("Sf_lfrac");
}

void import_atm_fields(const double *import_data, int ncols, int nfields,
                       const AtmCouplingIndices &idx, AtmFieldManager &fields) {
  for (int i = 0; i < ncols; ++i) {
    if (idx.Sx_t >= 0)
      fields.ts[i] = import_data[i * nfields + idx.Sx_t];
    if (idx.So_t >= 0)
      fields.sst[i] = import_data[i * nfields + idx.So_t];
    if (idx.Faxx_sen >= 0)
      fields.shf[i] = import_data[i * nfields + idx.Faxx_sen];
    if (idx.Faxx_lat >= 0)
      fields.lhf[i] = import_data[i * nfields + idx.Faxx_lat];
    if (idx.Faxx_taux >= 0)
      fields.wsx[i] = import_data[i * nfields + idx.Faxx_taux];
    if (idx.Faxx_tauy >= 0)
      fields.wsy[i] = import_data[i * nfields + idx.Faxx_tauy];
    if (idx.Faxx_lwup >= 0)
      fields.lwup[i] = import_data[i * nfields + idx.Faxx_lwup];
    if (idx.Sx_avsdr >= 0)
      fields.asdir[i] = import_data[i * nfields + idx.Sx_avsdr];
    if (idx.Sx_anidr >= 0)
      fields.aldir[i] = import_data[i * nfields + idx.Sx_anidr];
    if (idx.Sx_avsdf >= 0)
      fields.asdif[i] = import_data[i * nfields + idx.Sx_avsdf];
    if (idx.Sx_anidf >= 0)
      fields.aldif[i] = import_data[i * nfields + idx.Sx_anidf];
    if (idx.Sl_snowh >= 0)
      fields.snowhland[i] = import_data[i * nfields + idx.Sl_snowh];
    if (idx.Si_snowh >= 0)
      fields.snowhice[i] = import_data[i * nfields + idx.Si_snowh];
    if (idx.Sf_ifrac >= 0)
      fields.icefrac[i] = import_data[i * nfields + idx.Sf_ifrac];
    if (idx.Sf_ofrac >= 0)
      fields.ocnfrac[i] = import_data[i * nfields + idx.Sf_ofrac];
    if (idx.Sf_lfrac >= 0)
      fields.lndfrac[i] = import_data[i * nfields + idx.Sf_lfrac];
  }
}

void export_atm_fields(double *export_data, int ncols, int nfields,
                       const AtmCouplingIndices &idx,
                       const AtmFieldManager &fields) {
  for (int i = 0; i < ncols; ++i) {
    if (idx.Sa_z >= 0)
      export_data[i * nfields + idx.Sa_z] = fields.zbot[i];
    if (idx.Sa_u >= 0)
      export_data[i * nfields + idx.Sa_u] = fields.ubot[i];
    if (idx.Sa_v >= 0)
      export_data[i * nfields + idx.Sa_v] = fields.vbot[i];
    if (idx.Sa_tbot >= 0)
      export_data[i * nfields + idx.Sa_tbot] = fields.tbot[i];
    if (idx.Sa_ptem >= 0)
      export_data[i * nfields + idx.Sa_ptem] = fields.ptem[i];
    if (idx.Sa_shum >= 0)
      export_data[i * nfields + idx.Sa_shum] = fields.shum[i];
    if (idx.Sa_dens >= 0)
      export_data[i * nfields + idx.Sa_dens] = fields.dens[i];
    if (idx.Sa_pbot >= 0)
      export_data[i * nfields + idx.Sa_pbot] = fields.pbot[i];
    if (idx.Sa_pslv >= 0)
      export_data[i * nfields + idx.Sa_pslv] = fields.pslv[i];
    if (idx.Faxa_lwdn >= 0)
      export_data[i * nfields + idx.Faxa_lwdn] = fields.lwdn[i];
    if (idx.Faxa_rainc >= 0)
      export_data[i * nfields + idx.Faxa_rainc] = fields.rainc[i];
    if (idx.Faxa_rainl >= 0)
      export_data[i * nfields + idx.Faxa_rainl] = fields.rainl[i];
    if (idx.Faxa_snowc >= 0)
      export_data[i * nfields + idx.Faxa_snowc] = fields.snowc[i];
    if (idx.Faxa_snowl >= 0)
      export_data[i * nfields + idx.Faxa_snowl] = fields.snowl[i];
    if (idx.Faxa_swndr >= 0)
      export_data[i * nfields + idx.Faxa_swndr] = fields.swndr[i];
    if (idx.Faxa_swvdr >= 0)
      export_data[i * nfields + idx.Faxa_swvdr] = fields.swvdr[i];
    if (idx.Faxa_swndf >= 0)
      export_data[i * nfields + idx.Faxa_swndf] = fields.swndf[i];
    if (idx.Faxa_swvdf >= 0)
      export_data[i * nfields + idx.Faxa_swvdf] = fields.swvdf[i];
    if (idx.Faxa_swnet >= 0)
      export_data[i * nfields + idx.Faxa_swnet] = fields.swnet[i];
  }
}

} // namespace impl
} // namespace emulator
```


