#include "shoc_f90.hpp"
#include "shoc_functions_f90.hpp"
#include "physics/share/physics_constants.hpp"
#include "shoc_ic_cases.hpp"

#include "ekat/ekat_assert.hpp"

using scream::Real;
using scream::Int;
extern "C" {
  void shoc_init_c(int nlev, Real gravit, Real rair, Real rh2o, Real cpair,
                   Real zvir, Real latvap, Real latice, Real karman);
  void shoc_use_cxx_c(bool use_cxx);
  void shoc_main_c(int shcol, int nlev, int nlevi, Real dtime, int nadv,
                   Real* host_dx, Real* host_dy, Real* thv, Real* zt_grid,
                   Real* zi_grid, Real* pres, Real* presi, Real* pdel,
                   Real* wthl_sfc, Real* wqw_sfc, Real* uw_sfc, Real* vw_sfc,
                   Real* wtracer_sfc, int num_qtracers, Real* w_field,
                   Real* exner, Real* phis, Real* host_dse, Real* tke,
                   Real* thetal, Real* qw, Real* u_wind, Real* v_wind,
                   Real* qtracers, Real* wthv_sec, Real* tkh, Real* tk,
                   Real* shoc_ql, Real* shoc_cldfrac, Real* pblh,
                   Real* shoc_mix, Real* isotropy, Real* w_sec, Real* thl_sec,
                   Real* qw_sec, Real* qwthl_sec, Real* wthl_sec, Real* wqw_sec,
                   Real* wtke_sec, Real* uw_sec, Real* vw_sec, Real* w3,
                   Real* wqls_sec, Real* brunt, Real* shoc_ql2);
}

namespace scream {
namespace shoc {

FortranData::FortranData(Int shcol_, Int nlev_, Int nlevi_,
                         Int num_qtracers_)
  : shcol(shcol_), nlev(nlev_), nlevi(nlevi_), num_qtracers(num_qtracers_)
{
  dtime = -1; // model time step [s]; set to invalid -1
  nadv = -1;  // depends on timestep, so also invalid

  // In variables
  host_dx = Array1("grid spacing of host model in x direction [m]", shcol);
  host_dy = Array1("grid spacing of host model in y direction [m]", shcol);
  zt_grid = Array2("heights for thermo grid [m]", shcol, nlev);
  zi_grid = Array2("heights for interface grid [m]", shcol, nlevi);
  pres = Array2("Pressure levels on thermo grid [Pa]", shcol, nlev);
  presi = Array2("Pressure levels on interface grid [Pa]", shcol, nlevi);
  pdel = Array2("Differences in pressure levels [Pa]", shcol, nlev);
  thv = Array2("Virtual potential temperature [K]", shcol, nlev);
  w_field = Array2("Large-scale vertical velocity [m/s]", shcol, nlev);
  wthl_sfc = Array1("Surface sensible heat flux [K m/s]", shcol);
  wqw_sfc = Array1("Surface latent heat flux [kg/kg m/s]", shcol);
  uw_sfc = Array1("Surface momentum flux (u direction) [m2/s2]", shcol);
  vw_sfc = Array1("Surface momentum flux (v direction) [m2/s2]", shcol);
  wtracer_sfc = Array2("Surface tracer flux [various]", shcol, num_qtracers);
  exner = Array2("Exner function [-]", shcol, nlev);
  phis = Array1("Host model surface geopotential height [m?]", shcol);

  // In/out variables
  host_dse = Array2("Dry static energy [J/kg]", shcol, nlev);
  tke = Array2("Turbulent kinetic energy [m2/s2]", shcol, nlev);
  thetal = Array2("Liquid water potential temperature [K]", shcol, nlev);
  qw = Array2("Total water mixing ratio [kg/kg]", shcol, nlev);
  u_wind = Array2("U wind component [m/s]", shcol, nlev);
  v_wind = Array2("V wind component [m/s]", shcol, nlev);
  wthv_sec = Array2("Buoyancy flux [K m/s]", shcol, nlev);
  qtracers = Array3("Tracers [varies]", shcol, nlev, num_qtracers);
  tk = Array2("Eddy coefficient for momentum [m2/s]", shcol, nlev);
  tkh = Array2("Eddy coefficent for heat [m2/s]", shcol, nlev);

  // Out variables (including diagnostics)
  shoc_cldfrac = Array2("Cloud fraction [-]", shcol, nlev);
  shoc_ql = Array2("Cloud liquid mixing ratio [kg/kg]", shcol, nlev);
  pblh = Array1("Planetary boundary layer depth [m]", shcol);
  shoc_mix = Array2("Turbulent length scale [m]", shcol, nlev);
  w_sec = Array2("Vertical velocity variance [m2/s2]", shcol, nlev);
  thl_sec = Array2("Temperature variance [K^2]", shcol, nlevi);
  qw_sec = Array2("Moisture variance [kg2/kg2]", shcol, nlevi);
  qwthl_sec = Array2("Temp moisture covariance [K kg/kg]", shcol, nlevi);
  wthl_sec = Array2("Vertical heat flux [K m/s]", shcol, nlevi);
  wqw_sec = Array2("Vertical moisture flux [K m/s]", shcol, nlevi);
  wtke_sec = Array2("Vertical tke flux [m3/s3]", shcol, nlevi);
  uw_sec = Array2("Vertical zonal momentum flux [m2/s2]", shcol, nlevi);
  vw_sec = Array2("Vertical meridional momentum flux [m2/s2]", shcol, nlevi);
  w3 = Array2("Third moment vertical velocity [m3/s3]", shcol, nlevi);
  wqls_sec = Array2("Liquid water flux [kg/kg m/s]", shcol, nlev);
  brunt = Array2("Brunt-Vaisala frequency [s-1]", shcol, nlev);
  isotropy = Array2("Return to isotropic timescale [s]", shcol, nlev);
  shoc_ql2 = Array2("Variance in liquid water [kg2/kg2]", shcol,nlev);
}

FortranDataIterator::FortranDataIterator (const FortranData::Ptr& d) {
  init(d);
}

void FortranDataIterator::init (const FortranData::Ptr& dp) {
  d_ = dp;
#define fdipb(name)                                                     \
  fields_.push_back({#name,                                             \
        2,                                                              \
        {d_->name.extent_int(0), d_->name.extent_int(1), d_->name.extent_int(2)}, \
        d_->name.data(),                                                \
        d_->name.size()})
  fdipb(host_dx); fdipb(host_dy);
  fdipb(zt_grid); fdipb(zi_grid); fdipb(pres); fdipb(presi); fdipb(pdel);
  fdipb(thv); fdipb(w_field);
  fdipb(wthl_sfc); fdipb(wqw_sfc); fdipb(uw_sfc); fdipb(vw_sfc);
  fdipb(wtracer_sfc); fdipb(exner);
  fdipb(phis);

  fdipb(host_dse); fdipb(tke); fdipb(thetal); fdipb(qw);
  fdipb(u_wind); fdipb(v_wind); fdipb(wthv_sec);
  fdipb(qtracers);
  fdipb(tk); fdipb(tkh);

  fdipb(shoc_cldfrac); fdipb(shoc_ql);
  fdipb(pblh);
  fdipb(shoc_mix); fdipb(w_sec); fdipb(thl_sec); fdipb(qw_sec);
  fdipb(qwthl_sec); fdipb(wthl_sec); fdipb(wqw_sec); fdipb(wtke_sec);
  fdipb(uw_sec); fdipb(vw_sec); fdipb(w3); fdipb(wqls_sec); fdipb(isotropy);
  fdipb(brunt); fdipb(shoc_ql2);
#undef fdipb
}

const FortranDataIterator::RawArray&
FortranDataIterator::getfield (Int i) const {
  EKAT_ASSERT(i >= 0 || i < nfield());
  return fields_[i];
}

void shoc_init(Int nlev, bool use_fortran) {
  static bool is_init = false;
  if (!is_init) {
    using Scalar = Real;
    using C = scream::physics::Constants<Scalar>;

    shoc_init_c((int)nlev, C::gravit, C::Rair, C::RH2O, C::Cpair, C::ZVIR,
                C::LatVap, C::LatIce, C::Karman);
    is_init = true;
  }
  shoc_use_cxx_c(!use_fortran);
}

void shoc_main(FortranData& d) {
  shoc_main_c((int)d.shcol, (int)d.nlev, (int)d.nlevi, d.dtime, (int)d.nadv,
              d.host_dx.data(), d.host_dy.data(), d.thv.data(),
              d.zt_grid.data(), d.zi_grid.data(), d.pres.data(), d.presi.data(),
              d.pdel.data(), d.wthl_sfc.data(), d.wqw_sfc.data(), d.uw_sfc.data(),
              d.vw_sfc.data(), d.wtracer_sfc.data(), (int)d.num_qtracers,
              d.w_field.data(), d.exner.data(), d.phis.data(), d.host_dse.data(),
              d.tke.data(), d.thetal.data(), d.qw.data(), d.u_wind.data(),
              d.v_wind.data(), d.qtracers.data(), d.wthv_sec.data(), d.tkh.data(),
              d.tk.data(), d.shoc_ql.data(), d.shoc_cldfrac.data(), d.pblh.data(),
              d.shoc_mix.data(), d.isotropy.data(), d.w_sec.data(),
              d.thl_sec.data(), d.qw_sec.data(), d.qwthl_sec.data(),
              d.wthl_sec.data(), d.wqw_sec.data(), d.wtke_sec.data(),
              d.uw_sec.data(), d.vw_sec.data(), d.w3.data(), d.wqls_sec.data(),
              d.brunt.data(), d.shoc_ql2.data() );
}

int test_FortranData () {
  Int shcol = 1;
  Int nlev = 128, num_tracers = 1;
  FortranData d(shcol, nlev, nlev+1, num_tracers);
  return 0;
}

int test_shoc_init (bool use_fortran) {
  Int nz = 160;
  shoc_init(nz, use_fortran);
  return 0;
}

namespace {

using KT     = KokkosTypes<HostDevice>;
using Scalar = Real;
using Array2 = typename KT::template lview<Scalar**>;
using Array3 = typename KT::template lview<Scalar***>;

// Returns a string representation of the given 2D array.
std::string array_as_string(const Array2& array)
{
  std::stringstream s;
  s << "npy.array([";
  int n0 = array.extent(0);
  int n1 = array.extent(1);
  for (int i = 0; i < n0; ++i)
  {
    s << '[';
    for (int j = 0; j < n1; ++j)
      s << array(i, j) << ", ";
    s << "],";
  }
  s << "])" << std::ends;
  return s.str();
}

// Returns a string representation of the given 3D array.
std::string array_as_string(const Array3& array)
{
  std::stringstream s;
  s << "npy.array([";
  int n0 = array.extent(0);
  int n1 = array.extent(1);
  int n2 = array.extent(2);
  for (int i = 0; i < n0; ++i)
  {
    s << '[';
    for (int j = 0; j < n1; ++j)
    {
      s << '[';
      for (int k = 0; k < n2; ++k)
        s << array(i, j, k) << ", ";
      s << "],";
    }
    s << "],";
  }
  s << "])" << std::ends;
  return s.str();
}

// Generates a Python script to the given path that plots data to a PDF file
// with the given prefix. Ѕilently fails if a file cannot be opened for writing.
void gen_plot_script(const std::vector<std::shared_ptr<FortranData> >& data,
                     const char* script_path, const char* pdf_prefix)
{
  FILE* fp = fopen(script_path, "w");
  if (fp != NULL)
  {
    fprintf(stderr, "Generating %s...\n", script_path);

    // Plotting code, cribbed from scream-docs/shoc-port/shocintr.py.
    const char* plotting_code =
      "import numpy as npy\n"
      "import matplotlib.pyplot as pl\n\n"
      "def dispfig(fn_prefix=None, format='pdf', tight=True):\n"
      "    if tight: pl.tight_layout()\n"
      "    if not fn_prefix or len(fn_prefix) == 0:\n"
      "        return pl.show()\n"
      "    else:\n"
      "        pl.savefig(fn_prefix + '.' + format, format=format, bbox_inches='tight')\n\n"
      "def pad_lim(lim, pad=0.05, mult=False):\n"
      "    if mult:\n"
      "        v = lim[0] * (1 - pad), lim[1] * (1 + pad)\n"
      "    else:\n"
      "        d = lim[1] - lim[0]\n"
      "        delta = pad * d\n"
      "        v = lim[0] - delta, lim[1] + delta\n"
      "    return v\n\n"
      "def axis_tight_pad(pad=0.05, mult=False):\n"
      "    pl.axis('tight')\n"
      "    xl = pl.xlim()\n"
      "    yl = pl.ylim()\n"
      "    pl.xlim(pad_lim(xl, pad, mult))\n"
      "    return pl.ylim(pad_lim(yl, pad, mult))\n\n"
      "class pl_plot:\n"
      "    def __init__(self, figsize, filename, format=None, tight=True):\n"
      "        self.filename = filename\n"
      "        self.format = 'pdf' if not None else format\n"
      "        self.tight = tight\n"
      "        pl.close()\n"
      "        pl.figure(num=1, figsize=figsize)\n"
      "    def cleanup(self):\n"
      "        dispfig(self.filename, format=self.format, tight=self.tight)\n"
      "    def __enter__(self): return self\n"
      "    def __exit__(self, *args): pass\n"
      "    def __del__(self): return self.cleanup()\n\n"
      "def plot_basics(states, filename):\n"
      "    axs = []\n"
      "    def _plot_basics(s, first):\n"
      "        def vec(a): return a.reshape(a.size)\n"
      "        plotno = [1, 1, 2, 3, 3, 4, 4, 5, 5, 5, 6, 7, 8, 9, 10,\n"
      "                  11, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22 ]\n"
      "        grids  = [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,\n"
      "                  1,  1,  1,  1,  0,  1,  1 , 1,  0,  0,  0,  0 , 0  ]\n"
      "        zs = s['zt_grid'][0], s['zi_grid'][0]\n"
      "        fields = ['pres', 'presi', 'pdel', 'thv', 'thetal', 'shoc_ql',\n"
      "                  'qw', 'u_wind', 'v_wind', 'w_field', 'tke', 'tkh',\n"
      "                  'shoc_mix', 'isotropy', 'wtke_sec', 'uw_sec', 'vw_sec', 'wthl_sec',\n"
      "                  'wqw_sec', 'w_sec', 'thl_sec', 'qw_sec',\n"
      "                  'w3', 'wqls_sec', 'brunt', 'qtracers', 'host_dse', 'exner','shoc_ql2']\n"
      "        for i in range(len(plotno)):\n"
      "            if i == 0 or plotno[i] != plotno[i-1]:\n"
      "                if first: axs.append(pl.subplot(5, 5, plotno[i]))\n"
      "                else: pl.sca(axs[plotno[i]-1])\n"
      "            name = fields[i]\n"
      "            if name == 'qtracers':\n"
      "                y = s[name][0,:,0]\n"
      "            else:\n"
      "                y = s[name]\n"
      "            pl.plot(vec(y), 1e-3 * zs[grids[i]], '-',\n"
      "                    label=name if first else None)\n"
      "            if first: pl.legend(loc='best', fontsize=12)\n"
      "            axis_tight_pad()\n\n"
      "    with pl_plot((20, 20), filename):\n"
      "        for i, s in enumerate(states):\n"
      "            _plot_basics(s, i == 0)\n\n";
    fprintf(fp, "%s", plotting_code);

    // Write out state data (a list of dicts).

    fprintf(fp, "data = [\n");
    for (size_t i = 0; i < data.size(); ++i)
    {
#define WRITE_FIELD(field) \
    auto field = array_as_string(d.field); \
    fprintf(fp, "        '%s': %s,\n", #field, field.c_str())
      const FortranData& d = *(data[i]);
      fprintf(fp, "    {\n");
      WRITE_FIELD(zt_grid);
      WRITE_FIELD(zi_grid);
      WRITE_FIELD(pres);
      WRITE_FIELD(presi);
      WRITE_FIELD(pdel);
      WRITE_FIELD(thv);
      WRITE_FIELD(thetal);
      WRITE_FIELD(shoc_ql);
      WRITE_FIELD(qw);
      WRITE_FIELD(u_wind);
      WRITE_FIELD(v_wind);
      WRITE_FIELD(w_field);
      WRITE_FIELD(tke);
      WRITE_FIELD(tkh);
      WRITE_FIELD(shoc_mix);
      WRITE_FIELD(isotropy);
      WRITE_FIELD(wtke_sec);
      WRITE_FIELD(uw_sec);
      WRITE_FIELD(vw_sec);
      WRITE_FIELD(wthl_sec);
      WRITE_FIELD(wqw_sec);
      WRITE_FIELD(w_sec);
      WRITE_FIELD(thl_sec);
      WRITE_FIELD(qw_sec);
      WRITE_FIELD(w3);
      WRITE_FIELD(wqls_sec);
      WRITE_FIELD(brunt);
      WRITE_FIELD(qtracers);
      WRITE_FIELD(host_dse);
      WRITE_FIELD(exner);
      WRITE_FIELD(shoc_ql2);
      fprintf(fp, "    },\n");
#undef WRITE_FIELD
    }
    fprintf(fp, "]\n\n");

    // Write the plot call.
    fprintf(fp, "print('Writing %s.pdf...')\n", pdf_prefix);
    fprintf(fp, "plot_basics(data, \"%s\")\n", pdf_prefix);
    fclose(fp);
  }
  else
    fprintf(stderr, "Couldn't open %s for writing (skipping).", script_path);
}

} // end anonymous namespace

int test_shoc_ic (bool use_fortran, bool gen_plot_scripts) {
  Int nz = 160;
  shoc_init(nz, use_fortran);
  // Here we:
  // 1. Initialize a standard case with settings identical to
  //    scream-doc/ѕhoc_port/shocintr.py's example_run_case method
  auto d = ic::Factory::create(ic::Factory::standard, 1, nz, 1);

  // 4. Generate a Python script that plots the initial conditions.
  {
    std::vector<std::shared_ptr<FortranData> > ics;
    ics.push_back(d);
    if (gen_plot_scripts) {
      gen_plot_script(ics, "plot_ics.py", "ics");
    }
  }

  // 3. Run 100 steps, each of size dtime = 10 (as in that method)
  d->nadv = 100;
  d->dtime = 10;
  shoc_main(*d);

  // 4. Generate a Python script that plots the results.
  {
    std::vector<std::shared_ptr<FortranData> > data;
    data.push_back(d);
    if (gen_plot_scripts) {
      gen_plot_script(data, "plot_results.py", "results");
      fprintf(stderr, "To generate PDFs of initial and final conditions, run "
                      "the plot_*.py python scripts using a Python 3 interpreter "
                      "with numpy and Matplotlib installed.\n");
    }
  }

  return 0;
}

} // namespace shoc
} // namespace scream
