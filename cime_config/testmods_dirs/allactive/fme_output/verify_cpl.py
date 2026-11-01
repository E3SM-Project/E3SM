#!/usr/bin/env python
"""
verify_cpl.py -- readiness check for coupler-native FME output.

Validates the as-exchanged merged-forcing streams written by
driver-mct/main/cpl_fme_mod.F90:

    <case>.cpl.fme.x2o1D.YYYY-MM.nc   x2o ocean merged forcing, 1-day mean
    <case>.cpl.fme.x2o5D.YYYY-MM.nc   x2o, 5-day mean
    <case>.cpl.fme.xao1D/5D.*.nc      xao atm-ocean bulk fluxes (Faox_*)
    <case>.cpl.fme.x2i1D/5D.*.nc      x2i sea-ice merged import

Checks (mirroring verify_production.py conventions):
  1. File inventory per stream.
  2. CF time scheme: lowercase `time` (unlimited) + `time_bnds`(nbnd,time)
     + `date`/`datesec`, units 'days since ...', calendar present; time =
     END of window (== time_bnds[...,1]); no legacy xtime.
  3. Lat-lon grid: `lon`/`lat` dims + coord vars, plausible ranges.
  4. Expected merged-forcing fields present for the bundle.
  5. Fill-leak scan: every data cell is finite and either physical or the
     canonical SHR_FILL_VALUE=1e20 (land) -- nothing in the fractional-fill
     band [pollution_threshold, 1e20) (gotcha #36).
  6. 5D-vs-mean(1D) cross-check (p99-tolerant, gotcha #42) when both the
     5D and 1D files for a stream/bundle are present.
  7. Time layout (0301 spec check 3): with cpl_fme_state_mean on (global
     attr fme_avg_mode = 'mean', or --require-state-mean) no field of a 1D/5D
     stream may carry cell_methods 'time: point'.  The *Dinst companions are
     the exception: they must carry ONLY state (S*) fields, all 'time: point'.
     Flux (F*) fields should carry an `area_basis` attribute.
  8. *Dinst companions (cpl_fme_inst_enable): present whenever the parent
     names one (global attr fme_inst_companion, or --require-inst), on the
     parent's exact time / time_bnds, and 5Dinst == 1Dinst at shared window
     ends (same bundle, map and phase -> identical values).
  9. 0301 spec check 5, *5Dinst state vs the 6-hourly a2x/x2a sample at the
     window end.  Exact equality is only expected where the two files share a
     source bundle AND remap file; none of the ocean-mesh pairs do (e.g.
     x2i Sa_tbot is a2x Sa_tbot mapped atm->ocean mesh by the coupler, then
     remapped with the OCEAN map; a2x is remapped from the atm grid with the
     ATM map), so they differ by the two remaps.  For those the check reports
     rms(inst - 6h) next to rms(5D mean - 6h): an end-of-window instant must
     sit much closer to the 6h end sample than the window mean does (0201:
     0.22-0.25 K vs 1.85-2.29 K for Sa_tbot), else FAIL.

Exit code is nonzero if any FAIL (or a PRIMARY stream is missing), 0 otherwise
-- so it gates production like verify_production.py. Writes an HTML summary to
<outdir>/verify_cpl/index.html and, with --share-dir, copies it under
<share-dir>/<label>/ and additively links it from <share-dir>/index.html.

Usage:
  micromamba run -n xgns python verify_cpl.py --rundir $RUNDIR \
      --outdir /tmp/fme-verify --share-dir $PSCRATCH/share --label mycase
"""

import argparse
import glob
import os
import re
import sys

import numpy as np

try:
    import xarray as xr
except ImportError:
    sys.stderr.write("verify_cpl.py requires xarray (micromamba run -n xgns ...)\n")
    raise

SHR_FILL_VALUE = 1e20

# Default field lists must track cpl_fme_mod.F90 (x2o_default / xao_default /
# x2i_default). Absent fields are skipped at runtime, so we treat these as the
# *maximal* expected set and only FAIL if NONE of a bundle's fields appear.
EXPECTED = {
    "x2o": ["Foxx_taux", "Foxx_tauy", "Foxx_sen", "Foxx_lat", "Foxx_lwup",
            "Foxx_evap", "Foxx_swnet", "Faxa_lwdn", "Faxa_rain", "Faxa_snow",
            "Fioi_melth", "Fioi_meltw", "Fioi_salt", "Foxx_rofl", "Foxx_rofi",
            "Sa_pslv", "Si_ifrac"],
    "xao": ["Faox_taux", "Faox_tauy", "Faox_lat", "Faox_sen", "Faox_lwup",
            "Faox_evap", "Faox_swdn", "Faox_swup", "So_tref", "So_qref",
            "So_u10", "So_ustar", "So_duu10n"],
    "x2i": ["Sa_z", "Sa_u", "Sa_v", "Sa_tbot", "Sa_ptem", "Sa_shum", "Sa_pbot",
            "Sa_dens", "Faxa_rain", "Faxa_snow", "Faxa_lwdn", "Faxa_swndr",
            "Faxa_swvdr", "Faxa_swndf", "Faxa_swvdf", "So_t", "So_s", "So_u",
            "So_v", "So_dhdx", "Fioo_q"],
    "o2x": ["So_t", "So_s", "So_u", "So_v", "So_dhdx", "So_dhdy", "So_ssh",
            "So_bldepth", "So_fswpen", "Fioo_q"],
    "i2x": ["Si_ifrac", "Si_t", "Si_avsdr", "Si_anidr", "Si_avsdf", "Si_anidf",
            "Si_tref", "Si_qref", "Si_u10", "Si_snowh", "Faii_taux", "Faii_tauy",
            "Faii_lat", "Faii_sen", "Faii_lwup", "Faii_evap", "Fioi_melth",
            "Fioi_meltw", "Fioi_salt", "Fioi_swpen"],
    "x2a": ["Sf_lfrac", "Sf_ifrac", "Sf_ofrac", "Sx_avsdr", "Sx_anidr",
            "Sx_avsdf", "Sx_anidf", "Sx_tref", "Sx_qref", "So_t", "Sx_t",
            "Sx_u10", "Faxx_taux", "Faxx_tauy", "Faxx_lat", "Faxx_sen",
            "Faxx_lwup", "Faxx_evap", "So_ssq"],
    "a2x": ["Sa_z", "Sa_u", "Sa_v", "Sa_tbot", "Sa_ptem", "Sa_shum", "Sa_pbot",
            "Sa_pslv", "Sa_dens", "Faxa_lwdn", "Faxa_rainc", "Faxa_rainl",
            "Faxa_snowc", "Faxa_snowl", "Faxa_swndr", "Faxa_swvdr",
            "Faxa_swndf", "Faxa_swvdf"],
    "x2l": ["Sa_z", "Sa_topo", "Sa_u", "Sa_v", "Sa_tbot", "Sa_ptem", "Sa_shum",
            "Sa_pbot", "Faxa_rainc", "Faxa_rainl", "Faxa_snowc", "Faxa_snowl",
            "Faxa_lwdn", "Faxa_swndr", "Faxa_swvdr", "Faxa_swndf", "Faxa_swvdf",
            "Flrr_flood", "Flrr_volr", "Flrr_volrmch", "Flrr_supply",
            "Flrr_deficit"],
    "l2x": ["Sl_avsdr", "Sl_anidr", "Sl_avsdf", "Sl_anidf", "Sl_tref", "Sl_qref",
            "Sl_t", "Sl_fv", "Sl_ram1", "Sl_snowh", "Sl_u10", "Fall_taux",
            "Fall_tauy", "Fall_lat", "Fall_sen", "Fall_lwup", "Fall_evap",
            "Fall_swnet", "Flrl_rofsur", "Flrl_rofgwl", "Flrl_rofsub",
            "Flrl_rofdto", "Flrl_rofi"],
    "x2r": ["Flrl_rofsur", "Flrl_rofgwl", "Flrl_rofsub", "Flrl_rofdto",
            "Flrl_rofi", "Flrl_demand", "Flrl_Tqsur", "Flrl_Tqsub"],
    "r2x": ["Forr_rofl", "Forr_rofi", "Firr_rofi", "Flrr_flood", "Flrr_volr",
            "Flrr_volrmch", "Flrr_supply", "Flrr_deficit"],
}

# The *Dinst companions carry only the state (S*) subset of their bundle.
OCN_MESH_BUNDLES = ("x2o", "xao", "x2i", "o2x", "i2x")
for _b in OCN_MESH_BUNDLES:
    EXPECTED[_b + "_inst"] = [v for v in EXPECTED[_b] if v.startswith("S")]

# Spec check 5 pairs: (inst bundle, field) -> (6h stream, field).  Only state
# fields that are the same physical quantity on both sides.
INST_VS_6H = {
    ("x2i", "Sa_z"): ("a2x", "Sa_z"), ("x2i", "Sa_u"): ("a2x", "Sa_u"),
    ("x2i", "Sa_v"): ("a2x", "Sa_v"), ("x2i", "Sa_tbot"): ("a2x", "Sa_tbot"),
    ("x2i", "Sa_ptem"): ("a2x", "Sa_ptem"), ("x2i", "Sa_shum"): ("a2x", "Sa_shum"),
    ("x2i", "Sa_pbot"): ("a2x", "Sa_pbot"), ("x2i", "Sa_dens"): ("a2x", "Sa_dens"),
    ("x2o", "Sa_pslv"): ("a2x", "Sa_pslv"),
    ("o2x", "So_t"): ("x2a", "So_t"), ("x2i", "So_t"): ("x2a", "So_t"),
}

# stream label -> (bundle, cadence). PRIMARY (1D) streams must exist for a
# bundle that is expected on; 5D companions are optional.
STREAMS = [
    ("x2o1D", "x2o", "1D"), ("x2o5D", "x2o", "5D"),
    ("xao1D", "xao", "1D"), ("xao5D", "xao", "5D"),
    ("x2i1D", "x2i", "1D"), ("x2i5D", "x2i", "5D"),
    ("o2x1D", "o2x", "1D"), ("o2x5D", "o2x", "5D"),
    ("i2x1D", "i2x", "1D"), ("i2x5D", "i2x", "5D"),
    # window-end instant companions (cpl_fme_inst_enable); never PRIMARY
    ("x2o1Dinst", "x2o_inst", "1Dinst"), ("x2o5Dinst", "x2o_inst", "5Dinst"),
    ("xao1Dinst", "xao_inst", "1Dinst"), ("xao5Dinst", "xao_inst", "5Dinst"),
    ("x2i1Dinst", "x2i_inst", "1Dinst"), ("x2i5Dinst", "x2i_inst", "5Dinst"),
    ("o2x1Dinst", "o2x_inst", "1Dinst"), ("o2x5Dinst", "o2x_inst", "5Dinst"),
    ("i2x1Dinst", "i2x_inst", "1Dinst"), ("i2x5Dinst", "i2x_inst", "5Dinst"),
    # atm streams: single configurable sub-daily cadence (no 5D companion)
    ("x2a", "x2a", "1D"), ("a2x", "a2x", "1D"),
    # lnd / rof streams: single configurable cadence, same pattern as the atm.
    # The rof pair samples only on river coupling steps (ROF_NCPL), so its
    # records carry fewer samples per window than the land pair.
    ("x2l", "x2l", "1D"), ("l2x", "l2x", "1D"),
    ("x2r", "x2r", "1D"), ("r2x", "r2x", "1D"),
]


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rundir", required=True, help="case run directory")
    ap.add_argument("--outdir", required=True, help="output dir for the report")
    ap.add_argument("--year", default=None, help="restrict to files matching .YYYY-")
    ap.add_argument("--pollution-threshold", type=float, default=1e15,
                    help="|x| in [thresh, 1e20) is flagged as fractional-fill "
                         "pollution (gotcha #36). Default 1e15.")
    ap.add_argument("--cross-check-rel-tol", type=float, default=5e-2,
                    help="5D-vs-mean(1D) p99 FAIL tolerance (fraction). Default 0.05")
    ap.add_argument("--require-state-mean", action="store_true",
                    help="FAIL any 'time: point' field in a 1D/5D stream even "
                         "if the file does not declare fme_avg_mode='mean' "
                         "(0301 spec check 3)")
    ap.add_argument("--require-inst", action="store_true",
                    help="FAIL if an enabled 1D/5D ocean-mesh stream has no "
                         "*Dinst companion files")
    ap.add_argument("--inst-ratio-max", type=float, default=0.5,
                    help="spec check 5: FAIL when rms(inst-6h) exceeds this "
                         "fraction of rms(5Dmean-6h). Default 0.5")
    ap.add_argument("--label", default=None, help="case label for the share dir")
    ap.add_argument("--share-dir", default=None,
                    help="copy the report under <share-dir>/<label>/ and link it")
    return ap.parse_args()


DATED_RE = re.compile(r"\.\d{4}-\d{2}\.nc$")


def stream_files(rundir, stream, year):
    pat = os.path.join(rundir, f"*.cpl.fme.{stream}.*.nc")
    # dated monthly files only: the glob also matches the restart sidecar
    # <case>.cpl.fme.<stream>.accum.nc, which sorts after the dated files
    files = sorted(f for f in glob.glob(pat) if DATED_RE.search(f))
    if year:
        files = [f for f in files if f".{year}-" in os.path.basename(f)]
    return files


def _bnds_2col(arr):
    """Return time_bnds as (time, 2) regardless of stored axis order."""
    a = np.asarray(arr)
    if a.ndim != 2:
        return None
    if a.shape[1] == 2:
        return a
    if a.shape[0] == 2:
        return a.T
    return None


def check_cf_time(ds):
    out = []
    if "time" not in ds.dims:
        return [("FAIL", "no lowercase `time` dim (pre-CF-rewrite / wrong file)")]
    for v in ("time_bnds", "date", "datesec"):
        if v not in ds.variables:
            out.append(("FAIL", f"missing CF time variable `{v}`"))
    if "xtime" in ds.variables:
        out.append(("FAIL", "legacy `xtime` present (should be CF time only)"))
    if "time" in ds.variables:
        units = ds["time"].attrs.get("units", "")
        if not units.startswith("days since "):
            out.append(("FAIL", f"time:units = '{units}' (expected 'days since ...')"))
        if not ds["time"].attrs.get("calendar"):
            out.append(("WARN", "time:calendar attribute missing"))
        if "time_bnds" in ds.variables:
            b = _bnds_2col(ds["time_bnds"].values)
            t = np.asarray(ds["time"].values)
            if b is not None and b.shape[0] == t.shape[0]:
                # time should equal the END of the window (matches EAM, #39)
                if not np.allclose(t, b[:, 1], rtol=0, atol=1e-6):
                    out.append(("WARN", "time != time_bnds[:,1] (END-of-window "
                                        "convention not satisfied)"))
                if np.any(b[:, 1] < b[:, 0]):
                    out.append(("FAIL", "time_bnds with hi < lo"))
    if not out:
        out.append(("OK", "CF time scheme: time/time_bnds/date/datesec present, "
                          "END-of-window, no xtime"))
    return out


def check_grid(ds):
    out = []
    for d in ("lon", "lat"):
        if d not in ds.dims:
            out.append(("FAIL", f"missing `{d}` dim"))
        if d not in ds.variables:
            out.append(("FAIL", f"missing `{d}` coordinate variable"))
    if "lon" in ds.variables:
        lon = ds["lon"].values
        if lon.min() < -360.1 or lon.max() > 360.1:
            out.append(("FAIL", f"lon out of range [{lon.min()},{lon.max()}]"))
    if "lat" in ds.variables:
        lat = ds["lat"].values
        if lat.min() < -90.1 or lat.max() > 90.1:
            out.append(("FAIL", f"lat out of range [{lat.min()},{lat.max()}]"))
    if not out:
        out.append(("OK", f"grid: lon/lat present "
                          f"({ds.dims.get('lon')}x{ds.dims.get('lat')})"))
    return out


def check_fields(ds, bundle):
    expected = EXPECTED[bundle]
    present = [v for v in expected if v in ds.variables]
    missing = [v for v in expected if v not in ds.variables]
    if not present:
        return [("FAIL", f"none of the {bundle} merged-forcing fields present")]
    out = [("OK", f"{len(present)}/{len(expected)} expected {bundle} fields present")]
    if missing:
        out.append(("WARN", f"{len(missing)} expected field(s) absent "
                            f"(compset-dependent; ok if intentional): "
                            f"{', '.join(missing[:8])}"
                            f"{'...' if len(missing) > 8 else ''}"))
    return out


def check_fills(ds, thresh):
    out = []
    spatial = [v for v in ds.data_vars
               if {"lon", "lat"} <= set(ds[v].dims)]
    nbad = 0
    for v in spatial:
        a = np.asarray(ds[v].values)
        absv = np.abs(a)
        n_nan = int(np.isnan(a).sum())
        n_pol = int(((absv >= thresh) & (absv < SHR_FILL_VALUE)).sum())
        if n_nan:
            out.append(("FAIL", f"'{v}': {n_nan} NaN cell(s) (use SHR_FILL_VALUE, "
                                f"not NaN)"))
            nbad += 1
        if n_pol:
            sample = absv[(absv >= thresh) & (absv < SHR_FILL_VALUE)].ravel()[:3]
            out.append(("FAIL", f"'{v}': {n_pol} fractional-fill cell(s) "
                                f"[{thresh:g},1e20) e.g. {sample} (gotcha #36)"))
            nbad += 1
    if not nbad:
        out.append(("OK", f"fill-leak scan: {len(spatial)} var(s) clean "
                          f"(finite; physical or canonical 1e20)"))
    return out


def cross_check_5d_1d(ds5, ds1, rel_tol):
    """p99(|5D - mean(overlapping 1D)|) / p99(|value|) per shared field."""
    out = []
    if "time_bnds" not in ds5.variables or "time_bnds" not in ds1.variables:
        return [("WARN", "5D-vs-1D cross-check: missing time_bnds")]
    b5 = _bnds_2col(ds5["time_bnds"].values)
    b1 = _bnds_2col(ds1["time_bnds"].values)
    if b5 is None or b1 is None:
        return [("WARN", "5D-vs-1D cross-check: time_bnds shape unrecognized")]
    shared = [v for v in ds5.data_vars
              if v in ds1.data_vars and {"lon", "lat", "time"} <= set(ds5[v].dims)]
    worst = 0.0
    nchk = 0
    for v in shared:
        a5 = np.asarray(ds5[v].values)
        a1 = np.asarray(ds1[v].values)
        for j in range(b5.shape[0]):
            lo, hi = b5[j, 0], b5[j, 1]
            idx = [i for i in range(b1.shape[0])
                   if b1[i, 1] > lo + 1e-9 and b1[i, 1] <= hi + 1e-9]
            # Only windows whose 1D records are ALL in this file: the 5D
            # window straddling a month start (e.g. Oct 28 - Nov 2) has its
            # first days in the previous month's file, and a partial mean of
            # 2 of 5 days is not comparable to the 5D mean.
            w1 = float(np.median(b1[:, 1] - b1[:, 0])) if b1.shape[0] else 0.0
            if len(idx) < 2 or w1 <= 0 or len(idx) != int(round((hi - lo) / w1)):
                continue
            v5 = a5[j]
            v1 = a1[idx]
            good = (np.isfinite(v1) & (np.abs(v1) < SHR_FILL_VALUE * 0.1)).all(axis=0)
            good &= np.isfinite(v5) & (np.abs(v5) < SHR_FILL_VALUE * 0.1)
            if not good.any():
                continue
            m1 = v1.mean(axis=0)
            diff = np.abs(v5[good] - m1[good])
            scale = np.percentile(np.abs(v5[good]), 99) or 1.0
            r = float(np.percentile(diff, 99) / scale)
            worst = max(worst, r)
            nchk += 1
    if nchk == 0:
        return [("WARN", "5D-vs-1D cross-check: no comparable windows yet")]
    status = "FAIL" if worst > rel_tol else "OK"
    out.append((status, f"5D-vs-mean(1D): worst p99 rel diff {worst:.2%} over "
                        f"{nchk} window(s) (tol {rel_tol:.0%})"))
    return out


# ---------------------------------------------------------------------------
# 0301 layout checks: window means, *Dinst companions, area_basis
# ---------------------------------------------------------------------------

def _valid(a):
    a = np.asarray(a, dtype=np.float64)
    return np.isfinite(a) & (np.abs(a) < SHR_FILL_VALUE * 0.1)


def _time_vars(ds):
    return [v for v in ds.data_vars
            if {"lon", "lat", "time"} <= set(ds[v].dims)]


def check_time_layout(ds, bundle, cad, require_state_mean):
    """Spec check 3 + area_basis.  `bundle` is the STREAMS bundle key."""
    out = []
    base = bundle[:-5] if bundle.endswith("_inst") else bundle
    if base not in OCN_MESH_BUNDLES:
        return out                     # 6h atm / lnd / rof: layout unchanged
    tv = _time_vars(ds)
    cm = {v: str(ds[v].attrs.get("cell_methods", "")).strip() for v in tv}
    points = [v for v in tv if cm[v] == "time: point"]
    if cad.endswith("inst"):
        nonpoint = [v for v in tv if cm[v] != "time: point"]
        fluxes = [v for v in tv if not v.startswith("S")]
        if fluxes:
            out.append(("FAIL", f"inst companion carries non-state field(s): "
                                f"{', '.join(fluxes[:8])}"))
        if nonpoint:
            out.append(("FAIL", f"inst companion field(s) not 'time: point': "
                                f"{', '.join(nonpoint[:8])}"))
        if not fluxes and not nonpoint:
            out.append(("OK", f"inst companion: {len(tv)} state field(s), all "
                              f"'time: point'"))
        return out
    mode = str(ds.attrs.get("fme_avg_mode", ""))
    if points:
        msg = (f"{len(points)} field(s) are window-end instants "
               f"('time: point'): {', '.join(points[:8])}"
               f"{'...' if len(points) > 8 else ''}")
        if mode == "mean":
            out.append(("FAIL", msg + " although fme_avg_mode='mean'"))
        elif require_state_mean:
            out.append(("FAIL", msg + " (spec check 3; set cpl_fme_state_mean)"))
        else:
            out.append(("WARN", msg + " (mixed layout; spec check 3 needs "
                                      "cpl_fme_state_mean=.true.)"))
    else:
        out.append(("OK", f"all {len(tv)} field(s) are window means "
                          f"(fme_avg_mode={mode or 'unset'})"))
    noab = [v for v in tv if v.startswith("F") and "area_basis" not in ds[v].attrs]
    if noab:
        out.append(("WARN", f"{len(noab)} flux field(s) lack an area_basis "
                            f"attribute: {', '.join(noab[:6])}"))
    unv = [v for v in tv if str(ds[v].attrs.get("area_basis", "")) == "unverified"]
    if unv:
        out.append(("WARN", f"area_basis 'unverified' for: {', '.join(unv[:6])}"))
    return out


def _bnds(ds):
    return _bnds_2col(ds["time_bnds"].values) if "time_bnds" in ds.variables else None


def check_inst_alignment(ds_inst, ds_parent):
    """The companion must sit on its parent's exact windows."""
    ti = np.asarray(ds_inst["time"].values)
    tp = np.asarray(ds_parent["time"].values)
    if ti.shape != tp.shape or not np.allclose(ti, tp, rtol=0, atol=1e-9):
        return [("FAIL", f"inst time axis != parent's ({ti.size} vs {tp.size} "
                         f"records)")]
    bi, bp = _bnds(ds_inst), _bnds(ds_parent)
    if bi is None or bp is None or bi.shape != bp.shape or \
            not np.allclose(bi, bp, rtol=0, atol=1e-9):
        return [("FAIL", "inst time_bnds != parent's")]
    return [("OK", f"inst on the parent's windows ({ti.size} record(s))")]


def _match_times(t_a, t_b, atol=1e-6):
    """Index pairs (i, j) with t_a[i] == t_b[j]."""
    pairs = []
    for i, t in enumerate(np.asarray(t_a)):
        j = np.where(np.abs(np.asarray(t_b) - t) <= atol)[0]
        if j.size:
            pairs.append((i, int(j[0])))
    return pairs


def check_inst_5d_vs_1d(ds5i, ds1i):
    """Same bundle, map and phase -> 5Dinst == 1Dinst at shared window ends."""
    pairs = _match_times(ds5i["time"].values, ds1i["time"].values)
    if not pairs:
        return [("WARN", "5Dinst vs 1Dinst: no shared window ends yet")]
    worst, worst_v = 0.0, ""
    for v in _time_vars(ds5i):
        if v not in ds1i.data_vars:
            continue
        a5 = np.asarray(ds5i[v].values, dtype=np.float64)
        a1 = np.asarray(ds1i[v].values, dtype=np.float64)
        for i, j in pairs:
            g5, g1 = _valid(a5[i]), _valid(a1[j])
            if (g5 != g1).any():
                return [("FAIL", f"5Dinst vs 1Dinst: '{v}' valid mask differs "
                                 f"at t={ds5i['time'].values[i]}")]
            if not g5.any():
                continue
            scale = float(np.max(np.abs(a5[i][g5]))) or 1.0
            d = float(np.max(np.abs(a5[i][g5] - a1[j][g5]))) / scale
            if d > worst:
                worst, worst_v = d, v
    status = "OK" if worst <= 1e-6 else "FAIL"
    return [(status, f"5Dinst vs 1Dinst at {len(pairs)} shared window end(s): "
                     f"max rel diff {worst:.2e}"
                     f"{' (' + worst_v + ')' if worst_v else ''} (expect 0)")]


def _rms_diff(a, b, idx_pairs):
    num, n = 0.0, 0
    for i, j in idx_pairs:
        g = _valid(a[i]) & _valid(b[j])
        if g.any():
            d = np.asarray(a[i], dtype=np.float64)[g] - np.asarray(b[j], dtype=np.float64)[g]
            num += float(np.sum(d * d))
            n += int(g.sum())
    return (np.sqrt(num / n) if n else None), n


def check_inst_vs_6h(inst_bundle, ds_inst, ds_mean, six, ratio_max):
    """Spec check 5.  `six` maps '6h stream name' -> dataset of the same month."""
    out = []
    for (b, v), (s6, v6) in INST_VS_6H.items():
        if b != inst_bundle or v not in ds_inst.data_vars:
            continue
        ds6 = six.get(s6)
        if ds6 is None or v6 not in ds6.data_vars:
            continue
        pairs = _match_times(ds_inst["time"].values, ds6["time"].values)
        if not pairs:
            out.append(("WARN", f"{v} vs {s6}.{v6}: no 6h record at any window end"))
            continue
        ai = np.asarray(ds_inst[v].values)
        a6 = np.asarray(ds6[v6].values)
        # both attrs must be present: files without them (pre-0301) would
        # otherwise compare None == None and demand round-off equality
        bi, b6 = ds_inst.attrs.get("fme_bundle"), ds6.attrs.get("fme_bundle")
        mi, m6 = ds_inst.attrs.get("remap_file"), ds6.attrs.get("remap_file")
        same_map = bool(bi and mi) and bi == b6 and mi == m6
        rms_i, n = _rms_diff(ai, a6, pairs)
        if rms_i is None:
            continue
        if same_map:
            st = "OK" if rms_i <= 1e-6 * (np.nanmax(np.abs(a6)) or 1.0) else "FAIL"
            out.append((st, f"{v} vs {s6}.{v6} (same bundle+map): rms {rms_i:.3e} "
                            f"over {n} cell-records (expect round-off)"))
            continue
        msg = (f"{v} vs {s6}.{v6} at {len(pairs)} window end(s): rms(inst-6h) "
               f"{rms_i:.4g}")
        # The ratio test needs a WINDOW-MEAN parent: in the mixed layout the
        # parent's S* fields are the same window-end instants (ratio ~1).
        # It is also only a hard gate for the fast atm states (Sa_*): for the
        # ocean states (So_*) the 5-day variability is tiny next to the
        # coastal remap residual (x2a So_t is ocn->atm mapped, 0 over land),
        # so rms(mean-6h) ~ rms(inst-6h) and the number is reported, not gated.
        hard = v.startswith("Sa_")
        if ds_mean is not None and str(ds_mean.attrs.get("fme_avg_mode", "")) != "mean":
            ds_mean = None
        if ds_mean is not None and v in ds_mean.data_vars:
            pm = _match_times(ds_mean["time"].values, ds6["time"].values)
            rms_m, _ = _rms_diff(np.asarray(ds_mean[v].values), a6, pm)
            if rms_m:
                ratio = rms_i / rms_m
                if ratio <= ratio_max:
                    out.append(("OK", msg + f" vs rms(mean-6h) {rms_m:.4g} "
                                            f"(ratio {ratio:.2f}, max {ratio_max}); "
                                            f"the residual is the ocean-map vs "
                                            f"atm-map remap, not round-off"))
                else:
                    out.append(("FAIL" if hard else "WARN",
                                msg + f" vs rms(mean-6h) {rms_m:.4g} "
                                f"(ratio {ratio:.2f} > {ratio_max})"
                                + (": the inst stream is not the window-end "
                                   "instant" if hard else
                                   ": slow ocean state, remap residual "
                                   "dominates (reported, not gated)")))
                continue
        out.append(("OK", msg + " (no window-mean parent to compare; "
                                "differs by remap)"))
    return out


def check_inst_vs_mixed_parent(ds_inst, ds_parent):
    """Mixed layout (cpl_fme_state_mean off): the parent's S* fields are the
    same window-end instants as the companion's -> must be identical."""
    if str(ds_parent.attrs.get("fme_avg_mode", "")) == "mean":
        return []
    pairs = _match_times(ds_inst["time"].values, ds_parent["time"].values)
    worst, worst_v = 0.0, ""
    for v in _time_vars(ds_inst):
        if v not in ds_parent.data_vars:
            continue
        ai = np.asarray(ds_inst[v].values, dtype=np.float64)
        ap_ = np.asarray(ds_parent[v].values, dtype=np.float64)
        for i, j in pairs:
            g = _valid(ai[i]) & _valid(ap_[j])
            if g.any():
                d = float(np.max(np.abs(ai[i][g] - ap_[j][g])))
                if d > worst:
                    worst, worst_v = d, v
    if not pairs:
        return [("WARN", "inst vs mixed parent: no shared records")]
    return [("OK" if worst == 0.0 else "FAIL",
             f"inst == parent's window-end S* (mixed layout): max abs diff "
             f"{worst:.3e}{' (' + worst_v + ')' if worst_v else ''} (expect 0)")]


def month_of(path):
    m = re.search(r"\.(\d{4}-\d{2})\.nc$", path)
    return m.group(1) if m else None


def open_month(rundir, stream, ym):
    if ym is None:
        return None
    files = [f for f in stream_files(rundir, stream, None) if month_of(f) == ym]
    if not files:
        return None
    return xr.open_dataset(files[0], decode_times=False, mask_and_scale=False)



def main():
    args = parse_args()
    rundir = args.rundir
    results = []          # (stream, status, msg)
    n_fail = 0
    missing_primary = []

    # which bundles look "on" (any file for either cadence)?
    bundle_has_files = {}
    for stream, bundle, cad in STREAMS:
        if stream_files(rundir, stream, args.year):
            bundle_has_files[bundle] = True

    open_ds = {}  # stream -> first dataset (for the cross-check)
    for stream, bundle, cad in STREAMS:
        files = stream_files(rundir, stream, args.year)
        if not files:
            if bundle in bundle_has_files and cad == "1D":
                results.append((stream, "FAIL", "PRIMARY stream has no files"))
                missing_primary.append(stream)
                n_fail += 1
            else:
                results.append((stream, "OK", "no files (stream disabled) -- skipped"))
            continue
        results.append((stream, "OK", f"{len(files)} file(s)"))
        try:
            # raw values: check_fills needs to see the 1e20 fill, not NaN
            ds = xr.open_dataset(files[-1], decode_times=False,
                                 mask_and_scale=False)
        except Exception as e:                       # noqa: BLE001
            results.append((stream, "FAIL", f"cannot open {os.path.basename(files[-1])}: {e}"))
            n_fail += 1
            continue
        open_ds[stream] = ds
        for grp in (check_cf_time(ds), check_grid(ds),
                    check_fields(ds, bundle), check_fills(ds, args.pollution_threshold),
                    check_time_layout(ds, bundle, cad, args.require_state_mean)):
            for status, msg in grp:
                results.append((stream, status, msg))
                if status == "FAIL":
                    n_fail += 1

    # 5D-vs-1D cross-checks per bundle
    for bundle in ("x2o", "xao", "x2i", "o2x", "i2x"):
        s1, s5 = f"{bundle}1D", f"{bundle}5D"
        if s1 in open_ds and s5 in open_ds:
            for status, msg in cross_check_5d_1d(open_ds[s5], open_ds[s1],
                                                 args.cross_check_rel_tol):
                results.append((s5, status, msg))
                if status == "FAIL":
                    n_fail += 1

    # *Dinst companions: presence, alignment, 5D==1D, spec check 5
    def note(stream, grp):
        nonlocal n_fail
        for status, msg in grp:
            results.append((stream, status, msg))
            if status == "FAIL":
                n_fail += 1

    for bundle in OCN_MESH_BUNDLES:
        for cad in ("1D", "5D"):
            parent, inst = f"{bundle}{cad}", f"{bundle}{cad}inst"
            if parent not in open_ds:
                continue
            dsp = open_ds[parent]
            named = str(dsp.attrs.get("fme_inst_companion", ""))
            if inst not in open_ds:
                if named or args.require_inst:
                    note(inst, [("FAIL", f"no files, but {parent} "
                                         f"{'names it' if named else 'is enabled'}")])
                elif str(dsp.attrs.get("fme_avg_mode", "")) == "mean":
                    note(inst, [("WARN", f"{parent} writes window-mean states "
                                         f"but no window-end instants "
                                         f"(cpl_fme_inst_enable off)")])
                continue
            dsi = open_ds[inst]
            ym = month_of(stream_files(rundir, inst, args.year)[-1])
            dsp_m = open_month(rundir, parent, ym)
            if dsp_m is None:
                note(inst, [("WARN", f"no {parent} file for {ym}; alignment "
                                     f"not checked")])
            else:
                note(inst, check_inst_alignment(dsi, dsp_m))
                note(inst, check_inst_vs_mixed_parent(dsi, dsp_m))
            if cad == "5D" and f"{bundle}1Dinst" in open_ds:
                ds1i = open_month(rundir, f"{bundle}1Dinst", ym)
                if ds1i is not None:
                    note(inst, check_inst_5d_vs_1d(dsi, ds1i))
            if cad == "5D":
                six = {s6: open_month(rundir, s6, ym) for s6 in ("a2x", "x2a")}
                six = {k: v for k, v in six.items() if v is not None}
                if six:
                    note(inst, check_inst_vs_6h(bundle, dsi, dsp_m, six,
                                                args.inst_ratio_max))
                else:
                    note(inst, [("WARN", f"spec check 5: no a2x/x2a file for "
                                         f"{ym}")])

    # report
    os.makedirs(os.path.join(args.outdir, "verify_cpl"), exist_ok=True)
    html = os.path.join(args.outdir, "verify_cpl", "index.html")
    write_html(html, args.label or os.path.basename(rundir.rstrip("/")),
               results, n_fail)

    print(f"\n=== coupler-native FME readiness: "
          f"{'PASS' if n_fail == 0 else 'FAIL (%d)' % n_fail} ===")
    for stream, status, msg in results:
        if status != "OK":
            print(f"  [{status}] {stream}: {msg}")
    print(f"report: {html}")

    if args.share_dir and args.label:
        publish(html, args.share_dir, args.label)

    sys.exit(1 if (n_fail or missing_primary) else 0)


def write_html(path, label, results, n_fail):
    color = {"OK": "#2e7d32", "WARN": "#f9a825", "FAIL": "#c62828"}
    rows = "\n".join(
        f"<tr><td>{s}</td><td style='color:{color.get(st,'#000')}'><b>{st}</b></td>"
        f"<td>{m}</td></tr>"
        for s, st, m in results)
    status = "PASS" if n_fail == 0 else f"FAIL ({n_fail})"
    with open(path, "w") as f:
        f.write(f"""<!doctype html><html><head><meta charset=utf-8>
<title>cpl-FME readiness: {label}</title>
<style>body{{font-family:sans-serif;margin:2em}}table{{border-collapse:collapse}}
td,th{{border:1px solid #ccc;padding:4px 8px}}</style></head><body>
<h1>Coupler-native FME readiness: {label}</h1>
<h2>Overall: {status}</h2>
<table><tr><th>stream</th><th>status</th><th>detail</th></tr>
{rows}</table></body></html>""")


def publish(html, share_dir, label):
    import shutil
    dest = os.path.join(share_dir, label)
    os.makedirs(dest, exist_ok=True)
    out = os.path.join(dest, "verify_cpl.html")
    shutil.copy(html, out)
    # additively link from the top-level index (never touch existing entries)
    idx = os.path.join(share_dir, "index.html")
    link = f'<li><a href="{label}/verify_cpl.html">{label} (cpl-FME)</a></li>'
    if os.path.exists(idx):
        txt = open(idx).read()
        if link not in txt and "</ul>" in txt:
            txt = txt.replace("</ul>", link + "\n</ul>", 1)
            open(idx, "w").write(txt)
    else:
        open(idx, "w").write(f"<html><body><ul>\n{link}\n</ul></body></html>")
    print(f"published: {out}")


if __name__ == "__main__":
    main()
