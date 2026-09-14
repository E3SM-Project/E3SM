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
}

# stream label -> (bundle, cadence). PRIMARY (1D) streams must exist for a
# bundle that is expected on; 5D companions are optional.
STREAMS = [
    ("x2o1D", "x2o", "1D"), ("x2o5D", "x2o", "5D"),
    ("xao1D", "xao", "1D"), ("xao5D", "xao", "5D"),
    ("x2i1D", "x2i", "1D"), ("x2i5D", "x2i", "5D"),
    ("o2x1D", "o2x", "1D"), ("o2x5D", "o2x", "5D"),
    ("i2x1D", "i2x", "1D"), ("i2x5D", "i2x", "5D"),
    # atm streams: single configurable sub-daily cadence (no 5D companion)
    ("x2a", "x2a", "1D"), ("a2x", "a2x", "1D"),
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
    ap.add_argument("--label", default=None, help="case label for the share dir")
    ap.add_argument("--share-dir", default=None,
                    help="copy the report under <share-dir>/<label>/ and link it")
    return ap.parse_args()


def stream_files(rundir, stream, year):
    pat = os.path.join(rundir, f"*.cpl.fme.{stream}.*.nc")
    files = sorted(glob.glob(pat))
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
            if len(idx) < 2:
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
            ds = xr.open_dataset(files[-1], decode_times=False)
        except Exception as e:                       # noqa: BLE001
            results.append((stream, "FAIL", f"cannot open {os.path.basename(files[-1])}: {e}"))
            n_fail += 1
            continue
        open_ds[stream] = ds
        for grp in (check_cf_time(ds), check_grid(ds),
                    check_fields(ds, bundle), check_fills(ds, args.pollution_threshold)):
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
