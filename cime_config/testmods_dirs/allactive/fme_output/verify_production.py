#!/usr/bin/env python3
"""verify_production.py -- production-tape validator for FME output.

Per spatial variable in EAM h0 + the three MPAS .remapped.nc streams,
generates two figures:
  * area-weighted global-mean timeseries (line plot over time)
  * temporal-mean climatology (lat-lon map)

Plus a sanity-check report (CF conventions, time coord, dim names,
expected globals like time_reference_date) and an HTML viewer at
<outdir>/verify_production/index.html with a thumbnail grid linking
to full-resolution PNGs.

The current FME testmod produces ~290 spatial variables across the
four streams. By default we ingest every month present in the run
directory via xr.open_mfdataset, which gives proper multi-month
timeseries and a meaningful climatology mean. Use --year YYYY to
restrict to one year for a faster pass.

Targets the post-2026-05-04 CF time scheme strictly: dim is `time`
(lowercase), `time_bnds` and `date`/`datesec` present, no xtime.
Older files (legacy `Time`/xtime) are detected and skipped with a
note in the sanity report.

Usage:
  python verify_production.py --rundir <CASE_RUNDIR> --outdir <ROOT>
  python verify_production.py --rundir <CASE_RUNDIR> --outdir <ROOT> --year 0411
  python verify_production.py --rundir <CASE_RUNDIR> --outdir <ROOT> --include-inst
"""

import argparse
import glob
import html
import os
import re
import sys
import traceback
from pathlib import Path

import numpy as np
import xarray as xr
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Stream definitions
# ---------------------------------------------------------------------------

STREAMS = [
    {
        "key": "eam_h0",
        "label": "EAM h0",
        "pattern": "*eam.h0.*.nc",
        "exclude": [".rh0."],
    },
    {
        "key": "mpaso_depth",
        "label": "MPAS-O fmeDepthCoarsening",
        "pattern": "*mpaso.hist.am.fmeDepthCoarsening.*.remapped.nc",
        # Exclude the 5D companion files which also share this prefix
        # in some glob orderings depending on locale.
        "exclude": ["fmeDepthCoarsening5D"],
    },
    {
        "key": "mpaso_depth_5D",
        "label": "MPAS-O fmeDepthCoarsening (5-day mean companion)",
        "pattern": "*mpaso.hist.am.fmeDepthCoarsening5D.*.remapped.nc",
        "exclude": [],
        "companion_of": "mpaso_depth",
    },
    {
        "key": "mpaso_derived",
        "label": "MPAS-O fmeDerivedFields",
        "pattern": "*mpaso.hist.am.fmeDerivedFields.*.remapped.nc",
        "exclude": ["fmeDerivedFields5D"],
    },
    {
        "key": "mpaso_derived_5D",
        "label": "MPAS-O fmeDerivedFields (5-day mean companion)",
        "pattern": "*mpaso.hist.am.fmeDerivedFields5D.*.remapped.nc",
        "exclude": [],
        "companion_of": "mpaso_derived",
    },
    {
        "key": "mpassi_derived",
        "label": "MPAS-SI fmeSeaiceDerivedFields",
        "pattern": "*mpassi.hist.am.fmeSeaiceDerivedFields.*.remapped.nc",
        "exclude": ["fmeSeaiceDerivedFields5D"],
    },
    {
        "key": "mpassi_derived_5D",
        "label": "MPAS-SI fmeSeaiceDerivedFields (5-day mean companion)",
        "pattern": "*mpassi.hist.am.fmeSeaiceDerivedFields5D.*.remapped.nc",
        "exclude": [],
        "companion_of": "mpassi_derived",
    },
]


YYYYMM_RE = re.compile(r"\.(\d{4})-(\d{2})\.")


# ---------------------------------------------------------------------------
# CLI + file selection
# ---------------------------------------------------------------------------

def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rundir", required=True,
                    help="Case run directory containing eam.h0.*.nc and *.remapped.nc files")
    ap.add_argument("--outdir", required=True,
                    help="Root directory under which verify_production/ is created")
    ap.add_argument("--year", default=None,
                    help="YYYY year to restrict to (default: all months present)")
    ap.add_argument("--include-inst", action="store_true",
                    help="Also plot *_inst instantaneous companion vars (default: skip)")
    ap.add_argument("--pollution-threshold", type=float, default=1e15,
                    help="|x| above this magnitude (but != 1e20 canonical fill) is "
                         "flagged as fill-pollution (gotcha #36). Default 1e15 -- "
                         "well above any geophysical value (max ~6e11 J/m^2 for OHC) "
                         "and well below SHR_FILL_VALUE=1e20. Catches subtler leaks "
                         "than the original 1e18 detection threshold.")
    ap.add_argument("--dpi", type=int, default=85,
                    help="Figure DPI (default 85)")
    ap.add_argument("--label", default=None,
                    help="Display label for the case (default: case name from rundir)")
    ap.add_argument("--share-dir", default=None,
                    help="If set, copy the rendered verify_production tree to "
                         "<share-dir>/<case_label>/ (e.g. $PSCRATCH/share). "
                         "Useful on NERSC where $PSCRATCH/share is web-viewable "
                         "via the science portal.")
    ap.add_argument("--cross-check-rel-tol", type=float, default=1e-2,
                    help="Relative-tolerance threshold for 5D vs mean(1D) "
                         "cross-check (default 1e-2 = 1%% of variable's "
                         "p99-scale at 99%% of cells). Differences above this "
                         "are flagged WARN; >5%% are flagged FAIL. "
                         "5D and 1D are sample-weighted vs day-weighted "
                         "respectively; count-flicker cells (ice edge etc.) "
                         "legitimately diverge so the metric uses p99(diff) "
                         "across cells, not max(diff).")
    ap.add_argument("--cross-check-abs-tol", type=float, default=1e-9,
                    help="Absolute-tolerance floor used when normalizing the "
                         "relative-diff denominator (default 1e-9).")
    return ap.parse_args()


def stream_files(rundir, stream, year_filter):
    files = sorted(glob.glob(os.path.join(rundir, stream["pattern"])))
    files = [f for f in files if not any(ex in f for ex in stream["exclude"])]
    if year_filter:
        files = [f for f in files if f".{year_filter}-" in os.path.basename(f)]
    return files


def filter_new_format(files):
    """Keep only files that have lowercase `time` dim and `time_bnds` var
    (the post-2026-05-04 CF rewrite). Files that pre-date the rewrite
    (`Time`/xtime) are dropped with a note for the sanity report.
    Returns (kept, skipped_legacy)."""
    kept, legacy = [], []
    for f in files:
        try:
            with xr.open_dataset(f, decode_times=False) as ds:
                if "time" in ds.dims and "time_bnds" in ds.variables:
                    kept.append(f)
                else:
                    legacy.append(f)
        except Exception:
            legacy.append(f)
    return kept, legacy


# ---------------------------------------------------------------------------
# Sanity checks
# ---------------------------------------------------------------------------

SHR_FILL_VALUE = 1e20


def check_fills(ds, pollution_threshold):
    """Scan every spatial variable for fill-value leaks. Returns a list of
    (level, msg) tuples. Categories per cell:
      * `|x| < threshold`     -- physical
      * `threshold <= |x| < 1e20` -- POLLUTION (fractional fill, gotcha #36)
      * `|x| == 1e20`         -- canonical SHR_FILL_VALUE (expected on land)
    Also flags NaN/inf (Fortran should never write these -- it writes 1e20).

    The dataset must have been opened with mask_and_scale=False so xarray
    didn't auto-convert _FillValue cells to NaN.
    """
    out = []
    spatial = sorted(v for v in ds.data_vars if {"lat", "lon"} <= set(ds[v].dims))
    if not spatial:
        out.append(("OK", "fill-leak scan: no spatial vars to check"))
        return out

    polluted_vars = []
    nan_vars = []
    bad_mask_vars = []
    fill_counts = []  # (var, n_canonical_fill, n_total) for informational summary

    for v in spatial:
        arr = ds[v].values
        # Use a flat copy for cheap scans
        n_nan = int(np.isnan(arr).sum())
        n_inf = int(np.isinf(arr).sum())
        if n_nan or n_inf:
            nan_vars.append((v, n_nan, n_inf))
            continue
        absv = np.abs(arr)
        n_polluted = int(((absv >= pollution_threshold) & (absv < SHR_FILL_VALUE)).sum())
        if n_polluted > 0:
            samples = arr[(absv >= pollution_threshold) & (absv < SHR_FILL_VALUE)].ravel()
            polluted_vars.append((v, n_polluted, samples[:5]))
        n_canonical = int((absv == SHR_FILL_VALUE).sum())
        fill_counts.append((v, n_canonical, arr.size))
        # Mask vars must be exactly 0 or 1
        if v.startswith("mask_"):
            uniq = np.unique(arr)
            bad = uniq[~np.isin(uniq, [0.0, 1.0])]
            if bad.size > 0:
                bad_mask_vars.append((v, bad[:5].tolist()))

    if polluted_vars:
        for v, n, samples in polluted_vars:
            sample_str = ", ".join(f"{x:.4g}" for x in samples)
            out.append(("FAIL", f"fill leak in '{v}': {n} cell(s) with "
                                f"|x| in [{pollution_threshold:.0e}, 1e20). "
                                f"e.g. {sample_str}"))
    if nan_vars:
        for v, nn, ni in nan_vars:
            out.append(("FAIL", f"non-finite in '{v}': {nn} NaN, {ni} Inf "
                                f"(Fortran should write only 1e20, never NaN/Inf)"))
    if bad_mask_vars:
        for v, bad in bad_mask_vars:
            out.append(("FAIL", f"mask var '{v}' has non-{{0,1}} values: {bad}"))

    if not (polluted_vars or nan_vars or bad_mask_vars):
        # Summarize the canonical-fill landscape so the report is informative
        # even when everything is clean.
        ranges = [c[1] for c in fill_counts]
        sz = fill_counts[0][2] if fill_counts else 0
        out.append(("OK", f"fill-leak scan: {len(spatial)} var(s) clean. "
                          f"canonical-fill cells per var: "
                          f"{min(ranges)}..{max(ranges)} of {sz} (expected: land/bathymetry)"))
    return out


def _bnds_to_2col(arr):
    """time_bnds is sometimes (time, nbnd), sometimes (nbnd, time) depending
    on writer / xarray decode order. Normalize to (n, 2)."""
    if arr.ndim != 2:
        return None
    if arr.shape[1] == 2:
        return arr
    if arr.shape[0] == 2:
        return arr.T
    return None


def cross_check_5d_vs_1d(ds_5d, ds_1d, rel_tol=1e-4, abs_tol=1e-9):
    """For each (variable, 5D record), find the 1D records whose time_bnds are
    fully inside the 5D record's bnds and verify
        ds_5d[v][j5] ~= mean(ds_1d[v][matching])
    Both AMs accumulate the same per-step samples, so the 5D mean is
    exactly mean(per-step) over a 5-day window, and the 1D-mean over the
    same window collapses to the same expression *iff* both windows are
    perfectly aligned and contain a whole number of 1D windows.

    Returns a list of (level, msg) tuples for the sanity report.
    """
    out = []
    if "time_bnds" not in ds_5d.variables or "time_bnds" not in ds_1d.variables:
        out.append(("WARN", "5D-vs-1D cross-check: missing time_bnds on 5D or 1D dataset"))
        return out

    bnds_5d = _bnds_to_2col(ds_5d["time_bnds"].values)
    bnds_1d = _bnds_to_2col(ds_1d["time_bnds"].values)
    if bnds_5d is None or bnds_1d is None:
        out.append(("WARN", "5D-vs-1D cross-check: time_bnds shape unrecognized"))
        return out

    # Restrict to vars present in both with (time, lat, lon) signature.
    common = sorted(
        v for v in set(ds_5d.data_vars) & set(ds_1d.data_vars)
        if "time" in ds_5d[v].dims
        and {"lat", "lon"}.issubset(ds_5d[v].dims)
        and "time" in ds_1d[v].dims
        and {"lat", "lon"}.issubset(ds_1d[v].dims)
    )
    if not common:
        out.append(("OK", "5D-vs-1D cross-check: no common time-varying spatial vars"))
        return out

    # Per-window membership: 1D window i is "inside" 5D window j iff
    # bnds_1d[i,0] >= bnds_5d[j,0] and bnds_1d[i,1] <= bnds_5d[j,1] (with eps).
    eps = 1e-6
    n5 = bnds_5d.shape[0]
    membership = []  # list of (j5, idx_array_into_1d) for each 5D window with >=2 members
    for j5 in range(n5):
        t_lo, t_hi = bnds_5d[j5, 0], bnds_5d[j5, 1]
        idx = np.where(
            (bnds_1d[:, 0] >= t_lo - eps)
            & (bnds_1d[:, 1] <= t_hi + eps)
        )[0]
        if idx.size >= 2:
            membership.append((j5, idx))

    if not membership:
        out.append(("WARN", f"5D-vs-1D cross-check: no 5D window contains >=2 "
                            f"matching 1D windows ({n5} 5D records, "
                            f"{bnds_1d.shape[0]} 1D records); skipping"))
        return out

    out.append(("OK", f"5D-vs-1D cross-check: matched {len(membership)} of {n5} "
                      f"5D windows to >=2 1D records each; comparing "
                      f"{len(common)} variable(s)"))

    fail_vars = []
    warn_vars = []
    rel_summary = []
    for v in common:
        v5_arr = ds_5d[v].values
        v1_arr = ds_1d[v].values
        # Per-variable scale: p99 of |valid samples| across all 1D records.
        # This is the robust normalizer for "how much of the variable's
        # natural range does the 5D-vs-1D disagreement amount to". Using
        # per-cell |v| as the normalizer (the obvious choice) blows up at
        # zero-crossings (e.g., the 0 °C isotherm in SST, or null-velocity
        # cells in U/V) where the cell value ~ 0 but the disagreement is
        # still tiny in absolute terms.
        v1_flat = v1_arr.ravel()
        ok = np.isfinite(v1_flat) & (np.abs(v1_flat) < SHR_FILL_VALUE * 0.1)
        scale = float(np.percentile(np.abs(v1_flat[ok]), 99)) if ok.any() else 0.0
        scale = max(scale, abs_tol)

        # IMPORTANT: 5D[j] = sum(per-step samples)/count(per-step samples)
        # is *sample-weighted*. Python's mean(1D[i_in_window]) is
        # *day-weighted*. They coincide only when each 1D record has the
        # same per-cell sample count, which holds for cells whose
        # validity is stable across the window. At fill-edge cells (ice
        # edge, land/ocean transitions, river-runoff cells) the per-cell
        # count flips day-to-day and the two means legitimately disagree
        # by O(spike-magnitude * count_imbalance). The model output
        # doesn't carry per-cell counts so we can't reconstruct the
        # sample-weighted 1D mean exactly.
        #
        # We work around it with TWO complementary metrics:
        #   - p99-diff: 99th-percentile abs-diff across all valid cells.
        #     Robust to a handful of count-flicker cells; sub-percent
        #     means "everywhere except ice/land/runoff edges".
        #   - max-diff: shown for context but not used for FAIL.
        max_p99_rel = 0.0
        max_p99_abs = 0.0
        max_max_abs = 0.0
        worst_j5 = -1
        worst_n_match = 0
        for (j5, idx) in membership:
            v5 = v5_arr[j5, ...]
            v1_window = v1_arr[idx, ...]
            v1_mean = v1_window.mean(axis=0)
            v1_window_valid = (
                np.isfinite(v1_window)
                & (np.abs(v1_window) < SHR_FILL_VALUE * 0.1)
            )
            all_1d_valid = v1_window_valid.all(axis=0)
            valid = (
                all_1d_valid
                & np.isfinite(v5) & (np.abs(v5) < SHR_FILL_VALUE * 0.1)
            )
            if not valid.any():
                continue
            diff = np.abs(v5[valid] - v1_mean[valid])
            p99 = float(np.percentile(diff, 99))
            mx = float(diff.max())
            r = p99 / scale
            if r > max_p99_rel:
                max_p99_rel = r
                max_p99_abs = p99
                max_max_abs = mx
                worst_j5 = j5
                worst_n_match = idx.size
        # Names kept for downstream packing
        max_rel, max_abs = max_p99_rel, max_p99_abs
        rel_summary.append((v, max_rel, max_abs, worst_j5, worst_n_match, scale))
        if max_rel > 0.01:           # >1% of variable's p99 — broken
            fail_vars.append((v, max_rel, max_abs, worst_j5, worst_n_match, scale))
        elif max_rel > rel_tol:      # rel_tol .. 1% — suspicious
            warn_vars.append((v, max_rel, max_abs, worst_j5, worst_n_match, scale))

    if fail_vars:
        for v, mr, ma, j, n, sc in fail_vars:
            out.append(("FAIL", f"5D vs mean(1D): {v}: p99(diff)={ma:.3e} "
                                f"({mr*100:.2f}% of p99-scale {sc:.3e}) at "
                                f"5D-record={j} (matched {n} 1D records)"))
    if warn_vars:
        for v, mr, ma, j, n, sc in warn_vars:
            out.append(("WARN", f"5D vs mean(1D): {v}: p99(diff)={ma:.3e} "
                                f"({mr*100:.2f}% of p99-scale {sc:.3e}) at "
                                f"5D-record={j} (matched {n} 1D records)"))
    if not (fail_vars or warn_vars):
        if rel_summary:
            v, mr, ma, j, n, sc = max(rel_summary, key=lambda r: r[1])
            out.append(("OK", f"5D vs mean(1D): all {len(common)} var(s) "
                              f"within rel_tol={rel_tol:.0e} of p99-scale "
                              f"(99% of cells); worst: {v} p99(diff)={ma:.3e} "
                              f"({mr*100:.3f}% of p99-scale {sc:.3e})"))
    return out


def sanity_check(ds, stream_label, n_files, n_legacy_skipped, is_mpas,
                 pollution_threshold):
    """Return list of (level, msg). level in {OK, WARN, FAIL}."""
    out = []

    def add(lvl, msg):
        out.append((lvl, msg))

    add("OK", f"merged {n_files} file(s)")
    if n_legacy_skipped:
        add("WARN", f"skipped {n_legacy_skipped} legacy-format file(s) "
                    f"(no `time` dim or `time_bnds` -- pre-CF-rewrite)")

    for d in ("lat", "lon", "time"):
        if d not in ds.dims:
            add("FAIL", f"missing dimension '{d}'")

    if "time" in ds.variables:
        units = ds["time"].attrs.get("units", "")
        cal = ds["time"].attrs.get("calendar", "")
        bounds = ds["time"].attrs.get("bounds", "")
        if not units.startswith("days since "):
            add("WARN", f"time:units='{units}' does not start with 'days since '")
        if cal != "noleap":
            add("WARN", f"time:calendar={cal!r} (expected 'noleap')")
        if bounds != "time_bnds":
            add("WARN", f"time:bounds={bounds!r} (expected 'time_bnds')")
        # Monotonicity
        tvals = ds["time"].values
        if tvals.ndim == 1 and tvals.size > 1:
            if not np.all(np.diff(tvals) > 0):
                add("FAIL", "time axis not strictly monotonically increasing")
            add("OK", f"time: {tvals.size} records spanning [{tvals[0]:.2f}, {tvals[-1]:.2f}] {units.split(' ', 2)[-1]}")
        for v in ("time_bnds", "date", "datesec"):
            if v not in ds.variables:
                add("WARN", f"missing CF time companion var '{v}'")

    if "xtime" in ds.variables:
        add("WARN", "xtime present (legacy MPAS char timestamp; should be gone after CF rewrite)")

    if "Conventions" not in ds.attrs:
        add("WARN", "missing global :Conventions attribute")
    elif not str(ds.attrs["Conventions"]).startswith("CF"):
        add("WARN", f"Conventions={ds.attrs['Conventions']!r} (expected to start with 'CF')")

    if is_mpas:
        if "time_reference_date" not in ds.attrs:
            add("WARN", "MPAS file lacks global :time_reference_date attribute "
                        "(should be set by *_write_accum sidecar persistence)")
        else:
            add("OK", f"time_reference_date = {ds.attrs['time_reference_date']!r}")

    if "lat" in ds.variables:
        latv = ds["lat"].values
        if latv.min() < -90.001 or latv.max() > 90.001:
            add("FAIL", f"lat range out of bounds: [{latv.min()}, {latv.max()}]")
    if "lon" in ds.variables:
        lonv = ds["lon"].values
        if lonv.min() < -0.001 or lonv.max() > 360.001:
            add("WARN", f"lon range unusual: [{lonv.min()}, {lonv.max()}]")

    spatial = [v for v in ds.data_vars if {"lat", "lon"} <= set(ds[v].dims)]
    add("OK", f"{len(spatial)} variable(s) with lat+lon dims")

    # Fill-leak scan -- the headline correctness check.
    out.extend(check_fills(ds, pollution_threshold))
    return out


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def cos_lat_weights(lat):
    return np.clip(np.cos(np.deg2rad(lat)), 0.0, None)


def is_plottable(da):
    """Plot only purely-2D-spatial variables: lat & lon dims, plus optionally
    `time`. Skip 3D-in-vertical (no production var hits this after vcoarsen)."""
    dims = set(da.dims)
    if not {"lat", "lon"}.issubset(dims):
        return False
    extras = dims - {"lat", "lon", "time", "nbnd"}
    return len(extras) == 0


def is_static(da):
    return "time" not in da.dims


def _mask_invalid(da, pollution_threshold=1e15):
    """Return a copy with non-finite and fill-magnitude values masked to NaN.
    Threshold matches the sanity-check pollution boundary so what we plot
    aligns with what we flag."""
    arr = da.where(np.isfinite(da)).where(np.abs(da) < pollution_threshold)
    return arr


def plot_timeseries(da, ds, title, outpath, dpi, pollution_threshold=1e15):
    if is_static(da):
        return False
    arr = _mask_invalid(da, pollution_threshold)
    lat = ds["lat"].values
    w = xr.DataArray(cos_lat_weights(lat), coords={"lat": lat}, dims=["lat"])
    gmean = arr.weighted(w).mean(("lat", "lon")).values

    tvals = ds["time"].values
    if tvals.ndim != 1 or tvals.size != gmean.size:
        tvals = np.arange(gmean.size)
        tlabel = "record index"
    else:
        tlabel = ds["time"].attrs.get("units", "time")

    fig, ax = plt.subplots(figsize=(7, 3))
    ax.plot(tvals, gmean, "-", lw=1)
    if gmean.size <= 200:
        ax.plot(tvals, gmean, "o", ms=2.5)
    ax.set_xlabel(tlabel)
    units = da.attrs.get("units", "")
    ax.set_ylabel(f"{da.name} [{units}]" if units else da.name)
    ax.set_title(f"{title}\nglobal mean (cos-lat weighted)")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(outpath, dpi=dpi)
    plt.close(fig)
    return True


def plot_climatology(da, ds, title, outpath, dpi, pollution_threshold=1e15):
    arr = _mask_invalid(da, pollution_threshold)
    if not is_static(da):
        clim = arr.mean(dim="time")
    else:
        clim = arr
    clim_vals = clim.values

    lat = ds["lat"].values
    lon = ds["lon"].values

    fig, ax = plt.subplots(figsize=(7.5, 3.6))
    flat = clim_vals[np.isfinite(clim_vals)]
    if flat.size == 0:
        ax.text(0.5, 0.5, "no valid data", ha="center", va="center", transform=ax.transAxes)
    else:
        vmin, vmax = np.percentile(flat, [2.0, 98.0])
        if vmin == vmax:
            vmin, vmax = float(flat.min()), float(flat.max() + 1e-12)
        if vmin < 0 < vmax and abs(vmin + vmax) < 0.5 * (vmax - vmin):
            v = max(abs(vmin), abs(vmax))
            vmin, vmax = -v, v
            cmap = "RdBu_r"
        else:
            cmap = "viridis"
        im = ax.imshow(
            clim_vals,
            origin="lower",
            extent=[lon.min(), lon.max(), lat.min(), lat.max()],
            aspect="auto",
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            interpolation="nearest",
        )
        cbar = fig.colorbar(im, ax=ax, shrink=0.85)
        units = da.attrs.get("units", "")
        cbar.set_label(units)
    ax.set_xlabel("longitude")
    ax.set_ylabel("latitude")
    suffix = "static" if is_static(da) else f"time-mean over {ds.sizes.get('time', '?')} records"
    ax.set_title(f"{title}\n{suffix}")
    fig.tight_layout()
    fig.savefig(outpath, dpi=dpi)
    plt.close(fig)
    return True


# ---------------------------------------------------------------------------
# HTML viewer
# ---------------------------------------------------------------------------

CSS = """
* { box-sizing: border-box; }
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
       max-width: 1500px; margin: 1em auto; padding: 0 1em; color: #222; }
h1 { border-bottom: 2px solid #333; padding-bottom: 0.3em; }
h2 { border-bottom: 1px solid #999; padding-bottom: 0.2em; margin-top: 2em;
     display: flex; align-items: baseline; gap: 0.6em; flex-wrap: wrap; }
h2 .browse-btn { font-size: 0.65em; font-weight: normal; padding: 4px 10px;
                 background: #2c3e6b; color: white; border-radius: 4px;
                 cursor: pointer; border: none; }
h2 .browse-btn:hover { background: #1a2744; }
.meta { background: #f4f4f4; padding: 0.6em 1em; border-radius: 4px; font-size: 0.9em; }
.sanity { background: #fafafa; border-left: 4px solid #999; padding: 0.4em 1em; margin-bottom: 1em;
          font-family: monospace; font-size: 0.85em; }
.sanity .ok    { color: #277; }
.sanity .warn  { color: #b85; }
.sanity .fail  { color: #b22; font-weight: bold; }
.grid { display: grid; grid-template-columns: 1fr 1fr; gap: 1.2em; }
.var-card { border: 1px solid #ddd; border-radius: 4px; padding: 0.4em;
            cursor: pointer; transition: border-color 0.15s, box-shadow 0.15s; }
.var-card:hover { border-color: #2c3e6b; box-shadow: 0 2px 8px rgba(0,0,0,0.08); }
.var-name { font-family: monospace; font-size: 0.95em; font-weight: bold; }
.var-units { color: #666; font-size: 0.85em; }
.thumbs { display: flex; gap: 0.4em; margin-top: 0.4em; }
.thumbs img { max-width: 320px; height: auto; border: 1px solid #ccc; }
nav.toc { background: #fafafa; padding: 0.6em; border-radius: 4px; }
nav.toc a { margin-right: 1em; }
.help { font-size: 0.85em; color: #666; margin-top: 0.5em; }
kbd { background: #eee; border: 1px solid #ccc; border-radius: 3px;
      padding: 1px 5px; font-family: monospace; font-size: 0.85em; }

/* ----- Modal viewer ----- */
#viewer { position: fixed; inset: 0; background: rgba(15, 20, 30, 0.95);
          z-index: 1000; display: none; flex-direction: column; padding: 14px; gap: 10px; }
#viewer.open { display: flex; }
#viewer .vh { color: #ddd; display: flex; align-items: center; gap: 1.4em;
              font-family: -apple-system, BlinkMacSystemFont, sans-serif; }
#viewer .vh .stream { color: #88c; font-size: 0.9em; }
#viewer .vh .name { font-family: monospace; font-size: 1.15em; color: #fff; font-weight: bold; }
#viewer .vh .units { color: #aaa; font-family: monospace; }
#viewer .vh .pos { color: #aaa; font-size: 0.9em; margin-left: auto; }
#viewer .vh button { background: rgba(255,255,255,0.1); color: white; border: 1px solid #555;
                     border-radius: 4px; padding: 6px 14px; cursor: pointer; font-size: 0.9em; }
#viewer .vh button:hover { background: rgba(255,255,255,0.2); }
#viewer .vbody { flex: 1; display: grid; grid-template-columns: 1fr 1fr; gap: 12px; min-height: 0; }
#viewer .vbody .pane { display: flex; flex-direction: column; gap: 6px;
                       background: rgba(0,0,0,0.3); border-radius: 4px; padding: 8px; }
#viewer .vbody .pane h3 { color: #aaa; font-size: 0.85em; font-weight: normal;
                          font-family: monospace; margin: 0; }
#viewer .vbody img { width: 100%; max-height: calc(100vh - 160px); object-fit: contain;
                     background: white; border-radius: 3px; }
#viewer .vbody .empty { color: #777; font-style: italic; padding: 20px;
                        text-align: center; flex: 1; display: flex;
                        align-items: center; justify-content: center; }
#viewer .vfoot { color: #aaa; font-size: 0.8em; text-align: center; }
#viewer .vfoot kbd { background: #222; color: #ddd; border-color: #444; }
"""


def render_sanity(checks):
    rows = []
    for lvl, msg in checks:
        rows.append(f"<div class='{lvl.lower()}'>[{lvl}] {html.escape(msg)}</div>")
    return "<div class='sanity'>" + "".join(rows) + "</div>"


def render_var_card(stream_key, varname, units, ts_path_rel, clim_path_rel):
    """Card includes data attrs so JS can open the viewer at this var."""
    units_s = f"<span class='var-units'>[{html.escape(units)}]</span>" if units else ""
    thumbs = []
    if ts_path_rel:
        thumbs.append(f"<img src='{html.escape(ts_path_rel)}' alt='timeseries' "
                      f"title='global timeseries'>")
    if clim_path_rel:
        thumbs.append(f"<img src='{html.escape(clim_path_rel)}' alt='climatology' "
                      f"title='time-mean'>")
    return (f"<div class='var-card' data-stream='{html.escape(stream_key)}' "
            f"data-var='{html.escape(varname)}' "
            f"onclick=\"openViewer('{html.escape(stream_key)}','{html.escape(varname)}')\">"
            f"<div><span class='var-name'>{html.escape(varname)}</span> {units_s}</div>"
            f"<div class='thumbs'>{''.join(thumbs)}</div>"
            "</div>")


def render_section(stream_key, label, anchor, files, sanity_checks, var_cards):
    parts = [
        f"<h2 id='{anchor}'>{html.escape(label)}"
        + (f" <button class='browse-btn' "
           f"onclick=\"openViewer('{html.escape(stream_key)}', null)\">Browse "
           f"&rarr;</button>" if var_cards else "")
        + "</h2>"
    ]
    if files:
        files_html = "<br>".join(f"<code>{html.escape(os.path.basename(p))}</code>" for p in files[:3])
        more = f" (+{len(files)-3} more)" if len(files) > 3 else ""
        parts.append(f"<div class='meta'>{len(files)} file(s):<br>{files_html}{more}</div>")
    parts.append(render_sanity(sanity_checks))
    if not var_cards:
        parts.append("<p><em>no plottable variables</em></p>")
    else:
        parts.append("<div class='grid'>")
        parts.extend(var_cards)
        parts.append("</div>")
    return "\n".join(parts)


VIEWER_HTML = """
<div id='viewer' onclick='if(event.target.id===\"viewer\")closeViewer()'>
  <div class='vh'>
    <span><span id='v-stream' class='stream'></span> /
          <span id='v-name' class='name'></span>
          <span id='v-units' class='units'></span></span>
    <span id='v-pos' class='pos'></span>
    <button onclick='prevVar()'>&larr; Prev</button>
    <button onclick='nextVar()'>Next &rarr;</button>
    <button onclick='closeViewer()'>Close (Esc)</button>
  </div>
  <div class='vbody'>
    <div class='pane'>
      <h3>global timeseries (cos-lat weighted)</h3>
      <img id='v-ts' alt='timeseries'>
    </div>
    <div class='pane'>
      <h3>time-mean climatology</h3>
      <img id='v-clim' alt='climatology'>
    </div>
  </div>
  <div class='vfoot'>
    <kbd>&larr;</kbd> / <kbd>h</kbd> prev &nbsp;|&nbsp;
    <kbd>&rarr;</kbd> / <kbd>l</kbd> next &nbsp;|&nbsp;
    <kbd>Esc</kbd> close &nbsp;|&nbsp;
    <kbd>Home</kbd>/<kbd>End</kbd> first/last
  </div>
</div>
"""


VIEWER_JS = r"""
let _curStream = null;
let _curIdx = 0;

function _setImg(id, src) {
  const el = document.getElementById(id);
  if (src) {
    el.src = src;
    el.style.display = '';
    el.parentElement.querySelector('.empty')?.remove();
  } else {
    el.style.display = 'none';
    if (!el.parentElement.querySelector('.empty')) {
      const div = document.createElement('div');
      div.className = 'empty';
      div.textContent = 'no figure';
      el.parentElement.appendChild(div);
    }
  }
}

function _show() {
  const list = MANIFEST[_curStream] || [];
  if (list.length === 0) return;
  if (_curIdx < 0) _curIdx = list.length - 1;
  if (_curIdx >= list.length) _curIdx = 0;
  const v = list[_curIdx];
  document.getElementById('v-stream').textContent = _curStream;
  document.getElementById('v-name').textContent = v.name;
  document.getElementById('v-units').textContent = v.units ? '[' + v.units + ']' : '';
  document.getElementById('v-pos').textContent = (_curIdx + 1) + ' / ' + list.length;
  _setImg('v-ts', v.ts);
  _setImg('v-clim', v.clim);
  history.replaceState(null, '', '#' + _curStream + '/' + v.name);
}

function openViewer(stream, varname) {
  _curStream = stream;
  const list = MANIFEST[stream] || [];
  if (varname) {
    _curIdx = list.findIndex(v => v.name === varname);
    if (_curIdx < 0) _curIdx = 0;
  } else {
    _curIdx = 0;
  }
  document.getElementById('viewer').classList.add('open');
  _show();
}
function closeViewer() {
  document.getElementById('viewer').classList.remove('open');
  history.replaceState(null, '', window.location.pathname + window.location.search);
}
function nextVar() { _curIdx++; _show(); }
function prevVar() { _curIdx--; _show(); }
function firstVar() { _curIdx = 0; _show(); }
function lastVar() { _curIdx = (MANIFEST[_curStream] || []).length - 1; _show(); }

document.addEventListener('keydown', e => {
  if (!document.getElementById('viewer').classList.contains('open')) return;
  switch (e.key) {
    case 'ArrowRight': case 'l': case ' ': nextVar(); e.preventDefault(); break;
    case 'ArrowLeft':  case 'h': prevVar(); e.preventDefault(); break;
    case 'Home': firstVar(); e.preventDefault(); break;
    case 'End':  lastVar();  e.preventDefault(); break;
    case 'Escape': closeViewer(); e.preventDefault(); break;
  }
});

// Honor URL fragment on load: #stream/varname
window.addEventListener('DOMContentLoaded', () => {
  const h = window.location.hash.slice(1);
  if (!h) return;
  const [s, v] = h.split('/');
  if (s && MANIFEST[s]) openViewer(s, v || null);
});
"""


def write_html(outdir, title, sections, toc, manifest):
    import json
    body = [
        f"<h1>{html.escape(title)}</h1>",
        "<nav class='toc'>" + " | ".join(toc) + "</nav>",
        "<p class='help'>Click any variable card to open the side-by-side viewer; "
        "use <kbd>&larr;</kbd>/<kbd>&rarr;</kbd> to navigate within a stream, "
        "<kbd>Esc</kbd> to close.</p>",
        *sections,
        VIEWER_HTML,
    ]
    manifest_js = "const MANIFEST = " + json.dumps(manifest, separators=(",", ":")) + ";\n"
    html_doc = (
        "<!doctype html><html><head><meta charset='utf-8'>"
        f"<title>{html.escape(title)}</title>"
        f"<style>{CSS}</style></head><body>"
        + "\n".join(body)
        + f"<script>{manifest_js}{VIEWER_JS}</script>"
        + "</body></html>"
    )
    out = outdir / "index.html"
    out.write_text(html_doc)
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def safe_var_filename(name):
    return re.sub(r"[^A-Za-z0-9_\-.]", "_", name)


def case_label_from_rundir(rundir):
    parts = rundir.rstrip("/").split("/")
    # Common layouts: .../<casename>/run or .../<casename>
    if parts and parts[-1] == "run" and len(parts) > 1:
        return parts[-2]
    return parts[-1] if parts else "case"


def push_to_share(out_root, share_dir, case_label):
    """Sync the rendered tree to <share_dir>/<case_label>/ and refresh a
    top-level <share_dir>/index.html that lists every case present in
    share_dir. The top-level index is rebuilt from the filesystem on
    every push, so concurrent pushes from multiple cases stay in sync.

    Uses staging-then-swap so a viewer reading <case_label>/ during a
    push won't see a half-built tree.
    """
    import shutil
    from datetime import datetime, timezone

    share_dir.mkdir(parents=True, exist_ok=True)
    target = share_dir / case_label
    staging = share_dir / f".{case_label}.staging"
    if staging.exists():
        shutil.rmtree(staging)
    print(f"verify_production: copying to share staging {staging} ...")
    shutil.copytree(out_root, staging)
    if target.exists():
        backup = share_dir / f".{case_label}.prev"
        if backup.exists():
            shutil.rmtree(backup)
        target.rename(backup)
    staging.rename(target)
    backup = share_dir / f".{case_label}.prev"
    if backup.exists():
        shutil.rmtree(backup)
    # Expose group/world readability so a portal user can fetch.
    try:
        for p in target.rglob("*"):
            mode = 0o644 if p.is_file() else 0o755
            p.chmod(mode)
        target.chmod(0o755)
    except Exception as e:
        print(f"verify_production: WARN chmod failed: {e}")
    print(f"verify_production: pushed to {target}")
    print(f"verify_production: index.html -> {target / 'index.html'}")

    # Add an entry for this case to <share_dir>/index.html if one isn't
    # already present. NEVER rewrite existing entries -- the user's
    # curated descriptions and ordering are theirs to manage. Idempotent:
    # if a link to <case_label>/ already exists (with any text), this is
    # a no-op so subsequent pushes don't accumulate duplicates.
    top_index = share_dir / "index.html"
    href_canonical = f"{case_label}/"
    new_li = (
        f'    <li><a href="{html.escape(href_canonical)}">\n'
        f'      <div class="title">{html.escape(case_label)}</div>\n'
        f'      <div class="desc">verify_production tape '
        f'(pushed {datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")}). '
        f'Edit this description in {top_index.name} to taste.</div>\n'
        f'    </a></li>'
    )

    if top_index.is_file():
        text = top_index.read_text()
        # Detect any existing link pointing at this case (any of the
        # common forms: "<label>/", "<label>/index.html", "./<label>/",
        # leading slash, etc.). If found, leave the index untouched.
        href_re = re.compile(
            r'href\s*=\s*[\"\']\.?/?'
            + re.escape(case_label)
            + r'/?(?:index\.html)?[\"\']',
            re.IGNORECASE,
        )
        if href_re.search(text):
            print(f"verify_production: {top_index} already links to "
                  f"{case_label}/ -- not modifying it.")
        else:
            # Insert just before the LAST </ul> in the document, so we
            # land inside the main listing if there's only one <ul>, and
            # don't break any nested submenus.
            i = text.rfind("</ul>")
            if i < 0:
                print(f"verify_production: {top_index} has no </ul> -- "
                      f"not auto-inserting. Add a link manually:")
                print(f"  <li><a href='{href_canonical}'>{case_label}</a></li>")
            else:
                new_text = text[:i] + new_li + "\n  " + text[i:]
                top_index.write_text(new_text)
                try: top_index.chmod(0o644)
                except Exception: pass
                print(f"verify_production: inserted {case_label} entry into "
                      f"{top_index} (existing entries left untouched)")
    else:
        # No top-level index yet -- write a minimal one with this entry.
        # User is free to restyle / curate.
        minimal = (
            "<!doctype html><html><head><meta charset='utf-8'>"
            f"<title>{html.escape(share_dir.name)} -- FME verify cases</title>"
            "<style>body{font-family:sans-serif;max-width:720px;margin:2em auto;"
            "padding:0 1em}h1{border-bottom:2px solid #333;padding-bottom:0.3em}"
            "ul{list-style:none;padding:0}li{margin-bottom:12px}"
            "li a{display:block;background:#fff;border:1px solid #ddd;"
            "border-radius:8px;padding:14px 18px;text-decoration:none;color:#222}"
            "li a:hover{border-color:#2c3e6b}.title{font-weight:600;color:#1a2744}"
            ".desc{font-size:0.85em;color:#666;margin-top:2px}</style></head><body>"
            f"<h1>FME verify cases</h1>\n<ul>\n{new_li}\n</ul>\n</body></html>"
        )
        top_index.write_text(minimal)
        try: top_index.chmod(0o644); share_dir.chmod(0o755)
        except Exception: pass
        print(f"verify_production: wrote minimal top-level index {top_index}")

    # NERSC science-portal hint (only meaningful on perlmutter / $PSCRATCH).
    pscratch = os.environ.get("PSCRATCH", "")
    if pscratch and str(target).startswith(pscratch):
        rel = Path(share_dir).relative_to(pscratch).as_posix()
        user = os.environ.get("USER", "")
        if user:
            print(f"verify_production: NERSC portal URL (if web-share enabled):")
            print(f"  https://portal.nersc.gov/pscratch/sd/{user[0]}/{user}/{rel}/")
            print(f"  https://portal.nersc.gov/pscratch/sd/{user[0]}/{user}/{rel}/{case_label}/")


def open_concat(files):
    """Open and concat multiple netCDF files along `time` without using
    dask. Avoids xr.open_mfdataset's dask requirement, which isn't always
    available in the analysis env on NERSC. For our use case (monthly
    files, dozens at most) eager-loading per-file then concatenating is
    fast enough and pulls only one chunk-manager dep (numpy)."""
    datasets = [
        xr.open_dataset(f, decode_times=False, mask_and_scale=False)
        for f in files
    ]
    if len(datasets) == 1:
        return datasets[0]
    # Use first dataset's attrs; per-record vars get concat'd along time;
    # static vars (no time dim) are taken from the first file.
    ds = xr.concat(datasets, dim="time", data_vars="minimal",
                   coords="minimal", compat="override",
                   combine_attrs="override")
    return ds


def main():
    args = parse_args()
    rundir = os.path.abspath(args.rundir)
    if not os.path.isdir(rundir):
        print(f"ERROR: rundir does not exist: {rundir}", file=sys.stderr)
        sys.exit(2)

    case_label = args.label or case_label_from_rundir(rundir)
    out_root = Path(args.outdir).expanduser().resolve() / "verify_production"
    out_root.mkdir(parents=True, exist_ok=True)
    title = f"FME production verify -- {case_label}" + (f" -- year {args.year}" if args.year else "")

    print(f"verify_production: rundir = {rundir}")
    print(f"verify_production: outdir = {out_root}")
    if args.year:
        print(f"verify_production: year   = {args.year}")
    print(f"verify_production: include _inst = {args.include_inst}")

    toc = []
    manifest = {}  # {stream_key: [{name, units, ts, clim}, ...]} for the JS viewer
    # Per-stream collected state. Rendering is deferred until after the
    # 5D-vs-1D cross-checks so cross-check sanity entries can be appended
    # to the 5D stream's sanity report before HTML is built.
    results = {}  # stream_key -> dict(files, sanity, cards, ds_or_None)

    for stream in STREAMS:
        anchor = f"sec_{stream['key']}"
        toc.append(f"<a href='#{anchor}'>{html.escape(stream['label'])}</a>")
        manifest[stream["key"]] = []
        results[stream["key"]] = {
            "anchor": anchor,
            "files": [],
            "sanity": [],
            "cards": [],
            "ds": None,
            "no_files": False,
        }

        all_files = stream_files(rundir, stream, args.year)
        if not all_files:
            results[stream["key"]]["no_files"] = True
            print(f"  [{stream['key']}] no files; skipping")
            continue

        files, legacy = filter_new_format(all_files)
        n_legacy = len(legacy)
        if not files:
            results[stream["key"]]["sanity"] = [
                ("FAIL", f"all {len(all_files)} file(s) are legacy format "
                         f"(no `time` dim / `time_bnds`); rebuild with the post-2026-05-04 "
                         f"FME testmod and rerun")
            ]
            print(f"  [{stream['key']}] all {len(all_files)} files are legacy; skipping")
            continue

        print(f"  [{stream['key']}] opening {len(files)} new-format file(s) "
              f"(skipping {n_legacy} legacy)")

        # mask_and_scale=False -> we see RAW values (no auto-NaN of _FillValue
        # cells), which is what check_fills needs to spot the gotcha #36
        # fractional-fill pollution. Plot helpers mask via _mask_invalid
        # so they work in either mode. open_concat avoids dask.
        try:
            ds = open_concat(files)
        except Exception as e:
            msg = f"open failed: {e}"
            results[stream["key"]]["files"] = files
            results[stream["key"]]["sanity"] = [("FAIL", msg)]
            print(f"    [FAIL] {msg}")
            continue

        is_mpas = stream["key"].startswith("mpas")
        sanity = sanity_check(ds, stream["label"], len(files), n_legacy, is_mpas,
                              args.pollution_threshold)
        for lvl, msg in sanity:
            if lvl != "OK":
                print(f"    [{lvl}] {msg}")

        cards = []
        var_dir = out_root / stream["key"]
        var_dir.mkdir(parents=True, exist_ok=True)
        plottable = sorted([v for v in ds.data_vars if is_plottable(ds[v])])
        if not args.include_inst:
            plottable = [v for v in plottable if not v.endswith("_inst")]
        print(f"    {len(plottable)} plottable variable(s)")

        for varname in plottable:
            try:
                da = ds[varname]
                base = safe_var_filename(varname)
                ts_file = var_dir / f"{base}.timeseries.png"
                clim_file = var_dir / f"{base}.climatology.png"

                ts_made = plot_timeseries(da, ds, varname, ts_file, args.dpi,
                                          args.pollution_threshold)
                clim_made = plot_climatology(da, ds, varname, clim_file, args.dpi,
                                             args.pollution_threshold)

                ts_rel = (ts_file.relative_to(out_root)).as_posix() if ts_made else None
                clim_rel = (clim_file.relative_to(out_root)).as_posix() if clim_made else None
                units = da.attrs.get("units", "")
                cards.append(
                    render_var_card(stream["key"], varname, units, ts_rel, clim_rel)
                )
                manifest[stream["key"]].append({
                    "name": varname,
                    "units": units,
                    "ts": ts_rel,
                    "clim": clim_rel,
                })
            except Exception as e:
                msg = f"plot failed: {e}"
                print(f"    [WARN] {varname}: {msg}")
                cards.append(
                    f"<div class='var-card'><span class='var-name'>{html.escape(varname)}</span>"
                    f" <span class='sanity fail'>{html.escape(msg)}</span></div>"
                )

        results[stream["key"]]["files"] = files
        results[stream["key"]]["sanity"] = sanity
        results[stream["key"]]["cards"] = cards
        results[stream["key"]]["ds"] = ds

    # ----- Pass 2: 5D-vs-1D cross-checks -----------------------------------
    for stream in STREAMS:
        peer_key = stream.get("companion_of")
        if not peer_key:
            continue
        ds_5d = results[stream["key"]]["ds"]
        ds_1d = results.get(peer_key, {}).get("ds")
        if ds_5d is None or ds_1d is None:
            results[stream["key"]]["sanity"].append(
                ("WARN", f"5D-vs-1D cross-check: peer stream '{peer_key}' "
                         f"unavailable; skipping cross-check")
            )
            continue
        print(f"  [{stream['key']}] cross-checking against {peer_key}")
        try:
            cc = cross_check_5d_vs_1d(ds_5d, ds_1d,
                                      rel_tol=args.cross_check_rel_tol,
                                      abs_tol=args.cross_check_abs_tol)
        except Exception as e:
            cc = [("WARN", f"5D-vs-1D cross-check failed: {e}")]
            traceback.print_exc()
        for lvl, msg in cc:
            if lvl != "OK":
                print(f"    [{lvl}] {msg}")
        results[stream["key"]]["sanity"].extend(cc)

    # ----- Pass 3: render ---------------------------------------------------
    sections = []
    for stream in STREAMS:
        res = results[stream["key"]]
        anchor = res["anchor"]
        if res["no_files"]:
            sections.append(
                f"<h2 id='{anchor}'>{html.escape(stream['label'])}</h2>"
                f"<p><em>no files found</em></p>"
            )
            continue
        sections.append(render_section(stream["key"], stream["label"], anchor,
                                       res["files"], res["sanity"], res["cards"]))
        if res["ds"] is not None:
            res["ds"].close()

    out_html = write_html(out_root, title, sections, toc, manifest)
    print(f"\nverify_production: wrote {out_html}")
    print(f"verify_production: open with file://{out_html} or copy {out_root} to a web server.")

    # Optional push to a shared web-viewable location (e.g. NERSC's $PSCRATCH/share).
    if args.share_dir:
        push_to_share(out_root, Path(args.share_dir).expanduser().resolve(),
                      case_label)

    # Aggregate sanity results into a process exit code so CI / callers can
    # distinguish a clean tape from FAILs (previously main() always fell through
    # to exit 0 even with FAIL sanity entries or missing streams). A FAIL in any
    # stream's sanity report -- or a missing PRIMARY stream -- fails the run;
    # WARN does not, and a missing 5D companion is tolerated (may be disabled).
    # The HTML/share output above is still produced on failure.
    n_fail = 0
    fail_lines = []
    for stream in STREAMS:
        res = results[stream["key"]]
        if res["no_files"]:
            if stream.get("companion_of"):
                print(f"  NOTE [{stream['key']}]: no files (5D companion "
                      f"disabled?) -- not counted as a failure")
            else:
                n_fail += 1
                fail_lines.append(f"  FAIL [{stream['key']}]: no files found")
            continue
        for lvl, msg in res["sanity"]:
            if lvl == "FAIL":
                n_fail += 1
                fail_lines.append(f"  FAIL [{stream['key']}]: {msg}")

    if n_fail:
        print(f"\nverify_production: {n_fail} FAIL(s) -- exiting nonzero:")
        for line in fail_lines:
            print(line)
        sys.exit(1)
    print("\nverify_production: all streams PASS (no FAIL sanity entries)")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(130)
    except Exception:
        traceback.print_exc()
        sys.exit(1)
