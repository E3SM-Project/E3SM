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
        "exclude": [],
    },
    {
        "key": "mpaso_derived",
        "label": "MPAS-O fmeDerivedFields",
        "pattern": "*mpaso.hist.am.fmeDerivedFields.*.remapped.nc",
        "exclude": [],
    },
    {
        "key": "mpassi_derived",
        "label": "MPAS-SI fmeSeaiceDerivedFields",
        "pattern": "*mpassi.hist.am.fmeSeaiceDerivedFields.*.remapped.nc",
        "exclude": [],
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

    sections = []
    toc = []
    manifest = {}  # {stream_key: [{name, units, ts, clim}, ...]} for the JS viewer

    for stream in STREAMS:
        anchor = f"sec_{stream['key']}"
        toc.append(f"<a href='#{anchor}'>{html.escape(stream['label'])}</a>")
        manifest[stream["key"]] = []

        all_files = stream_files(rundir, stream, args.year)
        if not all_files:
            sections.append(
                f"<h2 id='{anchor}'>{html.escape(stream['label'])}</h2>"
                f"<p><em>no files found</em></p>"
            )
            print(f"  [{stream['key']}] no files; skipping")
            continue

        files, legacy = filter_new_format(all_files)
        n_legacy = len(legacy)
        if not files:
            sanity = [("FAIL", f"all {len(all_files)} file(s) are legacy format "
                               f"(no `time` dim / `time_bnds`); rebuild with the post-2026-05-04 "
                               f"FME testmod and rerun")]
            sections.append(render_section(stream["key"], stream["label"],
                                           anchor, [], sanity, []))
            print(f"  [{stream['key']}] all {len(all_files)} files are legacy; skipping")
            continue

        print(f"  [{stream['key']}] opening {len(files)} new-format file(s) "
              f"(skipping {n_legacy} legacy)")

        # mask_and_scale=False -> we see RAW values (no auto-NaN of _FillValue
        # cells), which is what check_fills needs to spot the gotcha #36
        # fractional-fill pollution. Plot helpers mask via _mask_invalid
        # so they work in either mode.
        try:
            ds = xr.open_mfdataset(files, decode_times=False, combine="by_coords",
                                   mask_and_scale=False)
        except Exception as e:
            try:
                ds = xr.open_mfdataset(files, decode_times=False, combine="nested",
                                       concat_dim="time", mask_and_scale=False)
            except Exception as e2:
                msg = f"open_mfdataset failed: {e2}"
                sections.append(render_section(stream["label"], anchor, files,
                                               [("FAIL", msg)], []))
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
        sections.append(render_section(stream["key"], stream["label"],
                                       anchor, files, sanity, cards))
        ds.close()

    out_html = write_html(out_root, title, sections, toc, manifest)
    print(f"\nverify_production: wrote {out_html}")
    print(f"verify_production: open with file://{out_html} or copy {out_root} to a web server.")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(130)
    except Exception:
        traceback.print_exc()
        sys.exit(1)
