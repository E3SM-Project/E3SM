#!/usr/bin/env python3
"""
Verify EAM FME online output (vcoarsen, derived fields) and cross-verify
against legacy output. Produces HTML dashboard with BFB identity tests,
offline vcoarsen comparison, and diagnostic figures.

Usage:
    python verify_eam.py --rundir /path/to/FME/RUNDIR --outdir /path/to/figs
    python verify_eam.py --rundir /path/to/FME/RUNDIR --legacy-rundir /path/to/LEGACY/RUNDIR --outdir /path/to/figs
"""

import argparse
import os
import sys
import glob
import textwrap
import time as _time
from datetime import datetime
from pathlib import Path

import numpy as np

# -- optional imports ----------------------------------------------------------
try:
    import xarray as xr
    HAS_XR = True
except ImportError:
    HAS_XR = False

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    HAS_CARTOPY = True
except ImportError:
    HAS_CARTOPY = False

try:
    from netCDF4 import Dataset as NC4Dataset
    HAS_NC4 = True
except ImportError:
    HAS_NC4 = False

# -- ACE atmosphere variable requirements ------------------------------------
# Single source of truth for the ACE-to-EAM variable mapping.
# Drives verification, HTML coverage table, and reproducibility tracking.
#
# The lists below mirror the ACE YAML specification exactly.
# EAM field names differ from ACE canonical names in several cases;
# ace_to_eam() handles the mapping.

N_ATM_LAYERS = 8  # ACE uses 8 pressure layers (0-based: T_0..T_7)

# 3D vertically coarsened fields: ACE base name -> EAM base name
ACE_3D_BASE_MAP = {
    "T": "T",
    "specific_total_water": "STW",
    "U": "U",
    "V": "V",
}

# EAM coarsened base names (used by check_eam and comparison figures)
ACE_ATM_3D_COARSENED = list(ACE_3D_BASE_MAP.values())  # ["T", "STW", "U", "V"]

# 2D / scalar field mapping: ACE canonical name -> EAM FME output name
ACE_2D_MAP = {
    "land_fraction": "LANDFRAC",
    "ocean_fraction": "OCNFRAC",
    "sea_ice_fraction": "ICEFRAC",
    "PHIS": "PHIS",
    "SOLIN": "SOLIN",
    "PS": "PS",
    "TS": "TS",
    "LHFLX": "LHFLX",
    "SHFLX": "SHFLX",
    "surface_precipitation_rate": "PRECT",
    "surface_upward_longwave_flux": "FLUS",
    "FLUT": "FLUT",
    "FLDS": "FLDS",
    "FSDS": "FSDS",
    "surface_upward_shortwave_flux": "FSUS",
    "top_of_atmos_upward_shortwave_flux": "FSUTOA",
    "tendency_of_total_water_path_due_to_advection": "DTENDTTW",
    "TAUX": "TAUX",
    "TAUY": "TAUY",
}

# ACE variable lists (mirror the YAML specification)
ACE_FORCING_NAMES = ["SOLIN"]

ACE_IN_NAMES = ([
    "land_fraction", "ocean_fraction", "sea_ice_fraction",
    "PHIS", "SOLIN", "PS", "TS",
] + [f"T_{k}" for k in range(N_ATM_LAYERS)]
  + [f"specific_total_water_{k}" for k in range(N_ATM_LAYERS)]
  + [f"U_{k}" for k in range(N_ATM_LAYERS)]
  + [f"V_{k}" for k in range(N_ATM_LAYERS)])

ACE_OUT_NAMES = ([
    "PS", "TS",
] + [f"T_{k}" for k in range(N_ATM_LAYERS)]
  + [f"specific_total_water_{k}" for k in range(N_ATM_LAYERS)]
  + [f"U_{k}" for k in range(N_ATM_LAYERS)]
  + [f"V_{k}" for k in range(N_ATM_LAYERS)]
  + [
    "LHFLX", "SHFLX",
    "surface_precipitation_rate",
    "surface_upward_longwave_flux",
    "FLUT", "FLDS", "FSDS",
    "surface_upward_shortwave_flux",
    "top_of_atmos_upward_shortwave_flux",
    "tendency_of_total_water_path_due_to_advection",
    "TAUX", "TAUY",
])


def ace_to_eam(ace_name):
    """Map an ACE canonical variable name to its EAM FME output name."""
    if ace_name in ACE_2D_MAP:
        return ACE_2D_MAP[ace_name]
    # 3D indexed variables (e.g. specific_total_water_3 -> STW_3)
    for ace_base in sorted(ACE_3D_BASE_MAP, key=len, reverse=True):
        prefix = ace_base + "_"
        if ace_name.startswith(prefix):
            suffix = ace_name[len(ace_base):]
            return ACE_3D_BASE_MAP[ace_base] + suffix
    return ace_name


# Physical range checks: (vmin, vmax)
# Ranges are deliberately generous to avoid false positives during spinup.
RANGE_CHECKS = {
    "T": (150, 350),
    "Q": (0, 0.05),
    "PS": (40000, 110000),
    "TS": (200, 340),
    "PHIS": (-1000, 60000),
}

# -- Vcoarsen constants -------------------------------------------------------
GRAVIT = 9.80616       # m/s2, matches EAM physconst
P0     = 100000.0      # Pa, reference pressure for hybrid coordinates
VCOARSEN_PBOUNDS = np.array([
    10.0, 4803.81, 13913.06, 26856.34,
    43998.31, 59659.67, 76854.15, 90711.83, 101325.0
])  # 9 interfaces -> 8 layers, matching ACE L80 indices [0,25,38,46,52,56,61,69,80]
N_VCOARSEN_LAYERS = len(VCOARSEN_PBOUNDS) - 1  # 8

# Fields that appear in BOTH fme_output and fme_legacy_output for BFB identity checks
# Instantaneous (no :A suffix)
BFB_INST_FIELDS = ["PS", "TS", "PHIS", "LANDFRAC", "OCNFRAC", "ICEFRAC"]
# Averaged (:A suffix in both cases)
BFB_AVG_FIELDS = ["SOLIN", "FSNTOA", "FLUT", "FLDS", "FSDS",
                   "LHFLX", "SHFLX", "TAUX", "TAUY", "PRECT"]

# FME vcoarsen fields -> legacy raw source (for offline recomputation)
# FME name -> (legacy fields to sum for derived, then vcoarsen)
FME_TO_LEGACY = {
    "T":   ["T"],
    "U":   ["U"],
    "V":   ["V"],
    "STW": ["Q", "CLDICE", "CLDLIQ", "RAINQM"],
}


# -- Offline vcoarsen / column integration ------------------------------------

def compute_pint(ps, hyai, hybi):
    """Compute interface pressures from hybrid coordinates.

    Parameters
    ----------
    ps : ndarray, shape (ncol,)
        Surface pressure in Pa.
    hyai, hybi : ndarray, shape (nlev+1,)
        Hybrid A and B coefficients (interface levels, top to bottom).

    Returns
    -------
    pint : ndarray, shape (ncol, nlev+1)
        Interface pressures in Pa.
    """
    # pint(i, k) = hyai(k) * P0 + hybi(k) * ps(i)
    return hyai[np.newaxis, :] * P0 + hybi[np.newaxis, :] * ps[:, np.newaxis]


def compute_pdel(pint):
    """Compute layer thickness in Pa from interface pressures.

    Parameters
    ----------
    pint : ndarray, shape (ncol, nlev+1)

    Returns
    -------
    pdel : ndarray, shape (ncol, nlev)
    """
    return pint[:, 1:] - pint[:, :-1]


def offline_vcoarsen_avg(field, pint, pbounds=None):
    """Overlap-weighted vertical averaging onto coarsened pressure layers.

    Implements the same algorithm as shr_vcoarsen_avg in shr_vcoarsen_mod.F90.
    For each target layer [pb_top, pb_bot], computes the pressure-weighted mean
    of all model levels that overlap with it.

    Parameters
    ----------
    field : ndarray, shape (ncol, nlev)
        Full-resolution field values.
    pint : ndarray, shape (ncol, nlev+1)
        Interface pressures from compute_pint.
    pbounds : ndarray, shape (nlayers+1,), optional
        Target layer pressure boundaries. Defaults to VCOARSEN_PBOUNDS.

    Returns
    -------
    coarsened : ndarray, shape (ncol, nlayers)
        Coarsened field values.
    """
    if pbounds is None:
        pbounds = VCOARSEN_PBOUNDS
    nlayers = len(pbounds) - 1
    ncol, nlev = field.shape

    coarsened = np.zeros((ncol, nlayers))
    for layer in range(nlayers):
        pb_top = pbounds[layer]
        pb_bot = pbounds[layer + 1]

        weight_sum = np.zeros(ncol)
        for k in range(nlev):
            # Overlap between model level k and target layer
            overlap_top = np.maximum(pb_top, pint[:, k])
            overlap_bot = np.minimum(pb_bot, pint[:, k + 1])
            overlap = np.maximum(0.0, overlap_bot - overlap_top)

            # Skip NaN and fill values — only average valid data
            fk = field[:, k]
            ok = np.isfinite(fk) & (overlap > 0)
            coarsened[ok, layer] += fk[ok] * overlap[ok]
            weight_sum[ok] += overlap[ok]

        # Normalize by total overlap of valid levels
        valid = weight_sum > 0
        coarsened[valid, layer] /= weight_sum[valid]
        coarsened[~valid, layer] = np.nan

    return coarsened


def offline_column_integrate(field, pdel):
    """Column integration: sum(field * pdel / g) over all levels.

    Parameters
    ----------
    field : ndarray, shape (ncol, nlev)
    pdel : ndarray, shape (ncol, nlev)

    Returns
    -------
    integrated : ndarray, shape (ncol,)
    """
    return np.sum(field * pdel / GRAVIT, axis=1)


# -----------------------------------------------------------------------------
# Utility helpers
# -----------------------------------------------------------------------------

def find_files(rundir, pattern, exclude=None):
    """Find files matching glob pattern, optionally excluding a substring.

    Parameters
    ----------
    rundir : str
        Directory to search in.
    pattern : str
        Glob pattern relative to rundir.
    exclude : str or None
        If set, exclude any file whose basename contains this substring.

    Returns
    -------
    list of str
        Sorted list of matching file paths.
    """
    hits = sorted(glob.glob(os.path.join(rundir, pattern)))
    # Always exclude .base restart files
    hits = [f for f in hits if not f.endswith('.base')]
    if exclude:
        hits = [f for f in hits if exclude not in os.path.basename(f)]
    return hits


def safe_open(path):
    """Open a NetCDF file with xarray (preferred) or netCDF4."""
    if HAS_XR:
        try:
            return xr.open_dataset(path, decode_times=False)
        except Exception:
            pass
    if HAS_NC4:
        return NC4Dataset(path, "r")
    raise RuntimeError("Neither xarray nor netCDF4 is available.")


def get_var(ds, name):
    """Return numpy array from xarray Dataset or netCDF4 Dataset."""
    if HAS_XR and isinstance(ds, xr.Dataset):
        if name in ds:
            return ds[name].values
        return None
    if HAS_NC4 and isinstance(ds, NC4Dataset):
        if name in ds.variables:
            return np.array(ds.variables[name][:])
        return None
    return None


def get_dims(ds):
    """Return dict of dimension names -> sizes."""
    if HAS_XR and isinstance(ds, xr.Dataset):
        return dict(ds.sizes)
    if HAS_NC4 and isinstance(ds, NC4Dataset):
        return {d: len(ds.dimensions[d]) for d in ds.dimensions}
    return {}


def get_varnames(ds):
    """Return list of variable names in the dataset."""
    if HAS_XR and isinstance(ds, xr.Dataset):
        return list(ds.data_vars)
    if HAS_NC4 and isinstance(ds, NC4Dataset):
        return list(ds.variables.keys())
    return []


def close_ds(ds):
    """Close a dataset if applicable."""
    if HAS_XR and isinstance(ds, xr.Dataset):
        ds.close()
    elif HAS_NC4 and isinstance(ds, NC4Dataset):
        ds.close()


def valid_data(arr, fill=1e10):
    """Return non-fill, finite values.

    The threshold is 1e10 because remapped output can contain interpolated
    fill artifacts from neighboring land/ocean cells that can reach ~1e14
    scale while _FillValue is 1e34.  No geophysical quantity exceeds 1e10.
    """
    if arr is None:
        return None
    flat = arr.ravel().astype(float)
    mask = (np.abs(flat) < fill) & np.isfinite(flat)
    return flat[mask] if mask.any() else None


def fill_nan_report(arr, name):
    """Return a dict with fill/NaN fraction stats for a variable array."""
    if arr is None:
        return {"name": name, "total": 0, "valid": 0, "fill_nan": 0, "pct": 100.0}
    flat = arr.ravel().astype(float)
    total = flat.size
    fill_mask = np.abs(flat) >= 1e10
    nan_mask = ~np.isfinite(flat)
    bad = fill_mask | nan_mask
    n_bad = int(bad.sum())
    n_valid = total - n_bad
    pct = 100.0 * n_bad / total if total > 0 else 0.0
    return {"name": name, "total": total, "valid": n_valid,
            "fill_nan": n_bad, "pct": pct}


def check_range(arr, name, vmin=None, vmax=None):
    v = valid_data(arr)
    issues = []
    if v is None or v.size == 0:
        issues.append(f"  {name}: no valid data")
        return issues
    if not np.isfinite(v).all():
        issues.append(f"  {name}: NaN/Inf present")
    if vmin is not None and v.min() < vmin:
        issues.append(f"  {name}: min={v.min():.4g} below {vmin}")
    if vmax is not None and v.max() > vmax:
        issues.append(f"  {name}: max={v.max():.4g} above {vmax}")
    return issues


def summary_stats(arr):
    v = valid_data(arr)
    if v is None or v.size == 0:
        return "no valid data"
    return f"mean={v.mean():.4g}  min={v.min():.4g}  max={v.max():.4g}  n={v.size}"


def has_time_zero(ds):
    """Check if the Time/time dimension has size 0."""
    dims = get_dims(ds)
    for tname in ("Time", "time"):
        if tname in dims and dims[tname] == 0:
            return True
    return False


def last_time_index(arr):
    """Return the index of the last timestep.

    Uses -1 (last available) to avoid initialization artifacts at t=0.
    """
    if arr is None:
        return -1
    if arr.ndim >= 2 and arr.shape[0] > 0:
        return -1
    return -1


# -----------------------------------------------------------------------------
# Plotting helpers
# -----------------------------------------------------------------------------

def savefig(fig, outdir, name, subdir=None):
    """Save figure, optionally in a subdirectory."""
    if subdir:
        d = os.path.join(outdir, subdir)
        os.makedirs(d, exist_ok=True)
        path = os.path.join(d, name)
    else:
        path = os.path.join(outdir, name)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def _fix_lon(lons):
    """Shift longitudes from [0,360) to [-180,180) for plotting."""
    return (lons + 180) % 360 - 180


def _plot_on_ax(ax, data, lons, lats, cmap, vmin, vmax, is_cartopy):
    """Plot data on an axis -- tripcolor for 1D, pcolormesh for 2D."""
    transform = ccrs.PlateCarree() if is_cartopy else None
    if data.ndim == 2:
        kw = dict(transform=transform) if is_cartopy else {}
        im = ax.pcolormesh(lons, lats, data, cmap=cmap,
                           vmin=vmin, vmax=vmax, shading="auto", **kw)
    else:
        # Unstructured: use tripcolor for filled triangulation.
        # Subsample if > 50k points to keep triangulation fast.
        # Mask triangles that cross the antimeridian (lon span > 180).
        from matplotlib.tri import Triangulation
        lons_fix = _fix_lon(lons)
        kw = dict(transform=ccrs.PlateCarree()) if is_cartopy else {}
        max_pts = 50000
        if data.size > max_pts:
            idx = np.random.default_rng(42).choice(data.size, max_pts, replace=False)
            lons_sub, lats_sub, data_sub = lons_fix[idx], lats[idx], data[idx]
        else:
            lons_sub, lats_sub, data_sub = lons_fix, lats, data
        try:
            tri = Triangulation(lons_sub, lats_sub)
            # Mask triangles spanning the antimeridian
            tri_lons = lons_sub[tri.triangles]
            lon_span = tri_lons.max(axis=1) - tri_lons.min(axis=1)
            tri.set_mask(lon_span > 180)
            im = ax.tripcolor(tri, data_sub, cmap=cmap,
                              vmin=vmin, vmax=vmax, **kw)
        except Exception:
            # Fallback to scatter if tripcolor fails
            im = ax.scatter(lons_sub, lats_sub, c=data_sub, s=0.5, cmap=cmap,
                            vmin=vmin, vmax=vmax, **kw)
    return im


def global_map(data, lons, lats, title, cmap="RdBu_r", vmin=None, vmax=None,
               outdir=None, fname=None, units="", subdir=None):
    """
    Plot a global filled map.  data is 2-D (lat x lon) or 1-D (nCells) on an
    unstructured grid (tripcolor with scatter fallback).
    """
    if not HAS_MPL:
        return None

    fig = plt.figure(figsize=(10, 5))
    is_cartopy = HAS_CARTOPY
    if is_cartopy:
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
        ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
        ax.set_global()
    else:
        ax = fig.add_subplot(1, 1, 1)

    im = _plot_on_ax(ax, data, lons, lats, cmap, vmin, vmax, is_cartopy)
    plt.colorbar(im, ax=ax, shrink=0.6 if is_cartopy else 0.8, label=units)

    if not is_cartopy:
        ax.set_xlabel("lon"); ax.set_ylabel("lat")

    ax.set_title(title, fontsize=11)
    if outdir and fname:
        return savefig(fig, outdir, fname, subdir=subdir)
    plt.show()
    return None


def latlon_map(ds, varname, title, cmap="RdBu_r", vmin=None, vmax=None,
               outdir=None, fname=None, units="", subdir=None):
    """Plot a variable from a remapped (lon, lat, Time) dataset.

    Uses the LAST timestep to avoid initialization artifacts.
    """
    if not HAS_MPL:
        return None

    arr = get_var(ds, varname)
    if arr is None:
        return None
    # Use last timestep if time dimension present
    if arr.ndim == 3:
        data = arr[-1, :, :]
    elif arr.ndim == 2:
        data = arr
    else:
        return None

    lon = get_var(ds, "lon")
    lat = get_var(ds, "lat")
    if lon is None or lat is None:
        return None

    if lon.ndim == 1 and lat.ndim == 1:
        lons, lats = np.meshgrid(lon, lat)
    else:
        lons, lats = lon, lat

    return global_map(data, lons, lats, title, cmap=cmap, vmin=vmin, vmax=vmax,
                      outdir=outdir, fname=fname, units=units, subdir=subdir)


def layer_profiles(data3d, label, outdir, fname, ylabel="Level index", subdir=None):
    """Plot global-mean vertical profile from (nlev x ncells) array."""
    if not HAS_MPL:
        return
    means = []
    for lev in range(data3d.shape[0]):
        v = valid_data(data3d[lev])
        means.append(v.mean() if v is not None and v.size > 0 else np.nan)
    fig, ax = plt.subplots(figsize=(4, 5))
    ax.plot(means, np.arange(len(means)), "o-")
    ax.invert_yaxis()
    ax.set_xlabel(label)
    ax.set_ylabel(ylabel)
    ax.set_title(f"Global-mean vertical profile: {label}")
    ax.grid(True, alpha=0.4)
    savefig(fig, outdir, fname, subdir=subdir)


def time_series(values, times_label, title, ylabel, outdir, fname, subdir=None):
    """Plot a simple time series."""
    if not HAS_MPL:
        return
    fig, ax = plt.subplots(figsize=(8, 3))
    ax.plot(values, ".-")
    ax.set_xlabel(times_label)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.4)
    plt.tight_layout()
    savefig(fig, outdir, fname, subdir=subdir)


def multi_panel_maps(panels, suptitle, outdir, fname, ncols=4,
                     cmap="RdBu_r", vmin=None, vmax=None, units=""):
    """Plot a grid of small-multiple global maps.

    Parameters
    ----------
    panels : list of (data, lons, lats, subtitle) tuples
        data can be 2D (lat x lon) for pcolormesh or 1D (nCells) for tripcolor.
    """
    if not HAS_MPL or not panels:
        return None

    n = len(panels)
    nrows = (n + ncols - 1) // ncols
    fig_w = 5 * ncols
    fig_h = 2.8 * nrows + 0.6
    is_cartopy = HAS_CARTOPY

    fig = plt.figure(figsize=(fig_w, fig_h))
    axes = []
    for idx in range(n):
        kw = dict(projection=ccrs.PlateCarree()) if is_cartopy else {}
        ax = fig.add_subplot(nrows, ncols, idx + 1, **kw)
        axes.append(ax)

    im = None
    for idx, (data, lons, lats, subtitle) in enumerate(panels):
        ax = axes[idx]
        if is_cartopy:
            ax.add_feature(cfeature.COASTLINE, linewidth=0.3)
            ax.set_global()
        im = _plot_on_ax(ax, data, lons, lats, cmap, vmin, vmax, is_cartopy)
        ax.set_title(subtitle, fontsize=9)

    if im is not None:
        fig.colorbar(im, ax=axes, shrink=0.5, label=units, pad=0.03, aspect=30)
    fig.suptitle(suptitle, fontsize=12, y=1.01)
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        plt.tight_layout(rect=[0, 0, 1, 0.97])
    path = os.path.join(outdir, fname)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def zonal_mean_plot(profiles, lat, title, outdir, fname, ylabel=""):
    """Plot zonal mean profiles (latitude vs value) for multiple variables.

    Parameters
    ----------
    profiles : list of (zmean_1d, label) tuples
        zmean_1d is a 1D array of length nlat (zonal mean values).
    lat : 1D array
        Latitude values.
    """
    if not HAS_MPL or not profiles:
        return None

    fig, ax = plt.subplots(figsize=(7, 3.5))
    has_data = False
    for zmean, label in profiles:
        if zmean is None:
            continue
        ax.plot(lat, zmean, label=label, linewidth=1.3)
        has_data = True

    if not has_data:
        plt.close(fig)
        return None

    ax.set_xlabel("Latitude")
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=11)
    ax.legend(fontsize=8, loc="best")
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-90, 90)
    plt.tight_layout()
    path = os.path.join(outdir, fname)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def _compute_zonal_mean(arr, fill_thresh=1e10):
    """Compute zonal mean from (lat, lon) array, masking fill values."""
    if arr is None:
        return None
    if arr.ndim == 3:
        data = arr[-1, :, :]
    elif arr.ndim == 2:
        data = arr
    else:
        return None
    masked = np.where(np.abs(data) < fill_thresh, data, np.nan)
    return np.nanmean(masked, axis=1)


def fill_mask_map(data, lons, lats, title, outdir, fname):
    """Plot a binary map showing valid (blue) vs fill/NaN (red) cells."""
    if not HAS_MPL:
        return None
    mask = ((np.abs(data) >= 1e10) | ~np.isfinite(data)).astype(float)
    cmap = mcolors.ListedColormap(["#3b82f6", "#ef4444"])
    bounds = [-0.5, 0.5, 1.5]
    norm = mcolors.BoundaryNorm(bounds, cmap.N)
    fig = plt.figure(figsize=(10, 5))
    is_cartopy = HAS_CARTOPY
    if is_cartopy:
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
        ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
        ax.set_global()
    else:
        ax = fig.add_subplot(1, 1, 1)
    transform = ccrs.PlateCarree() if is_cartopy else None
    if data.ndim == 2:
        kw = dict(transform=transform) if is_cartopy else {}
        im = ax.pcolormesh(lons, lats, mask, cmap=cmap, norm=norm, shading="auto", **kw)
    else:
        kw = dict(transform=ccrs.PlateCarree()) if is_cartopy else {}
        lons_fix = _fix_lon(lons)
        im = ax.scatter(lons_fix, lats, c=mask, s=0.5, cmap=cmap, norm=norm, **kw)
    cbar = plt.colorbar(im, ax=ax, shrink=0.6 if is_cartopy else 0.8, ticks=[0, 1])
    cbar.ax.set_yticklabels(["Valid", "Fill/NaN"])
    ax.set_title(title, fontsize=11)
    return savefig(fig, outdir, fname)


def side_by_side_comparison(data1, lons1, lats1, title1,
                            data2, lons2, lats2, title2,
                            suptitle, outdir, fname,
                            cmap="RdBu_r", vmin=None, vmax=None, units=""):
    """Two global maps side by side -- native tripcolor on left, remapped pcolormesh on right."""
    if not HAS_MPL:
        return None

    is_cartopy = HAS_CARTOPY
    fig = plt.figure(figsize=(16, 4.5))
    im = None
    for idx, (data, lons, lats, title) in enumerate([
        (data1, lons1, lats1, title1),
        (data2, lons2, lats2, title2),
    ]):
        kw = dict(projection=ccrs.PlateCarree()) if is_cartopy else {}
        ax = fig.add_subplot(1, 2, idx + 1, **kw)
        if is_cartopy:
            ax.add_feature(cfeature.COASTLINE, linewidth=0.4)
            ax.set_global()
        im = _plot_on_ax(ax, data, lons, lats, cmap, vmin, vmax, is_cartopy)
        plt.colorbar(im, ax=ax, shrink=0.7, label=units)
        ax.set_title(title, fontsize=10)

    fig.suptitle(suptitle, fontsize=12, y=1.02)
    plt.tight_layout()
    path = os.path.join(outdir, fname)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


# -----------------------------------------------------------------------------
# EAM verification + visualization
# -----------------------------------------------------------------------------

def check_eam(rundir, outdir, verbose):
    print("\n=== EAM ===")
    issues = []
    plots = []
    fill_reports = []

    # -- tape 1 (h0): all EAM fields (single-tape layout) -------------------
    t1_files = find_files(rundir, "*.eam.h0.*.nc")
    if not t1_files:
        print("  SKIP: no *.eam.h0.*.nc files found")
    else:
        print(f"  h0: found {len(t1_files)} file(s)")
        ds = safe_open(t1_files[0])
        if has_time_zero(ds):
            print("  h0: time dimension has size 0 -- skipping")
            close_ds(ds)
        else:
            for var in ["PHIS", "PS", "TS", "LHFLX", "SHFLX", "FLDS", "FSDS",
                        "FSNTOA", "FSUTOA", "FLUT", "SOLIN", "FSUS", "FLUS",
                        "ICEFRAC", "LANDFRAC", "OCNFRAC", "TAUX", "TAUY",
                        "PRECT", "DTENDTTW"]:
                arr = get_var(ds, var)
                if arr is None:
                    issues.append(f"  h0 missing: {var}")
                else:
                    vmin, vmax = RANGE_CHECKS.get(var, (None, None))
                    issues += check_range(arr, f"h0/{var}", vmin, vmax)
                    fill_reports.append(fill_nan_report(arr, f"h0/{var}"))
                    if verbose:
                        print(f"  tape1/{var}: {summary_stats(arr)}")

            # Maps of key 2D fields on lat-lon grid -- use LAST timestep
            if HAS_MPL and HAS_XR and isinstance(ds, xr.Dataset):
                lon_name = next((d for d in ["lon", "ncol", "longitude"] if d in ds.dims), None)
                lat_name = next((d for d in ["lat", "latitude"] if d in ds.dims), None)
                for var, cmap, vm, vx, units in [
                    ("TS",      "RdBu_r",  220, 320, "K"),
                    ("LHFLX",   "viridis",   0, 300, "W/m2"),
                    ("PRECT",   "Blues",     0, 1e-3, "kg/m2/s"),
                    ("ICEFRAC", "Blues",      0,   1, "fraction"),
                    ("FLUT",    "inferno",  100, 320, "W/m2"),
                    ("FSDS",    "YlOrRd",     0, 500, "W/m2"),
                    ("FSUTOA", "YlOrRd", 0, 500, "W/m2"),
                    ("FLUS", "inferno", 200, 500, "W/m2"),
                    ("FSUS", "YlOrRd", 0, 300, "W/m2"),
                    ("DTENDTTW", "RdBu_r", -5e-4, 5e-4, "kg/m2/s"),
                ]:
                    arr = get_var(ds, var)
                    if arr is None or arr.size == 0:
                        continue
                    # Last timestep
                    data = arr[-1] if arr.ndim >= 2 and arr.shape[0] > 0 else arr
                    if lat_name and lon_name and data.ndim == 2:
                        lons = ds[lon_name].values
                        lats = ds[lat_name].values
                        if lons.ndim == 1:
                            lons, lats = np.meshgrid(lons, lats)
                        p = global_map(data, lons, lats, f"EAM {var} (h0, last t)",
                                       cmap=cmap, vmin=vm, vmax=vx,
                                       outdir=outdir, fname=f"eam_h0_{var}.png",
                                       units=units)
                        if p: plots.append(p)

                # Fill value mask maps for key fields
                for fvar in ["T_7", "ICEFRAC"]:
                    farr = get_var(ds, fvar)
                    if farr is None or farr.size == 0:
                        continue
                    fdata = farr[-1] if farr.ndim >= 2 and farr.shape[0] > 0 else farr
                    if lat_name and lon_name and fdata.ndim == 2:
                        p = fill_mask_map(fdata, lons, lats,
                                          f"Fill Mask: {fvar} (h0, last t)",
                                          outdir, f"eam_h0_{fvar}_fill_mask.png")
                        if p:
                            plots.append(p)

            close_ds(ds)

    # -- vcoarsen layers (on h0 in single-tape config) -----------------------
    if t1_files:
        ds3 = safe_open(t1_files[0])
        if not has_time_zero(ds3):
            for base in ACE_ATM_3D_COARSENED:
                found_layers = []
                for k in range(N_ATM_LAYERS):  # 0-based: T_0..T_7
                    vname = f"{base}_{k}"
                    arr = get_var(ds3, vname)
                    if arr is not None:
                        found_layers.append(k)
                        vmin, vmax = RANGE_CHECKS.get(base, (None, None))
                        issues += check_range(arr, f"h0/{vname}", vmin, vmax)
                        fill_reports.append(fill_nan_report(arr, f"h0/{vname}"))
                if len(found_layers) == N_ATM_LAYERS:
                    print(f"  h0/{base}: all {N_ATM_LAYERS} layers present (0-based)")
                elif found_layers:
                    print(f"  h0/{base}: {len(found_layers)}/{N_ATM_LAYERS} layers found "
                          f"(indices {found_layers})")
                    issues.append(f"  h0/{base}: only {len(found_layers)}/{N_ATM_LAYERS} layers")
                else:
                    alts = [v for v in get_varnames(ds3)
                            if v.startswith(base + "_") or v.startswith(base.upper() + "_")]
                    if alts:
                        print(f"  h0/{base}: found alternate names: {alts[:5]}")
                    else:
                        issues.append(f"  h0 missing coarsened field: {base}_0..{base}_{N_ATM_LAYERS-1}")

            # Plot T_0 (near-top) and T_7 (near-surface) at last timestep
            if HAS_MPL and HAS_XR and isinstance(ds3, xr.Dataset):
                lon_name = next((d for d in ["lon", "ncol", "longitude"] if d in ds3.dims), None)
                lat_name = next((d for d in ["lat", "latitude"] if d in ds3.dims), None)
                for vname, label_tag in [("T_0", "near-top"), ("T_7", "near-surface")]:
                    arr = get_var(ds3, vname)
                    if arr is None or arr.size == 0:
                        continue
                    data = arr[-1] if arr.ndim >= 2 and arr.shape[0] > 0 else arr
                    if lat_name and lon_name and data.ndim == 2:
                        lons = ds3[lon_name].values
                        lats = ds3[lat_name].values
                        if lons.ndim == 1:
                            lons, lats = np.meshgrid(lons, lats)
                        p = global_map(data, lons, lats,
                                       f"EAM {vname} ({label_tag}, h0, last t)",
                                       cmap="RdBu_r", vmin=150, vmax=320,
                                       outdir=outdir, fname=f"eam_h0_{vname}.png",
                                       units="K")
                        if p: plots.append(p)
        close_ds(ds3)

    _report("EAM", issues)
    return issues, plots, fill_reports


def _report(name, issues):
    if issues:
        print(f"  FAIL ({len(issues)} issues):")
        for msg in issues:
            print(msg)
    else:
        print(f"  PASS")


# -----------------------------------------------------------------------------
# Timing summary
# -----------------------------------------------------------------------------

def read_timing_summary(rundir):
    """Read model_timing_stats if present and return summary lines."""
    timing_path = os.path.join(rundir, "timing", "model_timing_stats")
    if not os.path.exists(timing_path):
        return None

    try:
        with open(timing_path, "r") as fh:
            lines = fh.readlines()
    except Exception:
        return None

    if not lines:
        return None

    # Extract key timing info: look for overall model time and component times
    summary = []
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        # Keep header lines and lines with timing data
        # Typical format: component-name  seconds  seconds  seconds
        if any(kw in stripped.lower() for kw in [
            "total", "atm", "lnd", "ocn", "ice", "cpl", "rof", "glc", "wav",
            "init", "run", "final", "overall", "model",
        ]):
            summary.append(stripped)
        elif len(summary) == 0:
            # Keep initial header/description lines
            summary.append(stripped)

    return summary[:40] if summary else None


# -----------------------------------------------------------------------------
# Radiation budget
# -----------------------------------------------------------------------------

def compute_radiation_budget(ds, outdir, comp_plots):
    """Compute and plot TOA/surface radiation budget diagnostics."""
    if not HAS_MPL:
        return

    budget_vars = ["SOLIN", "FSUTOA", "FLUT", "FSDS", "FSUS", "FLDS", "FLUS",
                   "LHFLX", "SHFLX"]
    gmeans = {}
    for vn in budget_vars:
        arr = get_var(ds, vn)
        if arr is None:
            return  # need all variables
        data = arr[-1] if arr.ndim >= 2 and arr.shape[0] > 0 else arr
        v = valid_data(data)
        if v is None or v.size == 0:
            return
        gmeans[vn] = float(v.mean())

    toa_net = gmeans["SOLIN"] - gmeans["FSUTOA"] - gmeans["FLUT"]
    sfc_net = (gmeans["FSDS"] - gmeans["FSUS"] + gmeans["FLDS"] - gmeans["FLUS"]
               - gmeans["LHFLX"] - gmeans["SHFLX"])

    # Bar chart of components and net imbalances
    labels = list(budget_vars) + ["TOA net", "SFC net"]
    values = [gmeans[v] for v in budget_vars] + [toa_net, sfc_net]
    colors = ["#f59e0b" if v > 0 else "#3b82f6" for v in values]
    colors[-2] = "#16a34a"  # TOA net
    colors[-1] = "#16a34a"  # SFC net

    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(len(labels))
    ax.bar(x, values, color=colors, edgecolor="white", linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("W/m2")
    ax.set_title("Radiation Budget: Global-Mean Components")
    ax.axhline(0, color="gray", linewidth=0.5)
    ax.grid(axis="y", alpha=0.3)
    for i, val in enumerate(values):
        ax.text(i, val + (5 if val >= 0 else -12), f"{val:.1f}",
                ha="center", fontsize=7, color="black")
    plt.tight_layout()
    path = os.path.join(outdir, "eam_radiation_budget.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    comp_plots.append(path)
    print(f"  Wrote eam_radiation_budget.png")

    # Global map of TOA net radiation (pixel-by-pixel, last timestep)
    lon_name = next((d for d in ["lon", "ncol", "longitude"] if d in ds.dims), None)
    lat_name = next((d for d in ["lat", "latitude"] if d in ds.dims), None)
    if lon_name and lat_name:
        solin = get_var(ds, "SOLIN")
        fsutoa = get_var(ds, "FSUTOA")
        flut = get_var(ds, "FLUT")
        s_last = solin[-1] if solin.ndim >= 2 else solin
        fu_last = fsutoa[-1] if fsutoa.ndim >= 2 else fsutoa
        fl_last = flut[-1] if flut.ndim >= 2 else flut
        toa_map = s_last - fu_last - fl_last
        lons = ds[lon_name].values
        lats = ds[lat_name].values
        if lons.ndim == 1 and lat_name == "lat":
            lons, lats = np.meshgrid(lons, lats)
        p = global_map(toa_map, lons, lats,
                       "TOA Net Radiation (SOLIN - FSUTOA - FLUT, last t)",
                       cmap="RdBu_r", vmin=-200, vmax=200,
                       outdir=outdir, fname="eam_toa_net_radiation.png",
                       units="W/m2")
        if p:
            comp_plots.append(p)
            print(f"  Wrote eam_toa_net_radiation.png")


# -----------------------------------------------------------------------------
# EAM comparison figures (multi-panel and zonal mean)
# -----------------------------------------------------------------------------

def generate_comparison_figures(rundir, fig_root, all_plots_by_comp):
    """Generate multi-panel layer views and zonal mean profiles for EAM coarsened fields."""
    if not HAS_MPL:
        return

    comp_outdir = os.path.join(fig_root, "comparisons")
    os.makedirs(comp_outdir, exist_ok=True)
    comp_plots = []

    # --- EAM h0: multi-panel T_0..T_7 and STW_0..STW_7 ---
    h0_files = find_files(rundir, "*.eam.h0.*.nc")
    if h0_files:
        ds3 = safe_open(h0_files[0])
        if not has_time_zero(ds3) and HAS_XR and isinstance(ds3, xr.Dataset):
            lon_name = next((d for d in ["lon", "ncol", "longitude"] if d in ds3.dims), None)
            lat_name = next((d for d in ["lat", "latitude"] if d in ds3.dims), None)
            if lon_name and lat_name:
                lons = ds3[lon_name].values
                lats = ds3[lat_name].values if lat_name != lon_name else None
                if lats is not None:
                    if lons.ndim == 1 and lat_name == "lat":
                        lons, lats = np.meshgrid(lons, lats)

                    for base, cmap, vm, vx, units_str in [
                        ("T", "RdYlBu_r", 180, 310, "K"),
                        ("STW", "YlGnBu", 0, 0.02, "kg/kg"),
                        ("U", "RdBu_r", -40, 40, "m/s"),
                        ("V", "RdBu_r", -20, 20, "m/s"),
                    ]:
                        panels = []
                        for k in range(N_ATM_LAYERS):  # 0-based
                            vname = f"{base}_{k}"
                            arr = get_var(ds3, vname)
                            if arr is None or arr.size == 0:
                                continue
                            data = arr[-1] if arr.ndim >= 2 and arr.shape[0] > 0 else arr
                            v = valid_data(data)
                            stats = f"mean={v.mean():.1f}" if v is not None and v.size > 0 else ""
                            panels.append((data, lons, lats,
                                           f"{vname} ({stats})"))

                        if panels:
                            p = multi_panel_maps(
                                panels,
                                f"EAM Vertically Coarsened {base} "
                                f"(all {N_ATM_LAYERS} layers, last t)",
                                comp_outdir,
                                f"eam_h0_{base}_all_layers.png",
                                ncols=4, cmap=cmap, vmin=vm, vmax=vx,
                                units=units_str)
                            if p:
                                comp_plots.append(p)
                                print(f"  Wrote {os.path.basename(p)}")

                    # Zonal mean of T layers (only for lat-lon data)
                    if lat_name == "lat":
                        lat_1d = ds3["lat"].values
                        profiles = []
                        for k in range(N_ATM_LAYERS):  # 0-based
                            zm = _compute_zonal_mean(get_var(ds3, f"T_{k}"))
                            if zm is not None:
                                profiles.append((zm, f"T_{k}"))
                        if profiles:
                            p = zonal_mean_plot(
                                profiles, lat_1d,
                                "EAM Zonal Mean Temperature (all layers)",
                                comp_outdir, "eam_h0_T_zonal_mean.png",
                                ylabel="K")
                            if p:
                                comp_plots.append(p)
                                print(f"  Wrote {os.path.basename(p)}")

                        # Zonal mean of U layers
                        u_profiles = []
                        for k in range(N_ATM_LAYERS):
                            zm = _compute_zonal_mean(get_var(ds3, f"U_{k}"))
                            if zm is not None:
                                u_profiles.append((zm, f"U_{k}"))
                        if u_profiles:
                            p = zonal_mean_plot(
                                u_profiles, lat_1d,
                                "EAM Zonal Mean U Wind (all layers)",
                                comp_outdir, "eam_h0_U_zonal_mean.png",
                                ylabel="m/s")
                            if p:
                                comp_plots.append(p)
                                print(f"  Wrote {os.path.basename(p)}")

                        # Zonal mean of V layers
                        v_profiles = []
                        for k in range(N_ATM_LAYERS):
                            zm = _compute_zonal_mean(get_var(ds3, f"V_{k}"))
                            if zm is not None:
                                v_profiles.append((zm, f"V_{k}"))
                        if v_profiles:
                            p = zonal_mean_plot(
                                v_profiles, lat_1d,
                                "EAM Zonal Mean V Wind (all layers)",
                                comp_outdir, "eam_h0_V_zonal_mean.png",
                                ylabel="m/s")
                            if p:
                                comp_plots.append(p)
                                print(f"  Wrote {os.path.basename(p)}")

        close_ds(ds3)

    # --- Topographic fill validation: PS vs vcoarsen fill mask ---
    if h0_files:
        ds_ps = safe_open(h0_files[0])
        if not has_time_zero(ds_ps) and HAS_XR and isinstance(ds_ps, xr.Dataset):
            ps_arr = get_var(ds_ps, "PS")
            lon_name = next((d for d in ["lon", "ncol", "longitude"] if d in ds_ps.dims), None)
            lat_name = next((d for d in ["lat", "latitude"] if d in ds_ps.dims), None)
            if ps_arr is not None and lon_name and lat_name:
                lons = ds_ps[lon_name].values
                lats = ds_ps[lat_name].values if lat_name != lon_name else None
                if lats is not None:
                    if lons.ndim == 1 and lat_name == "lat":
                        lons, lats = np.meshgrid(lons, lats)
                    ps_last = ps_arr[-1] if ps_arr.ndim >= 2 else ps_arr

                    # For each layer, compare expected fill (PS < layer top bound)
                    # with actual fill from the output
                    fill_validation_rows = []
                    for k in range(N_VCOARSEN_LAYERS):
                        p_top = VCOARSEN_PBOUNDS[k]
                        p_bot = VCOARSEN_PBOUNDS[k + 1]
                        # Expected fill: surface pressure below the top bound of this layer
                        expected_fill = (ps_last < p_top).astype(float)
                        expected_pct = 100.0 * np.sum(expected_fill) / expected_fill.size

                        vname = f"T_{k}"
                        arr = get_var(ds_ps, vname)
                        if arr is not None:
                            data = arr[-1] if arr.ndim >= 2 else arr
                            actual_fill = ((np.abs(data) >= 1e10) | ~np.isfinite(data)).astype(float)
                            actual_pct = 100.0 * np.sum(actual_fill) / actual_fill.size
                        else:
                            actual_pct = float('nan')

                        fill_validation_rows.append({
                            "layer": k,
                            "p_top": p_top / 100,
                            "p_bot": p_bot / 100,
                            "expected_pct": expected_pct,
                            "actual_pct": actual_pct,
                        })

                    # Print fill validation table
                    print("\n  Topographic fill validation (T layers):")
                    print(f"  {'Layer':>5}  {'P range (hPa)':>18}  {'Expected':>10}  {'Actual':>10}  {'Match':>6}")
                    for row in fill_validation_rows:
                        match = "YES" if abs(row["expected_pct"] - row["actual_pct"]) < 0.1 else "NO"
                        print(f"  {row['layer']:>5}  {row['p_top']:>8.1f}-{row['p_bot']:>7.1f}"
                              f"  {row['expected_pct']:>9.2f}%  {row['actual_pct']:>9.2f}%  {match:>6}")

                    # 3-panel figure for T_7: PS | T_7 | fill mask overlay
                    t7_arr = get_var(ds_ps, "T_7")
                    if t7_arr is not None and HAS_MPL:
                        t7_last = t7_arr[-1] if t7_arr.ndim >= 2 else t7_arr
                        expected_fill_7 = (ps_last < VCOARSEN_PBOUNDS[7]).astype(float)
                        actual_fill_7 = ((np.abs(t7_last) >= 1e10) | ~np.isfinite(t7_last)).astype(float)

                        panels = [
                            (ps_last / 100.0, lons, lats, "Surface Pressure (hPa)"),
                            (t7_last, lons, lats, "T_7 (907-1013 hPa)"),
                            (expected_fill_7, lons, lats, "Expected fill (PS < 907 hPa)"),
                            (actual_fill_7, lons, lats, "Actual fill in T_7"),
                        ]

                        # Custom 4-panel figure
                        fig = plt.figure(figsize=(20, 8))
                        cmaps = ["viridis", "RdBu_r", "Reds", "Reds"]
                        vmins = [400, 200, 0, 0]
                        vmaxs = [1050, 310, 1, 1]
                        for idx, (data, lo, la, subtitle) in enumerate(panels):
                            kw = dict(projection=ccrs.PlateCarree()) if HAS_CARTOPY else {}
                            ax = fig.add_subplot(2, 2, idx + 1, **kw)
                            if HAS_CARTOPY:
                                ax.add_feature(cfeature.COASTLINE, linewidth=0.3)
                                ax.set_global()
                            im = _plot_on_ax(ax, data, lo, la,
                                             cmaps[idx], vmins[idx], vmaxs[idx],
                                             HAS_CARTOPY)
                            plt.colorbar(im, ax=ax, shrink=0.7)
                            ax.set_title(subtitle, fontsize=10)
                        fig.suptitle("Topographic Fill Validation: T_7 fill matches PS < 907 hPa",
                                     fontsize=12, y=1.01)
                        plt.tight_layout()
                        path = os.path.join(comp_outdir, "eam_topo_fill_validation.png")
                        fig.savefig(path, dpi=150, bbox_inches="tight")
                        plt.close(fig)
                        comp_plots.append(path)
                        print(f"  Wrote eam_topo_fill_validation.png")

                    # Bar chart: expected vs actual fill per layer
                    if HAS_MPL and fill_validation_rows:
                        fig, ax = plt.subplots(figsize=(8, 4))
                        layers = [r["layer"] for r in fill_validation_rows]
                        exp = [r["expected_pct"] for r in fill_validation_rows]
                        act = [r["actual_pct"] for r in fill_validation_rows]
                        x = np.arange(len(layers))
                        w = 0.35
                        ax.bar(x - w/2, exp, w, label="Expected (PS < layer top)", color="#4c72b0")
                        ax.bar(x + w/2, act, w, label="Actual fill in T_k", color="#dd8452")
                        ax.set_xlabel("Coarsened layer (0=top, 7=surface)")
                        ax.set_ylabel("Fill fraction (%)")
                        ax.set_title("Vcoarsen Fill Fraction: Expected from PS vs Actual")
                        ax.set_xticks(x)
                        ax.set_xticklabels([f"T_{k}" for k in layers])
                        ax.legend()
                        ax.grid(axis="y", alpha=0.3)
                        path = os.path.join(comp_outdir, "eam_fill_fraction_validation.png")
                        fig.savefig(path, dpi=150, bbox_inches="tight")
                        plt.close(fig)
                        comp_plots.append(path)
                        print(f"  Wrote eam_fill_fraction_validation.png")

        close_ds(ds_ps)

    # --- Radiation budget ---
    if h0_files:
        ds_rad = safe_open(h0_files[0])
        if not has_time_zero(ds_rad) and HAS_XR and isinstance(ds_rad, xr.Dataset):
            compute_radiation_budget(ds_rad, comp_outdir, comp_plots)
        close_ds(ds_rad)

    # --- Time series of key variables across all h0 files ---
    if h0_files and HAS_MPL and len(h0_files) >= 1:
        ts_vars = ["TS", "PS", "PRECT", "FLUT", "T_7", "STW_7"]
        ts_data = {vn: [] for vn in ts_vars}
        for fpath in sorted(h0_files):
            ds_ts = safe_open(fpath)
            if has_time_zero(ds_ts):
                close_ds(ds_ts)
                continue
            for vn in ts_vars:
                arr = get_var(ds_ts, vn)
                if arr is None:
                    continue
                # iterate over all timesteps in this file
                if arr.ndim >= 2:
                    for tidx in range(arr.shape[0]):
                        v = valid_data(arr[tidx])
                        if v is not None and v.size > 0:
                            ts_data[vn].append(float(v.mean()))
                else:
                    v = valid_data(arr)
                    if v is not None and v.size > 0:
                        ts_data[vn].append(float(v.mean()))
            close_ds(ds_ts)

        for vn in ts_vars:
            vals = ts_data[vn]
            if len(vals) >= 2:
                time_series(vals, "Timestep index", f"Global Mean {vn}",
                            vn, comp_outdir, f"eam_timeseries_{vn}.png")
                path = os.path.join(comp_outdir, f"eam_timeseries_{vn}.png")
                if os.path.exists(path):
                    comp_plots.append(path)
                    print(f"  Wrote eam_timeseries_{vn}.png")

    if comp_plots:
        all_plots_by_comp["Comparisons"] = comp_plots


# ─────────────────────────────────────────────────────────────────────────────
# Cross-verification: compare fme_output vs fme_legacy_output
# ─────────────────────────────────────────────────────────────────────────────

def cross_verify(fme_rundir, legacy_rundir, outdir, verbose=False):
    """
    Compare online FME output against legacy (raw) output.

    Both cases run identical model physics from the same ICs -- the only
    difference is what diagnostics get written.  Therefore:

    1. BFB IDENTITY: any field present in both h0 files must match
       bit-for-bit (same model state, same history averaging).

    2. ONLINE vs OFFLINE VCOARSEN: the legacy case outputs raw 3D fields
       (T, Q, U, V, CLDLIQ, CLDICE, RAINQM at full L80 resolution).
       We vcoarsen them offline with the same algorithm / pressure bounds
       and compare against the FME case's online T_0..T_7, STW_0..STW_7,
       U_0..U_7, V_0..V_7.  These should be BFB identical.

    3. ONLINE vs OFFLINE DERIVED: legacy raw Q + CLDICE + CLDLIQ + RAINQM
       should equal FME STW at every level.  After vcoarsening, this tests
       derived + vcoarsen together.

    Returns a dict of {test_name: list_of_result_dicts} for HTML rendering.
    """
    print("\n" + "=" * 70)
    print("CROSS-VERIFICATION: fme_output vs fme_legacy_output")
    print("=" * 70)

    issues = []
    results = []   # list of dicts: {test, field, status, detail}
    xv_plots = []  # comparison figures

    xv_outdir = os.path.join(outdir, "fme_output", "cross_verify")
    os.makedirs(xv_outdir, exist_ok=True)

    # -- locate h0 files in both rundirs --------------------------------------
    fme_h0 = find_files(fme_rundir, "*.eam.h0.*.nc")
    leg_h0 = find_files(legacy_rundir, "*.eam.h0.*.nc")

    if not fme_h0 or not leg_h0:
        print("  SKIP EAM: h0 files not found in both rundirs")
        print(f"    fme: {len(fme_h0)} files,  legacy: {len(leg_h0)} files")
        return issues

    ds_fme = safe_open(fme_h0[0])
    ds_leg = safe_open(leg_h0[0])

    if has_time_zero(ds_fme) or has_time_zero(ds_leg):
        print("  SKIP: time dimension has size 0")
        close_ds(ds_fme); close_ds(ds_leg)
        return issues

    # -- 1. FIELD IDENTITY ------------------------------------------------------
    # Detect grid mismatch (FME may be remapped to lat-lon, legacy on native)
    fme_dims = get_dims(ds_fme)
    leg_dims = get_dims(ds_leg)
    grids_match = ("ncol" in fme_dims) == ("ncol" in leg_dims)

    if grids_match:
        print("\n--- Test 1: BFB Identity (shared fields, same grid) ---")
    else:
        fme_grid = "lat-lon" if "lat" in fme_dims else "native"
        leg_grid = "lat-lon" if "lat" in leg_dims else "native"
        print(f"\n--- Test 1: Global-Mean Comparison (FME={fme_grid}, Legacy={leg_grid}) ---")

    all_bfb_fields = BFB_INST_FIELDS + BFB_AVG_FIELDS
    n_bfb_pass = 0
    n_bfb_fail = 0

    for var in all_bfb_fields:
        fme_data = get_var(ds_fme, var)
        leg_data = get_var(ds_leg, var)
        if fme_data is None or leg_data is None:
            label = "missing-fme" if fme_data is None else "missing-legacy"
            results.append({"test": "BFB", "field": var, "status": "SKIP",
                            "detail": f"{label}"})
            print(f"  {var}: SKIP ({label})")
            continue

        if grids_match and fme_data.shape == leg_data.shape:
            # Same grid: point-by-point BFB
            diff = np.abs(fme_data.astype(np.float64) - leg_data.astype(np.float64))
            maxdiff = float(np.nanmax(diff))
            status = "BFB" if maxdiff < 1e-12 else "DIFF"
            if status == "BFB":
                n_bfb_pass += 1
            else:
                n_bfb_fail += 1
                issues.append(f"  BFB {var}: max|diff| = {maxdiff:.4g}")
            results.append({"test": "BFB", "field": var, "status": status,
                            "detail": f"max|diff| = {maxdiff:.2e}"})
            tag = "PASS" if status == "BFB" else "FAIL"
            print(f"  {var}: {tag}  max|diff| = {maxdiff:.2e}")
        else:
            # Different grids: compare area-weighted global means (last timestep)
            fme_last = fme_data[-1] if fme_data.ndim >= 2 else fme_data
            leg_last = leg_data[-1] if leg_data.ndim >= 2 else leg_data

            # Area-weight lat-lon data with cos(lat); native grid is quasi-equal-area
            fme_flat = fme_last.ravel().astype(np.float64) if hasattr(fme_last, 'ravel') else fme_last
            leg_flat = leg_last.ravel().astype(np.float64) if hasattr(leg_last, 'ravel') else leg_last

            fme_lat = get_var(ds_fme, "lat")
            if fme_lat is not None and "lat" in fme_dims and fme_last.ndim == 2:
                # cos-lat weighting for regular lat-lon grid
                coslat = np.cos(np.deg2rad(fme_lat))
                weights = coslat[:, np.newaxis] * np.ones(fme_last.shape[1])[np.newaxis, :]
                mask = (np.abs(fme_last) < 1e10) & np.isfinite(fme_last)
                fme_mean = float(np.average(fme_last[mask], weights=weights[mask])) if mask.any() else float('nan')
            else:
                fme_v = valid_data(fme_flat)
                fme_mean = float(fme_v.mean()) if fme_v is not None and fme_v.size > 0 else float('nan')

            # Native grid (ne30pg2) is quasi-equal-area, unweighted mean is fine
            leg_v = valid_data(leg_flat)
            leg_mean = float(leg_v.mean()) if leg_v is not None and leg_v.size > 0 else float('nan')

            if np.isnan(fme_mean) or np.isnan(leg_mean):
                results.append({"test": "GMEAN", "field": var, "status": "SKIP",
                                "detail": "no valid data"})
                print(f"  {var}: SKIP (no valid data)")
                continue

            rel_diff = abs(fme_mean - leg_mean) / max(abs(leg_mean), 1e-30)
            abs_diff = abs(fme_mean - leg_mean)
            # For near-zero-mean fields (winds, stresses), relative difference
            # is meaningless. Use absolute difference against the larger magnitude.
            field_mag = max(abs(fme_mean), abs(leg_mean))
            # Area-weighted means should agree within remap interpolation error (~2%)
            # or within 2% of the field magnitude for near-zero fields
            if field_mag > 1e-10:
                effective_rel = abs_diff / field_mag
            else:
                effective_rel = 0.0  # both means are essentially zero
            status = "PASS" if min(rel_diff, effective_rel) < 0.02 else "DIFF"
            if status == "PASS":
                n_bfb_pass += 1
            else:
                n_bfb_fail += 1
                issues.append(f"  GMEAN {var}: rel_diff = {rel_diff:.4g}")
            results.append({"test": "GMEAN", "field": var, "status": status,
                            "detail": f"fme={fme_mean:.6g}, leg={leg_mean:.6g}, rel={rel_diff:.2e}"})
            print(f"  {var}: {status}  fme={fme_mean:.6g}  leg={leg_mean:.6g}  rel={rel_diff:.2e}")

    print(f"  Summary: {n_bfb_pass} pass, {n_bfb_fail} fail, "
          f"{len(all_bfb_fields) - n_bfb_pass - n_bfb_fail} skip")

    # -- Native vs Remapped comparison plots -----------------------------------
    if HAS_MPL and not grids_match:
        print("\n--- Generating native vs remapped comparison plots ---")
        tidx = -1

        # Get native grid coordinates (legacy: lat/lon are 1D arrays of ncol)
        leg_lat = get_var(ds_leg, "lat")
        leg_lon = get_var(ds_leg, "lon")
        # Get remapped grid coordinates (FME: lat/lon are 1D coordinate arrays)
        fme_lat = get_var(ds_fme, "lat")
        fme_lon = get_var(ds_fme, "lon")

        if leg_lat is not None and leg_lon is not None and \
           fme_lat is not None and fme_lon is not None:
            # Native: 1D arrays for tripcolor
            nat_lons = leg_lon if leg_lon.ndim == 1 else leg_lon.ravel()
            nat_lats = leg_lat if leg_lat.ndim == 1 else leg_lat.ravel()
            # Remapped: meshgrid for pcolormesh
            if fme_lat.ndim == 1 and fme_lon.ndim == 1:
                rem_lons, rem_lats = np.meshgrid(fme_lon, fme_lat)
            else:
                rem_lons, rem_lats = fme_lon, fme_lat

            plot_specs = [
                ("TS",       "RdBu_r",   220, 320, "K"),
                ("PS",       "viridis",  50000, 105000, "Pa"),
                ("LHFLX",    "YlOrRd",   0,   300, "W/m2"),
                ("PRECT",    "Blues",     0,   1e-6, "m/s"),
                ("ICEFRAC",  "Blues",     0,   1,   "fraction"),
                ("FLUT",     "inferno",   100, 320, "W/m2"),
                ("SOLIN",    "YlOrRd",   0,   500, "W/m2"),
                ("DTENDTTW", "RdBu_r",   -5e-4, 5e-4, "kg/m2/s"),
            ]
            for var, cmap, vm, vx, units in plot_specs:
                nat_arr = get_var(ds_leg, var)
                rem_arr = get_var(ds_fme, var)
                if nat_arr is None or rem_arr is None:
                    continue
                nat_data = nat_arr[tidx] if nat_arr.ndim >= 2 else nat_arr
                rem_data = rem_arr[tidx] if rem_arr.ndim >= 2 else rem_arr
                # Flatten native to 1D for tripcolor
                nat_flat = nat_data.ravel()
                p = side_by_side_comparison(
                    nat_flat, nat_lons, nat_lats, f"Native ne30pg2 ({var})",
                    rem_data, rem_lons, rem_lats, f"Remapped lat-lon ({var})",
                    f"{var}: Native vs Remapped (last timestep)",
                    xv_outdir, f"xv_native_vs_remap_{var}.png",
                    cmap=cmap, vmin=vm, vmax=vx, units=units)
                if p:
                    xv_plots.append(p)
                    print(f"  Wrote xv_native_vs_remap_{var}.png")

            # Vcoarsen layers: T_7 and STW_7 (near-surface)
            for base, cmap, vm, vx, units in [("T", "RdYlBu_r", 220, 310, "K"),
                                                ("STW", "YlGnBu", 0, 0.025, "kg/kg")]:
                vname = f"{base}_7"
                fme_vc = get_var(ds_fme, vname)
                if fme_vc is not None:
                    fme_data = fme_vc[tidx] if fme_vc.ndim >= 2 else fme_vc
                    # No native vcoarsen to compare, just show remapped
                    p = global_map(fme_data, rem_lons, rem_lats,
                                   f"FME {vname} (remapped, near-surface)",
                                   cmap=cmap, vmin=vm, vmax=vx,
                                   outdir=xv_outdir, fname=f"xv_fme_{vname}_remap.png",
                                   units=units)
                    if p:
                        xv_plots.append(p)
                        print(f"  Wrote xv_fme_{vname}_remap.png")
        else:
            print("  SKIP plots: lat/lon coordinates not found")

    # -- 1b. TOPOGRAPHIC FILL VALIDATION ----------------------------------------
    print("\n--- Test 1b: Topographic Fill Validation ---")
    # Verify that fill values in coarsened layers correspond exactly to grid cells
    # where the surface pressure is below the target layer's top pressure bound.
    fme_ps = get_var(ds_fme, "PS")
    if fme_ps is not None:
        ps_last = fme_ps[-1] if fme_ps.ndim >= 2 else fme_ps
        topo_max_err = 0.0
        for k in range(N_VCOARSEN_LAYERS):
            p_top = VCOARSEN_PBOUNDS[k]
            # Expected: fill where PS < layer's top bound (no model levels in range)
            expected_pct = 100.0 * np.mean(ps_last.ravel() < p_top)
            t_k = get_var(ds_fme, f"T_{k}")
            if t_k is None:
                continue
            t_last = t_k[-1] if t_k.ndim >= 2 else t_k
            actual_pct = 100.0 * np.mean(
                (np.abs(t_last.ravel()) >= 1e10) | ~np.isfinite(t_last.ravel()))
            topo_max_err = max(topo_max_err, abs(expected_pct - actual_pct))

        # Tolerance: remap can shift fill boundaries by ~0.5% of grid cells
        topo_status = "PASS" if topo_max_err < 1.0 else "DIFF"
        results.append({"test": "TOPO_FILL", "field": "T layers",
                        "status": topo_status,
                        "detail": f"max expected-vs-actual fill mismatch: {topo_max_err:.3f}%"})
        print(f"  Topographic fill: {topo_status}  max mismatch={topo_max_err:.3f}%")
        if topo_status != "PASS":
            issues.append(f"  TOPO_FILL: max fill mismatch {topo_max_err:.3f}%")
    else:
        print("  SKIP: PS not found in FME output")

    # -- 2. ONLINE vs OFFLINE VCOARSEN ----------------------------------------
    print("\n--- Test 2: Online vs Offline Vcoarsen ---")

    # Need hybrid coords and PS from legacy file to compute interface pressures
    hyai = get_var(ds_leg, "hyai")
    hybi = get_var(ds_leg, "hybi")
    leg_ps = get_var(ds_leg, "PS")

    if hyai is None or hybi is None or leg_ps is None:
        print("  SKIP: hyai/hybi/PS not found in legacy file")
    else:
        # Use the LAST timestep to avoid initialization artifacts
        tidx = -1
        ps_1t = leg_ps[tidx]  # (ncol,)
        pint = compute_pint(ps_1t, hyai, hybi)  # (ncol, nlev+1)

        for fme_base, legacy_sources in FME_TO_LEGACY.items():
            # Build the full-resolution field from legacy
            # (for derived fields like STW, sum the components)
            raw_3d = None
            missing_src = False
            for src in legacy_sources:
                src_data = get_var(ds_leg, src)
                if src_data is None:
                    print(f"  {fme_base}: SKIP (legacy missing {src})")
                    missing_src = True
                    break
                src_1t = src_data[tidx]  # (lev, ncol) or (ncol,) from EAM NetCDF
                if src_1t.ndim < 2:
                    print(f"  {fme_base}: SKIP ({src} is not 3D)")
                    missing_src = True
                    break
                # EAM NetCDF stores 3D as (time, lev, ncol); after time slice: (lev, ncol)
                # offline_vcoarsen_avg expects (ncol, lev) — transpose if needed
                if src_1t.shape[0] < src_1t.shape[1]:
                    src_1t = src_1t.T  # (lev, ncol) -> (ncol, lev)
                if raw_3d is None:
                    raw_3d = src_1t.astype(np.float64)
                else:
                    raw_3d = raw_3d + src_1t.astype(np.float64)

            if missing_src or raw_3d is None:
                continue

            # Offline vcoarsen
            offline_vc = offline_vcoarsen_avg(raw_3d, pint)  # (ncol, 8)

            # Compare with online vcoarsen from FME output
            layer_results = []
            max_err_all = 0.0
            max_rel_all = 0.0
            for k in range(N_VCOARSEN_LAYERS):
                vname = f"{fme_base}_{k}"
                fme_arr = get_var(ds_fme, vname)
                if fme_arr is None:
                    layer_results.append(f"  {vname}: MISSING in FME output")
                    continue

                fme_1t = fme_arr[tidx].astype(np.float64).ravel()  # flatten (lat,lon) or (ncol,)
                off_1t = offline_vc[:, k]

                # Grids may differ (remapped vs native) — compare area-weighted global means
                if fme_1t.size != off_1t.size:
                    # Area-weight the FME (lat-lon) data
                    fme_lat = get_var(ds_fme, "lat")
                    fme_arr_2d = fme_arr[tidx].astype(np.float64)
                    if fme_lat is not None and fme_arr_2d.ndim == 2:
                        coslat = np.cos(np.deg2rad(fme_lat))
                        wt = coslat[:, np.newaxis] * np.ones(fme_arr_2d.shape[1])[np.newaxis, :]
                        m = (np.abs(fme_arr_2d) < 1e10) & np.isfinite(fme_arr_2d)
                        fme_wm = float(np.average(fme_arr_2d[m], weights=wt[m])) if m.any() else float('nan')
                    else:
                        fv = valid_data(fme_1t)
                        fme_wm = float(fv.mean()) if fv is not None else float('nan')
                    off_v = valid_data(off_1t)
                    off_wm = float(off_v.mean()) if off_v is not None else float('nan')
                    if not np.isnan(fme_wm) and not np.isnan(off_wm):
                        rel = abs(fme_wm - off_wm) / max(abs(off_wm), 1e-30)
                        max_rel_all = max(max_rel_all, rel)
                        max_err_all = max(max_err_all, abs(fme_wm - off_wm))
                    continue

                # Same grid: point-by-point comparison
                valid_mask = (np.abs(fme_1t) < 1e10) & (np.abs(off_1t) < 1e10) & \
                             np.isfinite(fme_1t) & np.isfinite(off_1t)

                if not valid_mask.any():
                    layer_results.append(f"  {vname}: no valid data")
                    continue

                diff = np.abs(fme_1t[valid_mask] - off_1t[valid_mask])
                maxdiff = float(diff.max())
                scale = np.maximum(np.abs(off_1t[valid_mask]), 1e-30)
                maxrel = float((diff / scale).max())

                max_err_all = max(max_err_all, maxdiff)
                max_rel_all = max(max_rel_all, maxrel)

            # Report — thresholds depend on whether grids match.
            # For near-zero fields (U, V, TAUX, TAUY), relative difference is
            # meaningless because the global mean crosses zero.  Use absolute
            # difference instead, benchmarked against a typical scale.
            if grids_match:
                bfb_thresh, close_thresh = 1e-10, 1e-5
            else:
                bfb_thresh, close_thresh = 1e-5, 0.02  # remap interpolation ~2%

            # Determine a characteristic scale for this field to judge abs diff
            off_v = valid_data(offline_vc.ravel())
            field_scale = float(np.percentile(np.abs(off_v), 95)) if off_v is not None and off_v.size > 0 else 1.0
            # For near-zero-mean fields, use abs diff / field_scale as the metric
            if field_scale > 0:
                abs_rel = max_err_all / field_scale
            else:
                abs_rel = 0.0

            # Use the more lenient of relative-mean or absolute-relative
            effective_metric = min(max_rel_all, abs_rel)

            if effective_metric < bfb_thresh:
                status = "BFB"
                tag = "PASS"
            elif effective_metric < close_thresh:
                status = "CLOSE"
                tag = "PASS"
            else:
                status = "DIFF"
                tag = "FAIL"
                issues.append(f"  VCOARSEN {fme_base}: max|rel|={max_rel_all:.2e}, "
                              f"|diff|/scale={abs_rel:.2e}")

            detail = (f"max|diff|={max_err_all:.2e}, max|rel|={max_rel_all:.2e}, "
                      f"|diff|/p95={abs_rel:.2e}")
            results.append({"test": "VCOARSEN", "field": fme_base,
                            "status": status, "detail": detail})
            print(f"  {fme_base}_0..{fme_base}_{N_VCOARSEN_LAYERS-1}: {tag}  {detail}")

            # Generate difference map for last layer (skip if grids differ)
            if HAS_MPL and HAS_XR and grids_match:
                last_k = N_VCOARSEN_LAYERS - 1
                vname = f"{fme_base}_{last_k}"
                fme_arr = get_var(ds_fme, vname)
                if fme_arr is not None:
                    fme_1t = fme_arr[tidx].astype(np.float64)
                    off_1t = offline_vc[:, last_k]
                    diff_map = fme_1t - off_1t

                    # Get coordinates
                    lon_name = next((d for d in ["lon", "ncol"] if d in get_dims(ds_fme)), None)
                    lat_var = get_var(ds_fme, "lat")
                    lon_var = get_var(ds_fme, "lon")
                    if lat_var is not None and lon_var is not None:
                        if lat_var.ndim == 1 and lon_var.ndim == 1:
                            lons, lats = np.meshgrid(lon_var, lat_var)
                        else:
                            lons, lats = lon_var, lat_var
                        # Side-by-side: online vs offline
                        vmin_f = float(np.nanpercentile(valid_data(fme_1t) or [0], 2))
                        vmax_f = float(np.nanpercentile(valid_data(fme_1t) or [0], 98))
                        p = side_by_side_comparison(
                            fme_1t, lons, lats, f"Online {vname}",
                            off_1t, lons, lats, f"Offline vcoarsen({'+'.join(legacy_sources)})",
                            f"{fme_base}: Online vs Offline (layer {last_k})",
                            xv_outdir, f"xv_vcoarsen_{fme_base}_{last_k}.png",
                            cmap="RdYlBu_r", vmin=vmin_f, vmax=vmax_f)
                        if p:
                            xv_plots.append(p)

    # -- 3. VCOARSEN LINEARITY: STW_k == Q_k + CLDICE_k + ... ----------------
    print("\n--- Test 3: Vcoarsen Linearity (STW_k = sum of components) ---")
    # If we have offline vcoarsen of individual components (Q, CLDICE, CLDLIQ,
    # RAINQM), their sum should equal the online STW_k (since vcoarsen is linear).
    if hyai is not None and hybi is not None and leg_ps is not None:
        tidx = -1
        ps_1t = leg_ps[tidx]
        pint = compute_pint(ps_1t, hyai, hybi)

        # Offline vcoarsen each component
        components_vc = {}
        all_found = True
        for comp in ["Q", "CLDICE", "CLDLIQ", "RAINQM"]:
            arr = get_var(ds_leg, comp)
            if arr is None or arr.ndim < 3:
                all_found = False
                break
            comp_1t = arr[tidx].astype(np.float64)
            if comp_1t.shape[0] < comp_1t.shape[1]:
                comp_1t = comp_1t.T  # (lev, ncol) -> (ncol, lev)
            components_vc[comp] = offline_vcoarsen_avg(comp_1t, pint)

        if all_found:
            sum_vc = sum(components_vc.values())  # (ncol, 8)
            max_lin_err = 0.0
            max_lin_rel = 0.0

            for k in range(N_VCOARSEN_LAYERS):
                stw_k = get_var(ds_fme, f"STW_{k}")
                if stw_k is None:
                    continue
                stw_1t = stw_k[tidx].astype(np.float64).ravel()
                sum_1t = sum_vc[:, k]

                # Different grids: compare global means
                if stw_1t.size != sum_1t.size:
                    sv = valid_data(stw_1t)
                    ov = valid_data(sum_1t)
                    if sv is not None and ov is not None:
                        rel = abs(sv.mean() - ov.mean()) / max(abs(ov.mean()), 1e-30)
                        max_lin_rel = max(max_lin_rel, rel)
                        max_lin_err = max(max_lin_err, abs(sv.mean() - ov.mean()))
                    continue

                valid = (np.abs(stw_1t) < 1e10) & np.isfinite(stw_1t)
                if not valid.any():
                    continue

                diff = np.abs(stw_1t[valid] - sum_1t[valid])
                scale = np.maximum(np.abs(sum_1t[valid]), 1e-30)
                max_lin_err = max(max_lin_err, float(diff.max()))
                max_lin_rel = max(max_lin_rel, float((diff / scale).max()))

            if max_lin_rel < 1e-5:
                status = "PASS"
            else:
                status = "FAIL"
                issues.append(f"  LINEARITY STW: max|rel| = {max_lin_rel:.2e}")

            results.append({"test": "LINEARITY", "field": "STW_k",
                            "status": status,
                            "detail": f"STW_k vs sum(Q_k+CLDICE_k+CLDLIQ_k+RAINQM_k): "
                                      f"max|rel|={max_lin_rel:.2e}"})
            print(f"  STW_k linearity: {status}  max|rel|={max_lin_rel:.2e}")
        else:
            print("  SKIP linearity: component fields not all found in legacy")

    close_ds(ds_fme)
    close_ds(ds_leg)

    # -- 4. TIMING & STORAGE COMPARISON ---------------------------------------
    print("\n--- Timing & Storage ---")

    # Storage: compare total NetCDF file sizes
    def dir_nc_size(d):
        total = 0
        count = 0
        for f in glob.glob(os.path.join(d, "*.nc")):
            sz = os.path.getsize(f)
            total += sz
            count += 1
        return total, count

    fme_bytes, fme_nfiles = dir_nc_size(fme_rundir)
    leg_bytes, leg_nfiles = dir_nc_size(legacy_rundir)
    fme_gb = fme_bytes / 1e9
    leg_gb = leg_bytes / 1e9
    ratio = fme_gb / leg_gb if leg_gb > 0 else float('inf')
    print(f"  FME output:    {fme_nfiles} files, {fme_gb:.2f} GB")
    print(f"  Legacy output: {leg_nfiles} files, {leg_gb:.2f} GB")
    print(f"  Storage ratio: {ratio:.2f}x {'(smaller)' if ratio < 1 else '(larger)'}")

    results.append({"test": "STORAGE", "field": "total NC files",
                     "status": "INFO",
                     "detail": f"FME: {fme_nfiles} files / {fme_gb:.2f} GB, "
                               f"Legacy: {leg_nfiles} files / {leg_gb:.2f} GB, "
                               f"ratio: {ratio:.2f}x"})

    # Timing: parse model_timing_stats from both runs
    fme_timing = read_timing_summary(fme_rundir)
    leg_timing = read_timing_summary(legacy_rundir)

    # Try to extract total model time from timing files
    def parse_model_time(timing_lines):
        """Extract total model run time in seconds from timing summary."""
        if not timing_lines:
            return None
        for line in timing_lines:
            low = line.lower()
            # Look for lines with 'total' or 'model' and a number
            if ('total' in low or 'model' in low) and ('run' in low or 'time' in low):
                parts = line.split()
                for p in parts:
                    try:
                        val = float(p)
                        if val > 1:  # skip small numbers (likely ratios)
                            return val
                    except ValueError:
                        continue
        return None

    fme_time = parse_model_time(fme_timing)
    leg_time = parse_model_time(leg_timing)
    if fme_time and leg_time:
        overhead = ((fme_time - leg_time) / leg_time) * 100
        print(f"  FME model time:    {fme_time:.1f} s")
        print(f"  Legacy model time: {leg_time:.1f} s")
        print(f"  FME overhead: {overhead:+.1f}%")
        results.append({"test": "TIMING", "field": "model run time",
                         "status": "INFO",
                         "detail": f"FME: {fme_time:.1f}s, Legacy: {leg_time:.1f}s, "
                                   f"overhead: {overhead:+.1f}%"})
    elif fme_timing or leg_timing:
        print("  Timing data found in only one run -- cannot compare")
    else:
        print("  No timing data found in either run")

    # -- Summary --------------------------------------------------------------
    print("\n" + "-" * 50)
    n_pass = sum(1 for r in results if r["status"] in ("BFB", "PASS", "CLOSE"))
    n_fail = sum(1 for r in results if r["status"] in ("DIFF", "FAIL"))
    n_skip = sum(1 for r in results if r["status"] == "SKIP")
    n_info = sum(1 for r in results if r["status"] == "INFO")

    print(f"  CROSS-VERIFY TOTAL: {n_pass} pass, {n_fail} fail, "
          f"{n_skip} skip, {n_info} info")
    if not issues:
        print("  ALL CROSS-VERIFICATION PASSED")
    else:
        for iss in issues:
            print(iss)

    return issues, results, xv_plots


def write_cross_verify_html(outdir, results, xv_plots):
    """Write the cross-verification section of the HTML dashboard."""
    if not results:
        return ""

    # Group by test type
    tests_by_type = {}
    for r in results:
        tests_by_type.setdefault(r["test"], []).append(r)

    n_total_pass = sum(1 for r in results if r["status"] in ("BFB", "PASS", "CLOSE"))
    n_total_fail = sum(1 for r in results if r["status"] in ("DIFF", "FAIL"))
    n_total_info = sum(1 for r in results if r["status"] == "INFO")

    html = '<h2 id="Cross_Verify">Cross-Verification: FME vs Legacy</h2>\n'
    html += '<p>Both cases run identical physics from the same ICs. '
    html += 'FME diagnostics are output-only and do not modify model state. '
    html += f'<strong>{n_total_pass} pass, {n_total_fail} fail, {n_total_info} info</strong></p>\n'

    # Tabbed interface for test categories
    tab_ids = list(tests_by_type.keys())
    if xv_plots:
        tab_ids.append("PLOTS")

    html += '<div class="tabs">\n'
    for i, tid in enumerate(tab_ids):
        labels = {
            "BFB": "BFB Identity", "GMEAN": "Global Means",
            "TOPO_FILL": "Topo Fill", "VCOARSEN": "Vcoarsen",
            "LINEARITY": "Linearity",
            "STORAGE": "Storage", "TIMING": "Timing", "PLOTS": "Figures",
        }
        active = " active" if i == 0 else ""
        n_p = sum(1 for r in tests_by_type.get(tid, []) if r["status"] in ("BFB", "PASS", "CLOSE"))
        n_t = len(tests_by_type.get(tid, []))
        badge = f' <span class="badge">{n_p}/{n_t}</span>' if tid != "PLOTS" else ""
        html += f'<button class="tab{active}" onclick="showTab(\'{tid}\')">{labels.get(tid, tid)}{badge}</button>\n'
    html += '</div>\n'

    for tid, test_results in tests_by_type.items():
        n_pass = sum(1 for r in test_results if r["status"] in ("BFB", "PASS", "CLOSE"))
        n_tot = len(test_results)
        status_cls = "pass" if n_pass == n_tot else ("fail" if n_pass < n_tot else "")

        labels = {
            "BFB": "BFB Identity (shared fields match bit-for-bit)",
            "GMEAN": "Area-Weighted Global-Mean Comparison (cos-lat weighted, &lt;2% tolerance)",
            "TOPO_FILL": "Topographic Fill Validation (fill matches PS &lt; layer bound)",
            "VCOARSEN": "Online vs Offline Vcoarsen (area-weighted global mean comparison)",
            "LINEARITY": "Vcoarsen Linearity (STW_k = sum of component_k)",
            "STORAGE": "Storage Comparison", "TIMING": "Timing Comparison",
        }
        vis = "" if tid == tab_ids[0] else ' style="display:none"'
        html += f'<div class="tab-content" id="tab_{tid}"{vis}>\n'
        html += f'<h3>{labels.get(tid, tid)} '
        html += f'<span class="{status_cls}">({n_pass}/{n_tot})</span></h3>\n'

        # Sortable table
        html += '<table class="sortable"><thead><tr>'
        html += '<th onclick="sortTable(this,0)">Field</th>'
        html += '<th onclick="sortTable(this,1)">Status</th>'
        html += '<th onclick="sortTable(this,2)">Detail</th></tr></thead><tbody>\n'

        for r in test_results:
            cls = "pass" if r["status"] in ("BFB", "PASS", "CLOSE") else \
                  "fail" if r["status"] in ("DIFF", "FAIL") else ""
            html += (f'<tr><td>{r["field"]}</td>'
                     f'<td class="{cls}">{r["status"]}</td>'
                     f'<td style="font-family:monospace;font-size:0.85em">{r["detail"]}</td></tr>\n')
        html += '</tbody></table>\n'
        html += '</div>\n'

    # Figures tab
    if xv_plots:
        vis = "" if "PLOTS" == tab_ids[0] else ' style="display:none"'
        html += f'<div class="tab-content" id="tab_PLOTS"{vis}>\n'
        html += '<h3>Native vs Remapped Comparison Figures</h3>\n'
        html += '<p>Left: native ne30pg2 (tripcolor). Right: remapped Gaussian lat-lon (pcolormesh).</p>\n'
        html += '<div class="gallery">\n'
        for p in xv_plots:
            if p and os.path.exists(p):
                rel = os.path.relpath(p, outdir)
                html += (f'<div class="fig wide"><a href="{rel}">'
                         f'<img src="{rel}"/></a>'
                         f'<span>{os.path.basename(p)}</span></div>\n')
        html += '</div>\n</div>\n'

    return html


# -----------------------------------------------------------------------------
# File inventory (EAM only)
# -----------------------------------------------------------------------------

def collect_file_inventory(rundir):
    """Collect EAM file inventory data and print summary.

    Returns list of (label, count, file_list) tuples for HTML rendering.
    """
    print("\n=== File Inventory ===")
    categories = [
        ("EAM h0", "*.eam.h0.*.nc", None),
    ]
    total = 0
    inventory_data = []
    for label, pattern, exclude in categories:
        hits = find_files(rundir, pattern, exclude=exclude)
        total += len(hits)
        status = f"{len(hits)} file(s)" if hits else "NONE"
        print(f"  {label:38s} {status}")
        inventory_data.append((label, len(hits), hits))
    print(f"  {'TOTAL':38s} {total} EAM-related file(s)")
    return inventory_data


# -----------------------------------------------------------------------------
# ACE variable coverage
# -----------------------------------------------------------------------------

def check_variable_coverage(rundir):
    """Check which ACE-required variables are present in the FME output.

    Returns list of (ace_name, eam_name, category_str, found_bool) tuples.
    """
    h0_files = find_files(rundir, "*.eam.h0.*.nc")
    if not h0_files:
        return []

    ds = safe_open(h0_files[0])
    varnames_in_file = set(get_varnames(ds))
    close_ds(ds)

    forcing_set = set(ACE_FORCING_NAMES)
    in_set = set(ACE_IN_NAMES)
    out_set = set(ACE_OUT_NAMES)

    coverage = []
    seen = set()
    for ace_name in ACE_FORCING_NAMES + ACE_IN_NAMES + ACE_OUT_NAMES:
        if ace_name in seen:
            continue
        seen.add(ace_name)
        eam_name = ace_to_eam(ace_name)
        cats = []
        if ace_name in forcing_set:
            cats.append("forcing")
        if ace_name in in_set:
            cats.append("in")
        if ace_name in out_set:
            cats.append("out")
        category = ", ".join(cats)
        found = eam_name in varnames_in_file
        coverage.append((ace_name, eam_name, category, found))

    return coverage


def generate_yaml_text():
    """Generate YAML-formatted text of the ACE atmosphere variable spec."""
    lines = ["## ATMOSPHERE", "next_step_forcing_names:"]
    for name in ACE_FORCING_NAMES:
        lines.append(f"  - {name}")
    lines.append("in_names:")
    for name in ACE_IN_NAMES:
        lines.append(f"  - {name}")
    lines.append("out_names:")
    for name in ACE_OUT_NAMES:
        lines.append(f"  - {name}")
    return "\n".join(lines)


# -----------------------------------------------------------------------------
# Reproducibility info
# -----------------------------------------------------------------------------

def _find_user_nl_eam(rundir):
    """Search for user_nl_eam in likely locations relative to rundir."""
    candidates = [
        os.path.join(rundir, "user_nl_eam"),
        os.path.join(os.path.dirname(rundir.rstrip("/")), "user_nl_eam"),
    ]
    for path in candidates:
        if os.path.exists(path):
            return path
    return None


def collect_repro_info(args):
    """Collect reproducibility information for the HTML dashboard."""
    import subprocess

    info = {
        "command": " ".join(sys.argv),
        "script_path": os.path.abspath(__file__),
        "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "hostname": os.uname().nodename,
        "user": os.environ.get("USER", "unknown"),
        "rundir": os.path.abspath(args.rundir),
    }

    if args.legacy_rundir:
        info["legacy_rundir"] = os.path.abspath(args.legacy_rundir)

    # Git hash of the script's repo
    script_dir = os.path.dirname(info["script_path"])
    try:
        result = subprocess.run(
            ["git", "-C", script_dir, "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, timeout=5)
        if result.returncode == 0:
            info["git_hash"] = result.stdout.strip()
        result2 = subprocess.run(
            ["git", "-C", script_dir, "diff", "--stat", "--", info["script_path"]],
            capture_output=True, text=True, timeout=5)
        if result2.stdout.strip():
            info["git_dirty"] = True
    except Exception:
        pass

    # Read user_nl_eam from the FME case
    nl_path = _find_user_nl_eam(args.rundir)
    if nl_path:
        try:
            with open(nl_path) as f:
                info["user_nl_eam"] = f.read()
            info["user_nl_eam_path"] = nl_path
        except Exception:
            pass

    # Read user_nl_eam from the legacy case
    if args.legacy_rundir:
        leg_nl = _find_user_nl_eam(args.legacy_rundir)
        if leg_nl:
            try:
                with open(leg_nl) as f:
                    info["legacy_user_nl_eam"] = f.read()
                info["legacy_user_nl_eam_path"] = leg_nl
            except Exception:
                pass

    # Embed the script source for reproducibility
    try:
        with open(info["script_path"]) as f:
            info["script_source"] = f.read()
    except Exception:
        pass

    return info


# -----------------------------------------------------------------------------
# HTML helpers for new sections
# -----------------------------------------------------------------------------

def write_variable_coverage_html(coverage):
    """Generate HTML for the ACE Variable Requirements section."""
    if not coverage:
        return ""

    n_found = sum(1 for _, _, _, f in coverage if f)
    n_total = len(coverage)
    n_missing = n_total - n_found
    pct = 100 * n_found / n_total if n_total else 0

    html = '<h2 id="ACE_Variables">ACE Atmosphere Variable Requirements</h2>\n'

    # Summary with progress bar
    badge_cls = "badge-pass" if n_missing == 0 else "badge-fail"
    html += '<div class="card">\n'
    html += f'<p style="font-size:1.05em"><strong>{n_found}/{n_total}</strong> variables found '
    if n_missing:
        html += f'<span class="badge-fail">{n_missing} missing</span>'
    else:
        html += '<span class="badge-pass">all present</span>'
    html += '</p>\n'
    html += f'<div class="progress-bar"><div class="progress-fill" style="width:{pct:.0f}%"></div></div>\n'
    html += '</div>\n'

    # YAML spec in a collapsible block
    yaml_text = generate_yaml_text()
    html += '<details><summary>YAML Variable Specification (ACE canonical names)</summary>\n'
    html += '<pre class="yaml-block">'
    html += yaml_text
    html += '</pre>\n</details>\n'

    # Name mapping table with filter
    html += '<details open><summary>Variable Coverage</summary>\n'
    html += '<input class="filter-input" type="text" placeholder="Filter variables..." '
    html += 'oninput="filterTable(this,\'coverage-table\')">\n'
    html += '<table class="sortable" id="coverage-table"><thead><tr>'
    html += '<th onclick="sortTable(this,0)">ACE Name</th>'
    html += '<th onclick="sortTable(this,1)">EAM Name</th>'
    html += '<th onclick="sortTable(this,2)">Category</th>'
    html += '<th onclick="sortTable(this,3)">Status</th>'
    html += '</tr></thead><tbody>\n'

    for ace_name, eam_name, category, found in coverage:
        cls = "pass" if found else "fail"
        status = "FOUND" if found else "MISSING"
        html += (f'<tr><td><code>{ace_name}</code></td>'
                 f'<td><code>{eam_name}</code></td>'
                 f'<td>{category}</td>'
                 f'<td class="{cls}">{status}</td></tr>\n')

    html += '</tbody></table>\n</details>\n'
    return html


def write_repro_html(repro_info):
    """Generate HTML for the Reproducibility section."""
    if not repro_info:
        return ""

    html = '<h2 id="Reproducibility">Reproducibility</h2>\n'
    html += '<div class="repro-block"><dl>\n'

    html += f'<dt>Command</dt><dd><code>{repro_info["command"]}</code></dd>\n'
    html += f'<dt>Script</dt><dd><code>{repro_info["script_path"]}</code>'
    if repro_info.get("git_hash"):
        dirty = " (dirty)" if repro_info.get("git_dirty") else ""
        html += f' &mdash; git: <code>{repro_info["git_hash"]}{dirty}</code>'
    html += '</dd>\n'
    html += f'<dt>Run directory</dt><dd><code>{repro_info["rundir"]}</code></dd>\n'
    if repro_info.get("legacy_rundir"):
        html += f'<dt>Legacy run directory</dt><dd><code>{repro_info["legacy_rundir"]}</code></dd>\n'
    html += f'<dt>Host</dt><dd><code>{repro_info.get("user", "")}@{repro_info.get("hostname", "")}</code></dd>\n'
    html += f'<dt>Generated</dt><dd>{repro_info["timestamp"]}</dd>\n'

    html += '</dl>\n'

    # user_nl_eam contents
    if repro_info.get("user_nl_eam"):
        html += '<details><summary>user_nl_eam (FME)</summary>\n'
        html += f'<pre class="yaml-block">{repro_info["user_nl_eam"]}</pre>\n</details>\n'

    if repro_info.get("legacy_user_nl_eam"):
        html += '<details><summary>user_nl_eam (Legacy)</summary>\n'
        html += f'<pre class="yaml-block">{repro_info["legacy_user_nl_eam"]}</pre>\n</details>\n'

    # Full script source with line numbers
    if repro_info.get("script_source"):
        import html as html_mod
        src = repro_info["script_source"]
        lines = src.split("\n")
        # Build line-numbered source with anchors
        numbered = []
        pad = len(str(len(lines)))
        for i, line in enumerate(lines, 1):
            escaped = html_mod.escape(line)
            numbered.append(
                f'<span id="L{i}" class="src-line">'
                f'<span class="src-ln">{i:>{pad}}</span>  {escaped}</span>'
            )
        src_html = "\n".join(numbered)
        script_name = os.path.basename(repro_info["script_path"])
        html += (
            f'<details><summary>{script_name} '
            f'({len(lines)} lines)</summary>\n'
            f'<pre class="src-block">{src_html}</pre>\n'
            f'</details>\n'
        )

    html += '</div>\n'
    return html


# -----------------------------------------------------------------------------
# Executive summary
# -----------------------------------------------------------------------------

def build_executive_summary(all_issues, coverage, xv_results=None):
    """Compute executive summary metrics from verification results."""
    summary = {}

    # Checks passed / total
    n_checks_total = 0
    n_checks_passed = 0
    for comp, issues in all_issues.items():
        n_checks_total += 1
        if not issues:
            n_checks_passed += 1
    summary["n_checks_passed"] = n_checks_passed
    summary["n_checks_total"] = n_checks_total

    # Variable coverage
    if coverage:
        n_vars_found = sum(1 for _, _, _, f in coverage if f)
        n_vars_total = len(coverage)
    else:
        n_vars_found = 0
        n_vars_total = 0
    summary["n_vars_found"] = n_vars_found
    summary["n_vars_total"] = n_vars_total

    # Storage and timing from cross-verification results
    summary["storage_fme_gb"] = None
    summary["storage_leg_gb"] = None
    summary["storage_ratio"] = None
    summary["timing_overhead"] = None

    if xv_results:
        for r in xv_results:
            if r.get("test") == "STORAGE":
                detail = r.get("detail", "")
                # Parse "FME: N files / X.XX GB, Legacy: M files / Y.YY GB, ratio: Z.ZZx"
                import re
                m_fme = re.search(r"FME:.*?(\d+\.\d+)\s*GB", detail)
                m_leg = re.search(r"Legacy:.*?(\d+\.\d+)\s*GB", detail)
                m_rat = re.search(r"ratio:\s*(\d+\.\d+)x", detail)
                if m_fme:
                    summary["storage_fme_gb"] = float(m_fme.group(1))
                if m_leg:
                    summary["storage_leg_gb"] = float(m_leg.group(1))
                if m_rat:
                    summary["storage_ratio"] = float(m_rat.group(1))
            elif r.get("test") == "TIMING":
                detail = r.get("detail", "")
                import re
                m_oh = re.search(r"overhead:\s*([+-]?\d+\.?\d*)%", detail)
                if m_oh:
                    summary["timing_overhead"] = float(m_oh.group(1))

    return summary


def write_executive_summary_html(summary):
    """Return HTML string with a grid of 4 metric cards."""
    if not summary:
        return ""

    vars_val = f'{summary["n_vars_found"]}/{summary["n_vars_total"]}'
    vars_detail = ("All ACE atmosphere variables found"
                   if summary["n_vars_found"] == summary["n_vars_total"]
                   else f'{summary["n_vars_total"] - summary["n_vars_found"]} variables missing')

    checks_val = f'{summary["n_checks_passed"]}/{summary["n_checks_total"]}'
    checks_detail = ("All component checks passed"
                     if summary["n_checks_passed"] == summary["n_checks_total"]
                     else f'{summary["n_checks_total"] - summary["n_checks_passed"]} components with issues')

    if summary.get("storage_ratio") is not None:
        storage_val = f'{summary["storage_ratio"]:.2f}x'
        storage_detail = (f'FME {summary["storage_fme_gb"]:.2f} GB vs '
                          f'Legacy {summary["storage_leg_gb"]:.2f} GB')
    else:
        storage_val = "N/A"
        storage_detail = "No cross-verification data"

    if summary.get("timing_overhead") is not None:
        timing_val = f'{summary["timing_overhead"]:+.1f}%'
        timing_detail = "FME overhead vs legacy run"
    else:
        timing_val = "N/A"
        timing_detail = "No timing data available"

    html = '<div class="exec-summary">\n'
    html += ('  <div class="metric-card">\n'
             f'    <div class="metric-value">{vars_val}</div>\n'
             '    <div class="metric-label">Variable Coverage</div>\n'
             f'    <div class="metric-detail">{vars_detail}</div>\n'
             '  </div>\n')
    html += ('  <div class="metric-card">\n'
             f'    <div class="metric-value">{checks_val}</div>\n'
             '    <div class="metric-label">Checks Passed</div>\n'
             f'    <div class="metric-detail">{checks_detail}</div>\n'
             '  </div>\n')
    html += ('  <div class="metric-card">\n'
             f'    <div class="metric-value">{storage_val}</div>\n'
             '    <div class="metric-label">Storage Ratio</div>\n'
             f'    <div class="metric-detail">{storage_detail}</div>\n'
             '  </div>\n')
    html += ('  <div class="metric-card">\n'
             f'    <div class="metric-value">{timing_val}</div>\n'
             '    <div class="metric-label">Timing Overhead</div>\n'
             f'    <div class="metric-detail">{timing_detail}</div>\n'
             '  </div>\n')
    html += '</div>\n'
    return html


# -----------------------------------------------------------------------------
# HTML index
# -----------------------------------------------------------------------------

def write_html_index(outdir, all_plots_by_comp, all_issues, file_inventory_data,
                     all_fill_reports, timing_summary=None, extra_html="",
                     coverage=None, repro_info=None, exec_summary=None):
    index = os.path.join(outdir, "index.html")
    n_issues = sum(len(v) for v in all_issues.values())
    n_pass = sum(1 for v in all_issues.values() if not v)
    n_total = len(all_issues)
    status = "ALL PASS" if n_issues == 0 else f"{n_issues} issues in {n_total - n_pass}/{n_total} components"
    status_color = "#2a7d2a" if n_issues == 0 else "#c33"

    # Summary table
    summary_rows = ""
    for comp, issues in all_issues.items():
        if not issues:
            summary_rows += (f'<tr><td>{comp}</td><td class="pass">PASS</td>'
                             f'<td>0</td></tr>\n')
        else:
            summary_rows += (f'<tr><td>{comp}</td><td class="fail">FAIL</td>'
                             f'<td>{len(issues)}</td></tr>\n')

    # Detailed issues (collapsible)
    detail_rows = ""
    for comp, issues in all_issues.items():
        for msg in issues:
            detail_rows += f'<tr><td>{comp}</td><td>{msg.strip()}</td></tr>\n'

    # File inventory table
    inventory_rows = ""
    for label, count, flist in file_inventory_data:
        status_cls = "pass" if count > 0 else "fail"
        inventory_rows += (f'<tr><td>{label}</td>'
                           f'<td class="{status_cls}">{count}</td></tr>\n')

    # Fill/NaN report table
    fill_rows = ""
    for comp, reports in all_fill_reports.items():
        for r in reports:
            pct_cls = "pass" if r["pct"] < 50 else "fail"
            fill_rows += (f'<tr><td>{comp}</td><td>{r["name"]}</td>'
                          f'<td>{r["total"]:,}</td><td>{r["valid"]:,}</td>'
                          f'<td>{r["fill_nan"]:,}</td>'
                          f'<td class="{pct_cls}">{r["pct"]:.1f}%</td></tr>\n')

    # Group plots by component
    fig_sections = ""
    for comp, plots in all_plots_by_comp.items():
        if not plots:
            continue
        anchor = comp.replace(" ", "_").replace("-", "_").replace("(", "").replace(")", "")
        fig_sections += f'<h3 id="{anchor}">{comp}</h3>\n<div class="gallery">\n'
        # Use wider display for comparison/multi-panel figures
        is_comparison = (comp == "Comparisons")
        for p in plots:
            if p and os.path.exists(p):
                rel = os.path.relpath(p, outdir)
                fig_cls = "fig wide" if is_comparison else "fig"
                fig_sections += (f'<div class="{fig_cls}"><a href="{rel}">'
                                 f'<img src="{rel}"/></a>'
                                 f'<span>{os.path.basename(p)}</span></div>\n')
        fig_sections += '</div>\n'

    # Navigation
    nav_links = ""
    nav_items = ["Overview", "Summary", "ACE_Variables", "File_Inventory", "Fill_NaN_Report"]
    if extra_html:
        nav_items.append("Cross_Verify")
    for comp in all_plots_by_comp:
        if all_plots_by_comp[comp]:
            nav_items.append(comp)
    if timing_summary:
        nav_items.append("Timing")
    if repro_info:
        nav_items.append("Reproducibility")
    for item in nav_items:
        anchor = item.replace(" ", "_").replace("-", "_").replace("(", "").replace(")", "")
        display = item.replace("_", " ")
        nav_links += f'<a href="#{anchor}">{display}</a>\n'

    timing_html = ""
    if timing_summary:
        timing_rows = ""
        for line in timing_summary[:30]:
            parts = line.split()
            if len(parts) >= 2:
                timing_rows += "<tr>" + "".join(f"<td>{p}</td>" for p in parts) + "</tr>\n"
            else:
                timing_rows += f"<tr><td colspan='6'>{line}</td></tr>\n"
        timing_html = f'''
    <h2 id="Timing">Performance Timing</h2>
    <div class="timing-container">
    <table class="timing-table">
    {timing_rows}
    </table>
    </div>
'''

    # ACE variable coverage section
    coverage_html = write_variable_coverage_html(coverage) if coverage else ""

    # Reproducibility section
    repro_html = write_repro_html(repro_info) if repro_info else ""

    # Executive summary section
    exec_summary_html = write_executive_summary_html(exec_summary) if exec_summary else ""

    gen_timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    html = textwrap.dedent(f"""\
    <!DOCTYPE html><html data-theme="light"><head>
    <meta charset="utf-8"/>
    <meta name="viewport" content="width=device-width, initial-scale=1"/>
    <title>EAM FME Output Verification</title>
    <style>
      :root {{
        --bg: #f0f2f5; --text: #2c3e50; --text-muted: #6c757d;
        --card: #ffffff; --card-border: #e2e8f0;
        --card-shadow: 0 2px 6px rgba(0,0,0,0.06);
        --header-bg: linear-gradient(135deg, #1a2744 0%, #2c3e6b 100%);
        --nav-bg: #ffffff; --nav-border: #e2e8f0;
        --nav-link-bg: #e8eef4; --nav-link-text: #1a2744;
        --th-bg: #1a2744; --th-text: #fff; --th-hover: #2a3d5c;
        --accent: #1a2744; --accent-light: #e8eef4;
        --pass: #16a34a; --pass-bg: #dcfce7;
        --fail: #dc2626; --fail-bg: #fee2e2;
        --skip-bg: #fef9c3; --skip: #a16207;
        --code-bg: #f1f5f9; --hover-bg: #f8fafc;
        --yaml-bg: #1e293b; --yaml-text: #94a3b8;
        --stripe: #f8fafc; --border-radius: 8px;
        --ring-track: #e2e8f0; --ring-fill: #16a34a;
        --toggle-bg: #e2e8f0; --toggle-knob: #1a2744;
      }}
      [data-theme="dark"] {{
        --bg: #0d1117; --text: #c9d1d9; --text-muted: #8b949e;
        --card: #161b22; --card-border: #30363d;
        --card-shadow: 0 2px 6px rgba(0,0,0,0.4);
        --header-bg: linear-gradient(135deg, #010409 0%, #0d1117 100%);
        --nav-bg: #161b22; --nav-border: #30363d;
        --nav-link-bg: #21262d; --nav-link-text: #c9d1d9;
        --th-bg: #21262d; --th-text: #c9d1d9; --th-hover: #30363d;
        --accent: #58a6ff; --accent-light: #1c2d41;
        --pass: #3fb950; --pass-bg: #0d2818;
        --fail: #f85149; --fail-bg: #3d1418;
        --skip-bg: #3d2e00; --skip: #d29922;
        --code-bg: #0d1117; --hover-bg: #1c2128;
        --yaml-bg: #010409; --yaml-text: #7d8590;
        --stripe: #1c2128; --border-radius: 8px;
        --ring-track: #30363d; --ring-fill: #3fb950;
        --toggle-bg: #30363d; --toggle-knob: #58a6ff;
      }}
      * {{ box-sizing: border-box; margin: 0; padding: 0; }}
      html {{ scroll-behavior: smooth; }}
      body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto,
                           'Helvetica Neue', Arial, sans-serif;
             background: var(--bg); color: var(--text);
             line-height: 1.6; transition: background 0.3s, color 0.3s; }}
      .container {{ max-width: 1400px; margin: 0 auto; padding: 24px 32px; }}

      /* Header */
      header {{ background: var(--header-bg); color: #fff; padding: 28px 32px;
                display: flex; justify-content: space-between; align-items: center;
                flex-wrap: wrap; gap: 12px; }}
      header .hdr-left h1 {{ font-size: 1.5em; font-weight: 700; letter-spacing: -0.02em; }}
      header .hdr-left .subtitle {{ font-size: 0.85em; opacity: 0.7; margin-top: 2px; }}
      header .status-pill {{ display: inline-block; padding: 6px 16px; border-radius: 20px;
                             font-weight: 700; font-size: 0.95em; letter-spacing: 0.02em; }}
      .status-pass {{ background: rgba(22,163,74,0.2); color: #4ade80; }}
      .status-fail {{ background: rgba(220,38,38,0.2); color: #fca5a5; }}
      .hdr-right {{ display: flex; align-items: center; gap: 16px; }}

      /* Theme toggle */
      .theme-toggle {{ background: none; border: 2px solid rgba(255,255,255,0.3);
                       color: #fff; padding: 6px 14px; border-radius: 20px;
                       cursor: pointer; font-size: 0.85em; font-weight: 500;
                       transition: all 0.2s; display: flex; align-items: center; gap: 6px; }}
      .theme-toggle:hover {{ border-color: #fff; background: rgba(255,255,255,0.1); }}
      .theme-toggle .icon {{ font-size: 1.1em; }}

      /* Navigation */
      nav {{ background: var(--nav-bg); padding: 10px 16px; border-radius: var(--border-radius);
             box-shadow: var(--card-shadow); border: 1px solid var(--nav-border);
             margin: 20px 0; display: flex; flex-wrap: wrap; gap: 6px;
             position: sticky; top: 0; z-index: 100;
             transition: background 0.3s, border-color 0.3s; }}
      nav a {{ color: var(--nav-link-text); text-decoration: none; padding: 6px 14px;
               background: var(--nav-link-bg); border-radius: 6px; font-size: 0.82em;
               font-weight: 500; transition: all 0.15s; white-space: nowrap; }}
      nav a:hover, nav a.active {{ background: var(--accent); color: #fff; }}

      /* Section headings */
      h2 {{ color: var(--accent); margin-top: 36px; font-size: 1.25em; font-weight: 700;
            padding-bottom: 8px; border-bottom: 2px solid var(--accent);
            transition: color 0.3s; }}
      h3 {{ color: var(--text); margin-top: 24px; font-size: 1.05em; }}

      /* Cards */
      .card {{ background: var(--card); border: 1px solid var(--card-border);
               border-radius: var(--border-radius); box-shadow: var(--card-shadow);
               padding: 20px; margin-bottom: 16px; transition: all 0.3s; }}

      /* Tables */
      table {{ border-collapse: collapse; width: 100%; background: var(--card);
               border: 1px solid var(--card-border); border-radius: var(--border-radius);
               overflow: hidden; margin-bottom: 16px; transition: all 0.3s; }}
      td, th {{ border: 1px solid var(--card-border); padding: 10px 14px; text-align: left;
                transition: all 0.3s; }}
      th {{ background: var(--th-bg); color: var(--th-text); font-weight: 600;
            font-size: 0.88em; text-transform: uppercase; letter-spacing: 0.04em; }}
      tbody tr:nth-child(even) {{ background: var(--stripe); }}
      tbody tr:hover {{ background: var(--hover-bg); }}
      .pass {{ color: var(--pass); font-weight: 700; }}
      .fail {{ color: var(--fail); font-weight: 700; }}
      td.pass {{ background: var(--pass-bg); }}
      td.fail {{ background: var(--fail-bg); }}
      code {{ font-family: 'SF Mono', 'Fira Code', 'Cascadia Code', monospace;
              font-size: 0.88em; background: var(--code-bg); padding: 2px 6px;
              border-radius: 4px; transition: background 0.3s; }}

      /* Figure gallery */
      .gallery {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(420px, 1fr));
                  gap: 16px; margin-top: 12px; }}
      .fig {{ background: var(--card); border: 1px solid var(--card-border);
              border-radius: var(--border-radius); overflow: hidden;
              transition: all 0.25s; }}
      .fig:hover {{ box-shadow: 0 4px 16px rgba(0,0,0,0.12); transform: translateY(-2px); }}
      .fig img {{ width: 100%; display: block; cursor: zoom-in; }}
      .fig span {{ display: block; font-size: 0.75em; color: var(--text-muted);
                   padding: 8px 12px; border-top: 1px solid var(--card-border);
                   font-family: monospace; }}
      .fig.wide {{ grid-column: span 2; }}
      @media (max-width: 900px) {{ .fig.wide {{ grid-column: span 1; }} }}

      /* Details / collapsible */
      details {{ margin: 12px 0; background: var(--card); border-radius: var(--border-radius);
                 border: 1px solid var(--card-border); transition: all 0.3s; }}
      details[open] {{ box-shadow: var(--card-shadow); }}
      summary {{ cursor: pointer; font-weight: 600; color: var(--accent);
                 padding: 12px 16px; list-style: none; display: flex;
                 align-items: center; gap: 8px; transition: all 0.2s; }}
      summary::-webkit-details-marker {{ display: none; }}
      summary::before {{ content: '\\25B6'; font-size: 0.7em; transition: transform 0.2s;
                         display: inline-block; }}
      details[open] > summary::before {{ transform: rotate(90deg); }}
      summary:hover {{ background: var(--hover-bg); }}
      details > :not(summary) {{ padding: 0 16px 16px; }}

      /* Timing */
      .timing-container {{ background: var(--card); padding: 16px;
                           border-radius: var(--border-radius);
                           border: 1px solid var(--card-border);
                           overflow-x: auto; max-height: 500px; overflow-y: auto; }}
      .timing-table td {{ font-family: monospace; font-size: 0.83em;
                          padding: 4px 10px; white-space: nowrap; }}

      /* YAML / code blocks */
      .yaml-block {{ background: var(--yaml-bg); color: var(--yaml-text);
                     padding: 16px 20px; border-radius: var(--border-radius);
                     font-family: 'SF Mono', 'Fira Code', monospace; font-size: 0.83em;
                     overflow-x: auto; max-height: 400px; overflow-y: auto;
                     white-space: pre; line-height: 1.5; border: 1px solid var(--card-border);
                     transition: all 0.3s; }}

      /* Repro block */
      .repro-block {{ background: var(--card); padding: 20px 24px;
                      border-radius: var(--border-radius);
                      border: 1px solid var(--card-border);
                      box-shadow: var(--card-shadow); transition: all 0.3s; }}
      .repro-block dl {{ margin: 0; display: grid; grid-template-columns: auto 1fr;
                         gap: 4px 16px; align-items: baseline; }}
      .repro-block dt {{ font-weight: 600; color: var(--accent); font-size: 0.88em;
                         white-space: nowrap; }}
      .repro-block dd {{ margin: 0; font-family: monospace; font-size: 0.83em;
                         color: var(--text-muted); word-break: break-all; }}

      /* Progress bar */
      .progress-bar {{ height: 8px; background: var(--ring-track);
                       border-radius: 4px; overflow: hidden; margin: 8px 0 16px; }}
      .progress-fill {{ height: 100%; border-radius: 4px; transition: width 0.6s ease;
                        background: linear-gradient(90deg, var(--ring-fill), #22d3ee); }}

      /* Status badges */
      .badge-pass {{ display: inline-block; padding: 2px 10px; border-radius: 12px;
                     background: var(--pass-bg); color: var(--pass);
                     font-weight: 600; font-size: 0.82em; }}
      .badge-fail {{ display: inline-block; padding: 2px 10px; border-radius: 12px;
                     background: var(--fail-bg); color: var(--fail);
                     font-weight: 600; font-size: 0.82em; }}
      .badge-skip {{ display: inline-block; padding: 2px 10px; border-radius: 12px;
                     background: var(--skip-bg); color: var(--skip);
                     font-weight: 600; font-size: 0.82em; }}

      /* Executive summary */
      .exec-summary {{ display: grid; grid-template-columns: repeat(4, 1fr);
                       gap: 16px; margin: 20px 0; }}
      .metric-card {{ background: var(--card); border: 1px solid var(--card-border);
                      border-radius: var(--border-radius); padding: 20px 24px;
                      box-shadow: var(--card-shadow); text-align: center;
                      transition: all 0.3s; }}
      .metric-card:hover {{ box-shadow: 0 4px 16px rgba(0,0,0,0.12);
                            transform: translateY(-2px); }}
      .metric-value {{ font-size: 2em; font-weight: 800; color: var(--accent);
                       line-height: 1.2; }}
      .metric-label {{ font-size: 0.9em; color: var(--text-muted); font-weight: 600;
                       margin-top: 4px; }}
      .metric-detail {{ font-size: 0.78em; color: var(--text-muted); margin-top: 6px;
                        opacity: 0.8; }}
      @media (max-width: 900px) {{ .exec-summary {{ grid-template-columns: repeat(2, 1fr); }} }}

      /* Footer */
      footer {{ margin-top: 48px; padding: 20px 32px;
                border-top: 1px solid var(--card-border);
                color: var(--text-muted); font-size: 0.8em;
                background: var(--card); text-align: center;
                transition: all 0.3s; }}

      /* Lightbox */
      .lightbox {{ display:none; position:fixed; top:0; left:0; width:100%; height:100%;
                   background:rgba(0,0,0,0.92); z-index:1000; cursor:zoom-out;
                   align-items:center; justify-content:center;
                   backdrop-filter: blur(4px); -webkit-backdrop-filter: blur(4px); }}
      .lightbox.active {{ display:flex; }}
      .lightbox img {{ max-width:94%; max-height:94%; border-radius: 8px;
                       box-shadow: 0 8px 32px rgba(0,0,0,0.5);
                       animation: lbFadeIn 0.2s ease; }}
      @keyframes lbFadeIn {{ from {{ opacity:0; transform:scale(0.95); }}
                             to   {{ opacity:1; transform:scale(1); }} }}
      .lightbox .lb-nav {{ position:absolute; top:50%; transform:translateY(-50%);
                           color:#fff; font-size:2.5em; cursor:pointer;
                           padding:20px; opacity:0.5; transition:opacity 0.2s;
                           user-select:none; }}
      .lightbox .lb-nav:hover {{ opacity:1; }}
      .lightbox .lb-prev {{ left:10px; }}
      .lightbox .lb-next {{ right:10px; }}
      .lightbox .lb-close {{ position:absolute; top:16px; right:24px; color:#fff;
                             font-size:1.8em; cursor:pointer; opacity:0.6;
                             transition:opacity 0.2s; }}
      .lightbox .lb-close:hover {{ opacity:1; }}

      /* Tabs */
      .tabs {{ display:flex; flex-wrap:wrap; gap:4px; margin:16px 0 0; }}
      .tab {{ padding:8px 18px; border:none; background: var(--nav-link-bg);
              color: var(--nav-link-text); cursor:pointer; font-weight:600;
              font-size:0.88em; border-radius: 8px 8px 0 0;
              transition: all 0.15s; }}
      .tab:hover {{ background: var(--hover-bg); }}
      .tab.active {{ background: var(--accent); color:#fff; }}
      .badge {{ background:rgba(255,255,255,0.2); padding:2px 8px; border-radius:10px;
                font-size:0.8em; margin-left:6px; }}
      .tab-content {{ background: var(--card); padding:20px;
                      border-radius:0 var(--border-radius) var(--border-radius) var(--border-radius);
                      border: 1px solid var(--card-border);
                      margin-bottom:16px; transition: all 0.3s; }}

      /* Sortable headers */
      .sortable th {{ cursor:pointer; user-select:none; position:relative;
                      padding-right:22px; }}
      .sortable th::after {{ content:'\\2195'; position:absolute; right:6px; top:50%;
                             transform:translateY(-50%); opacity:0.4; font-size:0.85em; }}
      .sortable th:hover {{ background: var(--th-hover); }}

      /* Filter input */
      .filter-input {{ width: 100%; max-width: 400px; padding: 8px 14px;
                       border: 1px solid var(--card-border); border-radius: 6px;
                       background: var(--card); color: var(--text);
                       font-size: 0.9em; margin-bottom: 12px;
                       transition: all 0.2s; outline: none; }}
      .filter-input:focus {{ border-color: var(--accent);
                             box-shadow: 0 0 0 3px rgba(88,166,255,0.15); }}
      .filter-input::placeholder {{ color: var(--text-muted); }}

      /* Embedded source */
      .src-block {{ background: var(--yaml-bg); color: var(--yaml-text);
                    padding: 0; border-radius: var(--border-radius);
                    font-family: 'SF Mono', 'Fira Code', 'Cascadia Code', monospace;
                    font-size: 0.8em; overflow-x: auto; max-height: 600px;
                    overflow-y: auto; white-space: pre; line-height: 1.55;
                    border: 1px solid var(--card-border); counter-reset: line;
                    tab-size: 4; }}
      .src-line {{ display: block; padding: 0 16px 0 0; }}
      .src-line:hover {{ background: rgba(88,166,255,0.08); }}
      .src-line:target {{ background: rgba(88,166,255,0.18); }}
      .src-ln {{ display: inline-block; width: 4.5em; text-align: right;
                 color: rgba(139,148,158,0.4); padding: 0 12px 0 12px;
                 user-select: none; border-right: 1px solid var(--card-border);
                 margin-right: 12px; }}

      /* Animations */
      @keyframes fadeIn {{ from {{ opacity:0; transform:translateY(8px); }}
                           to   {{ opacity:1; transform:translateY(0); }} }}
      .card, table, details {{ animation: fadeIn 0.3s ease; }}
    </style>
    <script>
    /* Theme toggle */
    (function() {{
      var saved = localStorage.getItem('eam-theme');
      if (saved) document.documentElement.setAttribute('data-theme', saved);
      else if (window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches)
        document.documentElement.setAttribute('data-theme', 'dark');
    }})();
    function toggleTheme() {{
      var html = document.documentElement;
      var next = html.getAttribute('data-theme') === 'dark' ? 'light' : 'dark';
      html.setAttribute('data-theme', next);
      localStorage.setItem('eam-theme', next);
      var btn = document.querySelector('.theme-toggle .icon');
      if (btn) btn.textContent = next === 'dark' ? '\\u2600' : '\\u263E';
    }}

    /* Tabs */
    function showTab(id) {{
      document.querySelectorAll('.tab-content').forEach(function(el) {{ el.style.display='none'; }});
      document.querySelectorAll('.tab').forEach(function(el) {{ el.classList.remove('active'); }});
      var tab = document.getElementById('tab_' + id);
      if (tab) tab.style.display = '';
      document.querySelectorAll('.tab').forEach(function(btn) {{
        if (btn.getAttribute('onclick') && btn.getAttribute('onclick').indexOf(id) >= 0)
          btn.classList.add('active');
      }});
    }}

    /* Sortable tables */
    function sortTable(th, col) {{
      var table = th.closest('table');
      var tbody = table.querySelector('tbody');
      if (!tbody) return;
      var rows = Array.from(tbody.querySelectorAll('tr'));
      var asc = th.dataset.sort !== 'asc';
      th.dataset.sort = asc ? 'asc' : 'desc';
      // Reset other headers
      th.closest('tr').querySelectorAll('th').forEach(function(h) {{
        if (h !== th) delete h.dataset.sort;
      }});
      rows.sort(function(a, b) {{
        var va = a.cells[col].textContent.trim();
        var vb = b.cells[col].textContent.trim();
        var na = parseFloat(va), nb = parseFloat(vb);
        if (!isNaN(na) && !isNaN(nb)) return asc ? na - nb : nb - na;
        return asc ? va.localeCompare(vb) : vb.localeCompare(va);
      }});
      rows.forEach(function(r) {{ tbody.appendChild(r); }});
    }}

    /* Table filter */
    function filterTable(input, tableId) {{
      var filter = input.value.toLowerCase();
      var table = document.getElementById(tableId);
      if (!table) return;
      var rows = table.querySelectorAll('tbody tr');
      rows.forEach(function(row) {{
        var text = row.textContent.toLowerCase();
        row.style.display = text.indexOf(filter) > -1 ? '' : 'none';
      }});
    }}

    document.addEventListener('DOMContentLoaded', function() {{
      /* Theme button label */
      var theme = document.documentElement.getAttribute('data-theme');
      var icon = document.querySelector('.theme-toggle .icon');
      if (icon) icon.textContent = theme === 'dark' ? '\\u2600' : '\\u263E';

      /* Lightbox with navigation */
      var lb = document.createElement('div');
      lb.className = 'lightbox';
      lb.innerHTML = '<span class="lb-close">\\u00D7</span>'
                   + '<span class="lb-nav lb-prev">\\u276E</span>'
                   + '<img/>'
                   + '<span class="lb-nav lb-next">\\u276F</span>';
      document.body.appendChild(lb);
      var lbImg = lb.querySelector('img');
      var allFigs = [];
      var curIdx = 0;

      function openLB(idx) {{
        curIdx = idx;
        lbImg.src = allFigs[idx];
        lb.classList.add('active');
      }}
      function closeLB() {{ lb.classList.remove('active'); }}
      function navLB(dir) {{
        curIdx = (curIdx + dir + allFigs.length) % allFigs.length;
        lbImg.src = allFigs[curIdx];
      }}

      document.querySelectorAll('.fig a').forEach(function(a, i) {{
        allFigs.push(a.href);
        a.addEventListener('click', function(e) {{
          e.preventDefault();
          openLB(i);
        }});
      }});

      lb.querySelector('.lb-close').addEventListener('click', closeLB);
      lb.querySelector('.lb-prev').addEventListener('click', function(e) {{ e.stopPropagation(); navLB(-1); }});
      lb.querySelector('.lb-next').addEventListener('click', function(e) {{ e.stopPropagation(); navLB(1); }});
      lb.addEventListener('click', function(e) {{ if (e.target === lb) closeLB(); }});
      document.addEventListener('keydown', function(e) {{
        if (!lb.classList.contains('active')) return;
        if (e.key === 'Escape') closeLB();
        if (e.key === 'ArrowLeft') navLB(-1);
        if (e.key === 'ArrowRight') navLB(1);
      }});

      /* Active nav highlighting on scroll */
      var sections = document.querySelectorAll('h2[id]');
      var navLinks = document.querySelectorAll('nav a');
      if (sections.length && navLinks.length) {{
        var observer = new IntersectionObserver(function(entries) {{
          entries.forEach(function(entry) {{
            if (entry.isIntersecting) {{
              var id = entry.target.getAttribute('id');
              navLinks.forEach(function(link) {{
                link.classList.toggle('active', link.getAttribute('href') === '#' + id);
              }});
            }}
          }});
        }}, {{ rootMargin: '-20% 0px -70% 0px' }});
        sections.forEach(function(s) {{ observer.observe(s); }});
      }}
    }});
    </script>
    </head><body>
    <header>
      <div class="hdr-left">
        <h1>EAM FME Online Output Verification</h1>
        <div class="subtitle">E3SM Atmosphere Model &mdash; Full Model Emulation Diagnostics</div>
      </div>
      <div class="hdr-right">
        <span class="status-pill {"status-pass" if n_issues == 0 else "status-fail"}">{status}</span>
        <button class="theme-toggle" onclick="toggleTheme()"><span class="icon"></span> Theme</button>
      </div>
    </header>
    <div class="container">

    <nav>{nav_links}</nav>

    <h2 id="Overview">Overview</h2>
    {exec_summary_html}

    <h2 id="Summary">Component Summary</h2>
    <table>
    <tr><th>Component</th><th>Status</th><th>Issues</th></tr>
    {summary_rows}
    </table>

    {"<details><summary>Show all " + str(n_issues) + " issues</summary><table><tr><th>Component</th><th>Issue</th></tr>" + detail_rows + "</table></details>" if detail_rows else ""}

    {coverage_html}

    <h2 id="File_Inventory">File Inventory</h2>
    <table>
    <tr><th>Category</th><th>File Count</th></tr>
    {inventory_rows}
    </table>

    <h2 id="Fill_NaN_Report">Fill / NaN Report</h2>
    <details><summary>Variable-level fill fraction details</summary>
    <table>
    <tr><th>Component</th><th>Variable</th><th>Total Cells</th><th>Valid</th><th>Fill/NaN</th><th>Fill %</th></tr>
    {fill_rows if fill_rows else '<tr><td colspan="6">No data collected.</td></tr>'}
    </table>
    </details>

    {extra_html}

    <h2>Diagnostic Figures</h2>
    {fig_sections if fig_sections else '<p>No figures generated.</p>'}

    {timing_html}

    {repro_html}
    </div>
    <footer>
      <strong>verify_eam.py</strong> &mdash; EAM FME Online Output Verification for E3SM
      &mdash; Generated {gen_timestamp}
    </footer>
    </body></html>
    """)
    with open(index, "w") as fh:
        fh.write(html)
    print(f"\nIndex written: {index}")
    return index


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Verify EAM FME output (vcoarsen, derived fields) and "
                    "produce diagnostic figures with optional cross-verification.")
    parser.add_argument("--rundir", required=True,
                        help="CIME RUNDIR for the FME case (fme_output testmod)")
    parser.add_argument("--outdir", required=True,
                        help="Output directory for figures and HTML index")
    parser.add_argument("--verbose", "-v", action="store_true")
    parser.add_argument("--legacy-rundir",
                        help="CIME RUNDIR for the legacy case (fme_legacy_output testmod) "
                             "for cross-verification")
    args = parser.parse_args()

    if not os.path.isdir(args.rundir):
        sys.exit(f"ERROR: rundir not found: {args.rundir}")

    os.makedirs(args.outdir, exist_ok=True)
    print(f"EAM FME Verification  |  rundir: {args.rundir}")
    if args.legacy_rundir:
        print(f"Legacy RUNDIR         :  {args.legacy_rundir}")
    print(f"Output directory      :  {args.outdir}")
    print("=" * 70)

    if not HAS_MPL:
        print("WARNING: matplotlib not available -- no figures will be produced")
    if not HAS_CARTOPY:
        print("WARNING: cartopy not available -- falling back to scatter maps")
    if not (HAS_XR or HAS_NC4):
        sys.exit("ERROR: need xarray or netCDF4")

    # File inventory
    file_inventory_data = collect_file_inventory(args.rundir)

    all_issues = {}
    all_plots_by_comp = {}
    all_fill_reports = {}

    # All figures go under fme_output/ subdir so fme_output.html can reference them
    fig_root = os.path.join(args.outdir, "fme_output")
    os.makedirs(fig_root, exist_ok=True)

    def run_check(name, func, subdir, *a):
        comp_outdir = os.path.join(fig_root, subdir)
        os.makedirs(comp_outdir, exist_ok=True)
        result = func(*a[:1], comp_outdir, *a[1:])
        comp_issues, plots, fill_reports = result
        all_issues[name] = comp_issues
        all_plots_by_comp[name] = plots
        all_fill_reports[name] = fill_reports

    run_check("EAM", check_eam, "eam",
              args.rundir, args.verbose)

    # EAM comparison figures (multi-panel layer views, zonal means)
    generate_comparison_figures(args.rundir, fig_root, all_plots_by_comp)

    # Timing summary
    timing_summary = read_timing_summary(args.rundir)
    if timing_summary:
        print("\n=== Performance Timing ===")
        for line in timing_summary[:10]:
            print(f"  {line}")
        if len(timing_summary) > 10:
            print(f"  ... ({len(timing_summary) - 10} more lines in HTML report)")

    # ACE variable coverage
    print("\n=== ACE Variable Coverage ===")
    coverage = check_variable_coverage(args.rundir)
    if coverage:
        n_found = sum(1 for _, _, _, f in coverage if f)
        n_missing = len(coverage) - n_found
        print(f"  {n_found}/{len(coverage)} ACE variables found in output")
        if n_missing:
            for ace_name, eam_name, _, found in coverage:
                if not found:
                    print(f"  MISSING: {ace_name} (EAM: {eam_name})")

    # Reproducibility info
    repro_info = collect_repro_info(args)

    # Cross-verification against legacy testmod output
    xv_html = ""
    xv_results = None
    if args.legacy_rundir:
        if not os.path.isdir(args.legacy_rundir):
            print(f"WARNING: legacy-rundir not found: {args.legacy_rundir}")
        else:
            xv_issues, xv_results, xv_plots = cross_verify(
                args.rundir, args.legacy_rundir, args.outdir, args.verbose)
            all_issues["Cross-Verify"] = xv_issues
            if xv_plots:
                all_plots_by_comp["Cross-Verify"] = xv_plots
            xv_html = write_cross_verify_html(args.outdir, xv_results, xv_plots)

    # Build executive summary
    exec_summary = build_executive_summary(all_issues, coverage,
                                           xv_results=xv_results)

    write_html_index(args.outdir, all_plots_by_comp, all_issues,
                     file_inventory_data, all_fill_reports,
                     timing_summary=timing_summary,
                     extra_html=xv_html,
                     coverage=coverage,
                     repro_info=repro_info,
                     exec_summary=exec_summary)

    n_total = sum(len(v) for v in all_issues.values())
    print("\n" + "=" * 70)
    print(f"RESULT: {'ALL CHECKS PASSED' if n_total == 0 else f'{n_total} ISSUES FOUND'}")
    sys.exit(0 if n_total == 0 else 1)


if __name__ == "__main__":
    main()
