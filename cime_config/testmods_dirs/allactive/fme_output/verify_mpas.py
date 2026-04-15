#!/usr/bin/env python3
"""
Verify MPAS-Ocean and MPAS-SeaIce FME online output (depth coarsening,
derived fields, sea ice diagnostics). Produces HTML dashboard with field
presence checks, range checks, and diagnostic figures.

Usage:
    python verify_mpas.py --rundir /path/to/RUNDIR --outdir /path/to/figs
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

# -- ACE MPAS variable requirements -------------------------------------------
# Single source of truth for the ACE-to-MPAS variable mapping.
# Drives verification, HTML coverage table, and reproducibility tracking.

N_OCN_LAYERS = 19  # ACE uses 19 depth layers (0-based)

# 3D depth-coarsened fields: ACE base name -> MPAS FME output base name
ACE_OCN_3D_BASE_MAP = {
    "thetao": "temperatureCoarsened",
    "so": "salinityCoarsened",
    "uo": "velocityZonalCoarsened",
    "vo": "velocityMeridionalCoarsened",
}

# 2D / scalar field mapping: ACE canonical name -> MPAS FME output name
# Note: atmospheric forcings (from EAM) are listed with source="atm"
ACE_MPAS_2D_MAP = {
    "zos": "ssh",
    "sst": "sst",
    "ocean_sea_ice_fraction": "iceAreaTotal",
    "sea_ice_volume": "iceVolumeTotal",
    "land_fraction": "LANDFRAC",        # from EAM
    "TAUX": "TAUX",                     # from EAM
    "TAUY": "TAUY",                     # from EAM
    "surface_precipitation_rate": "PRECT",           # from EAM
    "surface_upward_longwave_flux": "FLUS",          # from EAM
    "surface_upward_shortwave_flux": "FSUS",         # from EAM
    "FLDS": "FLDS",                     # from EAM
    "FSDS": "FSDS",                     # from EAM
    "LHFLX": "LHFLX",                  # from EAM
    "SHFLX": "SHFLX",                  # from EAM
}

# Variables that come from EAM (not in MPAS output files)
ACE_MPAS_ATM_FORCINGS = {
    "TAUX", "TAUY", "surface_precipitation_rate",
    "surface_upward_longwave_flux", "surface_upward_shortwave_flux",
    "FLDS", "FSDS", "LHFLX", "SHFLX", "land_fraction",
}

# ACE variable lists (mirror the YAML specification)
ACE_MPAS_FORCING_NAMES = [
    "TAUX", "TAUY", "surface_precipitation_rate",
    "surface_upward_longwave_flux", "surface_upward_shortwave_flux",
    "FLDS", "FSDS", "LHFLX", "SHFLX",
]

ACE_MPAS_IN_NAMES = ([
    "zos",
    "TAUX", "TAUY", "surface_precipitation_rate",
    "surface_upward_longwave_flux", "surface_upward_shortwave_flux",
    "FLDS", "FSDS", "LHFLX", "SHFLX",
    "land_fraction", "ocean_sea_ice_fraction", "sea_ice_volume", "sst",
] + [f"so_{k}" for k in range(N_OCN_LAYERS)]
  + [f"thetao_{k}" for k in range(N_OCN_LAYERS)]
  + [f"uo_{k}" for k in range(N_OCN_LAYERS)]
  + [f"vo_{k}" for k in range(N_OCN_LAYERS)])

ACE_MPAS_OUT_NAMES = ([
    "zos", "sst",
] + [f"so_{k}" for k in range(N_OCN_LAYERS)]
  + [f"thetao_{k}" for k in range(N_OCN_LAYERS)]
  + [f"uo_{k}" for k in range(N_OCN_LAYERS)]
  + [f"vo_{k}" for k in range(N_OCN_LAYERS)]
  + ["ocean_sea_ice_fraction", "sea_ice_volume"])


def mpas_ace_to_fme(ace_name):
    """Map an ACE canonical variable name to its MPAS FME output name."""
    if ace_name in ACE_MPAS_2D_MAP:
        return ACE_MPAS_2D_MAP[ace_name]
    # 3D indexed variables (e.g. thetao_3 -> temperatureCoarsened_3)
    for ace_base in sorted(ACE_OCN_3D_BASE_MAP, key=len, reverse=True):
        prefix = ace_base + "_"
        if ace_name.startswith(prefix):
            suffix = ace_name[len(ace_base):]
            return ACE_OCN_3D_BASE_MAP[ace_base] + suffix
    return ace_name


def generate_mpas_yaml_text():
    """Generate YAML-formatted text of the ACE MPAS variable spec."""
    lines = ["## MPAS", "next_step_forcing_names:"]
    for name in ACE_MPAS_FORCING_NAMES:
        lines.append(f"  - {name}")
    lines.append("in_names:")
    for name in ACE_MPAS_IN_NAMES:
        lines.append(f"  - {name}")
    lines.append("out_names:")
    for name in ACE_MPAS_OUT_NAMES:
        lines.append(f"  - {name}")
    return "\n".join(lines)


# -- existing MPAS variable requirements (used by verification checks) --------

# Fields ACE ocean model needs from MPAS-O fmeDepthCoarsening (19 layers).
ACE_OCN_COARSENED = {
    "temperatureCoarsened": "thetao",
    "salinityCoarsened": "so",
    "velocityZonalCoarsened": "uo",
    "velocityMeridionalCoarsened": "vo",
    "layerThicknessCoarsened": "layerThickness",
}
N_OCN_LAYERS = None  # auto-detect from file (was 19 for ACE, 25 for default E3SM)

# Fields ACE ocean model needs from MPAS-O fmeDerivedFields.
ACE_OCN_DERIVED = {
    "sst": "sst",
    "sss": "sss",
    "surfaceHeatFluxTotal": "surfaceHeatFluxTotal",
}

# Fields ACE needs from MPAS-O fmeVerticalReduce.
ACE_OCN_VERTREDUCE = {
    "oceanHeatContent": "oceanHeatContent",
    "freshwaterContent": "freshwaterContent",
    "kineticEnergy": "kineticEnergy",
}

# Fields ACE needs from MPAS-SI fmeSeaiceDerivedFields.
ACE_ICE_DERIVED = {
    "iceAreaTotal": "ocean_sea_ice_fraction",
    "iceVolumeTotal": "sea_ice_volume",
    "snowVolumeTotal": "snowVolumeTotal",
    "iceThicknessMean": "iceThicknessMean",
    "surfaceTemperatureMean": "surfaceTemperatureMean",
    "airStressZonal": "airStressZonal",
    "airStressMeridional": "airStressMeridional",
}

# Expected variables in remapped files (lat-lon 360x180 output)
REMAPPED_DERIVED_VARS = [
    "sst", "sss", "surfaceHeatFluxTotal",
    "ssh", "shortWaveHeatFlux", "longWaveHeatFluxDown",
    "latentHeatFlux", "sensibleHeatFlux",
    "windStressZonal", "windStressMeridional",
]
REMAPPED_SEAICE_VARS = [
    "iceAreaTotal", "iceVolumeTotal", "snowVolumeTotal",
    "iceThicknessMean", "surfaceTemperatureMean",
    "airStressZonal", "airStressMeridional",
]
REMAPPED_DEPTH_COARSENED_BASES = [
    "temperatureCoarsened", "salinityCoarsened",
    "velocityZonalCoarsened", "velocityMeridionalCoarsened",
    "layerThicknessCoarsened",
]
REMAPPED_VERTREDUCE_VARS = [
    "oceanHeatContent", "freshwaterContent", "kineticEnergy",
]

# Physical range checks: (vmin, vmax)
# Ranges are deliberately generous to avoid false positives during spinup.
RANGE_CHECKS = {
    "temperatureCoarsened": (-5, 40),
    "salinityCoarsened": (0, 50),      # Red Sea/Persian Gulf can exceed 42 PSU
    "velocityZonalCoarsened": (-5, 5),
    "velocityMeridionalCoarsened": (-5, 5),
    "layerThicknessCoarsened": (0, None),
    "sst": (-5, 40),
    "sss": (0, 50),
    "surfaceHeatFluxTotal": (-2000, 2000),  # spinup transients can be large
    "iceAreaTotal": (0, 1),
    "iceVolumeTotal": (0, None),
    "snowVolumeTotal": (-1e-20, None),      # allow floating-point noise near zero
    "iceThicknessMean": (0, None),
    "surfaceTemperatureMean": (-50, 5),     # degC (not Kelvin!), Arctic ice can be -45C
    "oceanHeatContent": (None, None),
    "freshwaterContent": (-1000, 1000),
    "kineticEnergy": (0, None),
    "ssh": (-5, 5),
    "shortWaveHeatFlux": (0, 500),
    "longWaveHeatFluxDown": (0, 500),
    "latentHeatFlux": (-500, 500),
    "sensibleHeatFlux": (-300, 300),
    "windStressZonal": (-2, 2),
    "windStressMeridional": (-2, 2),
}

# Depth bounds for depth-latitude cross-section plots (metres)
DEPTH_BOUNDS = [0, 20, 30, 40, 50, 80, 110, 140, 170, 230, 410, 530,
                1020, 1080, 1720, 1980, 2820, 3380, 4620, 6380]

# Expected lat-lon grid dimensions for remapped output
EXPECTED_NLON = 360
EXPECTED_NLAT = 180

# -- Legacy vs FME variable mapping -------------------------------------------
# (legacy_var, fme_var, fme_stream, label, extract, cmap, vmin, vmax, units)
# extract: "last2d" = last timestep of (Time, nCells)
#          "surface" = last timestep, level 0 of (Time, nVertLevels, nCells)

LEGACY_FME_OCN_PAIRS = [
    ("timeCustom_avg_ssh", "ssh", "fmeDerivedFields",
     "SSH", "last2d", "RdBu_r", -2, 2, "m"),
    ("timeCustom_avg_activeTracers_temperature", "sst", "fmeDerivedFields",
     "SST", "surface", "RdYlBu_r", -2, 30, "degC"),
    ("timeCustom_avg_activeTracers_salinity", "sss", "fmeDerivedFields",
     "SSS", "surface", "viridis", 30, 40, "PSU"),
    ("timeCustom_avg_shortWaveHeatFlux", "shortWaveHeatFlux", "fmeDerivedFields",
     "SW Heat Flux", "last2d", "YlOrRd", 0, 400, "W/m2"),
    ("timeCustom_avg_longWaveHeatFluxDown", "longWaveHeatFluxDown", "fmeDerivedFields",
     "LW Heat Flux", "last2d", "inferno", 200, 450, "W/m2"),
    ("timeCustom_avg_sensibleHeatFlux", "sensibleHeatFlux", "fmeDerivedFields",
     "Sensible Heat", "last2d", "RdBu_r", -100, 100, "W/m2"),
    ("timeCustom_avg_windStressZonal", "windStressZonal", "fmeDerivedFields",
     "Wind Stress (zonal)", "last2d", "RdBu_r", -0.3, 0.3, "N/m2"),
    ("timeCustom_avg_windStressMeridional", "windStressMeridional", "fmeDerivedFields",
     "Wind Stress (merid)", "last2d", "RdBu_r", -0.2, 0.2, "N/m2"),
]

LEGACY_FME_ICE_PAIRS = [
    ("timeCustom_avg_iceAreaCell", "iceAreaTotal", "fmeSeaiceDerivedFields",
     "Ice Area", "last2d", "Blues", 0, 1, "fraction"),
    ("timeCustom_avg_iceVolumeCell", "iceVolumeTotal", "fmeSeaiceDerivedFields",
     "Ice Volume", "last2d", "viridis", 0, 5, "m"),
    ("timeCustom_avg_snowVolumeCell", "snowVolumeTotal", "fmeSeaiceDerivedFields",
     "Snow Volume", "last2d", "PuBu", 0, 2, "m"),
]


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


def _extract_field(arr, mode):
    """Extract a 1D (nCells) field from a multi-dimensional array.

    mode: 'last2d' = last timestep of (Time, nCells)
          'surface' = last timestep, level 0 of (Time, nVertLevels, nCells)
    """
    if arr is None:
        return None
    if mode == "surface":
        if arr.ndim == 3:      # (Time, nVertLevels, nCells)
            return arr[-1, 0, :]
        elif arr.ndim == 2:    # (nVertLevels, nCells) -- no time
            return arr[0, :]
    else:  # last2d
        if arr.ndim == 2:      # (Time, nCells)
            return arr[-1, :]
        elif arr.ndim == 1:    # (nCells)
            return arr[:]
    return None


def _diff_stats(leg, fme, fill_thresh=1e10):
    """Compute difference statistics between two 1D fields on the same grid.

    Returns dict with n_valid, legacy_mean, fme_mean, diff_mean, diff_rms,
    diff_max_abs, correlation -- or None if no valid overlap.
    """
    if leg is None or fme is None:
        return None
    a = leg.ravel().astype(float)
    b = fme.ravel().astype(float)
    if a.size != b.size:
        return None
    ok = ((np.abs(a) < fill_thresh) & np.isfinite(a) &
          (np.abs(b) < fill_thresh) & np.isfinite(b))
    if ok.sum() == 0:
        return None
    a, b = a[ok], b[ok]
    d = b - a
    s = {
        "n_valid": int(ok.sum()),
        "legacy_mean": float(a.mean()),
        "fme_mean": float(b.mean()),
        "diff_mean": float(d.mean()),
        "diff_rms": float(np.sqrt((d ** 2).mean())),
        "diff_max_abs": float(np.abs(d).max()),
    }
    if a.std() > 0 and b.std() > 0:
        s["correlation"] = float(np.corrcoef(a, b)[0, 1])
    else:
        s["correlation"] = float("nan")
    return s


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
    # Mask fill values so they render as transparent (not clipped to vmax)
    data = np.where((np.abs(data) < 1e10) & np.isfinite(data),
                    data.astype(float), np.nan)
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


def depth_latitude_section(ds, base_varname, n_levels, title, outdir, fname,
                           cmap="RdBu_r", vmin=None, vmax=None, units=""):
    """Plot a depth-latitude cross-section from per-level remapped variables.

    For each level k (0 to n_levels-1) reads {base_varname}_{k}, extracts the
    last timestep, computes zonal mean, and stacks into a (n_levels, nlat) array.
    Depth midpoints are computed from DEPTH_BOUNDS.
    """
    if not HAS_MPL:
        return None

    lat = get_var(ds, "lat")
    if lat is None:
        return None

    data = []
    for k in range(n_levels):
        vname = f"{base_varname}_{k}"
        arr = get_var(ds, vname)
        zm = _compute_zonal_mean(arr)
        if zm is None:
            return None
        data.append(zm)

    data = np.array(data)  # (n_levels, nlat)
    n_bounds = min(n_levels + 1, len(DEPTH_BOUNDS))
    depth_mid = np.array([(DEPTH_BOUNDS[k] + DEPTH_BOUNDS[k + 1]) / 2.0
                          for k in range(n_bounds - 1)])
    # Trim data to match available depth bounds
    n_plot = min(len(depth_mid), data.shape[0])
    depth_mid = depth_mid[:n_plot]
    data = data[:n_plot, :]

    fig, ax = plt.subplots(figsize=(10, 5))
    im = ax.pcolormesh(lat, depth_mid, data, cmap=cmap,
                       vmin=vmin, vmax=vmax, shading="auto")
    ax.invert_yaxis()
    ax.set_xlabel("Latitude")
    ax.set_ylabel("Depth (m)")
    ax.set_title(title, fontsize=11)
    plt.colorbar(im, ax=ax, label=units)
    plt.tight_layout()
    return savefig(fig, outdir, fname)


def fill_mask_map(data, lons, lats, title, outdir, fname):
    """Plot a binary valid/fill mask map.

    1 = valid (|data| < 1e10 and finite), 0 = fill/NaN.
    """
    if not HAS_MPL:
        return None

    mask = np.where((np.abs(data) < 1e10) & np.isfinite(data), 1.0, 0.0)
    cmap = mcolors.ListedColormap(["#d4e6f1", "#1b4f72"])
    return global_map(mask, lons, lats, title, cmap=cmap, vmin=0, vmax=1,
                      outdir=outdir, fname=fname, units="valid=1")


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


def trio_comparison(d_nat, lon_nat, lat_nat, d_rem, lon_rem, lat_rem,
                    title, outdir, fname, cmap="RdBu_r",
                    vmin=None, vmax=None, units=""):
    """Three-panel diagnostic: native tripcolor | remapped pcolormesh | zonal mean overlay.

    Produces a single figure with the most informative comparison layout.
    """
    if not HAS_MPL:
        return None
    from matplotlib.tri import Triangulation

    fig = plt.figure(figsize=(18, 4.5))
    is_cart = HAS_CARTOPY

    # --- Panel 1: native (tripcolor, subsampled) ---
    kw = dict(projection=ccrs.PlateCarree()) if is_cart else {}
    ax1 = fig.add_subplot(1, 3, 1, **kw)
    if is_cart:
        ax1.add_feature(cfeature.COASTLINE, linewidth=0.3)
        ax1.set_global()
    lon_fix = _fix_lon(lon_nat)
    max_pts = 50000
    rng = np.random.default_rng(42)
    if d_nat.size > max_pts:
        idx = rng.choice(d_nat.size, max_pts, replace=False)
    else:
        idx = np.arange(d_nat.size)
    # Mask fill for tripcolor data
    d_sub = np.where((np.abs(d_nat[idx]) < 1e10) & np.isfinite(d_nat[idx]),
                     d_nat[idx].astype(float), np.nan)
    try:
        tri = Triangulation(lon_fix[idx], lat_nat[idx])
        tri_lons = lon_fix[idx][tri.triangles]
        tri.set_mask((tri_lons.max(axis=1) - tri_lons.min(axis=1)) > 180)
        tkw = dict(transform=ccrs.PlateCarree()) if is_cart else {}
        im1 = ax1.tripcolor(tri, d_sub, cmap=cmap, vmin=vmin, vmax=vmax, **tkw)
    except Exception:
        tkw = dict(transform=ccrs.PlateCarree()) if is_cart else {}
        im1 = ax1.scatter(lon_fix[idx], lat_nat[idx], c=d_sub, s=0.3,
                          cmap=cmap, vmin=vmin, vmax=vmax, **tkw)
    ax1.set_title("Native (MPAS)", fontsize=10)
    plt.colorbar(im1, ax=ax1, shrink=0.6, label=units)

    # --- Panel 2: remapped (pcolormesh) ---
    kw = dict(projection=ccrs.PlateCarree()) if is_cart else {}
    ax2 = fig.add_subplot(1, 3, 2, **kw)
    if is_cart:
        ax2.add_feature(cfeature.COASTLINE, linewidth=0.3)
        ax2.set_global()
    d_rem_masked = np.where((np.abs(d_rem) < 1e10) & np.isfinite(d_rem),
                            d_rem.astype(float), np.nan)
    if lon_rem.ndim == 1:
        lons_g, lats_g = np.meshgrid(lon_rem, lat_rem)
    else:
        lons_g, lats_g = lon_rem, lat_rem
    tkw = dict(transform=ccrs.PlateCarree()) if is_cart else {}
    im2 = ax2.pcolormesh(lons_g, lats_g, d_rem_masked, cmap=cmap,
                         vmin=vmin, vmax=vmax, shading="auto", **tkw)
    ax2.set_title("Remapped (lat-lon)", fontsize=10)
    plt.colorbar(im2, ax=ax2, shrink=0.6, label=units)

    # --- Panel 3: zonal mean overlay ---
    ax3 = fig.add_subplot(1, 3, 3)
    # Remapped zonal mean
    zm_rem = _compute_zonal_mean(d_rem_masked[np.newaxis, :, :] if d_rem_masked.ndim == 2
                                 else d_rem_masked)
    lat_1d = lat_rem if lat_rem.ndim == 1 else lat_rem[:, 0]
    if zm_rem is not None:
        ax3.plot(lat_1d, zm_rem, "b-", linewidth=1.5, label="Remapped")
    # Native zonal mean (bin into latitude bands)
    lat_edges = np.linspace(-90, 90, 181)
    lat_centers = 0.5 * (lat_edges[:-1] + lat_edges[1:])
    ok_nat = (np.abs(d_nat) < 1e10) & np.isfinite(d_nat)
    zm_nat = np.full(180, np.nan)
    for b in range(180):
        mask = ok_nat & (lat_nat >= lat_edges[b]) & (lat_nat < lat_edges[b + 1])
        if mask.sum() > 0:
            zm_nat[b] = d_nat[mask].mean()
    ax3.plot(lat_centers, zm_nat, "r--", linewidth=1.2, alpha=0.8, label="Native")
    ax3.set_xlabel("Latitude")
    ax3.set_ylabel(units)
    ax3.set_title("Zonal Mean", fontsize=10)
    ax3.legend(fontsize=8)
    ax3.set_xlim(-90, 90)
    ax3.grid(True, alpha=0.3)

    fig.suptitle(title, fontsize=12, y=1.02)
    plt.tight_layout()
    path = os.path.join(outdir, fname)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


# All remapped fields with their plot parameters: (varname, cmap, vmin, vmax, units)
ALL_DERIVED_FIELD_PARAMS = [
    ("sst",                   "RdYlBu_r",  -2,   30,  "degC"),
    ("sss",                   "viridis",    30,   40,  "PSU"),
    ("ssh",                   "RdBu_r",     -2,    2,  "m"),
    ("surfaceHeatFluxTotal",  "RdBu_r",  -300,  300,  "W/m2"),
    ("shortWaveHeatFlux",     "YlOrRd",     0,  400,  "W/m2"),
    ("longWaveHeatFluxDown",  "inferno",  200,  450,  "W/m2"),
    ("latentHeatFlux",        "RdBu_r",  -300,  300,  "W/m2"),
    ("sensibleHeatFlux",      "RdBu_r",  -100,  100,  "W/m2"),
    ("windStressZonal",       "RdBu_r",  -0.3,  0.3,  "N/m2"),
    ("windStressMeridional",  "RdBu_r",  -0.2,  0.2,  "N/m2"),
]

ALL_SEAICE_FIELD_PARAMS = [
    ("iceAreaTotal",           "Blues",       0,    1,  "fraction"),
    ("iceVolumeTotal",         "viridis",     0,    5,  "m"),
    ("snowVolumeTotal",        "PuBu",        0,    2,  "m"),
    ("iceThicknessMean",       "viridis",     0,    5,  "m"),
    ("surfaceTemperatureMean", "RdBu_r",    210,  275,  "K"),
    ("airStressZonal",         "RdBu_r",     -1,    1,  "N/m2"),
    ("airStressMeridional",    "RdBu_r",     -1,    1,  "N/m2"),
]


def _remapped_with_zonal(d_rem, lon_rem, lat_rem, title, outdir, fname,
                         cmap="RdBu_r", vmin=None, vmax=None, units=""):
    """Two-panel diagnostic: remapped pcolormesh | zonal mean.

    Used when native-grid coordinates are not available.
    """
    if not HAS_MPL:
        return None

    fig = plt.figure(figsize=(14, 4.5))
    is_cart = HAS_CARTOPY

    d_masked = np.where((np.abs(d_rem) < 1e10) & np.isfinite(d_rem),
                        d_rem.astype(float), np.nan)
    lat_1d = lat_rem if lat_rem.ndim == 1 else lat_rem[:, 0]
    if lon_rem.ndim == 1:
        lons_g, lats_g = np.meshgrid(lon_rem, lat_rem)
    else:
        lons_g, lats_g = lon_rem, lat_rem

    # Panel 1: remapped map
    kw = dict(projection=ccrs.PlateCarree()) if is_cart else {}
    ax1 = fig.add_subplot(1, 2, 1, **kw)
    if is_cart:
        ax1.add_feature(cfeature.COASTLINE, linewidth=0.3)
        ax1.set_global()
    tkw = dict(transform=ccrs.PlateCarree()) if is_cart else {}
    im = ax1.pcolormesh(lons_g, lats_g, d_masked, cmap=cmap,
                        vmin=vmin, vmax=vmax, shading="auto", **tkw)
    ax1.set_title("Remapped (lat-lon)", fontsize=10)
    plt.colorbar(im, ax=ax1, shrink=0.6, label=units)

    # Panel 2: zonal mean
    ax2 = fig.add_subplot(1, 2, 2)
    zm = _compute_zonal_mean(d_masked[np.newaxis, :, :] if d_masked.ndim == 2
                             else d_masked)
    if zm is not None:
        ax2.plot(lat_1d, zm, "b-", linewidth=1.5)
    ax2.set_xlabel("Latitude")
    ax2.set_ylabel(units)
    ax2.set_title("Zonal Mean", fontsize=10)
    ax2.set_xlim(-90, 90)
    ax2.grid(True, alpha=0.3)

    fig.suptitle(title, fontsize=12, y=1.02)
    plt.tight_layout()
    path = os.path.join(outdir, fname)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def generate_trio_comparisons(rundir, fig_root, all_plots_by_comp):
    """Generate native-vs-remapped trio plots for all fields."""
    if not HAS_MPL:
        return

    trio_outdir = os.path.join(fig_root, "trio_comparisons")
    os.makedirs(trio_outdir, exist_ok=True)
    trio_plots = []

    def _process_component(nat_pattern, rem_pattern, field_params, component):
        nat_files = find_files(rundir, nat_pattern, exclude=".remapped.")
        rem_files = find_files(rundir, rem_pattern)
        if not rem_files:
            return
        ds_nat = safe_open(nat_files[0]) if nat_files else None
        ds_rem = safe_open(rem_files[-1])
        if has_time_zero(ds_rem):
            if ds_nat: close_ds(ds_nat)
            close_ds(ds_rem)
            return

        # Native coords may be missing in MPAS-O AM output
        lon_nat = lat_nat = None
        lon_nat_deg = lat_nat_deg = None
        if ds_nat is not None:
            lon_nat = get_var(ds_nat, "lonCell")
            lat_nat = get_var(ds_nat, "latCell")
            if lon_nat is not None and lat_nat is not None:
                lon_nat_deg = np.degrees(lon_nat)
                lat_nat_deg = np.degrees(lat_nat)

        lon_rem = get_var(ds_rem, "lon")
        lat_rem = get_var(ds_rem, "lat")
        if lon_rem is None or lat_rem is None:
            if ds_nat: close_ds(ds_nat)
            close_ds(ds_rem)
            return

        for vname, cmap, vm, vx, units in field_params:
            arr_rem = get_var(ds_rem, vname)
            if arr_rem is None:
                continue
            d_rem = arr_rem[-1] if arr_rem.ndim == 3 else arr_rem

            # Get native data if available
            d_nat = None
            if ds_nat is not None:
                arr_nat = get_var(ds_nat, vname)
                if arr_nat is not None:
                    d_nat = arr_nat[-1] if arr_nat.ndim > 1 else arr_nat

            fname = f"trio_{component.lower().replace('-','')}_{vname}.png"
            if d_nat is not None and lon_nat_deg is not None:
                # Full trio: native + remapped + zonal
                p = trio_comparison(
                    d_nat, lon_nat_deg, lat_nat_deg,
                    d_rem, lon_rem, lat_rem,
                    f"{component} {vname}",
                    trio_outdir, fname,
                    cmap=cmap, vmin=vm, vmax=vx, units=units)
            else:
                # Remapped-only: map + zonal mean (2 panels)
                p = _remapped_with_zonal(
                    d_rem, lon_rem, lat_rem,
                    f"{component} {vname}",
                    trio_outdir, fname,
                    cmap=cmap, vmin=vm, vmax=vx, units=units)
            if p:
                trio_plots.append(p)
                print(f"  Wrote {fname}")

        if ds_nat: close_ds(ds_nat)
        close_ds(ds_rem)

    print("\n=== Generating trio comparison plots ===")
    _process_component(
        "*.mpaso.hist.am.fmeDerivedFields.*.nc",
        "*.mpaso.hist.am.fmeDerivedFields.*.remapped.nc",
        ALL_DERIVED_FIELD_PARAMS, "MPAS-O")

    _process_component(
        "*.mpassi.hist.am.fmeSeaiceDerivedFields.*.nc",
        "*.mpassi.hist.am.fmeSeaiceDerivedFields.*.remapped.nc",
        ALL_SEAICE_FIELD_PARAMS, "MPAS-SI")

    if trio_plots:
        all_plots_by_comp["Native vs Remapped (all fields)"] = trio_plots


# -----------------------------------------------------------------------------
# Remapped file verification helpers
# -----------------------------------------------------------------------------

def check_remapped_grid(ds, label):
    """Verify that a remapped dataset has the expected lat-lon grid.

    Returns a list of issues (empty if OK).
    """
    issues = []
    dims = get_dims(ds)

    # Check for lon/lat dimensions
    has_lon = "lon" in dims
    has_lat = "lat" in dims
    if not has_lon:
        issues.append(f"  {label}: missing 'lon' dimension")
    if not has_lat:
        issues.append(f"  {label}: missing 'lat' dimension")
    if has_lon and dims["lon"] != EXPECTED_NLON:
        issues.append(f"  {label}: lon={dims['lon']}, expected {EXPECTED_NLON}")
    if has_lat and dims["lat"] != EXPECTED_NLAT:
        issues.append(f"  {label}: lat={dims['lat']}, expected {EXPECTED_NLAT}")

    # Check for Time dimension
    has_time = "Time" in dims or "time" in dims
    if not has_time:
        issues.append(f"  {label}: missing 'Time' dimension")

    if not issues:
        nlon = dims.get("lon", "?")
        nlat = dims.get("lat", "?")
        print(f"  {label}: lat-lon grid OK ({nlon}x{nlat})")

    return issues


def check_remapped_vars(ds, expected_vars, label):
    """Check that all expected variables are present in a remapped dataset."""
    issues = []
    varnames = get_varnames(ds)
    for var in expected_vars:
        if var not in varnames:
            issues.append(f"  {label}: missing variable '{var}'")
    present = [v for v in expected_vars if v in varnames]
    missing = [v for v in expected_vars if v not in varnames]
    if present:
        print(f"  {label}: {len(present)}/{len(expected_vars)} expected variables present")
    if missing:
        print(f"  {label}: MISSING variables: {missing}")
    return issues


# -----------------------------------------------------------------------------
# Per-component verification + visualization
# -----------------------------------------------------------------------------

def check_mpaso_depth_coarsening(rundir, outdir, verbose):
    print("\n=== MPAS-O Depth Coarsening (native) ===")
    issues = []
    plots = []
    fill_reports = []

    # Native files: *.mpaso.hist.am.fmeDepthCoarsening.*.nc, excluding .remapped.
    files = find_files(rundir, "*.mpaso.hist.am.fmeDepthCoarsening.*.nc",
                       exclude=".remapped.")
    if not files:
        print("  SKIP: no native fmeDepthCoarsening files found")
        return issues, plots, fill_reports

    print(f"  Found {len(files)} native file(s)")
    ds = safe_open(files[0])
    if has_time_zero(ds):
        print("  time dimension has size 0 -- skipping")
        close_ds(ds)
        return issues, plots, fill_reports

    # Auto-detect depth level count from dataset dimensions
    dims = get_dims(ds)
    n_ocn_levels = N_OCN_LAYERS  # None = auto-detect
    for dname in ["nFmeDepthLevels", "nFmeCoarsenLevels"]:
        if dname in dims:
            n_ocn_levels = dims[dname]
            break
    if n_ocn_levels is None:
        # Fallback: smallest non-Time, non-nCells dimension
        candidates = [v for k, v in dims.items()
                      if k not in ("Time", "time") and v < 100]
        n_ocn_levels = min(candidates) if candidates else 19
    print(f"  Depth levels detected: {n_ocn_levels}")

    for var, ace_name in ACE_OCN_COARSENED.items():
        arr = get_var(ds, var)
        if arr is None:
            issues.append(f"  native missing: {var} (ACE: {ace_name})")
            continue
        vmin, vmax = RANGE_CHECKS.get(var, (None, None))
        issues += check_range(arr, var, vmin, vmax)
        fill_reports.append(fill_nan_report(arr, var))
        if verbose:
            print(f"  {var} ({ace_name}): {summary_stats(arr)}")

        # Vertical profile plot (global mean across all cells, last timestep)
        if HAS_MPL and arr.ndim == 3:
            layer_profiles(arr[-1], f"{var} [{ace_name}]", outdir,
                           f"mpaso_depth_profile_{var}.png",
                           ylabel="Depth layer (0=surface)")

    # Surface maps -- last timestep
    if HAS_MPL:
        lon = get_var(ds, "lonCell")
        lat = get_var(ds, "latCell")
        if lon is not None and lat is not None:
            lon_deg = np.degrees(lon)
            lat_deg = np.degrees(lat)
            for var, cmap, vm, vx, units in [
                ("temperatureCoarsened", "RdBu_r", -2, 30, "degC"),
                ("salinityCoarsened", "viridis", 30, 40, "PSU"),
            ]:
                arr = get_var(ds, var)
                if arr is None or arr.ndim < 3:
                    continue
                data = arr[-1, 0, :]  # last time, surface layer, all cells
                p = global_map(data, lon_deg, lat_deg,
                               f"MPAS-O {var} layer 0 (native, last t)",
                               cmap=cmap, vmin=vm, vmax=vx,
                               outdir=outdir, fname=f"mpaso_depth_{var}_surf_native.png",
                               units=units)
                if p: plots.append(p)

    close_ds(ds)

    _report("MPAS-O depth coarsening (native)", issues)
    return issues, plots, fill_reports


def check_mpaso_depth_coarsening_remapped(rundir, outdir, verbose):
    print("\n=== MPAS-O Depth Coarsening (remapped) ===")
    issues = []
    plots = []
    fill_reports = []

    files = find_files(rundir, "*.mpaso.hist.am.fmeDepthCoarsening.*.remapped.nc")
    if not files:
        print("  SKIP: no remapped fmeDepthCoarsening files found")
        return issues, plots, fill_reports

    print(f"  Found {len(files)} remapped file(s)")
    ds = safe_open(files[0])
    if has_time_zero(ds):
        print("  time dimension has size 0 -- skipping")
        close_ds(ds)
        return issues, plots, fill_reports

    # Verify lat-lon grid
    issues += check_remapped_grid(ds, "depth coarsening remapped")

    # Auto-detect number of depth levels from variable names in the file
    varnames = get_varnames(ds)
    n_remap_levels = 0
    for k in range(100):
        if f"temperatureCoarsened_{k}" in varnames:
            n_remap_levels = k + 1
        else:
            break
    if n_remap_levels == 0:
        n_remap_levels = 19  # fallback
    print(f"  Depth levels detected: {n_remap_levels}")

    # Check for per-level variables
    expected_vars = []
    for base in REMAPPED_DEPTH_COARSENED_BASES:
        for k in range(n_remap_levels):
            expected_vars.append(f"{base}_{k}")
    issues += check_remapped_vars(ds, expected_vars, "depth coarsening remapped")

    # Range checks on per-level variables
    for base in REMAPPED_DEPTH_COARSENED_BASES:
        found_levels = []
        for k in range(n_remap_levels):
            vname = f"{base}_{k}"
            if vname in varnames:
                found_levels.append(k)
                arr = get_var(ds, vname)
                vmin, vmax = RANGE_CHECKS.get(base, (None, None))
                issues += check_range(arr, f"remapped/{vname}", vmin, vmax)
                fill_reports.append(fill_nan_report(arr, f"remapped/{vname}"))
                if verbose and k == 0:
                    print(f"  remapped/{vname}: {summary_stats(arr)}")
        if found_levels:
            print(f"  {base}: {len(found_levels)}/{n_remap_levels} levels in remapped file")

    # Maps of surface fields from remapped data
    if HAS_MPL:
        for base, cmap, vm, vx, units in [
            ("temperatureCoarsened", "RdBu_r", -2, 30, "degC"),
            ("salinityCoarsened", "viridis", 30, 40, "PSU"),
        ]:
            vname = f"{base}_0"
            p = latlon_map(ds, vname,
                           f"MPAS-O {base} layer 0 (remapped, last t)",
                           cmap=cmap, vmin=vm, vmax=vx,
                           outdir=outdir,
                           fname=f"mpaso_depth_{base}_0_remapped.png",
                           units=units)
            if p: plots.append(p)

        # Vertical profile from remapped per-level variables
        for base in ["temperatureCoarsened", "salinityCoarsened"]:
            means = []
            for k in range(n_remap_levels):
                vname = f"{base}_{k}"
                arr = get_var(ds, vname)
                v = valid_data(arr)
                means.append(v.mean() if v is not None and v.size > 0 else np.nan)
            if any(np.isfinite(m) for m in means):
                fig, ax = plt.subplots(figsize=(4, 5))
                ax.plot(means, np.arange(len(means)), "o-")
                ax.invert_yaxis()
                ax.set_xlabel(base)
                ax.set_ylabel("Depth layer (0=surface)")
                ax.set_title(f"Global-mean profile: {base} (remapped)")
                ax.grid(True, alpha=0.4)
                p = savefig(fig, outdir, f"mpaso_depth_profile_{base}_remapped.png")
                if p: plots.append(p)

        # Depth-latitude cross-sections
        for base, cmap, vm, vx, units in [
            ("temperatureCoarsened", "RdYlBu_r", -2, 30, "degC"),
            ("salinityCoarsened", "viridis", 33, 37, "PSU"),
            ("velocityZonalCoarsened", "RdBu_r", -0.5, 0.5, "m/s"),
            ("velocityMeridionalCoarsened", "RdBu_r", -0.2, 0.2, "m/s"),
        ]:
            p = depth_latitude_section(
                ds, base, n_remap_levels,
                f"Depth-Latitude: {base} (remapped)",
                outdir, f"mpaso_depth_lat_{base}.png",
                cmap=cmap, vmin=vm, vmax=vx, units=units)
            if p:
                plots.append(p)

    close_ds(ds)

    _report("MPAS-O depth coarsening (remapped)", issues)
    return issues, plots, fill_reports


def check_mpaso_derived(rundir, outdir, verbose):
    print("\n=== MPAS-O Derived Fields (native) ===")
    issues = []
    plots = []
    fill_reports = []

    files = find_files(rundir, "*.mpaso.hist.am.fmeDerivedFields.*.nc",
                       exclude=".remapped.")
    if not files:
        print("  SKIP: no native fmeDerivedFields files found")
        return issues, plots, fill_reports

    print(f"  Found {len(files)} native file(s)")
    ds = safe_open(files[0])
    if has_time_zero(ds):
        print("  time dimension has size 0 -- skipping")
        close_ds(ds)
        return issues, plots, fill_reports

    for var, ace_name in ACE_OCN_DERIVED.items():
        arr = get_var(ds, var)
        if arr is None:
            issues.append(f"  native missing: {var} (ACE: {ace_name})")
            continue
        vmin, vmax = RANGE_CHECKS.get(var, (None, None))
        issues += check_range(arr, var, vmin, vmax)
        fill_reports.append(fill_nan_report(arr, var))
        if verbose:
            print(f"  {var} ({ace_name}): {summary_stats(arr)}")

    # List all variables for reference
    varnames = get_varnames(ds)
    print(f"  Native file has {len(varnames)} variables")
    if verbose:
        print(f"  Variables: {varnames}")

    # Maps -- use LAST timestep
    if HAS_MPL:
        lon = get_var(ds, "lonCell")
        lat = get_var(ds, "latCell")
        if lon is not None and lat is not None:
            lon_deg = np.degrees(lon)
            lat_deg = np.degrees(lat)
            for var, cmap, vm, vx, units in [
                ("sst",  "RdBu_r",  -2, 30, "degC"),
                ("sss",  "viridis",  30, 40, "PSU"),
                ("surfaceHeatFluxTotal", "RdBu_r", -300, 300, "W/m2"),
            ]:
                arr = get_var(ds, var)
                if arr is None:
                    continue
                data = arr[-1] if arr.ndim > 1 else arr
                p = global_map(data, lon_deg, lat_deg,
                               f"MPAS-O {var} (native, last t)",
                               cmap=cmap, vmin=vm, vmax=vx,
                               outdir=outdir, fname=f"mpaso_derived_{var}_native.png",
                               units=units)
                if p: plots.append(p)

    close_ds(ds)

    _report("MPAS-O derived fields (native)", issues)
    return issues, plots, fill_reports


def check_mpaso_derived_remapped(rundir, outdir, verbose):
    print("\n=== MPAS-O Derived Fields (remapped) ===")
    issues = []
    plots = []
    fill_reports = []

    files = find_files(rundir, "*.mpaso.hist.am.fmeDerivedFields.*.remapped.nc")
    if not files:
        print("  SKIP: no remapped fmeDerivedFields files found")
        return issues, plots, fill_reports

    print(f"  Found {len(files)} remapped file(s)")
    ds = safe_open(files[0])
    if has_time_zero(ds):
        print("  time dimension has size 0 -- skipping")
        close_ds(ds)
        return issues, plots, fill_reports

    # Verify lat-lon grid
    issues += check_remapped_grid(ds, "derived fields remapped")

    # Check expected variables
    issues += check_remapped_vars(ds, REMAPPED_DERIVED_VARS, "derived fields remapped")

    # Range checks
    for var in REMAPPED_DERIVED_VARS:
        arr = get_var(ds, var)
        if arr is not None:
            vmin, vmax = RANGE_CHECKS.get(var, (None, None))
            issues += check_range(arr, f"remapped/{var}", vmin, vmax)
            fill_reports.append(fill_nan_report(arr, f"remapped/{var}"))
            if verbose:
                print(f"  remapped/{var}: {summary_stats(arr)}")

    # Maps from remapped data -- last timestep via latlon_map
    if HAS_MPL:
        for var, cmap, vm, vx, units in [
            ("sst", "RdBu_r", -2, 30, "degC"),
            ("sss", "viridis", 30, 40, "PSU"),
            ("surfaceHeatFluxTotal", "RdBu_r", -300, 300, "W/m2"),
            ("ssh", "RdBu_r", -2, 2, "m"),
            ("shortWaveHeatFlux", "YlOrRd", 0, 400, "W/m2"),
            ("longWaveHeatFluxDown", "inferno", 200, 450, "W/m2"),
            ("latentHeatFlux", "RdBu_r", -300, 300, "W/m2"),
            ("sensibleHeatFlux", "RdBu_r", -100, 100, "W/m2"),
            ("windStressZonal", "RdBu_r", -0.3, 0.3, "N/m2"),
            ("windStressMeridional", "RdBu_r", -0.2, 0.2, "N/m2"),
        ]:
            p = latlon_map(ds, var,
                           f"MPAS-O {var} (remapped, last t)",
                           cmap=cmap, vmin=vm, vmax=vx,
                           outdir=outdir,
                           fname=f"mpaso_derived_{var}_remapped.png",
                           units=units)
            if p: plots.append(p)

    close_ds(ds)

    _report("MPAS-O derived fields (remapped)", issues)
    return issues, plots, fill_reports


def check_mpaso_vertical_reduce(rundir, outdir, verbose):
    print("\n=== MPAS-O Vertical Reduction (native) ===")
    issues = []
    plots = []
    fill_reports = []

    files = find_files(rundir, "*.mpaso.hist.am.fmeVerticalReduce.*.nc",
                       exclude=".remapped.")
    if not files:
        print("  SKIP: no native fmeVerticalReduce files found")
        return issues, plots, fill_reports

    print(f"  Found {len(files)} native file(s)")
    ds = safe_open(files[0])
    if has_time_zero(ds):
        print("  time dimension has size 0 -- skipping")
        close_ds(ds)
        return issues, plots, fill_reports

    for var, (vmin, vmax) in [("oceanHeatContent", (None, None)),
                               ("freshwaterContent", (-1000, 1000)),
                               ("kineticEnergy", (0, None))]:
        arr = get_var(ds, var)
        if arr is None:
            issues.append(f"  native missing: {var}")
        else:
            issues += check_range(arr, var, vmin, vmax)
            fill_reports.append(fill_nan_report(arr, var))
            if verbose:
                print(f"  {var}: {summary_stats(arr)}")

    # Native scatter plots for vertical reduce fields
    if HAS_MPL:
        lon = get_var(ds, "lonCell")
        lat = get_var(ds, "latCell")
        if lon is not None and lat is not None:
            lon_deg = np.degrees(lon)
            lat_deg = np.degrees(lat)
            for var, cmap, vm, vx, units in [
                ("oceanHeatContent", "inferno", None, None, "J/m2"),
                ("freshwaterContent", "RdBu_r", -500, 500, "m"),
                ("kineticEnergy", "YlOrRd", 0, None, "J/m2"),
            ]:
                arr = get_var(ds, var)
                if arr is None:
                    continue
                data = arr[-1] if arr.ndim > 1 else arr
                p = global_map(data, lon_deg, lat_deg,
                               f"MPAS-O {var} (native, last t)",
                               cmap=cmap, vmin=vm, vmax=vx,
                               outdir=outdir,
                               fname=f"mpaso_vertreduce_{var}_native.png",
                               units=units)
                if p: plots.append(p)

        # Time series of global-mean OHC
        ohc = get_var(ds, "oceanHeatContent")
        if ohc is not None and ohc.ndim >= 2:
            area = get_var(ds, "areaCell")
            if area is not None:
                ts = [(ohc[t] * area).sum() / area.sum() for t in range(ohc.shape[0])]
                time_series(ts, "Timestep", "Global-mean OHC (native)",
                            "J/m2", outdir, "mpaso_ohc_timeseries_native.png")
                plots.append(os.path.join(outdir, "mpaso_ohc_timeseries_native.png"))

    close_ds(ds)

    _report("MPAS-O vertical reduction (native)", issues)
    return issues, plots, fill_reports


def check_mpaso_vertical_reduce_remapped(rundir, outdir, verbose):
    print("\n=== MPAS-O Vertical Reduction (remapped) ===")
    issues = []
    plots = []
    fill_reports = []

    files = find_files(rundir, "*.mpaso.hist.am.fmeVerticalReduce.*.remapped.nc")
    if not files:
        print("  SKIP: no remapped fmeVerticalReduce files found")
        return issues, plots, fill_reports

    print(f"  Found {len(files)} remapped file(s)")
    ds = safe_open(files[0])
    if has_time_zero(ds):
        print("  time dimension has size 0 -- skipping")
        close_ds(ds)
        return issues, plots, fill_reports

    # Verify lat-lon grid
    issues += check_remapped_grid(ds, "vertical reduce remapped")

    # Check expected variables
    issues += check_remapped_vars(ds, REMAPPED_VERTREDUCE_VARS, "vertical reduce remapped")

    # Range checks
    for var in REMAPPED_VERTREDUCE_VARS:
        arr = get_var(ds, var)
        if arr is not None:
            vmin, vmax = RANGE_CHECKS.get(var, (None, None))
            issues += check_range(arr, f"remapped/{var}", vmin, vmax)
            fill_reports.append(fill_nan_report(arr, f"remapped/{var}"))
            if verbose:
                print(f"  remapped/{var}: {summary_stats(arr)}")

    # Maps from remapped data -- last timestep via latlon_map
    if HAS_MPL:
        for var, cmap, vm, vx, units in [
            ("oceanHeatContent", "inferno", None, None, "J/m2"),
            ("freshwaterContent", "RdBu_r", -500, 500, "m"),
            ("kineticEnergy", "YlOrRd", 0, None, "J/m2"),
        ]:
            p = latlon_map(ds, var,
                           f"MPAS-O {var} (remapped, last t)",
                           cmap=cmap, vmin=vm, vmax=vx,
                           outdir=outdir,
                           fname=f"mpaso_vertreduce_{var}_remapped.png",
                           units=units)
            if p: plots.append(p)

    close_ds(ds)

    _report("MPAS-O vertical reduction (remapped)", issues)
    return issues, plots, fill_reports


def check_mpassi_derived(rundir, outdir, verbose):
    print("\n=== MPAS-SI Derived Fields (native) ===")
    issues = []
    plots = []
    fill_reports = []

    files = find_files(rundir, "*.mpassi.hist.am.fmeSeaiceDerivedFields.*.nc",
                       exclude=".remapped.")
    if not files:
        print("  SKIP: no native fmeSeaiceDerivedFields files found")
        return issues, plots, fill_reports

    print(f"  Found {len(files)} native file(s)")
    ds = safe_open(files[0])
    if has_time_zero(ds):
        print("  time dimension has size 0 -- skipping")
        close_ds(ds)
        return issues, plots, fill_reports

    for var, ace_name in ACE_ICE_DERIVED.items():
        arr = get_var(ds, var)
        if arr is None:
            issues.append(f"  native missing: {var} (ACE: {ace_name})")
            continue
        vmin, vmax = RANGE_CHECKS.get(var, (None, None))
        issues += check_range(arr, var, vmin, vmax)
        fill_reports.append(fill_nan_report(arr, var))
        if verbose:
            print(f"  {var} ({ace_name}): {summary_stats(arr)}")

    # Maps -- last timestep
    if HAS_MPL:
        lon = get_var(ds, "lonCell")
        lat = get_var(ds, "latCell")
        if lon is not None and lat is not None:
            lon_deg = np.degrees(lon)
            lat_deg = np.degrees(lat)
            for var, cmap, vm, vx, units in [
                ("iceAreaTotal",   "Blues",   0, 1, "fraction"),
                ("iceVolumeTotal", "viridis", 0, 5, "m"),
                ("surfaceTemperatureMean", "RdBu_r", 210, 275, "K"),
                ("airStressZonal", "RdBu_r", -1, 1, "N/m2"),
            ]:
                arr = get_var(ds, var)
                if arr is None:
                    continue
                data = arr[-1] if arr.ndim > 1 else arr
                p = global_map(data, lon_deg, lat_deg,
                               f"MPAS-SI {var} (native, last t)",
                               cmap=cmap, vmin=vm, vmax=vx,
                               outdir=outdir, fname=f"mpassi_{var}_native.png",
                               units=units)
                if p: plots.append(p)

    close_ds(ds)

    _report("MPAS-SI derived fields (native)", issues)
    return issues, plots, fill_reports


def check_mpassi_derived_remapped(rundir, outdir, verbose):
    print("\n=== MPAS-SI Derived Fields (remapped) ===")
    issues = []
    plots = []
    fill_reports = []

    files = find_files(rundir, "*.mpassi.hist.am.fmeSeaiceDerivedFields.*.remapped.nc")
    if not files:
        print("  SKIP: no remapped fmeSeaiceDerivedFields files found")
        return issues, plots, fill_reports

    print(f"  Found {len(files)} remapped file(s)")
    ds = safe_open(files[0])
    if has_time_zero(ds):
        print("  time dimension has size 0 -- skipping")
        close_ds(ds)
        return issues, plots, fill_reports

    # Verify lat-lon grid
    issues += check_remapped_grid(ds, "sea ice derived remapped")

    # Check expected variables
    issues += check_remapped_vars(ds, REMAPPED_SEAICE_VARS, "sea ice derived remapped")

    # Range checks
    for var in REMAPPED_SEAICE_VARS:
        arr = get_var(ds, var)
        if arr is not None:
            vmin, vmax = RANGE_CHECKS.get(var, (None, None))
            issues += check_range(arr, f"remapped/{var}", vmin, vmax)
            fill_reports.append(fill_nan_report(arr, f"remapped/{var}"))
            if verbose:
                print(f"  remapped/{var}: {summary_stats(arr)}")

    # Maps from remapped data -- last timestep via latlon_map
    if HAS_MPL:
        for var, cmap, vm, vx, units in [
            ("iceAreaTotal", "Blues", 0, 1, "fraction"),
            ("iceVolumeTotal", "viridis", 0, 5, "m"),
            ("surfaceTemperatureMean", "RdBu_r", 210, 275, "K"),
            ("snowVolumeTotal", "PuBu", 0, 2, "m"),
            ("airStressZonal", "RdBu_r", -1, 1, "N/m2"),
        ]:
            p = latlon_map(ds, var,
                           f"MPAS-SI {var} (remapped, last t)",
                           cmap=cmap, vmin=vm, vmax=vx,
                           outdir=outdir,
                           fname=f"mpassi_{var}_remapped.png",
                           units=units)
            if p: plots.append(p)

    close_ds(ds)

    _report("MPAS-SI derived fields (remapped)", issues)
    return issues, plots, fill_reports


def _area_weighted_mean(data_1d, area_1d, fill_thresh=1e10):
    """Area-weighted mean of a 1D field, excluding fill values."""
    ok = (np.abs(data_1d) < fill_thresh) & np.isfinite(data_1d) & np.isfinite(area_1d)
    if ok.sum() == 0:
        return float("nan")
    return float(np.sum(data_1d[ok] * area_1d[ok]) / np.sum(area_1d[ok]))


def _cosine_weighted_mean(data_2d, lat_1d, fill_thresh=1e10):
    """Cosine-of-latitude-weighted mean of a (lat, lon) field."""
    w = np.cos(np.deg2rad(lat_1d))
    w2d = np.broadcast_to(w[:, None], data_2d.shape)
    ok = (np.abs(data_2d) < fill_thresh) & np.isfinite(data_2d)
    if ok.sum() == 0:
        return float("nan")
    return float(np.sum(data_2d[ok] * w2d[ok]) / np.sum(w2d[ok]))


# -----------------------------------------------------------------------------
# Self-consistency checks
# -----------------------------------------------------------------------------

def check_remap_conservation(rundir, outdir, verbose):
    """Verify area-weighted global means match between native and remapped grids."""
    print("\n=== Remap Conservation Check ===")
    issues = []
    plots = []
    fill_reports = []

    nat_files = find_files(rundir, "*.mpaso.hist.am.fmeDerivedFields.*.nc",
                           exclude=".remapped.")
    rem_files = find_files(rundir, "*.mpaso.hist.am.fmeDerivedFields.*.remapped.nc")
    if not nat_files or not rem_files:
        print("  SKIP: need both native and remapped fmeDerivedFields files")
        return issues, plots, fill_reports

    ds_nat = safe_open(nat_files[0])
    ds_rem = safe_open(rem_files[0])
    if has_time_zero(ds_nat) or has_time_zero(ds_rem):
        print("  SKIP: time dimension empty")
        close_ds(ds_nat); close_ds(ds_rem)
        return issues, plots, fill_reports

    area = get_var(ds_nat, "areaCell")  # may be None in AM output
    lat = get_var(ds_rem, "lat")
    if lat is None:
        print("  SKIP: missing lat in remapped file")
        close_ds(ds_nat); close_ds(ds_rem)
        return issues, plots, fill_reports

    tol = 0.05  # 5% relative tolerance (remap + fill boundary effects)
    for vname in ["sst", "ssh"]:
        arr_nat = get_var(ds_nat, vname)
        arr_rem = get_var(ds_rem, vname)
        if arr_nat is None or arr_rem is None:
            continue
        d_nat = arr_nat[-1] if arr_nat.ndim > 1 else arr_nat
        d_rem = arr_rem[-1] if arr_rem.ndim == 3 else arr_rem

        # Native: area-weighted if areaCell available, else unweighted
        # (MPAS quasi-uniform mesh: unweighted is a good approximation)
        if area is not None:
            m_nat = _area_weighted_mean(d_nat, area)
        else:
            v = valid_data(d_nat)
            m_nat = float(v.mean()) if v is not None and v.size > 0 else float("nan")
        m_rem = _cosine_weighted_mean(d_rem, lat)
        denom = max(abs(m_nat), 1e-30)
        rel_diff = abs(m_nat - m_rem) / denom
        status = "PASS" if rel_diff < tol else "FAIL"
        if status == "FAIL":
            issues.append(f"  {vname}: native={m_nat:.4g} remap={m_rem:.4g} "
                          f"rel_diff={rel_diff:.2%}")
        print(f"  {vname:6s}  native={m_nat:+.4f}  remap={m_rem:+.4f}  "
              f"rel_diff={rel_diff:.2%}  {status}")

    close_ds(ds_nat); close_ds(ds_rem)
    _report("Remap conservation", issues)
    return issues, plots, fill_reports


def check_heat_flux_closure(rundir, outdir, verbose):
    """Verify surfaceHeatFluxTotal = SW + LW_down + latent + sensible."""
    print("\n=== Heat Flux Budget Closure ===")
    issues = []
    plots = []
    fill_reports = []

    files = find_files(rundir, "*.mpaso.hist.am.fmeDerivedFields.*.nc",
                       exclude=".remapped.")
    if not files:
        print("  SKIP: no native fmeDerivedFields files")
        return issues, plots, fill_reports

    ds = safe_open(files[0])
    if has_time_zero(ds):
        print("  SKIP: time dimension empty")
        close_ds(ds)
        return issues, plots, fill_reports

    total = get_var(ds, "surfaceHeatFluxTotal")
    sw = get_var(ds, "shortWaveHeatFlux")
    lw = get_var(ds, "longWaveHeatFluxDown")
    lat_hf = get_var(ds, "latentHeatFlux")
    sens = get_var(ds, "sensibleHeatFlux")

    if any(v is None for v in [total, sw, lw, lat_hf, sens]):
        print("  SKIP: missing flux fields")
        close_ds(ds)
        return issues, plots, fill_reports

    t = total[-1] if total.ndim > 1 else total
    computed = (sw[-1] if sw.ndim > 1 else sw) + \
               (lw[-1] if lw.ndim > 1 else lw) + \
               (lat_hf[-1] if lat_hf.ndim > 1 else lat_hf) + \
               (sens[-1] if sens.ndim > 1 else sens)

    ok = ((np.abs(t) < 1e10) & np.isfinite(t) &
          (np.abs(computed) < 1e10) & np.isfinite(computed))
    if ok.sum() == 0:
        print("  SKIP: no valid overlapping cells")
        close_ds(ds)
        return issues, plots, fill_reports

    residual = t[ok].astype(float) - computed[ok].astype(float)
    rms = float(np.sqrt((residual ** 2).mean()))
    max_abs = float(np.abs(residual).max())
    mean_res = float(residual.mean())

    # These should be BFB identical (same Fortran computation)
    status = "PASS" if max_abs < 1.0 else "FAIL"
    if status == "FAIL":
        issues.append(f"  Heat flux residual: rms={rms:.4g} max={max_abs:.4g}")
    print(f"  residual: mean={mean_res:.4g}  rms={rms:.4g}  max={max_abs:.4g}  "
          f"n_valid={ok.sum()}  {status}")

    # Residual map
    if HAS_MPL and max_abs > 1e-10:
        lon = get_var(ds, "lonCell")
        lat_c = get_var(ds, "latCell")
        if lon is not None and lat_c is not None:
            res_full = np.full_like(t, np.nan, dtype=float)
            res_full[ok] = residual
            p = global_map(res_full, np.degrees(lon), np.degrees(lat_c),
                           "Heat Flux Closure Residual (total - components)",
                           cmap="RdBu_r", vmin=-max_abs, vmax=max_abs,
                           outdir=outdir, fname="heat_flux_residual.png",
                           units="W/m2")
            if p:
                plots.append(p)

    close_ds(ds)
    _report("Heat flux closure", issues)
    return issues, plots, fill_reports


def check_temporal_consistency(rundir, outdir, verbose):
    """Check global-mean SST/SSH time series for discontinuities."""
    print("\n=== Temporal Consistency ===")
    issues = []
    plots = []
    fill_reports = []

    files = find_files(rundir, "*.mpaso.hist.am.fmeDerivedFields.*.nc",
                       exclude=".remapped.")
    if not files:
        print("  SKIP: no native fmeDerivedFields files")
        return issues, plots, fill_reports

    ts_data = {v: [] for v in ["sst", "ssh"]}
    for f in files:
        ds = safe_open(f)
        if has_time_zero(ds):
            close_ds(ds)
            continue
        area = get_var(ds, "areaCell")  # may be None in AM output
        dims = get_dims(ds)
        nt = dims.get("Time", dims.get("time", 1))
        for vname in ts_data:
            arr = get_var(ds, vname)
            if arr is None:
                continue
            for t in range(nt):
                d = arr[t] if arr.ndim > 1 else arr
                if area is not None:
                    m = _area_weighted_mean(d, area)
                else:
                    v = valid_data(d)
                    m = float(v.mean()) if v is not None and v.size > 0 else float("nan")
                ts_data[vname].append(m)
        close_ds(ds)

    for vname, series in ts_data.items():
        if len(series) < 2:
            print(f"  {vname}: {len(series)} timestep(s), skipping")
            continue
        s = np.array(series)
        deltas = np.abs(np.diff(s))
        thresh = 5.0 if vname == "sst" else 1.0
        jumps = np.where(deltas > thresh)[0]
        status = "PASS" if len(jumps) == 0 else "FAIL"
        if len(jumps) > 0:
            issues.append(f"  {vname}: {len(jumps)} jump(s) > {thresh}")
        print(f"  {vname:6s}  {len(series)} steps  range=[{s.min():.4f}, {s.max():.4f}]  "
              f"max_delta={deltas.max():.4g}  {status}")

        if HAS_MPL and len(series) > 1:
            time_series(series, "Output record", f"Global-mean {vname.upper()}",
                        vname.upper(), outdir, f"temporal_{vname}.png")
            plots.append(os.path.join(outdir, f"temporal_{vname}.png"))

    _report("Temporal consistency", issues)
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
# MPAS comparison figures (native vs remapped side-by-side, zonal means)
# -----------------------------------------------------------------------------

def generate_comparison_figures(rundir, fig_root, all_plots_by_comp):
    """Generate native-vs-remapped side-by-side comparison figures.

    Opens both native and remapped files for each MPAS component and
    produces tripcolor (native) vs pcolormesh (remapped) comparisons.
    Also generates zonal mean profiles for remapped lat-lon data.
    """
    if not HAS_MPL:
        return

    comp_outdir = os.path.join(fig_root, "comparisons")
    os.makedirs(comp_outdir, exist_ok=True)
    comp_plots = []

    # --- MPAS-O Derived: native vs remapped side-by-side ---
    nat_files = find_files(rundir, "*.mpaso.hist.am.fmeDerivedFields.*.nc",
                           exclude=".remapped.")
    rem_files = find_files(rundir, "*.mpaso.hist.am.fmeDerivedFields.*.remapped.nc")
    if nat_files and rem_files:
        print("\n=== Generating native vs remapped comparisons ===")
        ds_nat = safe_open(nat_files[0])
        ds_rem = safe_open(rem_files[-1])  # latest remapped
        if not has_time_zero(ds_nat) and not has_time_zero(ds_rem):
            lon_nat = get_var(ds_nat, "lonCell")
            lat_nat = get_var(ds_nat, "latCell")
            lon_rem = get_var(ds_rem, "lon")
            lat_rem = get_var(ds_rem, "lat")
            if lon_nat is not None and lat_nat is not None and \
               lon_rem is not None and lat_rem is not None:
                lon_nat_deg = np.degrees(lon_nat)
                lat_nat_deg = np.degrees(lat_nat)
                if lon_rem.ndim == 1 and lat_rem.ndim == 1:
                    lons_rem, lats_rem = np.meshgrid(lon_rem, lat_rem)
                else:
                    lons_rem, lats_rem = lon_rem, lat_rem

                for var, cmap, vm, vx, units_str in [
                    ("sst", "RdYlBu_r", -2, 30, "degC"),
                    ("sss", "viridis", 30, 40, "PSU"),
                    ("surfaceHeatFluxTotal", "RdBu_r", -300, 300, "W/m2"),
                ]:
                    arr_nat = get_var(ds_nat, var)
                    arr_rem = get_var(ds_rem, var)
                    if arr_nat is None or arr_rem is None:
                        continue
                    d_nat = arr_nat[-1] if arr_nat.ndim > 1 else arr_nat
                    d_rem = arr_rem[-1] if arr_rem.ndim == 3 else arr_rem
                    p = side_by_side_comparison(
                        d_nat, lon_nat_deg, lat_nat_deg, f"Native ({var})",
                        d_rem, lons_rem, lats_rem, f"Remapped ({var})",
                        f"MPAS-O {var}: Native vs Remapped",
                        comp_outdir, f"compare_mpaso_{var}.png",
                        cmap=cmap, vmin=vm, vmax=vx, units=units_str)
                    if p:
                        comp_plots.append(p)
                        print(f"  Wrote {os.path.basename(p)}")
        close_ds(ds_nat)
        close_ds(ds_rem)

    # --- MPAS-SI Derived: native vs remapped side-by-side ---
    nat_files = find_files(rundir, "*.mpassi.hist.am.fmeSeaiceDerivedFields.*.nc",
                           exclude=".remapped.")
    rem_files = find_files(rundir, "*.mpassi.hist.am.fmeSeaiceDerivedFields.*.remapped.nc")
    if nat_files and rem_files:
        ds_nat = safe_open(nat_files[0])
        ds_rem = safe_open(rem_files[-1])
        if not has_time_zero(ds_nat) and not has_time_zero(ds_rem):
            lon_nat = get_var(ds_nat, "lonCell")
            lat_nat = get_var(ds_nat, "latCell")
            lon_rem = get_var(ds_rem, "lon")
            lat_rem = get_var(ds_rem, "lat")
            if lon_nat is not None and lat_nat is not None and \
               lon_rem is not None and lat_rem is not None:
                lon_nat_deg = np.degrees(lon_nat)
                lat_nat_deg = np.degrees(lat_nat)
                if lon_rem.ndim == 1 and lat_rem.ndim == 1:
                    lons_rem, lats_rem = np.meshgrid(lon_rem, lat_rem)
                else:
                    lons_rem, lats_rem = lon_rem, lat_rem

                for var, cmap, vm, vx, units_str in [
                    ("iceAreaTotal", "Blues", 0, 1, "fraction"),
                    ("iceVolumeTotal", "viridis", 0, 5, "m"),
                    ("surfaceTemperatureMean", "RdBu_r", 210, 275, "K"),
                ]:
                    arr_nat = get_var(ds_nat, var)
                    arr_rem = get_var(ds_rem, var)
                    if arr_nat is None or arr_rem is None:
                        continue
                    d_nat = arr_nat[-1] if arr_nat.ndim > 1 else arr_nat
                    d_rem = arr_rem[-1] if arr_rem.ndim == 3 else arr_rem
                    p = side_by_side_comparison(
                        d_nat, lon_nat_deg, lat_nat_deg, f"Native ({var})",
                        d_rem, lons_rem, lats_rem, f"Remapped ({var})",
                        f"MPAS-SI {var}: Native vs Remapped",
                        comp_outdir, f"compare_mpassi_{var}.png",
                        cmap=cmap, vmin=vm, vmax=vx, units=units_str)
                    if p:
                        comp_plots.append(p)
                        print(f"  Wrote {os.path.basename(p)}")
        close_ds(ds_nat)
        close_ds(ds_rem)

    # --- Zonal means for remapped MPAS-O fields ---
    rem_derived = find_files(rundir, "*.mpaso.hist.am.fmeDerivedFields.*.remapped.nc")
    if rem_derived:
        ds = safe_open(rem_derived[-1])
        if not has_time_zero(ds):
            lat = get_var(ds, "lat")
            if lat is not None:
                profiles = []
                for var in ["sst", "sss"]:
                    zm = _compute_zonal_mean(get_var(ds, var))
                    if zm is not None:
                        profiles.append((zm, var))
                if profiles:
                    p = zonal_mean_plot(
                        profiles, lat,
                        "MPAS-O Zonal Mean SST/SSS (remapped)",
                        comp_outdir, "mpaso_zonal_mean_sst_sss.png",
                        ylabel="degC / PSU")
                    if p:
                        comp_plots.append(p)
                        print(f"  Wrote {os.path.basename(p)}")
        close_ds(ds)

    # --- Multi-panel heat flux summary ---
    rem_derived2 = find_files(rundir, "*.mpaso.hist.am.fmeDerivedFields.*.remapped.nc")
    if rem_derived2:
        ds = safe_open(rem_derived2[-1])
        if not has_time_zero(ds):
            lon = get_var(ds, "lon")
            lat = get_var(ds, "lat")
            if lon is not None and lat is not None:
                if lon.ndim == 1 and lat.ndim == 1:
                    lons_g, lats_g = np.meshgrid(lon, lat)
                else:
                    lons_g, lats_g = lon, lat
                flux_vars = [
                    ("shortWaveHeatFlux", "SW Down"),
                    ("longWaveHeatFluxDown", "LW Down"),
                    ("latentHeatFlux", "Latent"),
                    ("sensibleHeatFlux", "Sensible"),
                    ("surfaceHeatFluxTotal", "Total"),
                ]
                panels = []
                for vname, label in flux_vars:
                    arr = get_var(ds, vname)
                    if arr is None:
                        continue
                    data = arr[-1] if arr.ndim == 3 else arr
                    panels.append((data, lons_g, lats_g, label))
                if panels:
                    p = multi_panel_maps(
                        panels, "Surface Heat Flux Components (remapped)",
                        comp_outdir, "mpaso_heat_flux_multipanel.png",
                        ncols=3, cmap="RdBu_r", vmin=-400, vmax=400,
                        units="W/m2")
                    if p:
                        comp_plots.append(p)
                        print(f"  Wrote {os.path.basename(p)}")

                # --- Fill value mask maps ---
                for vname, label in [
                    ("sst", "SST ocean coverage"),
                    ("iceAreaTotal", "Ice extent mask"),
                    ("temperatureCoarsened_0", "Near-surface ocean mask"),
                ]:
                    arr = get_var(ds, vname)
                    if arr is None:
                        continue
                    data = arr[-1] if arr.ndim == 3 else arr
                    if data.ndim == 2:
                        p = fill_mask_map(data, lons_g, lats_g,
                                          f"Fill Mask: {label}",
                                          comp_outdir,
                                          f"fill_mask_{vname}.png")
                        if p:
                            comp_plots.append(p)
                            print(f"  Wrote {os.path.basename(p)}")
        close_ds(ds)

    # Fill mask for temperatureCoarsened_0 from depth coarsening file
    rem_depth = find_files(rundir, "*.mpaso.hist.am.fmeDepthCoarsening.*.remapped.nc")
    if rem_depth:
        ds = safe_open(rem_depth[-1])
        if not has_time_zero(ds):
            lon = get_var(ds, "lon")
            lat = get_var(ds, "lat")
            if lon is not None and lat is not None:
                if lon.ndim == 1 and lat.ndim == 1:
                    lons_g, lats_g = np.meshgrid(lon, lat)
                else:
                    lons_g, lats_g = lon, lat
                arr = get_var(ds, "temperatureCoarsened_0")
                if arr is not None:
                    data = arr[-1] if arr.ndim == 3 else arr
                    if data.ndim == 2:
                        p = fill_mask_map(data, lons_g, lats_g,
                                          "Fill Mask: Near-surface ocean (depth coarsening)",
                                          comp_outdir,
                                          "fill_mask_temperatureCoarsened_0.png")
                        if p:
                            comp_plots.append(p)
                            print(f"  Wrote {os.path.basename(p)}")
        close_ds(ds)

    # Fill mask for iceAreaTotal from sea ice file
    rem_ice = find_files(rundir, "*.mpassi.hist.am.fmeSeaiceDerivedFields.*.remapped.nc")
    if rem_ice:
        ds = safe_open(rem_ice[-1])
        if not has_time_zero(ds):
            lon = get_var(ds, "lon")
            lat = get_var(ds, "lat")
            if lon is not None and lat is not None:
                if lon.ndim == 1 and lat.ndim == 1:
                    lons_g, lats_g = np.meshgrid(lon, lat)
                else:
                    lons_g, lats_g = lon, lat
                arr = get_var(ds, "iceAreaTotal")
                if arr is not None:
                    data = arr[-1] if arr.ndim == 3 else arr
                    if data.ndim == 2:
                        p = fill_mask_map(data, lons_g, lats_g,
                                          "Fill Mask: Ice extent (sea ice)",
                                          comp_outdir,
                                          "fill_mask_iceAreaTotal_seaice.png")
                        if p:
                            comp_plots.append(p)
                            print(f"  Wrote {os.path.basename(p)}")
        close_ds(ds)

    if comp_plots:
        all_plots_by_comp["Comparisons"] = comp_plots


# -----------------------------------------------------------------------------
# Legacy vs FME cross-verification
# -----------------------------------------------------------------------------

def compare_legacy_vs_fme(legacy_rundir, fme_rundir, outdir, verbose):
    """Compare legacy timeSeriesStatsCustom output with FME online output.

    Both cases must use the same MPAS mesh (same nCells).  The legacy case
    produces time-averaged fields via the standard timeSeriesStatsCustom AM;
    the FME case produces instantaneous native-grid snapshots plus remapped
    time-averaged output.  This function compares the native-grid fields
    side-by-side and reports difference statistics.

    Returns (issues, plots, stats):
        issues: list of issue strings
        plots:  list of figure paths
        stats:  list of dicts (one per variable pair) for the HTML table
    """
    print("\n=== Legacy vs FME Cross-Verification ===")
    issues = []
    plots = []
    stats = []

    # --- locate files ---
    # Legacy: native-grid time-averaged (timeSeriesStatsCustom)
    leg_ocn = find_files(legacy_rundir,
                         "*.mpaso.hist.am.timeSeriesStatsCustom.*.nc")
    leg_ice = find_files(legacy_rundir,
                         "*.mpassi.hist.am.timeSeriesStatsCustom.*.nc")
    # FME: prefer remapped (time-averaged) for apples-to-apples temporal comparison;
    # fall back to native (instantaneous) if remapped not available
    fme_ocn = find_files(fme_rundir,
                         "*.mpaso.hist.am.fmeDerivedFields.*.remapped.nc")
    if not fme_ocn:
        fme_ocn = find_files(fme_rundir,
                             "*.mpaso.hist.am.fmeDerivedFields.*.nc",
                             exclude=".remapped.")
    fme_ice = find_files(fme_rundir,
                         "*.mpassi.hist.am.fmeSeaiceDerivedFields.*.remapped.nc")
    if not fme_ice:
        fme_ice = find_files(fme_rundir,
                             "*.mpassi.hist.am.fmeSeaiceDerivedFields.*.nc",
                             exclude=".remapped.")

    if not leg_ocn and not leg_ice:
        print("  SKIP: no legacy timeSeriesStatsCustom files found")
        return issues, plots, stats

    os.makedirs(outdir, exist_ok=True)

    # --- helper: process one pair list ---
    def _compare_pairs(ds_leg, ds_fme, pairs, component,
                       leg_lon, leg_lat, fme_lon, fme_lat):
        """Compare legacy (native, time-averaged) vs FME (remapped or native).

        Coordinates may be on different grids:
        - Legacy: 1D (nCells) from native MPAS
        - FME: 1D coordinate arrays from remapped lat-lon, or 1D (nCells) native
        """
        # Build FME meshgrid if it's a lat-lon grid
        if fme_lon is not None and fme_lat is not None:
            if fme_lon.ndim == 1 and fme_lat.ndim == 1 and \
               fme_lon.size != fme_lat.size:
                # lat-lon coordinate arrays (different sizes = regular grid)
                fme_lons_2d, fme_lats_2d = np.meshgrid(fme_lon, fme_lat)
            else:
                fme_lons_2d, fme_lats_2d = fme_lon, fme_lat
        else:
            fme_lons_2d = fme_lats_2d = None

        for (lvar, fvar, _stream, label, extract,
             cmap, vm, vx, units) in pairs:
            arr_leg = get_var(ds_leg, lvar)
            arr_fme = get_var(ds_fme, fvar)
            if arr_leg is None or arr_fme is None:
                if arr_leg is None and verbose:
                    print(f"  {component}: legacy missing {lvar}")
                if arr_fme is None and verbose:
                    print(f"  {component}: FME missing {fvar}")
                continue

            d_leg = _extract_field(arr_leg, extract)
            d_fme = _extract_field(arr_fme, "last2d")
            if d_leg is None or d_fme is None:
                continue

            # Statistics (only meaningful when grids match)
            if d_leg.size == d_fme.size:
                s = _diff_stats(d_leg, d_fme)
            else:
                # Different grids: compute global means for comparison
                v_leg = valid_data(d_leg)
                v_fme = valid_data(d_fme.ravel() if d_fme.ndim > 1 else d_fme)
                s = {
                    "n_valid": int(v_leg.size) if v_leg is not None else 0,
                    "legacy_mean": float(v_leg.mean()) if v_leg is not None and v_leg.size else float("nan"),
                    "fme_mean": float(v_fme.mean()) if v_fme is not None and v_fme.size else float("nan"),
                    "diff_mean": float("nan"),
                    "diff_rms": float("nan"),
                    "diff_max_abs": float("nan"),
                    "correlation": float("nan"),
                }
                s["diff_mean"] = s["fme_mean"] - s["legacy_mean"]

            if s is not None:
                s["component"] = component
                s["label"] = label
                s["legacy_var"] = lvar
                s["fme_var"] = fvar
                stats.append(s)
                corr_str = f"{s['correlation']:.6f}" if np.isfinite(s['correlation']) else "N/A (diff grids)"
                print(f"  {label:24s}  corr={corr_str}"
                      f"  rms={s['diff_rms']:.4g}"
                      f"  bias={s['diff_mean']:.4g}")

            # side-by-side plot
            if HAS_MPL and leg_lon is not None:
                fname = (f"xverify_{component.lower().replace('-','')}"
                         f"_{fvar}.png")
                fme_lon_plot = fme_lons_2d if fme_lons_2d is not None else leg_lon
                fme_lat_plot = fme_lats_2d if fme_lats_2d is not None else leg_lat
                p = side_by_side_comparison(
                    d_leg, leg_lon, leg_lat,
                    f"Legacy ({lvar.replace('timeCustom_avg_','')})",
                    d_fme, fme_lon_plot, fme_lat_plot,
                    f"FME ({fvar})",
                    f"{component} {label}: Legacy vs FME",
                    outdir, fname,
                    cmap=cmap, vmin=vm, vmax=vx, units=units)
                if p:
                    plots.append(p)

    # --- Ocean ---
    if leg_ocn and fme_ocn:
        ds_leg = safe_open(leg_ocn[0])
        ds_fme = safe_open(fme_ocn[0])
        if not has_time_zero(ds_leg) and not has_time_zero(ds_fme):
            # Legacy coordinates (native MPAS grid)
            lon_leg = get_var(ds_leg, "lonCell")
            lat_leg = get_var(ds_leg, "latCell")
            lon_leg_deg = np.degrees(lon_leg) if lon_leg is not None else None
            lat_leg_deg = np.degrees(lat_leg) if lat_leg is not None else None
            # FME coordinates (may be remapped lat-lon or native)
            fme_lon = get_var(ds_fme, "lon")
            fme_lat = get_var(ds_fme, "lat")
            if fme_lon is None:
                fme_lon = get_var(ds_fme, "lonCell")
                if fme_lon is not None:
                    fme_lon = np.degrees(fme_lon)
            if fme_lat is None:
                fme_lat = get_var(ds_fme, "latCell")
                if fme_lat is not None:
                    fme_lat = np.degrees(fme_lat)
            _compare_pairs(ds_leg, ds_fme, LEGACY_FME_OCN_PAIRS,
                           "MPAS-O", lon_leg_deg, lat_leg_deg,
                           fme_lon, fme_lat)
        else:
            print("  SKIP: time dimension empty in ocean files")
        close_ds(ds_leg)
        close_ds(ds_fme)
    elif not fme_ocn:
        print("  SKIP: no FME fmeDerivedFields files found")

    # --- Sea ice ---
    if leg_ice and fme_ice:
        ds_leg = safe_open(leg_ice[0])
        ds_fme = safe_open(fme_ice[0])
        if not has_time_zero(ds_leg) and not has_time_zero(ds_fme):
            lon_leg = get_var(ds_leg, "lonCell")
            lat_leg = get_var(ds_leg, "latCell")
            lon_leg_deg = np.degrees(lon_leg) if lon_leg is not None else None
            lat_leg_deg = np.degrees(lat_leg) if lat_leg is not None else None
            fme_lon = get_var(ds_fme, "lon")
            fme_lat = get_var(ds_fme, "lat")
            if fme_lon is None:
                fme_lon = get_var(ds_fme, "lonCell")
                if fme_lon is not None:
                    fme_lon = np.degrees(fme_lon)
            if fme_lat is None:
                fme_lat = get_var(ds_fme, "latCell")
                if fme_lat is not None:
                    fme_lat = np.degrees(fme_lat)
            _compare_pairs(ds_leg, ds_fme, LEGACY_FME_ICE_PAIRS,
                           "MPAS-SI", lon_leg_deg, lat_leg_deg,
                           fme_lon, fme_lat)
        else:
            print("  SKIP: time dimension empty in sea-ice files")
        close_ds(ds_leg)
        close_ds(ds_fme)
    elif not fme_ice:
        print("  SKIP: no FME fmeSeaiceDerivedFields files found")

    n_pairs = len(stats)
    if n_pairs == 0:
        issues.append("No overlapping variable pairs found for comparison")
    else:
        high_rms = [s for s in stats if s["diff_rms"] > 1e3]
        low_corr = [s for s in stats
                    if np.isfinite(s["correlation"]) and s["correlation"] < 0.9]
        if high_rms:
            for s in high_rms:
                issues.append(f"  {s['label']}: RMS diff = {s['diff_rms']:.4g}")
        if low_corr:
            for s in low_corr:
                issues.append(
                    f"  {s['label']}: correlation = {s['correlation']:.4f}")
        print(f"  Compared {n_pairs} variable pairs, "
              f"{len(issues)} issue(s)")

    return issues, plots, stats


def write_cross_verify_html(stats, legacy_rundir):
    """Generate HTML for the Legacy vs FME cross-verification section."""
    if not stats:
        return ""

    import html as html_mod

    html = '<h2 id="Cross_Verify">Legacy vs FME Cross-Verification</h2>\n'
    html += '<div class="card">\n'
    html += (f'<p>Comparing <code>{html_mod.escape(legacy_rundir)}</code> '
             f'(timeSeriesStatsCustom) with FME native output.</p>\n')
    html += ('<p style="font-size:0.85em;color:var(--text-muted)">'
             'Note: legacy output is time-averaged; FME native output is '
             'instantaneous. Differences are expected — focus on spatial '
             'patterns and correlation.</p>\n')
    html += '</div>\n'

    html += ('<table class="sortable" id="xverify-table"><thead><tr>'
             '<th onclick="sortTable(this,0)">Component</th>'
             '<th onclick="sortTable(this,1)">Field</th>'
             '<th onclick="sortTable(this,2)">Legacy Var</th>'
             '<th onclick="sortTable(this,3)">FME Var</th>'
             '<th onclick="sortTable(this,4)">Legacy Mean</th>'
             '<th onclick="sortTable(this,5)">FME Mean</th>'
             '<th onclick="sortTable(this,6)">Bias</th>'
             '<th onclick="sortTable(this,7)">RMS Diff</th>'
             '<th onclick="sortTable(this,8)">Correlation</th>'
             '</tr></thead><tbody>\n')

    for s in stats:
        corr = s["correlation"]
        corr_cls = "pass" if (np.isfinite(corr) and corr > 0.95) else (
            "fail" if (np.isfinite(corr) and corr < 0.8) else "")
        corr_str = f"{corr:.6f}" if np.isfinite(corr) else "N/A"
        lv = s["legacy_var"].replace("timeCustom_avg_", "")
        html += (f'<tr>'
                 f'<td>{s["component"]}</td>'
                 f'<td>{s["label"]}</td>'
                 f'<td><code>{lv}</code></td>'
                 f'<td><code>{s["fme_var"]}</code></td>'
                 f'<td>{s["legacy_mean"]:.4g}</td>'
                 f'<td>{s["fme_mean"]:.4g}</td>'
                 f'<td>{s["diff_mean"]:.4g}</td>'
                 f'<td>{s["diff_rms"]:.4g}</td>'
                 f'<td class="{corr_cls}">{corr_str}</td>'
                 f'</tr>\n')

    html += '</tbody></table>\n'
    return html


# -----------------------------------------------------------------------------
# File inventory
# -----------------------------------------------------------------------------

def collect_file_inventory(rundir):
    """Collect file inventory data and print summary.

    Returns list of (label, count, file_list) tuples for HTML rendering.
    """
    print("\n=== File Inventory ===")
    categories = [
        ("MPAS-O depth coarsening",    "*.mpaso.hist.am.fmeDepthCoarsening.*.nc",      ".remapped."),
        ("MPAS-O depth coarsening (R)","*.mpaso.hist.am.fmeDepthCoarsening.*.remapped.nc", None),
        ("MPAS-O derived fields",      "*.mpaso.hist.am.fmeDerivedFields.*.nc",        ".remapped."),
        ("MPAS-O derived fields (R)",  "*.mpaso.hist.am.fmeDerivedFields.*.remapped.nc", None),
        ("MPAS-O vertical reduce",     "*.mpaso.hist.am.fmeVerticalReduce.*.nc",       ".remapped."),
        ("MPAS-O vertical reduce (R)", "*.mpaso.hist.am.fmeVerticalReduce.*.remapped.nc", None),
        ("MPAS-SI derived fields",     "*.mpassi.hist.am.fmeSeaiceDerivedFields.*.nc", ".remapped."),
        ("MPAS-SI derived fields (R)", "*.mpassi.hist.am.fmeSeaiceDerivedFields.*.remapped.nc", None),
    ]
    total = 0
    inventory_data = []
    for label, pattern, exclude in categories:
        hits = find_files(rundir, pattern, exclude=exclude)
        total += len(hits)
        status = f"{len(hits)} file(s)" if hits else "NONE"
        print(f"  {label:38s} {status}")
        inventory_data.append((label, len(hits), hits))
    print(f"  {'TOTAL':38s} {total} FME-related file(s)")
    return inventory_data


# -----------------------------------------------------------------------------
# Reproducibility info
# -----------------------------------------------------------------------------

def _find_user_nl(rundir, component):
    """Search for user_nl_{component} in likely locations relative to rundir."""
    candidates = [
        os.path.join(rundir, f"user_nl_{component}"),
        os.path.join(os.path.dirname(rundir.rstrip("/")), f"user_nl_{component}"),
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

    # Read MPAS namelists from the case
    namelists = {}
    for comp in ["mpaso", "mpassi"]:
        nl_path = _find_user_nl(args.rundir, comp)
        if nl_path:
            try:
                with open(nl_path) as f:
                    namelists[comp] = f.read()
            except Exception:
                pass
    if namelists:
        info["namelists"] = namelists

    # Embed the script source for reproducibility
    try:
        with open(info["script_path"]) as f:
            info["script_source"] = f.read()
    except Exception:
        pass

    return info


def check_variable_coverage(rundir):
    """Check which ACE-required MPAS variables are present in the FME output.

    Scans all MPAS FME output files (depth coarsening, derived, sea ice).
    Atmospheric forcings are marked as 'from EAM' and not checked.
    Returns list of (ace_name, fme_name, category_str, found_bool, source) tuples.
    """
    # Collect all variable names across MPAS output files
    all_vars = set()
    patterns = [
        "*.mpaso.hist.am.fmeDepthCoarsening.*.remapped.nc",
        "*.mpaso.hist.am.fmeDerivedFields.*.remapped.nc",
        "*.mpassi.hist.am.fmeSeaiceDerivedFields.*.remapped.nc",
        # Also check native files in case remapped don't exist
        "*.mpaso.hist.am.fmeDepthCoarsening.*.nc",
        "*.mpaso.hist.am.fmeDerivedFields.*.nc",
        "*.mpassi.hist.am.fmeSeaiceDerivedFields.*.nc",
    ]
    for pattern in patterns:
        files = find_files(rundir, pattern)
        for f in files[:1]:  # only need first file per pattern
            ds = safe_open(f)
            all_vars.update(get_varnames(ds))
            close_ds(ds)

    forcing_set = set(ACE_MPAS_FORCING_NAMES)
    in_set = set(ACE_MPAS_IN_NAMES)
    out_set = set(ACE_MPAS_OUT_NAMES)

    coverage = []
    seen = set()
    for ace_name in ACE_MPAS_FORCING_NAMES + ACE_MPAS_IN_NAMES + ACE_MPAS_OUT_NAMES:
        if ace_name in seen:
            continue
        seen.add(ace_name)
        fme_name = mpas_ace_to_fme(ace_name)
        cats = []
        if ace_name in forcing_set:
            cats.append("forcing")
        if ace_name in in_set:
            cats.append("in")
        if ace_name in out_set:
            cats.append("out")
        category = ", ".join(cats)

        is_atm = ace_name in ACE_MPAS_ATM_FORCINGS
        if is_atm:
            source = "EAM"
            found = True  # don't check MPAS files for atmospheric vars
        else:
            source = "MPAS"
            found = fme_name in all_vars

        coverage.append((ace_name, fme_name, category, found, source))

    return coverage


def write_variable_coverage_html(coverage):
    """Generate HTML for the ACE MPAS Variable Requirements section."""
    if not coverage:
        return ""

    # Count only MPAS-sourced variables for the progress bar
    mpas_vars = [(a, f, c, found, s) for a, f, c, found, s in coverage if s == "MPAS"]
    n_found = sum(1 for _, _, _, f, _ in mpas_vars if f)
    n_total = len(mpas_vars)
    n_missing = n_total - n_found
    pct = 100 * n_found / n_total if n_total else 0

    html = '<h2 id="ACE_Variables">ACE MPAS Variable Requirements</h2>\n'

    html += '<div class="card">\n'
    html += f'<p style="font-size:1.05em"><strong>{n_found}/{n_total}</strong> '
    html += 'MPAS-sourced variables found in output '
    if n_missing:
        html += f'<span class="badge-fail">{n_missing} missing</span>'
    else:
        html += '<span class="badge-pass">all present</span>'
    html += f'</p>\n<p style="font-size:0.85em;color:var(--text-muted)">'
    html += f'+ {sum(1 for _,_,_,_,s in coverage if s=="EAM")} atmospheric forcings (from EAM)</p>\n'
    html += f'<div class="progress-bar"><div class="progress-fill" style="width:{pct:.0f}%"></div></div>\n'
    html += '</div>\n'

    # YAML spec
    yaml_text = generate_mpas_yaml_text()
    html += '<details><summary>YAML Variable Specification (ACE canonical names)</summary>\n'
    html += f'<pre class="yaml-block">{yaml_text}</pre>\n</details>\n'

    # Coverage table with filter
    html += '<details open><summary>Variable Coverage</summary>\n'
    html += '<input class="filter-input" type="text" placeholder="Filter variables..." '
    html += 'oninput="filterTable(this,\'mpas-coverage-table\')">\n'
    html += '<table class="sortable" id="mpas-coverage-table"><thead><tr>'
    html += '<th onclick="sortTable(this,0)">ACE Name</th>'
    html += '<th onclick="sortTable(this,1)">MPAS FME Name</th>'
    html += '<th onclick="sortTable(this,2)">Category</th>'
    html += '<th onclick="sortTable(this,3)">Source</th>'
    html += '<th onclick="sortTable(this,4)">Status</th>'
    html += '</tr></thead><tbody>\n'

    for ace_name, fme_name, category, found, source in coverage:
        if source == "EAM":
            cls = ""
            status = "from EAM"
        elif found:
            cls = "pass"
            status = "FOUND"
        else:
            cls = "fail"
            status = "MISSING"
        html += (f'<tr><td><code>{ace_name}</code></td>'
                 f'<td><code>{fme_name}</code></td>'
                 f'<td>{category}</td>'
                 f'<td>{source}</td>'
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
    html += f'<dt>Host</dt><dd><code>{repro_info.get("user", "")}@{repro_info.get("hostname", "")}</code></dd>\n'
    html += f'<dt>Generated</dt><dd>{repro_info["timestamp"]}</dd>\n'

    html += '</dl>\n'

    # Namelist contents
    for comp, content in repro_info.get("namelists", {}).items():
        html += f'<details><summary>user_nl_{comp}</summary>\n'
        html += f'<pre class="yaml-block">{content}</pre>\n</details>\n'

    # Full script source with line numbers
    if repro_info.get("script_source"):
        import html as html_mod
        src = repro_info["script_source"]
        lines = src.split("\n")
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

def build_executive_summary(all_issues, coverage, file_inventory_data):
    """Compute executive summary metrics for the HTML dashboard.

    Returns a dict with n_checks_passed, n_checks_total, n_vars_found,
    n_vars_total, total_files, total_mb, n_components.
    """
    n_checks_total = len(all_issues)
    n_checks_passed = sum(1 for v in all_issues.values() if not v)

    # Count only MPAS-sourced variables (not EAM forcings)
    n_vars_found = 0
    n_vars_total = 0
    if coverage:
        mpas_vars = [(a, f, c, found, s) for a, f, c, found, s in coverage
                     if s == "MPAS"]
        n_vars_found = sum(1 for _, _, _, f, _ in mpas_vars if f)
        n_vars_total = len(mpas_vars)

    total_files = 0
    total_bytes = 0
    for label, count, flist in file_inventory_data:
        total_files += count
        for fpath in flist:
            try:
                total_bytes += os.path.getsize(fpath)
            except OSError:
                pass
    total_mb = total_bytes / (1024 * 1024)

    n_components = n_checks_total

    return {
        "n_checks_passed": n_checks_passed,
        "n_checks_total": n_checks_total,
        "n_vars_found": n_vars_found,
        "n_vars_total": n_vars_total,
        "total_files": total_files,
        "total_mb": total_mb,
        "n_components": n_components,
    }


def write_executive_summary_html(summary):
    """Return HTML string with a grid of 4 metric cards."""
    if summary is None:
        return ""

    vars_detail = "All ACE MPAS variables found" if summary["n_vars_found"] == summary["n_vars_total"] \
        else f'{summary["n_vars_total"] - summary["n_vars_found"]} variable(s) missing'
    checks_detail = "All checks passed" if summary["n_checks_passed"] == summary["n_checks_total"] \
        else f'{summary["n_checks_total"] - summary["n_checks_passed"]} component(s) with issues'
    storage_detail = f'{summary["total_files"]} FME output files'
    comp_detail = "Components verified"

    html = '<div class="exec-summary">\n'

    html += ('  <div class="metric-card">\n'
             f'    <div class="metric-value">{summary["n_vars_found"]}/{summary["n_vars_total"]}</div>\n'
             f'    <div class="metric-label">Variable Coverage</div>\n'
             f'    <div class="metric-detail">{vars_detail}</div>\n'
             '  </div>\n')

    html += ('  <div class="metric-card">\n'
             f'    <div class="metric-value">{summary["n_checks_passed"]}/{summary["n_checks_total"]}</div>\n'
             f'    <div class="metric-label">Checks Passed</div>\n'
             f'    <div class="metric-detail">{checks_detail}</div>\n'
             '  </div>\n')

    html += ('  <div class="metric-card">\n'
             f'    <div class="metric-value">{summary["total_mb"]:.0f} MB</div>\n'
             f'    <div class="metric-label">Storage</div>\n'
             f'    <div class="metric-detail">{storage_detail}</div>\n'
             '  </div>\n')

    html += ('  <div class="metric-card">\n'
             f'    <div class="metric-value">{summary["n_components"]}</div>\n'
             f'    <div class="metric-label">Components</div>\n'
             f'    <div class="metric-detail">{comp_detail}</div>\n'
             '  </div>\n')

    html += '</div>\n'
    return html


# -----------------------------------------------------------------------------
# HTML index
# -----------------------------------------------------------------------------

def write_html_index(outdir, all_plots_by_comp, all_issues, file_inventory_data,
                     all_fill_reports, timing_summary=None, extra_html="",
                     repro_info=None, coverage=None, exec_summary=None):
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

    # Group plots: pair native/remapped side by side where possible
    fig_sections = ""
    # Identify paired components
    paired = [
        ("MPAS-O Depth Coarsening", "native", "remapped"),
        ("MPAS-O Derived Fields", "native", "remapped"),
        ("MPAS-O Vertical Reduction", "native", "remapped"),
        ("MPAS-SI Derived Fields", "native", "remapped"),
    ]
    shown_comps = set()

    for base_name, native_suffix, remap_suffix in paired:
        native_key = f"{base_name} ({native_suffix})"
        remap_key = f"{base_name} ({remap_suffix})"
        native_plots = all_plots_by_comp.get(native_key, [])
        remap_plots = all_plots_by_comp.get(remap_key, [])
        if not native_plots and not remap_plots:
            continue
        anchor = base_name.replace(" ", "_").replace("-", "_")
        fig_sections += f'<h3 id="{anchor}">{base_name}</h3>\n'
        fig_sections += '<div class="comparison">\n'
        if native_plots:
            fig_sections += '<div class="col"><h4>Native Grid</h4>\n'
            for p in native_plots:
                if p and os.path.exists(p):
                    rel = os.path.relpath(p, outdir)
                    fig_sections += (f'<div class="fig"><a href="{rel}">'
                                     f'<img src="{rel}"/></a>'
                                     f'<span>{os.path.basename(p)}</span></div>\n')
            fig_sections += '</div>\n'
        if remap_plots:
            fig_sections += '<div class="col"><h4>Remapped (lat-lon)</h4>\n'
            for p in remap_plots:
                if p and os.path.exists(p):
                    rel = os.path.relpath(p, outdir)
                    fig_sections += (f'<div class="fig"><a href="{rel}">'
                                     f'<img src="{rel}"/></a>'
                                     f'<span>{os.path.basename(p)}</span></div>\n')
            fig_sections += '</div>\n'
        fig_sections += '</div>\n'
        shown_comps.add(native_key)
        shown_comps.add(remap_key)

    # Show remaining (unpaired) components
    for comp, plots in all_plots_by_comp.items():
        if comp in shown_comps or not plots:
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
    for base_name, _, _ in paired:
        nav_items.append(base_name)
    for comp in all_plots_by_comp:
        if comp not in shown_comps and all_plots_by_comp[comp]:
            nav_items.append(comp)
    if timing_summary:
        nav_items.append("Timing")
    if repro_info:
        nav_items.append("Reproducibility")
    for item in nav_items:
        anchor = item.replace(" ", "_").replace("-", "_").replace("(", "").replace(")", "")
        display = item.replace("_", " ")
        nav_links += f'<a href="#{anchor}">{display}</a>\n'

    # Reproducibility section
    repro_html = write_repro_html(repro_info) if repro_info else ""

    # ACE variable coverage section
    coverage_html = write_variable_coverage_html(coverage) if coverage else ""

    # Executive summary
    exec_summary_html = write_executive_summary_html(exec_summary) if exec_summary else ""

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

    gen_timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    html = textwrap.dedent(f"""\
    <!DOCTYPE html><html data-theme="light"><head>
    <meta charset="utf-8"/>
    <meta name="viewport" content="width=device-width, initial-scale=1"/>
    <title>FME MPAS Output Verification</title>
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
        --code-bg: #f1f5f9; --hover-bg: #f8fafc;
        --yaml-bg: #1e293b; --yaml-text: #94a3b8;
        --stripe: #f8fafc; --border-radius: 8px;
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
        --code-bg: #0d1117; --hover-bg: #1c2128;
        --yaml-bg: #010409; --yaml-text: #7d8590;
        --stripe: #1c2128; --border-radius: 8px;
        --toggle-bg: #30363d; --toggle-knob: #58a6ff;
      }}
      * {{ box-sizing: border-box; margin: 0; padding: 0; }}
      html {{ scroll-behavior: smooth; }}
      body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto,
                           'Helvetica Neue', Arial, sans-serif;
             background: var(--bg); color: var(--text);
             line-height: 1.6; transition: background 0.3s, color 0.3s; }}
      .container {{ max-width: 1400px; margin: 0 auto; padding: 24px 32px; }}
      header {{ background: var(--header-bg); color: #fff; padding: 28px 32px;
                display: flex; justify-content: space-between; align-items: center;
                flex-wrap: wrap; gap: 12px; }}
      header .hdr-left h1 {{ font-size: 1.5em; font-weight: 700; letter-spacing: -0.02em; }}
      header .hdr-left .subtitle {{ font-size: 0.85em; opacity: 0.7; margin-top: 2px; }}
      header .status-pill {{ display: inline-block; padding: 6px 16px; border-radius: 20px;
                             font-weight: 700; font-size: 0.95em; }}
      .status-pass {{ background: rgba(22,163,74,0.2); color: #4ade80; }}
      .status-fail {{ background: rgba(220,38,38,0.2); color: #fca5a5; }}
      .hdr-right {{ display: flex; align-items: center; gap: 16px; }}
      .theme-toggle {{ background: none; border: 2px solid rgba(255,255,255,0.3);
                       color: #fff; padding: 6px 14px; border-radius: 20px;
                       cursor: pointer; font-size: 0.85em; font-weight: 500;
                       transition: all 0.2s; display: flex; align-items: center; gap: 6px; }}
      .theme-toggle:hover {{ border-color: #fff; background: rgba(255,255,255,0.1); }}
      .theme-toggle .icon {{ font-size: 1.1em; }}
      nav {{ background: var(--nav-bg); padding: 10px 16px; border-radius: var(--border-radius);
             box-shadow: var(--card-shadow); border: 1px solid var(--nav-border);
             margin: 20px 0; display: flex; flex-wrap: wrap; gap: 6px;
             position: sticky; top: 0; z-index: 100; transition: background 0.3s; }}
      nav a {{ color: var(--nav-link-text); text-decoration: none; padding: 6px 14px;
               background: var(--nav-link-bg); border-radius: 6px; font-size: 0.82em;
               font-weight: 500; transition: all 0.15s; white-space: nowrap; }}
      nav a:hover, nav a.active {{ background: var(--accent); color: #fff; }}
      h2 {{ color: var(--accent); margin-top: 36px; font-size: 1.25em; font-weight: 700;
            padding-bottom: 8px; border-bottom: 2px solid var(--accent); transition: color 0.3s; }}
      h3 {{ color: var(--text); margin-top: 24px; font-size: 1.05em; }}
      h4 {{ color: var(--text-muted); margin: 8px 0 4px; font-size: 0.95em; }}
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
      code {{ font-family: 'SF Mono', 'Fira Code', monospace; font-size: 0.88em;
              background: var(--code-bg); padding: 2px 6px; border-radius: 4px; }}
      .comparison {{ display: flex; gap: 24px; flex-wrap: wrap; }}
      .col {{ flex: 1; min-width: 420px; }}
      .gallery {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(420px, 1fr));
                  gap: 16px; margin-top: 12px; }}
      .fig {{ background: var(--card); border: 1px solid var(--card-border);
              border-radius: var(--border-radius); overflow: hidden; transition: all 0.25s; }}
      .fig:hover {{ box-shadow: 0 4px 16px rgba(0,0,0,0.12); transform: translateY(-2px); }}
      .fig img {{ width: 100%; display: block; cursor: zoom-in; }}
      .fig span {{ display: block; font-size: 0.75em; color: var(--text-muted);
                   padding: 8px 12px; border-top: 1px solid var(--card-border); font-family: monospace; }}
      .fig.wide {{ grid-column: span 2; }}
      @media (max-width: 900px) {{ .fig.wide {{ grid-column: span 1; }} }}
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
      .timing-container {{ background: var(--card); padding: 16px;
                           border-radius: var(--border-radius);
                           border: 1px solid var(--card-border);
                           overflow-x: auto; max-height: 500px; overflow-y: auto; }}
      .timing-table td {{ font-family: monospace; font-size: 0.83em;
                          padding: 4px 10px; white-space: nowrap; }}
      .yaml-block {{ background: var(--yaml-bg); color: var(--yaml-text);
                     padding: 16px 20px; border-radius: var(--border-radius);
                     font-family: 'SF Mono', 'Fira Code', monospace; font-size: 0.83em;
                     overflow-x: auto; max-height: 400px; overflow-y: auto;
                     white-space: pre; line-height: 1.5; border: 1px solid var(--card-border); }}
      .repro-block {{ background: var(--card); padding: 20px 24px;
                      border-radius: var(--border-radius);
                      border: 1px solid var(--card-border); box-shadow: var(--card-shadow); }}
      .repro-block dl {{ margin: 0; display: grid; grid-template-columns: auto 1fr;
                         gap: 4px 16px; align-items: baseline; }}
      .repro-block dt {{ font-weight: 600; color: var(--accent); font-size: 0.88em;
                         white-space: nowrap; }}
      .repro-block dd {{ margin: 0; font-family: monospace; font-size: 0.83em;
                         color: var(--text-muted); word-break: break-all; }}
      .src-block {{ background: var(--yaml-bg); color: var(--yaml-text);
                    padding: 0; border-radius: var(--border-radius);
                    font-family: 'SF Mono', 'Fira Code', monospace; font-size: 0.8em;
                    overflow-x: auto; max-height: 600px; overflow-y: auto;
                    white-space: pre; line-height: 1.55;
                    border: 1px solid var(--card-border); tab-size: 4; }}
      .src-line {{ display: block; padding: 0 16px 0 0; }}
      .src-line:hover {{ background: rgba(88,166,255,0.08); }}
      .src-line:target {{ background: rgba(88,166,255,0.18); }}
      .src-ln {{ display: inline-block; width: 4.5em; text-align: right;
                 color: rgba(139,148,158,0.4); padding: 0 12px 0 12px;
                 user-select: none; border-right: 1px solid var(--card-border);
                 margin-right: 12px; }}
      footer {{ margin-top: 48px; padding: 20px 32px;
                border-top: 1px solid var(--card-border);
                color: var(--text-muted); font-size: 0.8em;
                background: var(--card); text-align: center; transition: all 0.3s; }}
      .lightbox {{ display:none; position:fixed; top:0; left:0; width:100%; height:100%;
                   background:rgba(0,0,0,0.92); z-index:1000; cursor:zoom-out;
                   align-items:center; justify-content:center;
                   backdrop-filter: blur(4px); }}
      .lightbox.active {{ display:flex; }}
      .lightbox img {{ max-width:94%; max-height:94%; border-radius:8px;
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
      .card {{ background: var(--card); border: 1px solid var(--card-border);
               border-radius: var(--border-radius); box-shadow: var(--card-shadow);
               padding: 20px; margin-bottom: 16px; transition: all 0.3s; }}
      .progress-bar {{ height: 8px; background: #e2e8f0; border-radius: 4px;
                       overflow: hidden; margin: 8px 0 16px; }}
      [data-theme="dark"] .progress-bar {{ background: #30363d; }}
      .progress-fill {{ height: 100%; border-radius: 4px; transition: width 0.6s ease;
                        background: linear-gradient(90deg, #16a34a, #22d3ee); }}
      .badge-pass {{ display: inline-block; padding: 2px 10px; border-radius: 12px;
                     background: var(--pass-bg); color: var(--pass);
                     font-weight: 600; font-size: 0.82em; }}
      .badge-fail {{ display: inline-block; padding: 2px 10px; border-radius: 12px;
                     background: var(--fail-bg); color: var(--fail);
                     font-weight: 600; font-size: 0.82em; }}
      .filter-input {{ width: 100%; max-width: 400px; padding: 8px 14px;
                       border: 1px solid var(--card-border); border-radius: 6px;
                       background: var(--card); color: var(--text);
                       font-size: 0.9em; margin-bottom: 12px;
                       transition: all 0.2s; outline: none; }}
      .filter-input:focus {{ border-color: var(--accent);
                             box-shadow: 0 0 0 3px rgba(88,166,255,0.15); }}
      .filter-input::placeholder {{ color: var(--text-muted); }}
      .sortable th {{ cursor:pointer; user-select:none; position:relative; padding-right:22px; }}
      .sortable th::after {{ content:'\\2195'; position:absolute; right:6px; top:50%;
                             transform:translateY(-50%); opacity:0.4; font-size:0.85em; }}
      .sortable th:hover {{ background: var(--th-hover); }}
      .exec-summary {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
                       gap: 16px; margin: 20px 0; }}
      .metric-card {{ background: var(--card); border: 1px solid var(--card-border);
                      border-radius: var(--border-radius); padding: 20px; text-align: center;
                      box-shadow: var(--card-shadow); transition: all 0.3s; }}
      .metric-value {{ font-size: 2em; font-weight: 800; color: var(--accent); }}
      .metric-label {{ font-size: 0.9em; font-weight: 600; color: var(--text); margin-top: 4px; }}
      .metric-detail {{ font-size: 0.78em; color: var(--text-muted); margin-top: 4px; }}
      @keyframes fadeIn {{ from {{ opacity:0; transform:translateY(8px); }}
                           to   {{ opacity:1; transform:translateY(0); }} }}
      table, details {{ animation: fadeIn 0.3s ease; }}
    </style>
    <script>
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
    document.addEventListener('DOMContentLoaded', function() {{
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
      function openLB(idx) {{ curIdx = idx; lbImg.src = allFigs[idx]; lb.classList.add('active'); }}
      function closeLB() {{ lb.classList.remove('active'); }}
      function navLB(dir) {{ curIdx = (curIdx + dir + allFigs.length) % allFigs.length; lbImg.src = allFigs[curIdx]; }}
      document.querySelectorAll('.fig a').forEach(function(a, i) {{
        allFigs.push(a.href);
        a.addEventListener('click', function(e) {{ e.preventDefault(); openLB(i); }});
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

      /* Sortable tables */
      function sortTable(th, col) {{
        var table = th.closest('table');
        var tbody = table.querySelector('tbody');
        if (!tbody) return;
        var rows = Array.from(tbody.querySelectorAll('tr'));
        var asc = th.dataset.sort !== 'asc';
        th.dataset.sort = asc ? 'asc' : 'desc';
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
        table.querySelectorAll('tbody tr').forEach(function(row) {{
          row.style.display = row.textContent.toLowerCase().indexOf(filter) > -1 ? '' : 'none';
        }});
      }}

      /* Active nav highlighting */
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
        <h1>FME MPAS Online Output Verification</h1>
        <div class="subtitle">MPAS-Ocean &amp; MPAS-Sea Ice &mdash; Full Model Emulation Diagnostics</div>
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
      <strong>verify_mpas.py</strong> &mdash; FME MPAS Online Output Verification for E3SM
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
        description="Verify MPAS-Ocean and MPAS-SeaIce FME online output "
                    "and produce diagnostic figures.")
    parser.add_argument("--rundir", required=True,
                        help="CIME RUNDIR containing MPAS FME output files")
    parser.add_argument("--outdir", required=True,
                        help="Output directory for figures and HTML index")
    parser.add_argument("--verbose", "-v", action="store_true")
    parser.add_argument("--legacy-rundir", default=None,
                        help="Legacy FME output rundir (timeSeriesStatsCustom) "
                             "for cross-verification against FME online output")
    args = parser.parse_args()

    if not os.path.isdir(args.rundir):
        sys.exit(f"ERROR: rundir not found: {args.rundir}")

    os.makedirs(args.outdir, exist_ok=True)
    print(f"MPAS FME Verification  |  rundir: {args.rundir}")
    print(f"Output directory       :  {args.outdir}")
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

    run_check("MPAS-O Depth Coarsening (native)", check_mpaso_depth_coarsening, "mpaso_native",
              args.rundir, args.verbose)
    run_check("MPAS-O Depth Coarsening (remapped)", check_mpaso_depth_coarsening_remapped, "mpaso_remapped",
              args.rundir, args.verbose)
    run_check("MPAS-O Derived Fields (native)", check_mpaso_derived, "mpaso_native",
              args.rundir, args.verbose)
    run_check("MPAS-O Derived Fields (remapped)", check_mpaso_derived_remapped, "mpaso_remapped",
              args.rundir, args.verbose)
    run_check("MPAS-O Vertical Reduction (native)", check_mpaso_vertical_reduce, "mpaso_native",
              args.rundir, args.verbose)
    run_check("MPAS-O Vertical Reduction (remapped)", check_mpaso_vertical_reduce_remapped, "mpaso_remapped",
              args.rundir, args.verbose)
    run_check("MPAS-SI Derived Fields (native)", check_mpassi_derived, "mpassi_native",
              args.rundir, args.verbose)
    run_check("MPAS-SI Derived Fields (remapped)", check_mpassi_derived_remapped, "mpassi_remapped",
              args.rundir, args.verbose)

    # Self-consistency checks
    run_check("Remap Conservation", check_remap_conservation, "self_checks",
              args.rundir, args.verbose)
    run_check("Heat Flux Closure", check_heat_flux_closure, "self_checks",
              args.rundir, args.verbose)
    run_check("Temporal Consistency", check_temporal_consistency, "self_checks",
              args.rundir, args.verbose)

    # Cross-component comparison figures (native vs remapped side-by-side,
    # zonal means)
    generate_comparison_figures(args.rundir, fig_root, all_plots_by_comp)

    # Trio comparison plots: native | remapped | zonal mean for all fields
    generate_trio_comparisons(args.rundir, fig_root, all_plots_by_comp)

    # Timing summary
    timing_summary = read_timing_summary(args.rundir)
    if timing_summary:
        print("\n=== Performance Timing ===")
        for line in timing_summary[:10]:
            print(f"  {line}")
        if len(timing_summary) > 10:
            print(f"  ... ({len(timing_summary) - 10} more lines in HTML report)")

    # ACE variable coverage
    print("\n=== ACE MPAS Variable Coverage ===")
    coverage = check_variable_coverage(args.rundir)
    if coverage:
        mpas_vars = [(a, f, c, found, s) for a, f, c, found, s in coverage if s == "MPAS"]
        n_found = sum(1 for _, _, _, f, _ in mpas_vars if f)
        n_missing = len(mpas_vars) - n_found
        print(f"  {n_found}/{len(mpas_vars)} MPAS-sourced variables found in output")
        if n_missing:
            for ace_name, fme_name, _, found, src in coverage:
                if src == "MPAS" and not found:
                    print(f"  MISSING: {ace_name} (MPAS: {fme_name})")

    # Reproducibility info
    repro_info = collect_repro_info(args)

    # Legacy vs FME cross-verification
    extra_html = ""
    if args.legacy_rundir:
        if not os.path.isdir(args.legacy_rundir):
            print(f"WARNING: legacy-rundir not found: {args.legacy_rundir}")
        else:
            xv_outdir = os.path.join(fig_root, "cross_verify")
            os.makedirs(xv_outdir, exist_ok=True)
            xv_issues, xv_plots, xv_stats = compare_legacy_vs_fme(
                args.legacy_rundir, args.rundir, xv_outdir, args.verbose)
            all_issues["Legacy vs FME"] = xv_issues
            all_plots_by_comp["Cross-Verification"] = xv_plots
            extra_html = write_cross_verify_html(xv_stats, args.legacy_rundir)

    # Executive summary
    exec_summary = build_executive_summary(all_issues, coverage,
                                           file_inventory_data)

    write_html_index(args.outdir, all_plots_by_comp, all_issues,
                     file_inventory_data, all_fill_reports,
                     timing_summary=timing_summary,
                     extra_html=extra_html,
                     repro_info=repro_info,
                     coverage=coverage,
                     exec_summary=exec_summary)

    n_total = sum(len(v) for v in all_issues.values())
    print("\n" + "=" * 70)
    print(f"RESULT: {'ALL CHECKS PASSED' if n_total == 0 else f'{n_total} ISSUES FOUND'}")
    sys.exit(0 if n_total == 0 else 1)


if __name__ == "__main__":
    main()
