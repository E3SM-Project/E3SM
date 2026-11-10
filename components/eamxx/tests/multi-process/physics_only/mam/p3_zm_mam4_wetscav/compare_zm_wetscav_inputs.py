#!/usr/bin/env python3
"""
Compare the inputs that the EAMxx MAM4 wetscav interface passes to
mam4::wetdep::aero_model_wetdep with the inputs of the mam4xx validation
test aero_model_wetdep_convproc.

The EAMxx side comes from the netCDF output of the p3_zm_mam4_wetscav test
(see output.yaml). The script applies the same unit conversions as
eamxx_mam_wetscav_process_interface.cpp, so both sides are compared in the
units mam4xx expects (mb, mb/s, 0-based indices with exclusive kbot).

When the validation file is the Fortran data file (mam_<test>.py), the
script also compares the wet deposition outputs aerdepwetis and aerdepwetcw
per species with the Fortran values (and, optionally, with the mam4xx
standalone result).

The two sides are different columns and states, so values are not expected
to match. The script checks orders of magnitude, signs, zero patterns, value
ranges and the ktop/kbot index convention.

Usage:
  compare_zm_wetscav_inputs.py EAMXX_OUTPUT.nc VALIDATION_FILE
                               [--mam4xx MAM4XX_OUTPUT.py]
                               [--time IDX] [--report FILE.md]

VALIDATION_FILE is one of (from mam_x_validation/aero_model/):
  mam_aero_model_wetdep_convproc_ts_3.py   inputs + Fortran outputs (preferred)
  aero_model_wetdep_convproc_ts_3.yaml     inputs only
MAM4XX_OUTPUT.py is the file the mam4xx validation test writes
(mam4xx_aero_model_wetdep_convproc_ts_3.py in its build directory).

Using output from a full EAMxx case (e.g. ne30pg2):
  - The case must run ZM with mam4_wetscav convproc_do_aer: true.
  - Request the fields listed in output.yaml in the case's output stream,
    with averaging_type: instant, on the physics grid (no horizontal or
    vertical remap). Averaged ktop/kbot/zm_activity are meaningless; the
    script warns if the output looks averaged.
  - Pick a time slice with --time. Only that slice is read, so files with
    many slices are fine.
  - ZM, P3 and wetscav output fields (zm_*, evap_dq_out, nevapr,
    precip_total_tend, aerdepwetis, aerdepwetcw) are exactly what wetscav
    used or produced. In a full case, processes after wetscav (radiation,
    dynamics, ...) also change T_mid, p_mid, pseudo_density and qv/qc/qi/nc/ni
    before output, so those are close to, not exactly, the wetscav inputs.
  - Different level counts (e.g. 128 vs 72) are handled: layer thicknesses
    are rescaled before comparing.
"""

import argparse
import ast
import math
import sys
import warnings

import numpy as np
import yaml

GRAVITY = 9.80616   # mam4::Constants::gravity
PA_TO_MB = 0.01
KGM2S_TO_MBS = GRAVITY * PA_TO_MB   # same as the wetscav interface
DP1, DP2 = 0.018, 500.0             # dp_frac formula in the wetscav interface

# Known unit factors, used to label a magnitude mismatch
KNOWN_FACTORS = {
    "100": 100.0,
    "1/100": 0.01,
    "g": GRAVITY,
    "1/g": 1.0 / GRAVITY,
    "100/g": 100.0 / GRAVITY,
    "g/100": GRAVITY / 100.0,
    "1000": 1000.0,
    "1/1000": 0.001,
}

# key: (units mam4xx expects, expected sign: '+', '-', or None, is_fraction)
VARIABLES = {
    "temperature": ("K", "+", False),
    "pmid": ("Pa", "+", False),
    "pdel": ("Pa", "+", False),
    "qv": ("kg/kg", "+", False),
    "qc": ("kg/kg", "+", False),
    "qi": ("kg/kg", "+", False),
    "nc": ("#/kg", "+", False),
    "ni": ("#/kg", "+", False),
    "cldt": ("fraction", "+", True),
    "prain": ("kg/kg/s", None, False),
    "evapr": ("kg/kg/s", "+", False),
    "rprddp": ("kg/kg/s", None, False),
    "evapcdp": ("kg/kg/s", "+", False),
    "icwmrdp": ("kg/kg", "+", False),
    "dp_frac": ("fraction", "+", True),
    "dlf": ("kg/kg/s", "+", False),
    "mu": ("mb/s", "+", False),
    "md": ("mb/s", "-", False),
    "eu": ("1/s", "+", False),
    "du": ("1/s", "+", False),
    "ed": ("1/s", None, False),
    "dp": ("mb", "+", False),
    "dpdry": ("mb", "+", False),
}

# Mapping from EAMxx variable to the netCDF field(s) and the conversion the
# wetscav interface applies. Shown in the report.
EAMXX_SOURCE = {
    "temperature": "T_mid",
    "pmid": "p_mid",
    "pdel": "pseudo_density",
    "qv": "qv (wet; interface passes dry)",
    "qc": "qc (wet; interface passes dry)",
    "qi": "qi (wet; interface passes dry)",
    "nc": "nc (wet; interface passes dry)",
    "ni": "ni (wet; interface passes dry)",
    "cldt": "cldfrac_liq",
    "prain": "precip_total_tend",
    "evapr": "nevapr",
    "rprddp": "zm_rain_prod",
    "evapcdp": "evap_dq_out",
    "icwmrdp": "zm_icw_liq",
    "dp_frac": "f(zm_mass_flux_int[k+1])",
    "dlf": "zm_detr_qc + zm_detr_qi",
    "mu": "zm_mflx_up * g/100",
    "md": "zm_mflx_dn * g/100",
    "eu": "zm_entr_up",
    "du": "zm_detr_up",
    "ed": "zm_entr_dn",
    "dp": "pseudo_density * 0.01",
    "dpdry": "pseudo_density_dry * 0.01",
}

VALIDATION_SOURCE = {
    "temperature": "temperature",
    "pmid": "pmid",
    "pdel": "pdel",
    "qv": "state_q[:,0]",
    "qc": "state_q[:,1]",
    "qi": "state_q[:,2]",
    "nc": "state_q[:,3]",
    "ni": "state_q[:,4]",
    "cldt": "inputs_cldt",
    "prain": "inputs_prain",
    "evapr": "inputs_evapr",
    "rprddp": "rprddp",
    "evapcdp": "evapcdp",
    "icwmrdp": "p_icwmrdp",
    "dp_frac": "p_dp_frac",
    "dlf": "dlf",
    "mu": "mu",
    "md": "md",
    "eu": "eu",
    "du": "du",
    "ed": "ed",
    "dp": "dp",
    "dpdry": "pdeldry / 100",
}


# Tracer order of the pcnst = 40 constituents (same in EAM, mam4xx and EAMxx):
# 9 non-chemistry tracers, then the 31 gas-phase chemistry species
# (mam4xx mo_setsox.hpp). Aerosols start at index 15.
PCNST_NAMES = [
    "Q", "CLDLIQ", "CLDICE", "NUMLIQ", "NUMICE", "RAINQM", "SNOWQM", "NUMRAI",
    "NUMSNO", "O3", "H2O2", "H2SO4", "SO2", "DMS", "SOAG",
    "so4_a1", "pom_a1", "soa_a1", "bc_a1", "dst_a1", "ncl_a1", "mom_a1",
    "num_a1", "so4_a2", "soa_a2", "ncl_a2", "mom_a2", "num_a2", "dst_a3",
    "ncl_a3", "so4_a3", "bc_a3", "pom_a3", "soa_a3", "mom_a3", "num_a3",
    "pom_a4", "bc_a4", "mom_a4", "num_a4",
]
AERO_START = 15

# Fortran outputs compared with the EAMxx output fields of the same name.
DEPOSITION_OUTPUTS = ("aerdepwetis", "aerdepwetcw")


# -----------------------------------------------------------------------------
# Readers
# -----------------------------------------------------------------------------
def _literal(node):
    """ast.literal_eval that also accepts the names nan and inf."""
    if isinstance(node, ast.Name) and node.id in ("nan", "inf"):
        return float(node.id)
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        return -_literal(node.operand)
    if isinstance(node, (ast.List, ast.Tuple)):
        return [_literal(e) for e in node.elts]
    return ast.literal_eval(node)


def read_skywalker_py(path):
    """Read a Skywalker Python data file (e.g. mam_<test>.py) without running it.

    Returns (inputs, outputs, header). Each value is a flat float array.
    """
    with open(path) as f:
        text = f.read()
    header = text.splitlines()[0].strip() if text else ""
    groups = {"input": {}, "output": {}}
    for node in ast.parse(text).body:
        if not isinstance(node, ast.Assign) or len(node.targets) != 1:
            continue
        t = node.targets[0]
        if (isinstance(t, ast.Attribute) and isinstance(t.value, ast.Name)
                and t.value.id in groups):
            vals = np.asarray(_literal(node.value), dtype=float)
            # Skywalker stores one row per ensemble member; take the first.
            groups[t.value.id][t.attr] = vals[0].ravel() if vals.ndim > 1 \
                else vals.ravel()
    return groups["input"], groups["output"], header


def read_validation(path):
    """Read the validation inputs (and Fortran outputs, for a .py file).

    Accepts the Skywalker input YAML (aero_model_wetdep_convproc_ts_3.yaml)
    or the Fortran data file (mam_aero_model_wetdep_convproc_ts_3.py), which
    holds the same inputs plus the Fortran outputs.
    """
    if path.endswith(".py"):
        a, outputs, header = read_skywalker_py(path)
    else:
        with open(path) as f:
            header = f.readline().strip()
            f.seek(0)
            data = yaml.safe_load(f)["input"]["fixed"]
        a = {k: np.asarray(v, dtype=float) for k, v in data.items()}
        outputs = {}
    return _validation_inputs(a) + (outputs, header)


def _validation_inputs(a):
    """Inputs as the mam4xx test passes them to aero_model_wetdep."""
    nlev = int(a["pver"][0])
    pcnst = int(a["pcnst"][0])
    # state_q is stored column-major (Fortran order), as read by
    # convert_1d_vector_to_2d_view_device in the mam4xx test.
    state_q = a["state_q"].reshape(pcnst, nlev).T

    v = {
        "temperature": a["temperature"],
        "pmid": a["pmid"],
        "pdel": a["pdel"],
        "qv": state_q[:, 0],
        "qc": state_q[:, 1],
        "qi": state_q[:, 2],
        "nc": state_q[:, 3],
        "ni": state_q[:, 4],
        "cldt": a["inputs_cldt"],
        "prain": a["inputs_prain"],
        "evapr": a["inputs_evapr"],
        "rprddp": a["rprddp"],
        "evapcdp": a["evapcdp"],
        "icwmrdp": a["p_icwmrdp"],
        "dp_frac": a["p_dp_frac"],
        "dlf": a["dlf"],
        "mu": a["mu"],
        "md": a["md"],
        "eu": a["eu"],
        "du": a["du"],
        "ed": a["ed"],
        "dp": a["dp"],
        "dpdry": a["pdeldry"] / 100.0,
    }
    # same conversions as aero_model_wetdep_convproc.cpp
    scalars = {
        "ktop": int(round(a["jt"][0] - 1)),
        "kbot": int(round(a["maxg"][0])),
        "zm_activity": 1,
        "dt": float(a["dt"][0]),
        "convproc_do_aer": bool(a["convproc_do_aer"][0]),
        "convproc_do_gas": bool(a["convproc_do_gas"][0]),
    }
    shallow = {k: a[k] for k in ("rprdsh", "evapcsh", "p_sh_frac", "p_icwmrsh")
               if k in a}
    # Wrap 1-column profiles as (1, nlev) so both sides look the same.
    v = {k: x.reshape(1, -1) for k, x in v.items()}
    s = {k: np.array([x]) for k, x in scalars.items()}
    return v, s, shallow


def _get(nc, name, tidx):
    if name not in nc.variables:
        raise KeyError(f"field '{name}' not found in EAMxx output; "
                       f"add it to output.yaml")
    var = nc.variables[name]
    # Read only the requested time slice: a ne30pg2 file with many slices
    # would not fit in memory otherwise.
    idx = tuple(tidx if d == "time" else slice(None) for d in var.dimensions)
    x = var[idx]
    return np.ma.filled(np.ma.asarray(x, dtype=float), np.nan)


def read_eamxx(path, tidx):
    """Return inputs as the wetscav interface passes them, for all columns."""
    from netCDF4 import Dataset

    with Dataset(path) as nc:
        g = lambda n: _get(nc, n, tidx)  # noqa: E731
        ntime = len(nc.dimensions["time"]) if "time" in nc.dimensions else 1
        if not -ntime <= tidx < ntime:
            sys.exit(f"error: --time {tidx} is out of range; {path} has "
                     f"{ntime} time slice(s) (valid: 0..{ntime - 1} or -1)")

        mflx_int = g("zm_mass_flux_int")              # (ncol, nlev+1)
        cmfmc = np.maximum(0.0, mflx_int[:, 1:])      # interface below mid-level k
        dp_frac = np.clip(DP1 * np.log(1.0 + DP2 * cmfmc), 0.0, 0.60)

        v = {
            "temperature": g("T_mid"),
            "pmid": g("p_mid"),
            "pdel": g("pseudo_density"),
            "qv": g("qv"),
            "qc": g("qc"),
            "qi": g("qi"),
            "nc": g("nc"),
            "ni": g("ni"),
            "cldt": g("cldfrac_liq"),
            "prain": g("precip_total_tend"),
            "evapr": g("nevapr"),
            "rprddp": g("zm_rain_prod"),
            "evapcdp": g("evap_dq_out"),
            "icwmrdp": g("zm_icw_liq"),
            "dp_frac": dp_frac,
            "dlf": g("zm_detr_qc") + g("zm_detr_qi"),
            "mu": g("zm_mflx_up") * KGM2S_TO_MBS,
            "md": g("zm_mflx_dn") * KGM2S_TO_MBS,
            "eu": g("zm_entr_up"),
            "du": g("zm_detr_up"),
            "ed": g("zm_entr_dn"),
            "dp": g("pseudo_density") * PA_TO_MB,
            "dpdry": g("pseudo_density_dry") * PA_TO_MB,
        }
        s = {
            "ktop": np.rint(np.nan_to_num(g("zm_cldtop_klev"))).astype(int),
            # mam4xx loops over ktop <= k < kbot, so the interface passes the
            # 0-based ZM cloud-base level + 1 (same as Fortran maxg).
            "kbot": np.rint(np.nan_to_num(g("zm_msemax_klev"))).astype(int) + 1,
            "zm_activity": np.where(np.isnan(g("zm_activity")), 0,
                                    np.rint(np.nan_to_num(g("zm_activity")))
                                    ).astype(int),
            # raw values, to detect averaged (non-instant) output
            "zm_activity_raw": np.nan_to_num(g("zm_activity")),
            "ktop_raw": np.nan_to_num(g("zm_cldtop_klev")),
        }
        lat = g("lat") if "lat" in nc.variables else None
        lon = g("lon") if "lon" in nc.variables else None
        dep = {}
        for name in DEPOSITION_OUTPUTS + ("aerdepwetis_convproc",):
            if name in nc.variables:
                x = g(name)
                # make it (ncol, pcnst)
                if x.ndim == 2 and x.shape[0] == len(PCNST_NAMES) \
                        and x.shape[1] != len(PCNST_NAMES):
                    x = x.T
                dep[name] = x
    return v, s, lat, lon, ntime, dep


# -----------------------------------------------------------------------------
# Checks
# -----------------------------------------------------------------------------
def stats(x):
    x = np.asarray(x).ravel()
    x = x[~np.isnan(x)]          # fill values were turned into nan
    nz = np.abs(x[np.abs(x) > TINY * np.max(np.abs(x))]) if x.size else x
    return {
        "min": float(x.min()) if x.size else math.nan,
        "max": float(x.max()) if x.size else math.nan,
        "medabs": float(np.median(nz)) if nz.size else 0.0,
        "fneg": float(np.mean(x < 0)) if x.size else 0.0,
        "fnz": float(np.mean(x != 0)) if x.size else 0.0,
    }


def nearest_factor(ratio):
    best, best_d = None, None
    for name, f in KNOWN_FACTORS.items():
        d = abs(math.log10(ratio) - math.log10(f))
        if best_d is None or d < best_d:
            best, best_d = name, d
    return best if best_d is not None and best_d < math.log10(1.5) else None


THICKNESS_VARS = ("pdel", "dp", "dpdry")
TINY = 1e-8          # values below TINY * max|x| are ignored in medians
FRAC_FLAG = 0.25     # flag when more than this fraction of columns is off
MAX_INDEX_ROWS = 20  # EAMxx rows shown in the ktop/kbot table
RATIO_FLAG = 5.0   # flag a column whose median|x| is more than 5x from validation


def spread_verdict(ratios, cols, ref="validation"):
    """Judge per-column ratios (EAMxx / reference; nan = not defined).

    Returns (flag or None, note or None). A flag means a systematic factor
    (the median over columns is more than RATIO_FLAG off, which is what a unit
    error looks like) or more than FRAC_FLAG of the columns off. A few columns
    off is normal column-to-column variation and only gets a note.
    """
    r = np.asarray(ratios, dtype=float)
    ok = ~np.isnan(r) & (r > 0)
    if not np.any(ok):
        return None, None
    lr = np.log10(r[ok])
    off = np.abs(lr) > math.log10(RATIO_FLAG)
    med = float(10 ** np.median(lr))
    cols_off = np.asarray(cols)[ok][off]
    worst = float(10 ** lr[np.argmax(np.abs(lr))])
    p5, p95 = (float(10 ** q) for q in np.percentile(lr, [5, 95]))
    spread = f"p5-p95 {p5:.3g} to {p95:.3g}x"
    if abs(math.log10(med)) > math.log10(RATIO_FLAG):
        fac = nearest_factor(med)
        hint = f", close to {fac}" if fac else ""
        if abs(math.log10(med)) > 8:
            hint = ", mass and number mixed up: species index shift?"
        return (f"systematic: median over columns {med:.3g}x {ref}{hint} "
                f"({spread})"), None
    if off.mean() > FRAC_FLAG:
        fac = nearest_factor(worst)
        hint = f", close to {fac}" if fac else ""
        if abs(math.log10(worst)) > 8:
            hint = ", mass and number mixed up: species index shift?"
        return (f"more than {RATIO_FLAG:g}x from {ref} in "
                + _cols_list(cols_off, np.asarray(cols)[ok])
                + f", worst {worst:.3g}x{hint}"), None
    if off.any():
        return None, (f"{off.sum()} of {ok.sum()} columns more than "
                      f"{RATIO_FLAG:g}x off ({spread}), likely column-to-"
                      "column variation")
    return None, None


def check_variable(key, e, v, cols):
    """e: (ncol_selected, nlev) EAMxx values; v: (1, nlev) validation values."""
    units, sign, is_frac = VARIABLES[key]
    se, sv = stats(e), stats(v)
    flags = []

    if se["medabs"] > 0 and sv["medabs"] > 0:
        ratio = se["medabs"] / sv["medabs"]
    else:
        ratio = math.nan
        if (se["medabs"] > 0) != (sv["medabs"] > 0):
            flags.append("all zero on one side only")

    # Per-column ratios: one bad column must not hide in a pooled median.
    notes = []
    if sv["medabs"] > 0:
        a = np.abs(np.asarray(e, dtype=float))
        amax = np.nanmax(a) if np.any(~np.isnan(a)) else 0.0
        a[a <= TINY * amax] = np.nan
        with warnings.catch_warnings():   # all-nan columns are expected
            warnings.simplefilter("ignore", RuntimeWarning)
            med = np.nanmedian(a, axis=1)          # per-column median|x|
        flag, note = spread_verdict(med / sv["medabs"], cols)
        flags += [flag] if flag else []
        notes += [note] if note else []

    if sign == "+" and se["fneg"] > 0:
        flags.append(f"EAMxx has {se['fneg']:.1%} negative values")
    if sign == "-" and se["min"] < 0 and se["max"] > 0:
        flags.append("EAMxx has mixed signs (expected <= 0)")
    if sign is not None and (se["fneg"] > 0) != (sv["fneg"] > 0):
        flags.append(f"sign pattern differs (neg: EAMxx {se['fneg']:.1%}, "
                     f"validation {sv['fneg']:.1%})")
    if is_frac and (se["min"] < 0 or se["max"] > 1):
        flags.append("EAMxx fraction outside [0, 1]")
    if key in ("dp", "dpdry") and se["max"] > 100:
        flags.append("EAMxx dp > 100 mb: looks like Pa")
    if key in ("dp", "dpdry") and sv["max"] > 100:
        flags.append("validation dp > 100 mb: looks like Pa")

    return {"key": key, "units": units, "e": se, "v": sv, "ratio": ratio,
            "flags": flags, "notes": notes}


def conv_layers(eu, du):
    """0-based first and last level with convective entrainment/detrainment."""
    idx = np.nonzero((eu != 0) | (du != 0))[0]
    if idx.size == 0:
        return None, None
    return int(idx[0]), int(idx[-1])


def check_indices(name, vals, scal, cols):
    """Compare ktop/kbot with the layers where eu/du are nonzero.

    mam4xx loops over ktop <= k < kbot, so a correct kbot is
    (last convective layer) + 1, and ktop is the first convective layer.
    """
    rows = []
    for c in cols:
        first, last = conv_layers(vals["eu"][c], vals["du"][c])
        ktop, kbot = int(scal["ktop"][c]), int(scal["kbot"][c])
        if first is None:
            rows.append((name, c, ktop, kbot, "-", "-", "no eu/du > 0"))
            continue
        notes = []
        if ktop != first:
            notes.append(f"ktop != first conv layer ({first})")
        if kbot == last:
            notes.append("kbot == last conv layer: kbot is 1 too small "
                         "(needs +1)")
        elif kbot != last + 1:
            notes.append(f"kbot != last conv layer + 1 ({last + 1})")
        rows.append((name, c, ktop, kbot, first, last,
                     "; ".join(notes) if notes else "OK"))
    return rows


def _cols_list(bad, cols, limit=8):
    shown = ", ".join(str(int(c)) for c in bad[:limit])
    more = f" +{len(bad) - limit} more" if len(bad) > limit else ""
    return f"{len(bad)} of {len(cols)} columns (cols {shown}{more})"


def check_deposition(fortran, eamxx, cols, mam4xx=None):
    """Per-species comparison of a wet deposition flux [kg/m2/s or #/m2/s].

    fortran, mam4xx: (pcnst,) single-column values; eamxx: (ncol, pcnst).
    Every check is done per column, so one bad column is not hidden by a
    median over columns. Returns table rows, pooled per-species ratios, and
    per-column ratios (ncol_sel, pcnst; nan where not defined).
    """
    e = eamxx[cols]
    with np.errstate(divide="ignore", invalid="ignore"):
        rcol = np.where((fortran != 0) & (e != 0),
                        np.abs(e) / np.abs(fortran), np.nan)
    rows, ratios = [], {}
    for i, sp in enumerate(PCNST_NAMES):
        f = float(fortran[i])
        ei = e[:, i]
        flags, notes = [], []
        if f == 0.0:
            bad = cols[ei != 0]
            if bad.size:
                flags.append("nonzero where Fortran is zero in "
                             + _cols_list(bad, cols))
        else:
            if not np.any(ei != 0):
                flags.append("zero in all EAMxx columns, nonzero in Fortran")
            wrong = cols[np.sign(ei) == -np.sign(f)]
            if wrong.size:
                flags.append("opposite sign to Fortran in "
                             + _cols_list(wrong, cols))
            ri = rcol[:, i]
            if np.any(~np.isnan(ri)):
                ratios[sp] = float(np.nanmedian(ri))
            flag, note = spread_verdict(ri, cols, "Fortran")
            flags += [flag] if flag else []
            notes += [note] if note else []
        m = "" if mam4xx is None else fmt(float(mam4xx[i]))
        rows.append((i, sp, fmt(f), m, fmt(float(ei.min())),
                     fmt(float(np.median(ei))), fmt(float(ei.max())),
                     f"{ratios[sp]:.3g}" if sp in ratios else "",
                     "⚠️ " + "; ".join(flags) if flags else
                     "OK" + (f" ({'; '.join(notes)})" if notes else "")))
    return rows, ratios, rcol


def column_group_ratios(rcol, number):
    """Per-column median ratio over aerosol mass (or number) species."""
    idx = [i for i, sp in enumerate(PCNST_NAMES)
           if i >= AERO_START and sp.startswith("num_") == number]
    sub = rcol[:, idx]
    with np.errstate(all="ignore"):
        return np.array([np.nanmedian(r) if np.any(~np.isnan(r)) else np.nan
                         for r in sub])


# -----------------------------------------------------------------------------
# Report
# -----------------------------------------------------------------------------
def fmt(x):
    return "nan" if x is None or (isinstance(x, float) and math.isnan(x)) \
        else f"{x:.3e}"


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("eamxx_nc")
    p.add_argument("validation_file",
                   help="mam_<test>.py (inputs + Fortran outputs) or <test>.yaml")
    p.add_argument("--mam4xx",
                   help="mam4xx validation output (mam4xx_<test>.py), to show "
                        "the standalone mam4xx deposition next to Fortran")
    p.add_argument("--time", type=int, default=-1,
                   help="time index in the EAMxx file (default: last)")
    p.add_argument("--report", help="write the Markdown report to this file")
    args = p.parse_args()

    vv, vs, vshallow, vout, vheader = read_validation(args.validation_file)
    ev, es, lat, lon, ntime, edep = read_eamxx(args.eamxx_nc, args.time)
    mout = read_skywalker_py(args.mam4xx)[1] if args.mam4xx else {}

    ncol = ev["temperature"].shape[0]
    nlev_e, nlev_v = ev["temperature"].shape[1], vv["temperature"].shape[1]
    raw_act = es["zm_activity_raw"]
    averaged = bool(np.any(np.abs(raw_act - np.round(raw_act)) > 1e-6)
                    or np.any(np.abs(es["ktop_raw"] - np.round(es["ktop_raw"]))
                              > 1e-6))
    active = np.nonzero(es["zm_activity"] == 1)[0]

    out = []
    w = out.append
    w("# ZM → wetscav input check: EAMxx vs mam4xx validation\n")
    w(f"- EAMxx file: `{args.eamxx_nc}` (time index {args.time} of {ntime}, "
      f"{ncol} columns, {nlev_e} levels)")
    w(f"- Validation file: `{args.validation_file}` ({nlev_v} levels)")
    if args.mam4xx:
        w(f"- mam4xx standalone output: `{args.mam4xx}`")
    w(f"- Validation header: `{vheader.lstrip('# ')}`")
    w(f"- Columns with ZM active: {active.size} of {ncol}")
    w(f"- Validation dt = {vs['dt'][0]:g} s, convproc_do_aer = "
      f"{vs['convproc_do_aer'][0]}, convproc_do_gas = {vs['convproc_do_gas'][0]}")
    if averaged:
        w("- **⚠️ `zm_activity` or `zm_cldtop_klev` has non-integer values: the "
          "output looks averaged. Use `averaging_type: instant` for this "
          "check; averaged indices, activity flags and profiles are not what "
          "wetscav received.**")
    if nlev_e != nlev_v:
        w(f"- Level count differs: EAMxx {nlev_e}, validation {nlev_v}. "
          f"EAMxx layer thicknesses ({', '.join(THICKNESS_VARS)}) are "
          f"multiplied by {nlev_e}/{nlev_v} before comparing, so thinner "
          "layers are not taken for a unit error. Level indices (ktop, kbot) "
          "are only checked against each side's own profiles.")
    w("")

    if active.size == 0:
        w("**No active ZM column at this time.** Convective fields are all zero; "
          "try another `--time`. Only non-convective fields are compared, over "
          "all columns.\n")
        cols = np.arange(ncol)
    else:
        cols = active

    # --- variable table
    w("## Values as passed to `aero_model_wetdep`\n")
    w("EAMxx statistics are over the active columns. `median|x|` uses nonzero "
      "values only. `ratio` = EAMxx median|x| / validation median|x|, pooled "
      f"over columns. A column is flagged when its own ratio is more than "
      f"{RATIO_FLAG:g}x from 1. Different columns can differ that much "
      "physically, so a flag means: look at it. A unit error shows up as the "
      "same factor (100, g, 100/g, ...) in most columns.\n")
    w("| Arg | Units | EAMxx source | Validation source | EAMxx min | EAMxx max "
      "| EAMxx median\\|x\\| | Val min | Val max | Val median\\|x\\| | ratio "
      "| Status |")
    w("|---|---|---|---|---|---|---|---|---|---|---|---|")
    nflag = 0
    # Layer thicknesses scale with the number of levels (72 vs e.g. 128).
    thick_scale = nlev_e / nlev_v
    for key in VARIABLES:
        e_key = ev[key][cols]
        if key in THICKNESS_VARS and nlev_e != nlev_v:
            e_key = e_key * thick_scale
        r = check_variable(key, e_key, vv[key], cols)
        status = ("⚠️ " + "; ".join(r["flags"]) if r["flags"] else
                  "OK" + (f" ({'; '.join(r['notes'])})" if r["notes"] else ""))
        nflag += bool(r["flags"])
        ratio = "nan" if math.isnan(r["ratio"]) else f"{r['ratio']:.3g}"
        w(f"| `{key}` | {r['units']} | {EAMXX_SOURCE[key]} | "
          f"{VALIDATION_SOURCE[key]} | {fmt(r['e']['min'])} | "
          f"{fmt(r['e']['max'])} | {fmt(r['e']['medabs'])} | "
          f"{fmt(r['v']['min'])} | {fmt(r['v']['max'])} | "
          f"{fmt(r['v']['medabs'])} | {ratio} | {status} |")
    w("")

    # --- shallow convection
    w("## Shallow convection (EAMxx passes zeros)\n")
    for k, x in vshallow.items():
        mx = float(np.max(np.abs(x)))
        w(f"- validation `{k}`: max|x| = {mx:.3e} "
          f"{'OK' if mx == 0 else '⚠️ nonzero, EAMxx uses 0'}")
    w("")

    # --- index check
    w("## ktop / kbot convention\n")
    w("mam4xx loops over `ktop <= k < kbot`, so a correct `ktop` is the first "
      "layer with `eu` or `du` > 0, and a correct `kbot` is the last such layer "
      "+ 1. The validation test uses `ktop = jt - 1`, `kbot = maxg`. The EAMxx "
      "interface uses `ktop = round(zm_cldtop_klev)`, "
      "`kbot = round(zm_msemax_klev) + 1`.\n")
    w("| Side | col | ktop | kbot | first conv layer | last conv layer | Result |")
    w("|---|---|---|---|---|---|---|")
    rows = check_indices("validation", vv, vs, [0])
    if active.size:
        rows += check_indices("EAMxx", ev, es, active)
    n_off = sum(1 for r in rows if r[0] == "EAMxx" and "1 too small" in r[6])
    eam_rows = [r for r in rows if r[0] == "EAMxx"]
    eam_bad = [r for r in eam_rows if r[6] != "OK"]
    # Show the validation column, then problem columns first, at most
    # MAX_INDEX_ROWS EAMxx rows (a ne30pg2 run has thousands of active columns).
    shown = [r for r in rows if r[0] == "validation"] + \
        (eam_bad + [r for r in eam_rows if r[6] == "OK"])[:MAX_INDEX_ROWS]
    for side, c, kt, kb, f, l, note in shown:
        loc = ""
        if side == "EAMxx" and lat is not None and lon is not None:
            loc = f" ({float(lat[c]):.1f}, {float(lon[c]):.1f})"
        w(f"| {side} | {c}{loc} | {kt} | {kb} | {f} | {l} | {note} |")
    if len(eam_rows) > MAX_INDEX_ROWS:
        w(f"\n{len(eam_rows) - MAX_INDEX_ROWS} more EAMxx columns not shown.")
    if eam_rows:
        w(f"\nEAMxx: {len(eam_rows) - len(eam_bad)} of {len(eam_rows)} active "
          f"columns OK, {len(eam_bad)} with a problem.")
    w("")

    # --- wet deposition outputs
    dep_summary = []
    w("## Wet deposition outputs\n")
    if not vout:
        w("Skipped: pass the Fortran data file (`mam_<test>.py`) to compare "
          "`aerdepwetis` and `aerdepwetcw`.\n")
    else:
        w("Fortran is the single validation column. EAMxx min / median / max "
          "are over the columns used above. Units: kg/m2/s for mass species, "
          "#/m2/s for `num_*` (the EAMxx field metadata says kg/m2/s for all). "
          "Negative = removal. Deposition depends on the column's precipitation, "
          "so a single species more than "
          f"{RATIO_FLAG:g}x off can be physical; a unit or index error shows up "
          "as the same factor across the mass (or number) species, or as a "
          "shifted zero pattern.\n")
        for name in DEPOSITION_OUTPUTS:
            if name not in vout:
                w(f"- `{name}`: not in the validation file\n")
                continue
            if name not in edep:
                w(f"- `{name}`: not in the EAMxx output; add it to output.yaml\n")
                continue
            if edep[name].shape[1] != len(PCNST_NAMES):
                w(f"- `{name}`: EAMxx shape {edep[name].shape} is not "
                  f"(ncol, {len(PCNST_NAMES)})\n")
                continue
            rows, ratios, rcol = check_deposition(
                vout[name], edep[name], cols,
                mout.get(name) if mout else None)
            w(f"### `{name}`\n")
            mh = " mam4xx |" if mout else ""
            ms = "---|" if mout else ""
            w(f"| idx | species | Fortran |{mh} EAMxx min | EAMxx median | "
              "EAMxx max | ratio | Status |")
            w(f"|---|---|---|{ms}---|---|---|---|---|")
            nbad = 0
            for i, sp, f, m, emin, emed, emax, r, st in rows:
                mcol = f" {m} |" if mout else ""
                w(f"| {i} | {sp} | {f} |{mcol} {emin} | {emed} | {emax} | "
                  f"{r} | {st} |")
                nbad += st.startswith("⚠️")
            w("")
            w("`ratio` is the median over columns of \\|EAMxx\\| / "
              "\\|Fortran\\|.\n")
            for label, number in (("mass", False), ("number", True)):
                g = column_group_ratios(rcol, number)
                ok = ~np.isnan(g)
                if not np.any(ok):
                    continue
                w(f"- {label} species: per-column median ratio ranges "
                  f"{np.nanmin(g):.3g} to {np.nanmax(g):.3g} "
                  f"(median {np.nanmedian(g):.3g})")
                flag, note = spread_verdict(g, cols, "Fortran")
                if flag:
                    w(f"  - ⚠️ all {label} species together: {flag}")
                elif note:
                    w(f"  - {note}")
            if mout and name in mout:
                d = np.abs(mout[name] - vout[name])
                sc = np.maximum(np.abs(vout[name]), 1e-300)
                rel = float(np.max(np.where(vout[name] != 0, d / sc, 0.0)))
                w(f"- mam4xx standalone vs Fortran: max relative difference "
                  f"{rel:.3g} (this is the error of the validation test "
                  "itself)")
            w("")
            dep_summary.append(f"- `{name}`: {nbad} of {len(PCNST_NAMES)} "
                               "species flagged")
        if "aerdepwetis_convproc" in edep and "aerdepwetis" in edep:
            tot = np.abs(edep["aerdepwetis"][cols][:, AERO_START:]).sum(axis=1)
            cnv = np.abs(edep["aerdepwetis_convproc"][cols][:, AERO_START:]
                         ).sum(axis=1)
            share = cnv[tot > 0] / tot[tot > 0]
            if share.size:
                w(f"- EAMxx share of `aerdepwetis` from convective processing "
                  f"(`aerdepwetis_convproc`, aerosol species, summed |flux|): "
                  f"median {np.median(share):.3g}, max {share.max():.3g} over "
                  f"{share.size} columns. The Fortran file has no separate "
                  "convproc output.\n")

    # --- summary
    w("## Summary\n")
    w(f"- Input variables flagged: {nflag} of {len(VARIABLES)}")
    for line in dep_summary:
        w(line)
    if active.size:
        w(f"- EAMxx active columns where kbot is 1 too small: "
          f"{n_off} of {active.size}")
    w("- Not checked here: aerosol mixing ratios, dgnumwet, dgnum, qaerwat, "
      "wetdens. The output is written after wetscav updates them, so it does "
      "not show what wetscav received.")
    w("- EAMxx passes dry mixing ratios (qv, qc, qi, nc, ni, aerosols); the "
      "Fortran data are wet. The difference is a factor of about 1/(1-qv), "
      "well below the flag threshold.")

    text = "\n".join(out)
    print(text)
    if args.report:
        with open(args.report, "w") as f:
            f.write(text + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
