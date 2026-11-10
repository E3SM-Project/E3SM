#!/usr/bin/env python3
"""
Compare the inputs that the EAMxx MAM4 wetscav interface passes to
mam4::wetdep::aero_model_wetdep with the inputs of the mam4xx validation
test aero_model_wetdep_convproc.

The EAMxx side comes from the netCDF output of the p3_zm_mam4_wetscav test
(see output.yaml). The script applies the same unit conversions as
eamxx_mam_wetscav_process_interface.cpp, so both sides are compared in the
units mam4xx expects (mb, mb/s, 0-based indices with exclusive kbot).

The two sides are different columns and states, so values are not expected
to match. The script checks orders of magnitude, signs, value ranges and the
ktop/kbot index convention.

Usage:
  compare_zm_wetscav_inputs.py EAMXX_OUTPUT.nc VALIDATION_INPUT.yaml
                               [--time IDX] [--report FILE.md]

VALIDATION_INPUT.yaml is e.g. mam_x_validation/aero_model/aero_model_wetdep_convproc_ts_3.yaml
"""

import argparse
import math
import sys

import numpy as np
import yaml

GRAVITY = 9.80616   # mam4::Constants::gravity
PA_TO_MB = 0.01
KGM2S_TO_MBS = GRAVITY * PA_TO_MB   # same as the wetscav interface
DP1, DP2 = 0.018, 500.0             # dp_frac formula in the wetscav interface

# Known unit factors, used to label a magnitude mismatch
KNOWN_FACTORS = {
    "100 (Pa vs mb)": 100.0,
    "1/100 (mb vs Pa)": 0.01,
    "g": GRAVITY,
    "1/g": 1.0 / GRAVITY,
    "100/g (mb/s -> kg/m2/s)": 100.0 / GRAVITY,
    "g/100 (kg/m2/s -> mb/s)": GRAVITY / 100.0,
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


# -----------------------------------------------------------------------------
# Readers
# -----------------------------------------------------------------------------
def read_validation(path):
    """Return (inputs dict as passed to aero_model_wetdep, header string)."""
    with open(path) as f:
        header = f.readline().strip()
        f.seek(0)
        data = yaml.safe_load(f)["input"]["fixed"]

    a = {k: np.asarray(v, dtype=float) for k, v in data.items()}
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
    return v, s, shallow, header


def _get(nc, name, tidx):
    if name not in nc.variables:
        raise KeyError(f"field '{name}' not found in EAMxx output; "
                       f"add it to output.yaml")
    var = nc.variables[name]
    x = np.asarray(var[:], dtype=float)
    if "time" in var.dimensions:
        x = np.take(x, tidx, axis=var.dimensions.index("time"))
    return x


def read_eamxx(path, tidx):
    """Return inputs as the wetscav interface passes them, for all columns."""
    from netCDF4 import Dataset

    with Dataset(path) as nc:
        g = lambda n: _get(nc, n, tidx)  # noqa: E731
        ntime = len(nc.dimensions["time"]) if "time" in nc.dimensions else 1

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
            "ktop": np.rint(g("zm_cldtop_klev")).astype(int),
            # mam4xx loops over ktop <= k < kbot, so the interface passes the
            # 0-based ZM cloud-base level + 1 (same as Fortran maxg).
            "kbot": np.rint(g("zm_msemax_klev")).astype(int) + 1,
            "zm_activity": np.rint(g("zm_activity")).astype(int),
        }
        lat = g("lat") if "lat" in nc.variables else None
        lon = g("lon") if "lon" in nc.variables else None
    return v, s, lat, lon, ntime


# -----------------------------------------------------------------------------
# Checks
# -----------------------------------------------------------------------------
def stats(x):
    x = np.asarray(x).ravel()
    nz = np.abs(x[x != 0])
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


RATIO_FLAG = 5.0   # flag a column whose median|x| is more than 5x from validation


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
    if sv["medabs"] > 0:
        bad = []
        for i, c in enumerate(cols):
            m = stats(e[i])["medabs"]
            if m > 0:
                rc = m / sv["medabs"]
                if abs(math.log10(rc)) > math.log10(RATIO_FLAG):
                    bad.append((int(c), rc))
        if bad:
            worst = max(bad, key=lambda t: abs(math.log10(t[1])))
            fac = nearest_factor(worst[1])
            shown = ", ".join(f"{c}" for c, _ in bad[:8])
            more = f" (+{len(bad) - 8} more)" if len(bad) > 8 else ""
            flags.append(f"{len(bad)} of {len(cols)} columns more than "
                         f"{RATIO_FLAG:g}x from validation (cols {shown}{more}; "
                         f"worst {worst[1]:.3g}x"
                         + (f", close to {fac}" if fac else "") + ")")

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
            "flags": flags}


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
    p.add_argument("validation_yaml")
    p.add_argument("--time", type=int, default=-1,
                   help="time index in the EAMxx file (default: last)")
    p.add_argument("--report", help="write the Markdown report to this file")
    args = p.parse_args()

    vv, vs, vshallow, vheader = read_validation(args.validation_yaml)
    ev, es, lat, lon, ntime = read_eamxx(args.eamxx_nc, args.time)

    ncol = ev["temperature"].shape[0]
    nlev_e, nlev_v = ev["temperature"].shape[1], vv["temperature"].shape[1]
    active = np.nonzero(es["zm_activity"] == 1)[0]

    out = []
    w = out.append
    w("# ZM → wetscav input check: EAMxx vs mam4xx validation\n")
    w(f"- EAMxx file: `{args.eamxx_nc}` (time index {args.time} of {ntime}, "
      f"{ncol} columns, {nlev_e} levels)")
    w(f"- Validation file: `{args.validation_yaml}` ({nlev_v} levels)")
    w(f"- Validation header: `{vheader.lstrip('# ')}`")
    w(f"- Columns with ZM active: {active.size} of {ncol}")
    w(f"- Validation dt = {vs['dt'][0]:g} s, convproc_do_aer = "
      f"{vs['convproc_do_aer'][0]}, convproc_do_gas = {vs['convproc_do_gas'][0]}")
    if nlev_e != nlev_v:
        w(f"- **Level count differs: EAMxx {nlev_e}, validation {nlev_v}**")
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
    for key in VARIABLES:
        r = check_variable(key, ev[key][cols], vv[key], cols)
        status = "OK" if not r["flags"] else "⚠️ " + "; ".join(r["flags"])
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
    n_off = 0
    for side, c, kt, kb, f, l, note in rows:
        loc = ""
        if side == "EAMxx" and lat is not None and lon is not None:
            loc = f" ({float(lat[c]):.1f}, {float(lon[c]):.1f})"
        w(f"| {side} | {c}{loc} | {kt} | {kb} | {f} | {l} | {note} |")
        if side == "EAMxx" and "1 too small" in note:
            n_off += 1
    w("")

    # --- summary
    w("## Summary\n")
    w(f"- Variables flagged: {nflag} of {len(VARIABLES)}")
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
