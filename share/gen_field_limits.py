#!/usr/bin/env python3
"""Generate share/util/shr_field_limits_mod.F90 from share/field_limits.yaml.

Checks each entry on its own (keys, limits, names); it does not read the
driver's field lists. Whether every field a component sends has an entry is
checked by the model itself when it starts (shr_bounds_require, called from
seq_flds_set).

Usage:
    python share/gen_field_limits.py           # validate and write the Fortran module
    python share/gen_field_limits.py --check   # validate; fail if the module is out of date
    python share/gen_field_limits.py --check --base BASE.yaml
                                               # also fail if a field is unreviewed here
                                               # but not in BASE.yaml (e.g. from master)

Requires PyYAML.
"""

import argparse
import math
import sys
from collections import Counter
from pathlib import Path

import yaml

SHARE = Path(__file__).resolve().parent
YAML_FILE = SHARE / "field_limits.yaml"
FORTRAN_FILE = SHARE / "util" / "shr_field_limits_mod.F90"
COMPONENTS = {"atm", "lnd", "ice", "ocn", "rof", "glc", "wav", "iac"}
NAMELEN = 32
UNITSLEN = 24
STATUSLEN = 10
CHUNK = 100  # rows per Fortran array constructor, to stay within continuation-line limits
KEYS = {"component", "units", "min", "max", "min_inclusive", "max_inclusive",
        "where_positive", "owner", "reason", "unbounded", "status"}


def is_pattern(name):
    return name.endswith("*")


def kind(entry):
    if entry.get("status") == "unreviewed":
        return "unreviewed"
    if "unbounded" in entry:
        return "unbounded"
    return "limited"


def validate(fields):
    errors = []
    for name, entry in fields.items():
        where = f"{name}:"
        if len(name) > NAMELEN:
            errors.append(f"{where} name longer than {NAMELEN} characters")
        if is_pattern(name) and (len(name) < 2 or "*" in name[:-1]):
            errors.append(f"{where} a pattern is a base name followed by a single '*'")
        unknown = set(entry) - KEYS
        if unknown:
            errors.append(f"{where} unknown keys {sorted(unknown)}")
        for key in ("component", "units", "owner"):
            if key not in entry:
                errors.append(f"{where} missing '{key}'")
        comp = entry.get("component")
        if comp not in COMPONENTS:
            errors.append(f"{where} unknown component '{comp}'")

        has_limit = "min" in entry or "max" in entry
        if "status" in entry and entry["status"] != "unreviewed":
            errors.append(f"{where} status must be 'unreviewed'")
        if sum([has_limit, "unbounded" in entry, "status" in entry]) != 1:
            errors.append(f"{where} needs exactly one of: 'min'/'max', "
                          f"'unbounded: <reason>', 'status: unreviewed'")
        if "unbounded" in entry and not str(entry["unbounded"] or "").strip():
            errors.append(f"{where} 'unbounded' needs a reason")
        if not has_limit:
            for key in ("min_inclusive", "max_inclusive", "where_positive"):
                if key in entry:
                    errors.append(f"{where} '{key}' only applies to an entry with limits")

        if len(str(entry.get("units", ""))) > UNITSLEN:
            errors.append(f"{where} units longer than {UNITSLEN} characters")
        nums = {}
        for key in ("min", "max"):
            if key not in entry:
                continue
            v = entry[key]
            if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v):
                errors.append(f"{where} '{key}' is not a finite number")
            else:
                nums[key] = v
        if len(nums) == 2 and nums["min"] > nums["max"]:
            errors.append(f"{where} min > max")

        wp = entry.get("where_positive")
        if wp is not None and (fields.get(wp) or {}).get("component") != comp:
            errors.append(f"{where} where_positive '{wp}' is not an entry of {comp}")
    return errors


def load(path):
    with open(path) as fh:
        return yaml.safe_load(fh)["field_limits"]


def new_unreviewed(fields, base_file):
    """Fields unreviewed in fields but not in the base dictionary."""
    old = {n for n, e in load(base_file)["fields"].items() if kind(e) == "unreviewed"}
    return sorted(n for n, e in fields.items() if kind(e) == "unreviewed" and n not in old)


def fortran_real(value):
    # repr gives the shortest round-trip form, e.g. 0.0 or 1e-05
    return repr(float(value)) + "_r8"


def fortran_logical(value):
    return ".true." if value else ".false."


def generate(meta, fields):
    rows = []
    for name, e in fields.items():
        args = [f"'{name}'", f"'{e['component']}'",
                f"'{e['units']}'", f"'{kind(e)}'",
                fortran_logical("min" in e), fortran_real(e.get("min", 0.0)),
                fortran_logical(e.get("min_inclusive", True)),
                fortran_logical("max" in e), fortran_real(e.get("max", 0.0)),
                fortran_logical(e.get("max_inclusive", True)),
                f"'{e.get('where_positive', '')}'"]
        rows.append("       shr_field_limit_type(" + ", ".join(args) + ")")
    chunks = [rows[i:i + CHUNK] for i in range(0, len(rows), CHUNK)]
    decls = []
    for n, chunk in enumerate(chunks, 1):
        decls += [f"  type(shr_field_limit_type), parameter :: limits{n}({len(chunk)}) = [ &",
                  ", &\n".join(chunk) + " ]", ""]
    parts = ", ".join(f"limits{n}" for n in range(1, len(chunks) + 1))
    lines = [
        "! WARNING! DO NOT EDIT THIS FILE!",
        "! This file was generated automatically from share/field_limits.yaml",
        "! by share/gen_field_limits.py. Edit the YAML file and rerun the generator.",
        "",
        "module shr_field_limits_mod",
        "",
        "  ! Physical limits on fields sent from components to the coupler.",
        "",
        "  use shr_kind_mod, only: r8 => shr_kind_r8",
        "",
        "  implicit none",
        "  private",
        "",
        f"  character(len=*), parameter, public :: shr_field_limits_version = \"{meta['version_number']}\"",
        f"  integer, parameter, public :: shr_field_limits_nflds = {len(rows)}",
        "",
        "  type, public :: shr_field_limit_type",
        f"     character(len={NAMELEN}) :: name            ! coupler field name; 'base*' matches base followed by digits",
        "     character(len=3)  :: component       ! component that sends it",
        f"     character(len={UNITSLEN}) :: units",
        f"     character(len={STATUSLEN}) :: status          ! 'limited', 'unbounded' or 'unreviewed'",
        "     logical           :: has_min         ! false: no lower limit",
        "     real(r8)          :: min_value",
        "     logical           :: min_inclusive",
        "     logical           :: has_max         ! false: no upper limit",
        "     real(r8)          :: max_value",
        "     logical           :: max_inclusive",
        f"     character(len={NAMELEN}) :: where_positive  ! check only where this field > 0",
        "  end type shr_field_limit_type",
        "",
        *decls,
        "  type(shr_field_limit_type), parameter, public :: shr_field_limits(shr_field_limits_nflds) = &",
        f"       [ {parts} ]",
        "",
        "end module shr_field_limits_mod",
        "",
    ]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--check", action="store_true",
                        help="validate and fail if the Fortran module is out of date")
    parser.add_argument("--base", metavar="BASE_YAML",
                        help="fail if a field is unreviewed but was not unreviewed in BASE_YAML")
    args = parser.parse_args()

    data = load(YAML_FILE)
    fields = data["fields"]
    errors = validate(fields)
    if args.base:
        errors += [f"{n}: new 'status: unreviewed' entry; a new field needs min/max or "
                   f"'unbounded: <reason>'" for n in new_unreviewed(fields, args.base)]
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1

    counts = Counter(kind(e) for e in fields.values())
    summary = (f"{len(fields)} entries valid (" +
               ", ".join(f"{counts[k]} {k}" for k in ("limited", "unbounded", "unreviewed")) + ")")
    text = generate(data, fields)
    if args.check:
        if not FORTRAN_FILE.exists() or FORTRAN_FILE.read_text() != text:
            print(f"{FORTRAN_FILE} is out of date; run share/gen_field_limits.py", file=sys.stderr)
            return 1
        print(f"{summary}; {FORTRAN_FILE.name} is up to date")
        return 0
    FORTRAN_FILE.write_text(text)
    print(f"{summary}; wrote {FORTRAN_FILE}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
