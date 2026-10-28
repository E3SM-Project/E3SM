#!/usr/bin/env python3
"""Audit explicitly inventoried inputs for the experimental Gregorian input policy."""

import argparse
import json
from pathlib import Path
import re
import sys
import xml.etree.ElementTree as ET


ROOT = Path(__file__).resolve().parents[2]
ALIASES = {alias: "NO_LEAP" for alias in ("noleap", "no_leap", "365_day", "365day")}
ALIASES.update({alias: "GREGORIAN" for alias in ("gregorian", "standard", "proleptic_gregorian")})
ADAPTED_READERS = {"stream", "chemistry", "eamxx"}
UNSUPPORTED_READERS = {"mpas", "ww3", "flbc", "tracer", "eamxx-iop"}


def normalize_calendar(value):
    token = str(value).lower().translate(str.maketrans({"\t": " ", "\n": " ", "\r": " ", "\0": " "}))
    return ALIASES.get(token.strip(" "), "")


def calendar_flavor(value):
    canonical = normalize_calendar(value)
    if canonical != "GREGORIAN":
        return canonical
    token = str(value).lower().replace("\0", " ").strip()
    return "PROLEPTIC_GREGORIAN" if token == "proleptic_gregorian" else "STANDARD"


def make_rules(selectors, calendars):
    if len(selectors) != len(calendars) or len(selectors) > 128:
        raise ValueError("NOLEAP_INPUT: invalid override array sizes")
    rules = {}
    for selector, calendar in zip(selectors, calendars):
        if len(selector.rstrip(" ")) > 1024:
            raise ValueError("NOLEAP_INPUT: invalid exact override selector=" + selector)
        selector, calendar = selector.strip(" "), calendar.strip(" ")
        if not selector and not calendar:
            continue
        prefix, separator, identifier = selector.partition(":")
        if (not separator or prefix not in {"file", "stream"} or not identifier or
                len(selector) > 1024 or "\n" in selector or "\r" in selector):
            raise ValueError("NOLEAP_INPUT: invalid exact override selector=" + selector)
        canonical = calendar_flavor(calendar)
        if not canonical:
            raise ValueError("NOLEAP_INPUT: unsupported override calendar=" + calendar)
        if selector in rules:
            raise ValueError("NOLEAP_INPUT: duplicate override selector=" + selector)
        rules[selector] = canonical
    return rules


def resolve_calendar(declared, metadata_present, filename, stream, rules, include_flavor=False):
    matched = [(key, rules[key]) for key in ("file:" + filename, "stream:" + stream) if key in rules]
    overrides = {calendar for _, calendar in matched}
    if len(overrides) > 1:
        raise ValueError("NOLEAP_INPUT: conflicting file and stream calendar overrides")
    override = next(iter(overrides), "")
    canonical = normalize_calendar(declared) if metadata_present else ""
    empty = not str(declared).translate(str.maketrans({"\t": " ", "\n": " ", "\r": " ", "\0": " "})).strip(" ")
    if metadata_present and not canonical and not empty:
        raise ValueError("NOLEAP_INPUT: unsupported calendar=" + str(declared))
    if canonical and override and calendar_flavor(declared) != calendar_flavor(override):
        raise ValueError("NOLEAP_INPUT: calendar override conflicts with file metadata")
    if not canonical and not override:
        message = "empty calendar metadata" if metadata_present else "explicit time calendar metadata is required"
        raise ValueError("NOLEAP_INPUT: " + message)
    provenance = "metadata+override" if canonical and override else "metadata" if canonical else "override"
    result = (canonical or normalize_calendar(override), provenance, [key for key, _ in matched])
    flavor = calendar_flavor(declared) if canonical else calendar_flavor(override)
    return result + (flavor,) if include_flavor else result


def check_mixed_calendar_range(flavor, units, values):
    import cftime
    import netCDF4
    import numpy as np

    if flavor != "STANDARD":
        return
    cutoff = netCDF4.date2num(cftime.DatetimeProlepticGregorian(1582, 10, 15), units,
                             calendar="proleptic_gregorian")
    if cutoff > 0 or np.min(values) < cutoff:
        raise ValueError("INPUT_CALENDAR: unsupported historical calendar=" + flavor +
                         "; units=" + str(units) + "; range=" + str((float(np.min(values)), float(np.max(values)))) +
                         "; require epoch and range >= 1582-10-15")


def read_namelist(path, groupless=False):
    sys.path.insert(0, str(ROOT / "cime"))
    from CIME import namelist

    configuration = namelist.parse(in_file=str(path), groupless=groupless)
    if groupless:
        return configuration
    values = {}
    relevant = {"streams", "fudge_noleap_for_gregorian", "input_calendar_override_selectors",
                "input_calendar_override_calendars", "fudge_gregorian_to_noleap", "calendar"}
    for group in configuration.get_group_names():
        for variable in configuration.get_variable_names(group):
            if variable not in relevant:
                continue
            if variable in values:
                raise ValueError("Ambiguous namelist variable=" + variable)
            values[variable] = configuration.get_variable_value(group, variable)
    return values


def decoded_strings(values):
    from CIME.namelist import character_literal_to_string

    return [character_literal_to_string(value) for value in values]


def read_policy(path, groupless=False, model_calendar=None):
    values = read_namelist(path, groupless)
    flag = values.get("fudge_noleap_for_gregorian", [".false."])
    reverse = values.get("fudge_gregorian_to_noleap", [".false."])
    if len(flag) != 1 or flag[0].lower() not in {".true.", ".false."}:
        raise ValueError("Invalid fudge_noleap_for_gregorian namelist value")
    if len(reverse) != 1 or reverse[0].lower() not in {".true.", ".false."}:
        raise ValueError("Invalid fudge_gregorian_to_noleap namelist value")
    enabled = flag[0].lower() == ".true."
    reverse_enabled = reverse[0].lower() == ".true."
    if enabled and reverse_enabled:
        raise ValueError("INPUT_CALENDAR: input calendar policies are mutually exclusive")
    selectors = decoded_strings(values.get("input_calendar_override_selectors", ["''"]))
    calendars = decoded_strings(values.get("input_calendar_override_calendars", ["''"]))
    rules = make_rules(selectors, calendars)
    if not (enabled or reverse_enabled):
        raise ValueError("INPUT_CALENDAR: preflight requires an enabled input calendar policy")
    declared_model = decoded_strings(values["calendar"])[0] if "calendar" in values else None
    expected = "NO_LEAP" if reverse_enabled else "GREGORIAN"
    for candidate in (declared_model, model_calendar):
        if candidate is not None and normalize_calendar(candidate) != expected:
            raise ValueError("INPUT_CALENDAR: selected input policy requires CALENDAR=" + expected)
    if reverse_enabled and declared_model is None and model_calendar is None:
        raise ValueError("GREGORIAN_INPUT: supply --model-calendar NO_LEAP or a generated calendar declaration")
    return {"rules": rules, "direction": "gregorian_to_noleap" if reverse_enabled else "noleap_to_gregorian",
            "model_calendar": expected, "model_calendar_verified": declared_model is not None or model_calendar is not None}


def stream_inventory(identifier, rundir):
    path = Path(identifier)
    if not path.is_absolute():
        path = rundir / path
    text = path.read_text()
    try:
        tree = ET.fromstring(text)
    except ET.ParseError:
        tree = ET.fromstring("<stream>" + text + "</stream>")
    fields = tree if tree.tag == "fieldInfo" else tree.find(".//fieldInfo")
    if fields is None:
        raise ValueError("No fieldInfo in stream=" + identifier)
    folder = fields.findtext("filePath", "").strip()
    if ":" in folder:
        raise ValueError("Remote staging needs an explicit local input inventory: " + identifier)
    folder = folder.rstrip("/") + "/" if folder else "./"
    names = fields.findtext("fileNames", "").split()
    if not names:
        raise ValueError("Empty stream inventory=" + identifier)
    if any("$" in name or "%" in name for name in names) or "$" in folder or "%" in folder:
        raise ValueError("Unexpanded stream paths need generated namelists: " + identifier)
    return [{"reader": "stream", "file": folder + name, "stream": identifier} for name in names]


def inspect_input(entry, rules, rundir, direction="noleap_to_gregorian"):
    import netCDF4
    import numpy as np

    report = dict(entry, status="error", declared=None, resolved=None, provenance=None, rules=[])
    try:
        if entry["reader"] not in ADAPTED_READERS:
            raise ValueError("NOLEAP_INPUT: unsupported reader=" + entry["reader"])
        path = Path(entry["file"])
        if not path.is_absolute():
            path = rundir / path
        with netCDF4.Dataset(path) as dataset:
            if "time" not in dataset.variables:
                raise ValueError("NOLEAP_INPUT: input has no time coordinate")
            time = dataset.variables["time"]
            present = "calendar" in time.ncattrs()
            declared = time.getncattr("calendar") if present else ""
            report["declared"] = str(declared) if present else None
            stream = entry.get("stream", "") if entry["reader"] == "stream" else ""
            calendar, origin, matched, flavor = resolve_calendar(declared, present, entry["file"], stream, rules, True)
            report.update(resolved=calendar, provenance=origin, rules=matched, source_flavor=flavor)
            adapted = (direction == "noleap_to_gregorian" and calendar == "NO_LEAP" or
                       direction == "gregorian_to_noleap" and calendar == "GREGORIAN")
            report["interpolation_calendar"] = "GREGORIAN" if adapted else calendar
            values = np.ma.asarray(time[:], dtype="f8")
            if (values.ndim != 1 or not values.size or np.ma.getmaskarray(values).any() or
                    not np.isfinite(values).all() or not (np.diff(values) > 0).all()):
                raise ValueError("NOLEAP_INPUT: time coordinate must be nonempty, finite and strictly increasing")
            if entry["reader"] == "eamxx" and adapted and (values < 0).any():
                raise ValueError("NOLEAP_INPUT: EAMxx requires nonnegative no-leap offsets")
            if "units" not in time.ncattrs():
                raise ValueError("NOLEAP_INPUT: explicit CF time units are required by preflight")
            units = time.getncattr("units")
            source_calendar = "noleap" if calendar == "NO_LEAP" else "proleptic_gregorian"
            check_mixed_calendar_range(flavor, units, values)
            bounds_name = getattr(time, "bounds", "time_bnds" if entry["reader"] == "chemistry" else "")
            if bounds_name and (bounds_name in dataset.variables or "bounds" in time.ncattrs()):
                if bounds_name not in dataset.variables:
                    raise ValueError("INPUT_CALENDAR: missing declared time bounds")
                bounds = np.ma.asarray(dataset.variables[bounds_name][:], dtype="f8")
                if not bounds.size or np.ma.getmaskarray(bounds).any() or not np.isfinite(bounds).all():
                    raise ValueError("INPUT_CALENDAR: invalid time bounds")
                check_mixed_calendar_range(flavor, units, bounds)
            decoded = netCDF4.num2date(values, units, calendar=source_calendar, only_use_cftime_datetimes=True)
            if "date" in dataset.variables:
                dates = np.ma.asarray(dataset.variables["date"][:])
                seconds = np.ma.asarray(dataset.variables["datesec"][:]) if "datesec" in dataset.variables else np.zeros(values.size)
                if (dates.shape != values.shape or seconds.shape != values.shape or
                        np.ma.getmaskarray(dates).any() or np.ma.getmaskarray(seconds).any()):
                    raise ValueError("NOLEAP_INPUT: invalid explicit date coordinate shape or missing values")
                for date, second, numeric in zip(dates, seconds, decoded):
                    if not np.isfinite(date) or date < 0 or date != int(date) or not np.isfinite(second) or not 0 <= second < 86400:
                        raise ValueError("NOLEAP_INPUT: invalid explicit date/datesec")
                    year, monthday = divmod(int(date), 10000)
                    month, day = divmod(monthday, 100)
                    explicit = type(numeric)(year, month, day)
                    difference = (numeric - explicit).total_seconds() - float(second)
                    if abs(difference) > 1.e-4:
                        raise ValueError("NOLEAP_INPUT: explicit dates contradict numeric time coordinate")
            if adapted:
                have_point_method = False
                have_mean_method = False
                allowed = {"point", "instantaneous", "mean"} if entry["reader"] == "chemistry" else {"point", "instantaneous"}
                for variable in dataset.variables.values():
                    methods = str(getattr(variable, "cell_methods", "")).lower()
                    declarations = re.findall(r"\btime\s*:\s*([a-z_]*)", methods)
                    if len(declarations) > 1 or any(method not in allowed for method in declarations):
                        raise ValueError("NOLEAP_INPUT: unsupported or ambiguous time cell_methods")
                    have_point_method = have_point_method or any(method in {"point", "instantaneous"} for method in declarations)
                    have_mean_method = have_mean_method or "mean" in declarations
                if entry["reader"] == "chemistry":
                    bounds_name = getattr(time, "bounds", "time_bnds")
                    if have_mean_method and bounds_name not in dataset.variables:
                        raise ValueError("NOLEAP_INPUT: chemistry time means require mapped time bounds")
                    if bounds_name in dataset.variables:
                        bounds_variable = dataset.variables[bounds_name]
                        if time.dimensions[0] not in bounds_variable.dimensions:
                            raise ValueError("NOLEAP_INPUT: chemistry bounds lack the time dimension")
                        bounds = np.moveaxis(np.ma.asarray(bounds_variable[:], dtype="f8"),
                                             bounds_variable.dimensions.index(time.dimensions[0]), 0)
                        if (bounds.shape != (values.size, 2) or np.ma.getmaskarray(bounds).any() or
                                not np.isfinite(bounds).all() or not (bounds[:, 0] < values).all() or
                                not (values < bounds[:, 1]).all()):
                            raise ValueError("NOLEAP_INPUT: invalid chemistry time bounds")
                        netCDF4.num2date(bounds.ravel(), units, calendar=source_calendar)
                elif "bounds" in time.ncattrs() and not have_point_method:
                    raise ValueError("NOLEAP_INPUT: ambiguous bounded stream/EAMxx quantities need explicit point metadata")
            report.update(status="ok", records=int(values.size), first=str(decoded[0]), last=str(decoded[-1]), units=units,
                          leap_day_records=sum(date.month == 2 and date.day == 29 for date in decoded))
    except (OSError, ValueError, TypeError, IndexError, OverflowError) as error:
        report["error"] = str(error) + "; file=" + str(entry["file"])
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy-namelist", type=Path, required=True)
    parser.add_argument("--groupless", action="store_true", help="Read user_nl_cpl rather than generated drv_in")
    parser.add_argument("--model-calendar", choices=("NO_LEAP", "GREGORIAN"))
    parser.add_argument("--rundir", type=Path, default=Path.cwd())
    parser.add_argument("--namelist", type=Path, action="append", default=[], help="Generated component namelist containing streams")
    parser.add_argument("--stream", action="append", default=[], help="Exact stream identifier used by the reader")
    parser.add_argument("--input", nargs=2, action="append", default=[], metavar=("READER", "FILE"),
                        help="Explicit native input; readers: " + ", ".join(sorted(ADAPTED_READERS | UNSUPPORTED_READERS)))
    options = parser.parse_args(argv)
    result = {"inventory_complete": False, "notice": "Only explicitly supplied inputs and streams are audited; native inventory gaps remain.",
              "inputs": [], "errors": []}
    try:
        policy = read_policy(options.policy_namelist, options.groupless, options.model_calendar)
        rules = policy.pop("rules")
        result.update(policy)
        if policy["direction"] == "gregorian_to_noleap":
            result["leap_day_policy"] = "Retain Gregorian source brackets; the no-leap model never visits February 29."
        streams = list(options.stream)
        for path in options.namelist:
            values = read_namelist(path)
            for descriptor in decoded_strings(values.get("streams", [])):
                arguments = descriptor.split()
                if len(arguments) != 4:
                    raise ValueError("Invalid configured stream descriptor=" + descriptor)
                streams.append(arguments[0])
        entries = [{"reader": reader, "file": filename, "stream": ""} for reader, filename in options.input]
        for stream in dict.fromkeys(streams):
            try:
                entries.extend(stream_inventory(stream, options.rundir))
            except (OSError, ValueError, ET.ParseError) as error:
                result["errors"].append(str(error))
        if not entries:
            result["errors"].append("No configured inputs were inventoried")
        result["inputs"] = [inspect_input(entry, rules, options.rundir, policy["direction"]) for entry in entries]
        for stream in dict.fromkeys(entry.get("stream", "") for entry in entries):
            if not stream:
                continue
            calendars = {entry["resolved"] for entry in result["inputs"] if entry.get("stream") == stream and entry["status"] == "ok"}
            if len(calendars) > 1:
                result["errors"].append("NOLEAP_INPUT: mixed source calendars in stream=" + stream)
        used = {rule for entry in result["inputs"] for rule in entry["rules"]}
        result["unmatched_rules"] = sorted(set(rules) - used)
        if result["unmatched_rules"]:
            result["errors"].append("Explicit calendar rules matched no audited input")
    except (OSError, ValueError, ImportError) as error:
        result["errors"].append(str(error))
    print(json.dumps(result, indent=2))
    return int(bool(result["errors"]) or any(entry["status"] != "ok" for entry in result["inputs"]))


if __name__ == "__main__":
    sys.exit(main())