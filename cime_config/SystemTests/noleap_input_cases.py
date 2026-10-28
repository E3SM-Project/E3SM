"""Deterministic fixtures and independent checks for no-leap input tests."""

from bisect import bisect_right
from datetime import datetime
import gzip
import json
from pathlib import Path
import xml.etree.ElementTree as ET


CASES = {
    "daily": {},
    "subdaily": {"step": 0.25},
    "monthly": {"monthly": True},
    "gregorian": {"calendar": "gregorian"},
    "off": {"enabled": False},
    "default": {"enabled": False},
    "cycling": {"cycle": True},
    "missing_calendar": {"calendar": None,
                         "error": "NOLEAP_INPUT: explicit time calendar metadata is required"},
    "empty_calendar": {"calendar": "", "error": "NOLEAP_INPUT: empty calendar metadata"},
    "bad_method": {"method": "time: sum", "error": "NOLEAP_INPUT: unsupported time method=sum"},
    "mean": {"method": "time: mean", "error": "NOLEAP_INPUT: unsupported time method=mean"},
    "ambiguous_bounds": {"method": None, "bounds": True,
                         "error": "NOLEAP_INPUT: ambiguous bounded stream data require an interval adapter"},
    "bad_algorithm": {"algorithm": "nearest",
                      "error": "NOLEAP_INPUT: only linear point samples and coszen means are supported"},
    "bad_lower": {"algorithm": "lower",
                  "error": "NOLEAP_INPUT: only linear point samples and coszen means are supported"},
    "bad_upper": {"algorithm": "upper",
                  "error": "NOLEAP_INPUT: only linear point samples and coszen means are supported"},
    "bad_model_calendar": {"error": "NOLEAP_INPUT: fudge_noleap_for_gregorian requires CALENDAR=GREGORIAN"},
    "restart_on": {"enabled": False, "restart_enabled": True,
                   "error": "NOLEAP_INPUT: restart policy differs"},
    "restart_off": {"restart_enabled": False,
                    "error": "NOLEAP_INPUT: restart policy differs"},
}

GREGORIAN_CASES = {name: dict(settings) for name, settings in CASES.items() if name != "gregorian"}
GREGORIAN_CASES.update({
    "noleap": {"calendar": "noleap"},
    "unsupported_calendar": {"calendar": "360_day", "error": "NOLEAP_INPUT: unsupported calendar=360_day"},
    "mixed_calendar": {"mixed_calendar": True, "error": "NOLEAP_INPUT: mixed-calendar stream"},
    "both_flags": {"error": "INPUT_CALENDAR: input calendar policies are mutually exclusive"},
    "override_conflict": {"override": "noleap", "error": "NOLEAP_INPUT: calendar override conflicts with file metadata"},
    "restart_override": {"restart_override": True, "error": "NOLEAP_INPUT: restart calendar overrides differ"},
    "restart_missing": {"restart_marker": "missing", "error": "GREGORIAN_INPUT: restart policy differs"},
    "restart_corrupt": {"restart_marker": "corrupt", "error": "GREGORIAN_INPUT: invalid restart policy"},
})
for _name in ("bad_algorithm", "bad_lower", "bad_upper"):
    GREGORIAN_CASES[_name]["error"] = "GREGORIAN_INPUT: only linear point samples and coszen means are supported"
GREGORIAN_CASES["bad_model_calendar"]["error"] = "GREGORIAN_INPUT: fudge_gregorian_to_noleap requires CALENDAR=NO_LEAP"
for _name in ("restart_on", "restart_off"):
    GREGORIAN_CASES[_name]["error"] = "GREGORIAN_INPUT: restart policy differs"


def selected_case(test_name):
    reverse = "allactive-calendar_gregorian_to_noleap_" in test_name
    cases = GREGORIAN_CASES if reverse else CASES
    prefix = "allactive-calendar_gregorian_to_noleap_" if reverse else "allactive-calendar_noleap_"
    matches = [name for name in cases if prefix + name in test_name.split("--")[-1].split(".")]
    if len(matches) != 1:
        raise ValueError("Exactly one calendar_noleap case modifier is required: " + test_name)
    return matches[0], {"enabled": True, "calendar": "gregorian" if reverse else "noleap", "step": 1.0,
                       "direction": "gregorian_to_noleap" if reverse else "noleap_to_gregorian",
                       "model_calendar": "NO_LEAP" if reverse else "GREGORIAN",
                       "method": "time: point", "algorithm": "linear", **cases[matches[0]]}


def expected_value(records, query, enabled=True):
    """Interpolate civil labels using Python's independent Gregorian clock."""
    if not enabled and query.month == 2 and query.day == 29:
        query = query.replace(day=28)
    def coordinate(stamp):
        if enabled:
            return (stamp - datetime(1, 1, 1)).total_seconds()
        return stamp.year*365*86400 + (stamp.replace(year=2001) - datetime(2001, 1, 1)).total_seconds()

    times = [coordinate(datetime.fromisoformat(record[0])) for record in records]
    query = coordinate(query)
    index = bisect_right(times, query) - 1
    if index < 0 or index >= len(times):
        raise ValueError("Query precedes fixture coverage")
    if times[index] == query:
        return records[index][1]
    if index + 1 == len(times):
        raise ValueError("Query exceeds fixture coverage")
    fraction = (query - times[index]) / (times[index+1] - times[index])
    return records[index][1] + fraction * (records[index+1][1] - records[index][1])


def log_snapshot(rundir):
    snapshot = {}
    for path in Path(rundir).glob("*.log.*"):
        if not path.is_file():
            continue
        stat = path.stat()
        length = stat.st_size
        if path.suffix == ".gz":
            with gzip.open(path, "rb") as stream:
                stream.seek(0, 2)
                length = stream.tell()
        snapshot[path] = (stat.st_mtime_ns, stat.st_size, length)
    return snapshot


def fresh_logs(rundir, before):
    contents = []
    for path, state in log_snapshot(rundir).items():
        previous = before.get(path)
        if previous == state:
            continue
        opener = gzip.open if path.suffix == ".gz" else open
        if previous is None and opener is gzip.open:
            previous = before.get(path.with_suffix(""))
        with opener(path, "rb") as stream:
            if previous is not None:
                stream.seek(previous[2] if state[2] >= previous[2] else 0)
            contents.append(stream.read().decode(errors="replace"))
    return "\n".join(contents)


def require_expected_failure(error, diagnostic, logs):
    if error is None:
        raise AssertionError("Expected rejection did not occur: " + diagnostic)
    if diagnostic not in logs:
        raise AssertionError("Run failed without the expected diagnostic: " + diagnostic) from error


def validate_activation(manifest, logs):
    reverse = manifest.get("direction") == "gregorian_to_noleap"
    banner = "fudge_gregorian_to_noleap is ENABLED" if reverse else "fudge_noleap_for_gregorian is ENABLED"
    if manifest["enabled"] != (banner in logs):
        raise AssertionError("Current-run policy banner does not match the requested policy")
    if manifest["enabled"] and manifest["calendar"] == ("gregorian" if reverse else "noleap"):
        activation = "source=GREGORIAN model=NO_LEAP" if reverse else "source=NO_LEAP model=GREGORIAN"
        if activation not in logs:
            raise AssertionError("No adapted shared-stream activation was logged")
        prefix = "WARNING GREGORIAN_INPUT:" if reverse else "WARNING NOLEAP_INPUT:"
        descriptions = [line for line in logs.splitlines() if prefix + " description=" in line]
        for stream in manifest["streams"]:
            if not any(Path(stream).name in line for line in descriptions):
                raise AssertionError("No adapted shared-stream activation for " + stream)


def _read_stream(path):
    text = path.read_text()
    try:
        return ET.fromstring(text)
    except ET.ParseError:
        return ET.fromstring("<stream>" + text + "</stream>")


def _write_fixture(source, destination, mappings, offsets, epoch, calendar, method, bounds=False, reverse=False):
    import netCDF4

    source_calendar = calendar if calendar in {"noleap", "gregorian"} else ("gregorian" if reverse else "noleap")
    units = "days since " + epoch.strftime("%Y-%m-%d 00:00:00")
    records = []
    with netCDF4.Dataset(source) as original, netCDF4.Dataset(destination, "w") as output:
        source_time = original.variables["time"]
        time_dimension = source_time.dimensions[0]
        for name, dimension in original.dimensions.items():
            output.createDimension(name, None if name == time_dimension else len(dimension))
        time = output.createVariable("time", "f8", (time_dimension,))
        time.units = units
        if calendar is not None:
            time.calendar = calendar
        time[:] = offsets
        if bounds:
            output.createDimension("calendar_noleap_bounds", 2)
            endpoints = output.createVariable("time_bnds", "f8", (time_dimension, "calendar_noleap_bounds"))
            endpoints[:, 0] = offsets - 0.125
            endpoints[:, 1] = offsets + 0.125
            time.bounds = "time_bnds"
        for name, variable in original.variables.items():
            if name == "time" or time_dimension in variable.dimensions:
                continue
            copied = output.createVariable(name, variable.datatype, variable.dimensions)
            copied.setncatts({key: variable.getncattr(key) for key in variable.ncattrs()
                             if key != "_FillValue"})
            copied[:] = variable[:]
        for source_name, target_name in mappings:
            variable = original.variables[source_name]
            if time_dimension not in variable.dimensions:
                raise ValueError("Fixture field is not time dependent: " + source_name)
            axis = variable.dimensions.index(time_dimension)
            copied = output.createVariable(source_name, "f8",
                                           variable.dimensions, zlib=True)
            copied.setncatts({key: variable.getncattr(key) for key in variable.ncattrs()
                             if key not in ("_FillValue", "cell_methods", "scale_factor", "add_offset")})
            if method is not None:
                copied.cell_methods = method
            selection = [slice(None)] * variable.ndim
            selection[axis] = 0
            baseline = variable[tuple(selection)]
            probe = target_name.lower().endswith("tbot")
            if probe:
                copied.units = "K"
            for index, offset in enumerate(offsets):
                selection[axis] = index
                value = 280.0 + 4.0 * ((offset % 365.0) % 7.0)
                stamp = netCDF4.num2date(offset, units, calendar=source_calendar)
                if reverse and stamp.month == 2 and stamp.day == 29:
                    value = 320.0
                copied[tuple(selection)] = value if probe else baseline
                if probe:
                    label = datetime(stamp.year, stamp.month, stamp.day,
                                     stamp.hour, stamp.minute, stamp.second)
                    records.append((label.isoformat(), value))
    return records


def stage_fixtures(case, specification, stop_days):
    """Replace generated stream inputs, keeping their spatial schemas/domains."""
    import netCDF4
    import numpy as np
    from CIME import namelist

    rundir = Path(case.get_value("RUNDIR"))
    start = datetime.fromisoformat(case.get_value("RUN_STARTDATE"))
    reverse = specification.get("direction") == "gregorian_to_noleap"
    source_year = 2000 if reverse and specification.get("cycle") else start.year
    epoch = datetime(source_year, 1, 1)
    units = "days since " + epoch.strftime("%Y-%m-%d 00:00:00")
    calendar = specification["calendar"]
    if calendar not in {"noleap", "gregorian"}:
        calendar = "gregorian" if reverse else "noleap"
    first = netCDF4.date2num(start.replace(year=source_year), units, calendar=calendar)
    offsets = np.arange(max(0, first-2), first+stop_days+4, specification["step"])
    if specification.get("cycle"):
        cycle_days = netCDF4.date2num(datetime(source_year+1, 1, 1), units, calendar=calendar)
        offsets = np.arange(0, cycle_days+specification["step"], specification["step"])
    if specification.get("monthly"):
        offsets = np.array([netCDF4.date2num(datetime(year, month, 1), units, calendar=calendar)
                            for year in range(start.year, start.year+stop_days//365+2)
                            for month in range(1, 13)])
    fixture_dir = rundir / "calendar_noleap_fixtures"
    fixture_dir.mkdir(exist_ok=True)
    manifest = {"streams": [], "solar_streams": [], "records": [], "enabled": specification["enabled"],
                "calendar": specification["calendar"], "datamode": None, "files": [],
                "direction": specification.get("direction", "noleap_to_gregorian"),
                "model_calendar": specification.get("model_calendar", "GREGORIAN"),
                "source_year": source_year, "cycle": specification.get("cycle", False)}
    for path in sorted(rundir.glob("*_in*")):
        if not path.is_file() or path.name.startswith("drv_in"):
            continue
        try:
            configuration = namelist.parse(in_file=str(path))
        except (ValueError, AssertionError):
            continue
        changed = False
        for group in configuration.get_group_names():
            if path.name.startswith("datm_in") and "datamode" in configuration.get_variable_names(group):
                mode = configuration.get_variable_value(group, "datamode")[0].strip("\"'").upper()
                if mode not in ("COPYALL", "CLMNCEP", "CORE2_NYF"):
                    raise ValueError("Unverified temperature-export transformation for DATM mode=" + mode)
                manifest["datamode"] = mode
            if "streams" not in configuration.get_variable_names(group):
                continue
            streams = configuration.get_variable_value(group, "streams")
            rewritten = []
            algorithms = []
            for stream_index, item in enumerate(streams):
                arguments = item.strip("\"'").split()
                stream_path = Path(arguments[0])
                if not stream_path.is_absolute():
                    stream_path = rundir / stream_path
                tree = _read_stream(stream_path)
                fields = tree.find("fieldInfo")
                if fields is None:
                    raise ValueError("No fieldInfo in " + str(stream_path))
                mappings = [line.split() for line in fields.findtext("variableNames").splitlines()
                            if line.strip()]
                solar_fields = [target for _, target in mappings
                                if target.lower().rsplit("_", 1)[-1] in {"swdn", "swndr", "swndf", "swvdr", "swvdf"}]
                solar_only = bool(solar_fields) and len(solar_fields) == len(mappings)
                algorithm = specification["algorithm"]
                if solar_only and algorithm == "linear":
                    algorithm = "coszen"
                    manifest["solar_streams"].append(str(stream_path))
                algorithms.append("'" + algorithm + "'")
                if any(target.lower().endswith("tbot_af") for _, target in mappings):
                    raise ValueError("Temperature anomaly forcing needs a separate export oracle")
                names = fields.findtext("fileNames").split()
                source = Path(fields.findtext("filePath").strip()) / names[0]
                if not source.is_absolute():
                    source = rundir / source
                destination = fixture_dir / (path.name + "_" + str(stream_index) + ".nc")
                file_offsets = offsets[:1] if specification.get("mixed_calendar") else offsets
                records = _write_fixture(source, destination, mappings, file_offsets, epoch,
                                         specification["calendar"], specification["method"],
                                         specification.get("bounds", False), reverse)
                names = [destination.name]
                manifest["files"].append(str(destination))
                if specification.get("mixed_calendar"):
                    second = destination.with_name(destination.stem + "_noleap.nc")
                    records.extend(_write_fixture(source, second, mappings, offsets[1:], epoch, "noleap",
                                                  specification["method"], reverse=reverse))
                    names.append(second.name)
                    manifest["files"].append(str(second))
                if records:
                    if manifest["records"]:
                        raise ValueError("More than one atmospheric temperature probe stream")
                    if specification.get("cycle") and not reverse:
                        mapped = {}
                        for shift in range(stop_days//365+3):
                            for label, value in records:
                                date = datetime.fromisoformat(label)
                                mapped[date.replace(year=date.year+shift).isoformat()] = value
                        records = sorted(mapped.items())
                    manifest["records"] = records
                fields.find("filePath").text = "\n" + str(fixture_dir) + "\n"
                fields.find("fileNames").text = "\n" + "\n".join(names) + "\n"
                offset = fields.find("offset")
                if offset is not None:
                    offset.text = "\n0\n"
                ET.indent(tree)
                ET.ElementTree(tree).write(stream_path, encoding="unicode")
                last_year = netCDF4.num2date(float(offsets[-1]), units, calendar=calendar).year
                if specification.get("cycle"):
                    last_year = source_year
                rewritten.append("'{} {} {} {}'".format(arguments[0], source_year, last_year, start.year))
                manifest["streams"].append(str(stream_path))
            configuration.set_variable_value(group, "streams", rewritten, var_size=len(streams))
            configuration.set_variable_value(group, "tintalgo",
                                             algorithms,
                                             var_size=len(streams))
            mode = "'cycle'" if specification.get("cycle") else "'limit'"
            configuration.set_variable_value(group, "taxmode", [mode] * len(streams),
                                             var_size=len(streams))
            configuration.set_variable_value(group, "offset", ["0"] * len(streams),
                                             var_size=len(streams))
            configuration.set_variable_value(group, "dtlimit", ["1.0e30"] * len(streams),
                                             var_size=len(streams))
            changed = True
        if changed:
            configuration.write(str(path))
    if not manifest["streams"] or not manifest["records"]:
        raise AssertionError("No active prescribed streams/temperature probe were staged")
    if "override" in specification:
        selectors = ["file:" + filename for filename in manifest["files"]]
        found = False
        for path in rundir.glob("drv_in*"):
            if not path.is_file():
                continue
            configuration = namelist.parse(in_file=str(path))
            if "seq_infodata_inparm" not in configuration.get_group_names():
                continue
            for name, entries in (("input_calendar_override_selectors", selectors),
                                  ("input_calendar_override_calendars", [specification["override"]]*len(selectors))):
                configuration.set_variable_value("seq_infodata_inparm", name,
                                                 ["'"+entry+"'" for entry in entries], var_size=len(entries))
            configuration.write(str(path))
            found = True
        if not found:
            raise AssertionError("No generated driver namelist was found for calendar override staging")
    (rundir / "calendar_noleap_manifest.json").write_text(json.dumps(manifest, indent=2))
    return manifest


def mutate_restart_marker(case, action):
    import netCDF4
    from CIME import namelist

    rundir = Path(case.get_value("RUNDIR"))
    candidates = []
    for path in rundir.glob("drv_in*"):
        if not path.is_file():
            continue
        configuration = namelist.parse(in_file=str(path))
        if "seq_infodata_inparm" in configuration.get_group_names():
            pointer = Path(configuration.get_value("restart_pfile")[0].strip("\"'"))
            if not pointer.is_absolute():
                pointer = rundir / pointer
            candidates.append(pointer)
    if len(candidates) != 1:
        raise AssertionError("Expected exactly one generated coupler restart pointer")
    path = Path(candidates[0].read_text().strip())
    if not path.is_absolute():
        path = rundir / path
    if not path.resolve().is_relative_to(rundir.resolve()):
        raise AssertionError("Refuse to mutate a restart outside the test run directory")
    with netCDF4.Dataset(path, "a") as dataset:
        name = "seq_infodata_gregorian_input"
        if name not in dataset.variables:
            raise AssertionError("Reverse restart marker is absent before intentional mutation")
        if action == "missing":
            dataset.renameVariable(name, name+"_removed")
        elif action == "corrupt":
            dataset.variables[name][...] = 99
        else:
            raise ValueError("Unknown restart marker mutation: " + action)


def export_snapshot(rundir):
    return {path: (path.stat().st_mtime_ns, path.stat().st_size)
            for path in Path(rundir).glob("*.cpl*.hi.*.nc") if path.is_file()}


def validate_solar_snapshot(dataset, index, flux_limit=2000.0):
    import numpy as np

    names = ("xaoo_Faox_swdn", "xaoo_Faox_swup", "xaoo_So_anidr", "domo_mask")
    if any(name not in dataset.variables for name in names):
        raise AssertionError("Solar qualification requires instantaneous ocean SW/albedo and domain-mask history")
    values = []
    for name in names:
        variable = dataset.variables[name]
        value = np.ma.asarray(variable[:])
        if "time" in variable.dimensions:
            value = np.take(value, index, axis=variable.dimensions.index("time"))
        values.append(value)
    down, up, albedo, mask = values
    if np.ma.getmaskarray(mask).any() or not np.isfinite(mask.data).all():
        raise AssertionError("Invalid ocean domain mask in solar history")
    active = mask.data > 0
    if not active.any():
        raise AssertionError("No active ocean points in solar history")
    for name, value in zip(names, values):
        if value.shape != mask.shape or np.ma.getmaskarray(value)[active].any() or not np.isfinite(value.data[active]).all():
            raise AssertionError("Missing or non-finite active solar history: " + name)
    tolerance = 1.e-6
    if ((down.data[active] < -tolerance).any() or (down.data[active] > flux_limit).any() or
            (up.data[active] > tolerance).any() or (-up.data[active] > down.data[active]+tolerance).any() or
            (albedo.data[active] < 0).any() or (albedo.data[active] > 1).any()):
        raise AssertionError("Out-of-range ocean shortwave/albedo history")


def validate_exports(rundir, manifest, before=None):
    import netCDF4
    import numpy as np

    samples = 0
    leap_samples = 0
    for path, state in export_snapshot(rundir).items():
        if before is not None and before.get(path) == state:
            continue
        with netCDF4.Dataset(path) as dataset:
            probes = [name for name in dataset.variables if name == "a2x_Sa_tbot"]
            if not probes:
                continue
            time = dataset.variables["time"]
            reverse = manifest.get("direction") == "gregorian_to_noleap"
            if reverse and time.calendar.lower() not in {"noleap", "no_leap", "365_day", "365day"}:
                raise AssertionError("Reverse forcing test history must declare a NO_LEAP calendar")
            dates = netCDF4.num2date(time[:], time.units, calendar=time.calendar)
            for index, stamp in enumerate(dates):
                if manifest.get("solar_streams"):
                    validate_solar_snapshot(dataset,index,manifest.get("solar_flux_limit",2000.0))
                query = datetime(stamp.year, stamp.month, stamp.day,
                                 stamp.hour, stamp.minute, stamp.second)
                enabled = manifest["enabled"] or manifest["calendar"] == "gregorian"
                source_query = query
                if reverse:
                    enabled = manifest["calendar"] == "gregorian"
                    if manifest.get("cycle"):
                        source_query = query.replace(year=manifest["source_year"])
                expected = expected_value(manifest["records"], source_query, enabled)
                field = dataset.variables[probes[0]]
                actual = np.ma.asarray(np.take(field[:], index, axis=field.dimensions.index("time")))
                active = np.ones(actual.shape, dtype=bool)
                if manifest.get("datamode") == "CORE2_NYF":
                    latitude = dataset.variables["doma_lat"]
                    latitudes = latitude[:]
                    if "time" in latitude.dimensions:
                        latitudes = np.take(latitudes, index, axis=latitude.dimensions.index("time"))
                    active &= np.abs(latitudes) <= 60
                if (not active.any() or np.ma.getmaskarray(actual)[active].any() or
                        not np.isfinite(actual.data[active]).all() or
                        not np.allclose(actual.data[active], expected,
                                        rtol=0, atol=1e-4)):
                    raise AssertionError("Wrong forcing export at {} in {}: expected {}".format(
                        query, path, expected))
                samples += 1
                leap_samples += query.month == 2 and query.day == 29
    if samples == 0:
        raise AssertionError("No instantaneous a2x_Sa_tbot forcing exports were found")
    return samples, leap_samples


def set_user_policy(case, enabled, direction="noleap_to_gregorian"):
    from CIME import namelist

    path = Path(case.get_value("CASEROOT")) / "user_nl_cpl"
    values = namelist.parse(in_file=str(path), groupless=True)
    name = "fudge_gregorian_to_noleap" if direction == "gregorian_to_noleap" else "fudge_noleap_for_gregorian"
    values[name] = [".true." if enabled else ".false."]
    configuration = namelist.Namelist({"seq_infodata_inparm": values})
    configuration.write(str(path), format_="nmlcontents")