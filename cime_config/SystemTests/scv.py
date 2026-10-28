"""Shortwave TOA/heating, consumption, transfer and restart qualification."""

from datetime import datetime
import gzip
import hashlib
import json
import math
from pathlib import Path

from CIME.SystemTests.ers import ERS


def exchange_records(logs):
    records = []
    for line in logs.splitlines():
        if "SOLAR_EXCHANGE_V1" not in line:
            continue
        fields = line.split("SOLAR_EXCHANGE_V1", 1)[1].split()
        if len(fields) != 7:
            raise AssertionError("Malformed SOLAR_EXCHANGE_V1 record")
        records.append((int(fields[0]), float(fields[1]),
                        [float(value.replace("D", "E")) for value in fields[2:]]))
    return records


def schedule_records(logs):
    records = []
    for line in logs.splitlines():
        if "SOLAR_SCHEDULE_V1" not in line:
            continue
        fields = line.split("SOLAR_SCHEDULE_V1", 1)[1].split()
        if len(fields) != 8 or fields[2] not in ("T", "F"):
            raise AssertionError("Malformed SOLAR_SCHEDULE_V1 record")
        records.append((int(fields[0]), float(fields[1]), fields[2] == "T",
                        int(fields[3]), *[float(value.replace("D", "E")) for value in fields[4:]]))
    return records


def check_schedule(records, calendar, exchanges, interval, width, restart_version=0):
    """Check EAM's two-step driver-clock lead, separately from averaging width."""
    if len(records) != len(exchanges) or not records:
        raise AssertionError("Missing solar schedule records")
    valid_events = 0
    for record, exchange in zip(records, exchanges):
        ymd, tod, valid, next_ymd, next_tod, cday, dt_avg, version = record
        if (ymd, tod) != exchange[:2] or dt_avg != width or version != restart_version:
            raise AssertionError("Solar schedule interval or provenance differs")
        if not valid:
            continue
        current = calendar_seconds(ymd, tod, calendar)
        scheduled = calendar_seconds(next_ymd, next_tod, calendar)
        year_start = calendar_seconds(next_ymd // 10000 * 10000 + 101, 0, calendar)
        expected_cday = 1 + (scheduled - year_start) / 86400
        if not math.isfinite(cday) or abs(cday - expected_cday) > 1e-10:
            raise AssertionError("Solar schedule civil date and calendar day disagree")
        if scheduled - current != 2 * interval:
            raise AssertionError("Solar albedo validity differs from EAM's future event")
        valid_events += 1
    if not valid_events:
        raise AssertionError("No valid future solar events were checked")
    if width == interval and valid_events != len(records):
        raise AssertionError("Every-step radiation requires a valid schedule at every exchange")
    if width == 2 * interval:
        flags = [record[2] for record in records[4:]]
        if any(previous == current for previous, current in zip(flags, flags[1:])):
            raise AssertionError("Two-step radiation schedule does not alternate after startup")
    return valid_events


def diagnostic_fields(logs, tag, count):
    for line in logs.splitlines():
        if tag not in line:
            continue
        fields = line.split(tag, 1)[1].split()
        if len(fields) != count:
            raise AssertionError("Malformed " + tag + " record")
        yield fields


def column_records(logs):
    return [(int(row[0]), float(row[1]), int(row[2]), int(row[3]),
             *map(float, row[4:])) for row in diagnostic_fields(logs, "SOLAR_COLUMN_V1", 10)]


def consumption_records(logs):
    return [(int(row[0]), float(row[1]), int(row[2]), int(row[3]))
            for row in diagnostic_fields(logs, "SOLAR_CONSUME_V1", 4)]


def buffer_records(logs):
    return [(int(row[0]), float(row[1]), int(row[2]), row[3], int(row[4]),
             *map(float, row[5:])) for row in diagnostic_fields(logs, "SOLAR_BUFFER_V1", 9)]


def require_close(residual, scale, absolute=1e-7):
    if not math.isfinite(residual) or not math.isfinite(scale) or abs(residual) > absolute + 1e-10 * scale:
        raise AssertionError("Solar diagnostic closure failed: {} (scale {})".format(residual, scale))


def check_columns(columns, consumed, schedules, exchanges, calendar, interval, width):
    groups = {}
    seen = {}
    for row in columns:
        ymd, tod, chunk, count, residual, toa, surface, heating, scale, area = row
        stamp = calendar_seconds(ymd, tod, calendar)
        key = (stamp, chunk)
        if key in seen or chunk < 1 or count < 1:
            raise AssertionError("Duplicate or invalid atmosphere chunk")
        if not all(math.isfinite(value) for value in row[4:]) or area <= 0 or residual < 0 or scale < 0:
            raise AssertionError("Invalid atmosphere column diagnostic")
        require_close(residual, scale)
        require_close(toa - surface - heating, abs(toa) + abs(surface) + abs(heating), 1e-7 * area)
        seen[key] = count
        groups.setdefault(stamp, []).append(row)
    expected_times = {calendar_seconds(row[0], row[1], calendar) + interval for row in exchanges}
    if not expected_times or set(groups) != expected_times:
        raise AssertionError("Atmosphere heating does not cover every exchange interval")
    energy = [0., 0., 0.]
    chunk_sets = []
    for exchange in exchanges:
        stamp = calendar_seconds(exchange[0], exchange[1], calendar) + interval
        rows = groups[stamp]
        require_close(math.fsum(row[9] for row in rows) - 1., 1., 1e-10)
        chunk_sets.append({row[2] for row in rows})
        totals = [math.fsum(row[index] for row in rows) for index in (5, 6, 7)]
        require_close(totals[1] + exchange[2][0], abs(totals[1]) + abs(exchange[2][0]))
        energy = [total + value * interval for total, value in zip(energy, totals)]
    if any(chunks != chunk_sets[0] for chunks in chunk_sets):
        raise AssertionError("Atmosphere chunk coverage changed within a run")
    consumed_keys = set()
    consumed_by_time = {}
    for ymd, tod, chunk, count in consumed:
        key = (calendar_seconds(ymd, tod, calendar), chunk)
        if key in consumed_keys or seen.get(key) != count:
            raise AssertionError("Missing, duplicate or unmatched albedo consumption")
        consumed_keys.add(key)
        consumed_by_time.setdefault(key[0], set()).add(chunk)
    expected_events = {calendar_seconds(row[3], row[4], calendar)
                       for row in schedules if row[2]} & expected_times
    observed_events = set(consumed_by_time)
    if width == interval:
        if observed_events != expected_times:
            raise AssertionError("Every-step radiation lacks an albedo consumption event")
    elif (observed_events - expected_events) - {min(expected_times)} or expected_events - observed_events:
        raise AssertionError("Consumed albedo does not match the published future events")
    for stamp in observed_events:
        if consumed_by_time[stamp] != chunk_sets[0]:
            raise AssertionError("Albedo consumption does not cover every atmosphere chunk")
    return {"column_records": len(columns), "albedo_consumption_records": len(consumed),
            "toa_energy_j_m2": energy[0], "atmosphere_heating_j_m2": energy[2],
            "atmosphere_surface_j_m2": energy[1]}


def check_buffers(records, exchanges, calendar, interval, ocean_interval):
    if not records or ocean_interval < interval or ocean_interval % interval:
        raise AssertionError("Missing or unsupported ocean buffer diagnostics")
    ratio = int(ocean_interval / interval)
    incoming_by_time = {}
    previous_after = None
    previous_count = None
    previous_stamp = None
    first_before = records[0][5]
    incoming_total = 0.
    delivered_total = 0.
    deliveries = 0
    full_deliveries = 0
    for ymd, tod, instance, operation, count, before, incoming, delivered, after in records:
        stamp = calendar_seconds(ymd, tod, calendar)
        if instance != 1 or count < 0 or count > ratio or operation not in ("ADD", "SEND"):
            raise AssertionError("Invalid ocean buffer instance, operation or sample count")
        if not all(math.isfinite(value) for value in (before, incoming, delivered, after)):
            raise AssertionError("Nonfinite ocean buffer inventory")
        if previous_stamp is not None and stamp < previous_stamp:
            raise AssertionError("Ocean buffer records are not chronological")
        scale = abs(before) + abs(incoming) + abs(delivered) + abs(after)
        require_close(before + incoming - delivered - after, scale)
        if previous_after is not None:
            require_close(before - previous_after, abs(before) + abs(previous_after))
        if operation == "ADD":
            if count == 0 or stamp in incoming_by_time or delivered != 0:
                raise AssertionError("Invalid or duplicated ocean accumulation")
            if previous_count is not None and count != previous_count + 1:
                raise AssertionError("Ocean accumulation counter skipped a sample")
            incoming_by_time[stamp] = incoming
            previous_count = count
        else:
            if incoming != 0 or after != 0 or (previous_count is not None and count != previous_count):
                raise AssertionError("Ocean delivery counter or reset differs")
            if deliveries and count != ratio:
                raise AssertionError("Ocean delivery has the wrong accumulated duration")
            deliveries += 1
            full_deliveries += count == ratio
            previous_count = 0
        incoming_total += incoming
        delivered_total += delivered
        previous_after = after
        previous_stamp = stamp
    expected = {calendar_seconds(row[0], row[1], calendar): row[2][2] * interval for row in exchanges}
    if set(expected) != set(incoming_by_time) or not full_deliveries:
        raise AssertionError("Ocean buffer does not cover the exchange intervals and a full delivery")
    for stamp, value in expected.items():
        require_close(incoming_by_time[stamp] - value, abs(value) + abs(incoming_by_time[stamp]))
    change = previous_after - first_before
    require_close(incoming_total - delivered_total - change,
                  abs(incoming_total) + abs(delivered_total) + abs(change), 1e-7 * len(exchanges) * interval)
    return {"ocean_incoming_j_m2": incoming_total, "ocean_delivered_j_m2": delivered_total,
            "ocean_pending_start_j_m2": first_before, "ocean_pending_end_j_m2": previous_after,
            "ocean_pending_change_j_m2": change, "ocean_delivery_records": deliveries}


def log_files(rundir):
    files = {}
    for path in sorted(Path(rundir).glob("*.log.*")):
        key = path.with_suffix("") if path.suffix == ".gz" else path
        files.setdefault(key, []).append(path)
    return files


def log_fingerprint(path):
    digest = hashlib.sha256()
    length = 0
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", errors="replace") as stream:
        for block in iter(lambda: stream.read(65536), ""):
            length += len(block)
            digest.update(block.encode("utf-8"))
    return length, digest.digest()


def coupler_logs(rundir):
    """Snapshot decoded lengths and hashes, not whole rank-log contents."""
    snapshots = {}
    for key, paths in log_files(rundir).items():
        snapshots[key] = log_fingerprint(paths[0])
        if any(log_fingerprint(path) != snapshots[key] for path in paths[1:]):
            raise AssertionError("Conflicting compressed and uncompressed coupler logs")
    return snapshots


def fresh_log_lines(rundir, before):
    """Verify each saved prefix before streaming new records, including gzip logs."""
    for key, paths in log_files(rundir).items():
        old_length, old_digest = before.get(key, (0, hashlib.sha256().digest()))
        digest = hashlib.sha256()
        length = old_length
        opener = gzip.open if paths[0].suffix == ".gz" else open
        with opener(paths[0], "rt", errors="replace") as stream:
            remaining = old_length
            while remaining:
                block = stream.read(min(remaining, 65536))
                if not block:
                    raise AssertionError("Coupler log was replaced during SCV: " + str(key))
                digest.update(block.encode("utf-8"))
                remaining -= len(block)
            if digest.digest() != old_digest:
                raise AssertionError("Coupler log was replaced during SCV: " + str(key))
            for line in stream:
                digest.update(line.encode("utf-8"))
                length += len(line)
                yield line
        if any(log_fingerprint(path) != (length, digest.digest()) for path in paths[1:]):
            raise AssertionError("Conflicting compressed and uncompressed coupler logs")


def read_solar_records(lines):
    parsers = (("SOLAR_EXCHANGE_V1", exchange_records), ("SOLAR_SCHEDULE_V1", schedule_records),
               ("SOLAR_COLUMN_V1", column_records), ("SOLAR_CONSUME_V1", consumption_records),
               ("SOLAR_BUFFER_V1", buffer_records))
    records = [[] for parser in parsers]
    for line in lines:
        for rows, (tag, parser) in zip(records, parsers):
            if tag in line:
                rows.extend(parser(line))
    return records


def calendar_seconds(ymd, tod, calendar):
    """Use a civil-date oracle independent of the model's calendar routines."""
    stamp = datetime.strptime(str(ymd), "%Y%m%d")
    if not math.isfinite(tod) or not 0 <= tod < 86400:
        raise ValueError("Invalid solar budget time of day")
    if calendar == "GREGORIAN":
        days = stamp.toordinal() - 1
    elif calendar == "NO_LEAP":
        if stamp.month == 2 and stamp.day == 29:
            raise ValueError("February 29 in a NO_LEAP solar budget")
        days = (stamp.year - 1) * 365 + stamp.replace(year=2001).timetuple().tm_yday - 1
    else:
        raise ValueError("Unsupported solar budget calendar: " + calendar)
    return days * 86400 + tod


def check_exchange(records, calendar, start, stop, interval, atol=1e-7, rtol=1e-10):
    """Integrate signed, already area-weighted exchange terms in J/m2."""
    if not all(math.isfinite(value) for value in (start, stop, interval, atol, rtol)):
        raise ValueError("Non-finite solar budget configuration")
    if interval <= 0 or stop <= start or atol < 0 or rtol < 0:
        raise ValueError("Invalid solar budget interval or tolerance")
    previous = start
    energies = [0.0] * 5
    samples = 0
    for ymd, tod, terms in records:
        current = calendar_seconds(ymd, tod, calendar)
        if current - previous != interval or current > stop:
            raise AssertionError("Missing, duplicate or shifted solar budget interval")
        if len(terms) != 5 or not all(math.isfinite(value) for value in terms):
            raise AssertionError("Invalid solar exchange terms")
        if terms[0] > atol or any(value < -atol for value in terms[1:]):
            raise AssertionError("Invalid solar exchange sign convention")
        residual = math.fsum(terms)
        scale = math.fsum(abs(value) for value in terms)
        if abs(residual) > atol + rtol * scale:
            raise AssertionError("Solar exchange residual at {}:{}: {} W/m2".format(ymd, tod, residual))
        energies = [total + value * interval for total, value in zip(energies, terms)]
        previous = current
        samples += 1
    if previous != stop or samples < 2:
        raise AssertionError("Incomplete solar budget coverage")
    if energies[0] >= 0:
        raise AssertionError("No solar energy crossed the global surface interface")
    if abs(math.fsum(energies)) > atol * (stop - start) + rtol * math.fsum(map(abs, energies)):
        raise AssertionError("Integrated solar exchange does not close")
    return {"samples": samples, "duration_seconds": stop - start,
            "component_energy_j_m2": energies, "residual_j_m2": math.fsum(energies)}


class SCV(ERS):
    """Run the same conservation assertions on both legs of an ERS test."""

    def run_indv(self, *args, **kwargs):
        case = self._case
        frequency = case.get_value("ATM_NCPL")
        if (case.get_value("COMP_INTERFACE") != "mct" or
                case.get_value("COMPSET") != "WCYCL1850NS" or
                case.get_value("NINST") != 1 or
                case.get_value("NCPL_BASE_PERIOD") != "day" or
                frequency != 48 or case.get_value("STOP_OPTION") != "ndays"):
            raise ValueError("SCV requires single-instance MCT WCYCL1850NS, day-based 48/day coupling and ndays")
        for component in ("LND", "ICE"):
            if case.get_value(component + "_NCPL") != frequency:
                raise ValueError("SCV requires equal atmosphere/land/ice exchange intervals")
        for component in case.get_values("COMP_CLASSES"):
            if component != "CPL":
                component_frequency = case.get_value(component + "_NCPL")
                if component_frequency and component_frequency > frequency:
                    raise ValueError("SCV does not support a driver interval faster than the atmosphere")
        ocean_frequency = case.get_value("OCN_NCPL")
        if not ocean_frequency or ocean_frequency <= 0 or frequency % ocean_frequency:
            raise ValueError("SCV requires ocean exchange frequency to divide the atmosphere frequency")
        if not case.get_value("BUDGETS"):
            raise ValueError("SCV requires BUDGETS=TRUE")
        calendar = case.get_value("CALENDAR")
        interval = 86400 / frequency
        if case.get_value("CONTINUE_RUN"):
            date, seconds = self._rest_time.lstrip(".").rsplit("-", 1)
            start = calendar_seconds(int(date.replace("-", "")), float(seconds), calendar)
        else:
            start = calendar_seconds(int(case.get_value("RUN_STARTDATE").replace("-", "")),
                                     case.get_value("START_TOD"), calendar)
        stop = start + case.get_value("STOP_N") * 86400
        rundir = case.get_value("RUNDIR")
        before = coupler_logs(rundir)
        super().run_indv(*args, **kwargs)
        exchanges, schedules, columns, consumed, buffers = read_solar_records(fresh_log_lines(rundir, before))
        report = check_exchange(exchanges, calendar, start, stop, interval)
        width = interval * (2 if "solar_conservation_rad2" in case.get_value("CASE") else 1)
        report["valid_solar_events"] = check_schedule(
            schedules, calendar, exchanges, interval, width,
            restart_version=1 if case.get_value("CONTINUE_RUN") else 0)
        report.update(check_columns(columns, consumed, schedules, exchanges, calendar, interval, width))
        report.update(check_buffers(buffers, exchanges, calendar, interval, 86400 / ocean_frequency))
        component_energy = report["component_energy_j_m2"]
        full_residual = (report["toa_energy_j_m2"] - report["atmosphere_heating_j_m2"] -
                         math.fsum(component_energy[index] for index in (1, 3, 4)) -
                         report["ocean_delivered_j_m2"] - report["ocean_pending_change_j_m2"])
        require_close(full_residual, math.fsum(abs(value) for value in component_energy) +
                      abs(report["toa_energy_j_m2"]) + abs(report["atmosphere_heating_j_m2"]),
                      1e-7 * (stop - start))
        report["toa_to_delivery_residual_j_m2"] = full_residual
        comparable = {"exchange": exchanges, "schedule": [record[:-1] for record in schedules],
                      "column": sorted(columns), "consumption": sorted(consumed), "buffer": buffers}
        if kwargs.get("suffix", "base") == "rest":
            reference = {name: [record for record in rows
                               if calendar_seconds(record[0], record[1], calendar) >
                               start + (interval if name in ("column", "consumption") else 0)]
                         for name, rows in self._solar_reference.items()}
            if comparable != reference:
                raise AssertionError("Restart solar energy, consumption or buffer records differ from the continuous run")
        else:
            self._solar_reference = comparable
        report["calendar"] = calendar
        report["scope"] = "TOA, atmospheric heating, albedo consumption, surface exchange and ocean buffer energy"
        suffix = kwargs.get("suffix", "base")
        (Path(rundir) / ("solar-conservation." + suffix + ".json")).write_text(
            json.dumps(report, indent=2, allow_nan=False) + "\n")