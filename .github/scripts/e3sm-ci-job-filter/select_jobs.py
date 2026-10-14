#!/usr/bin/env python3
"""
Decide which jobs of a GitHub workflow must run, based on the triggering event.
Used by .github/workflows/e3sm-ci-job-filter.yml (see there for the rules).

Inputs are read from the environment (never interpolated in a shell command), and the
result is written to $GITHUB_OUTPUT (jobs=<json array>, any=<true|false>) and to the
step summary.
"""

import json
import os
import re
import sys
import urllib.request

MAX_FILES = 3000   # The GitHub API does not list more files than this for a PR
PER_PAGE = 100

# Files that are always relevant to a workflow using the filter
SELF_PATHS = [
    ".github/workflows/e3sm-ci-job-filter.yml",
    ".github/scripts/e3sm-ci-job-filter/**",
]

###############################################################################
def parse_names(text):
    """Split on whitespace/commas, preserving order, rejecting duplicates"""
    names = [n for n in re.split(r"[\s,]+", text.strip()) if n]
    if len(set(names)) != len(names):
        raise ValueError(f"Duplicate job names in '{text}'")
    return names

###############################################################################
def glob_to_regex(pattern):
    """Translate a GitHub path pattern to a regex ('**' crosses '/', '*' and '?' do not)"""
    i, out = 0, []
    while i < len(pattern):
        if pattern.startswith("**/", i):
            out.append("(?:.*/)?")
            i += 3
        elif pattern.startswith("**", i):
            out.append(".*")
            i += 2
        elif pattern[i] == "*":
            out.append("[^/]*")
            i += 1
        elif pattern[i] == "?":
            out.append("[^/]")
            i += 1
        else:
            out.append(re.escape(pattern[i]))
            i += 1
    return re.compile("".join(out) + r"\Z")

###############################################################################
def path_matches(path, patterns):
    """GitHub semantics: patterns are scanned in order, '!' excludes, last match wins"""
    included = False
    for pat in patterns:
        negated = pat.startswith("!")
        if glob_to_regex(pat[1:] if negated else pat).match(path):
            included = not negated
    return included

###############################################################################
def list_pr_files(event, token):
    """Names of the files changed by the PR (None if the list may be truncated)"""
    api = os.environ.get("GITHUB_API_URL", "https://api.github.com")
    repo = os.environ["GITHUB_REPOSITORY"]
    number = event["pull_request"]["number"]
    files = []
    for page in range(1, MAX_FILES // PER_PAGE + 1):
        req = urllib.request.Request(
            f"{api}/repos/{repo}/pulls/{number}/files?per_page={PER_PAGE}&page={page}",
            headers={"Authorization": f"Bearer {token}",
                     "Accept": "application/vnd.github+json"})
        with urllib.request.urlopen(req) as resp:
            batch = json.load(resp)
        # A renamed file is relevant through both its old and new name
        for f in batch:
            files.append(f["filename"])
            if f.get("previous_filename"):
                files.append(f["previous_filename"])
        if len(batch) < PER_PAGE:
            return files
    return None

###############################################################################
def select_jobs(all_jobs, event_name, event, skip_drafts=True, dispatch_input="",
                patterns=(), caller_path=None, get_files=None):
    """Returns (selected jobs, reason)"""
    if event_name == "pull_request":
        if skip_drafts and event["pull_request"].get("draft"):
            return [], "draft pull request"
        if patterns:
            files = get_files()
            if files is None:
                return list(all_jobs), "too many changed files to filter; running all jobs"
            relevant = list(SELF_PATHS) + ([caller_path] if caller_path else [])
            if not any(path_matches(f, relevant) or path_matches(f, patterns) for f in files):
                return [], "no changed file matches the paths filter"
            return list(all_jobs), "changed files match the paths filter"
        return list(all_jobs), "pull request"

    if event_name == "workflow_dispatch" and dispatch_input:
        choice = ((event.get("inputs") or {}).get(dispatch_input) or "all").strip()
        if choice == "all":
            return list(all_jobs), f"dispatch input '{dispatch_input}' is 'all'"
        if choice not in all_jobs:
            raise ValueError(f"Dispatch input '{dispatch_input}'='{choice}' is not one of: "
                             f"{', '.join(all_jobs)}, all")
        return [choice], f"dispatch input '{dispatch_input}' selects '{choice}'"

    return list(all_jobs), f"{event_name} event"

###############################################################################
def main():
    env = os.environ
    all_jobs = parse_names(env["INPUT_JOBS"])
    patterns = [p.strip() for p in env.get("INPUT_PATHS", "").splitlines() if p.strip()]
    skip_drafts = env.get("INPUT_SKIP_DRAFTS", "true").lower() == "true"

    # CALLER_WORKFLOW_REF is 'owner/repo/<workflow path>@<ref>'
    caller_path = None
    ref = env.get("CALLER_WORKFLOW_REF", "")
    if ref:
        caller_path = ref.split("@")[0].split("/", 2)[-1]

    event_name = env["GITHUB_EVENT_NAME"]
    with open(env["GITHUB_EVENT_PATH"]) as f:
        event = json.load(f)

    try:
        selected, reason = select_jobs(
            all_jobs, event_name, event, skip_drafts=skip_drafts,
            dispatch_input=env.get("INPUT_DISPATCH_INPUT", ""), patterns=patterns,
            caller_path=caller_path,
            get_files=lambda: list_pr_files(event, env["GH_TOKEN"]))
    except ValueError as e:
        print(f"::error::{e}")
        return 1

    print(f"Selected jobs ({reason}): {selected}")
    with open(env["GITHUB_OUTPUT"], "a") as f:
        f.write(f"jobs={json.dumps(selected)}\n")
        f.write(f"any={'true' if selected else 'false'}\n")

    summary = env.get("GITHUB_STEP_SUMMARY")
    if summary:
        with open(summary, "a") as f:
            f.write(f"### Job selection\n\nReason: {reason}\n\n| Job | Run |\n|---|---|\n")
            for j in all_jobs:
                f.write(f"| {j} | {'yes' if j in selected else 'no'} |\n")
    return 0

###############################################################################
if __name__ == "__main__":
    sys.exit(main())
