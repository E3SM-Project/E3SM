#!/usr/bin/env bash
#
# Make sure the branch of every PR recently merged into master is also in next.
#
# Meant to be run (by .github/workflows/merge-to-next.yml) from a clone of the repo, with
# the git credentials of the clone allowed to push to next, and GH_TOKEN set for `gh`.
#
# For each PR merged into master during the last LOOKBACK_DAYS days (plus PR_NUM, if set),
# oldest first:
#   - if the merged branch head is already in next, do nothing;
#   - otherwise merge it into next (--no-ff), with the same message body as the merge
#     commit on master and "into next" in the first line;
#   - on conflicts, skip that PR and open an issue (once per PR).
# All the merges are pushed to next in one go; if the push is rejected (e.g., someone
# pushed to next by hand in the meantime), everything is redone on top of the new next.
#
# Because every run looks at *all* recent PRs (not just the one that triggered it), a run
# that gets dropped from the concurrency queue, or a missed event, is harmless: the next
# run catches up.
#
# Environment:
#   GITHUB_REPOSITORY  owner/name                                     (required)
#   PR_NUM             PR number that triggered the run, if any
#   LOOKBACK_DAYS      how far back to look for merged PRs            (default: 3)
#   DRY_RUN            1: merge locally, but do not push or open issues (default: 0)
#   MAX_PUSH_ATTEMPTS  attempts at pushing to next                    (default: 3)
#   RUN_URL            link to the workflow run, used in issues
#
# Exit codes: 0 all good; 3 some PRs failed to merge, and issues were opened for them;
#             anything else: unexpected error.

set -euo pipefail

REPO=${GITHUB_REPOSITORY:?GITHUB_REPOSITORY is not set}
PR_NUM=${PR_NUM:-}
LOOKBACK_DAYS=${LOOKBACK_DAYS:-3}
DRY_RUN=${DRY_RUN:-0}
MAX_PUSH_ATTEMPTS=${MAX_PUSH_ATTEMPTS:-3}
RUN_URL=${RUN_URL:-}
ISSUE_LABEL="merge-to-next"
EXIT_REPORTED=3

TMP=$(mktemp -d "${RUNNER_TEMP:-${TMPDIR:-/tmp}}/merge-to-next.XXXXXX")
WT="${TMP}/next"       # worktree where next is built
PRS="${TMP}/prs.json"  # merged PRs to look at, oldest first

log() { echo "[merge-to-next] $*"; }

git config user.name  "github-actions[bot]"
git config user.email "41898282+github-actions[bot]@users.noreply.github.com"

# ---------------------------------------------------------------------------------------
# Collect the PRs to look at
# ---------------------------------------------------------------------------------------
collect_prs() {
  local fields="number,headRefName,headRefOid,mergeCommit,mergedAt,mergedBy,author,body"
  local since
  since=$(date -u -d "${LOOKBACK_DAYS} days ago" +%F)

  gh pr list --repo "${REPO}" --state merged --base master \
     --search "merged:>=${since}" --limit 200 --json "${fields}" > "${TMP}/list.json"

  # The search index can lag a bit behind: always include the PR that triggered the run
  if [ -n "${PR_NUM}" ]; then
    gh pr view "${PR_NUM}" --repo "${REPO}" --json "${fields}" > "${TMP}/one.json"
  else
    echo '{}' > "${TMP}/one.json"
  fi

  jq -s '(.[0] + (if (.[1] | has("number")) then [.[1]] else [] end))
         | unique_by(.number) | map(select(.mergedAt != null)) | sort_by(.mergedAt)' \
     "${TMP}/list.json" "${TMP}/one.json" > "${PRS}"
}

# ---------------------------------------------------------------------------------------
# Merge a single PR into the (worktree of) next
# ---------------------------------------------------------------------------------------
declare -a MERGED_PRS=()
declare -a FAILED_PRS=()
declare -A FAIL_KIND FAIL_DETAIL FAIL_SRC FAIL_SUBJECT FAIL_BRANCH FAIL_WHO FAIL_AUTHOR

record_failure() { # n kind detail
  local n=$1
  FAILED_PRS+=("${n}")
  FAIL_KIND[$n]=$2
  FAIL_DETAIL[$n]=$3
  FAIL_SRC[$n]=${src:-}
  FAIL_SUBJECT[$n]=${subject:-}
  FAIL_BRANCH[$n]=${branch:-}
  FAIL_WHO[$n]=${integrator:-}
  FAIL_AUTHOR[$n]=${author:-}
}

process_pr() {
  local pr=$1
  local n head mc p2 use_mc=0 body
  n=$(jq -r '.number' <<<"${pr}")
  head=$(jq -r '.headRefOid' <<<"${pr}")
  mc=$(jq -r '.mergeCommit.oid // empty' <<<"${pr}")
  # These are also read by record_failure
  branch=$(jq -r '.headRefName' <<<"${pr}")
  integrator=$(jq -r '.mergedBy.login // empty' <<<"${pr}")
  author=$(jq -r '.author.login // empty' <<<"${pr}")
  src=""
  subject="Merge branch '${branch}' into next (PR #${n})"

  # The PR head may be missing (e.g., squash/rebase merge of a PR from a fork)
  if ! git cat-file -e "${head}^{commit}" 2>/dev/null; then
    git fetch -q origin "refs/pull/${n}/head" 2>/dev/null || true
  fi
  if ! git cat-file -e "${head}^{commit}" 2>/dev/null; then
    record_failure "${n}" error "Could not find the head commit (${head}) of the PR branch."
    return
  fi

  # What did master actually get? If the PR was merged with a true merge commit, take the
  # second parent of that commit. Otherwise (squash, rebase, odd merges), use the PR head.
  # Beware that the PR head can itself be a merge commit (e.g., master merged into the
  # branch), hence the extra sanity checks on the merge commit.
  src=${head}
  if [ -n "${mc}" ] && [ "${mc}" != "${head}" ] \
     && git cat-file -e "${mc}^{commit}" 2>/dev/null \
     && [ "$(git rev-list --parents -n1 "${mc}" | wc -w)" = 3 ] \
     && git merge-base --is-ancestor "${mc}" origin/master \
     && p2=$(git rev-parse "${mc}^2") \
     && git merge-base --is-ancestor "${p2}" "${head}"; then
    src=${p2}
    use_mc=1
  fi

  if git -C "${WT}" merge-base --is-ancestor "${src}" HEAD; then
    log "PR #${n}: already in next"
    return
  fi

  # Commit message: new first line, same body as the master merge commit. If there is no
  # such commit, use the part of the PR description above the first '---' line.
  if [ "${use_mc}" = 1 ]; then
    body=$(git log -1 --format=%b "${mc}")
  else
    body=$(jq -r '.body // ""' <<<"${pr}" | tr -d '\r' | awk '/^---[[:space:]]*$/ {exit} {print}')
  fi
  {
    echo "${subject}"
    if [ -n "${body//[[:space:]]/}" ]; then
      echo
      echo "${body}"
    fi
  } > "${TMP}/msg-${n}.txt"

  # Keep the authorship of the master merge (the integrator), if there is one
  local -a authorship=()
  if [ "${use_mc}" = 1 ]; then
    authorship=("GIT_AUTHOR_NAME=$(git log -1 --format=%an "${mc}")"
                "GIT_AUTHOR_EMAIL=$(git log -1 --format=%ae "${mc}")"
                "GIT_AUTHOR_DATE=$(git log -1 --format=%aI "${mc}")")
  fi

  local out conflicts
  if out=$(env "${authorship[@]}" git -C "${WT}" merge --no-ff -F "${TMP}/msg-${n}.txt" "${src}" 2>&1); then
    log "PR #${n}: merged into next"
    MERGED_PRS+=("${n}")
  else
    conflicts=$(git -C "${WT}" diff --name-only --diff-filter=U || true)
    git -C "${WT}" merge --abort 2>/dev/null || git -C "${WT}" reset -q --hard HEAD
    if [ -n "${conflicts}" ]; then
      log "PR #${n}: conflicts in next"
      record_failure "${n}" conflict "${conflicts}"
    else
      log "PR #${n}: merge into next failed"
      record_failure "${n}" error "${out}"
    fi
  fi
}

# Build everything on top of the current origin/next
sweep() {
  MERGED_PRS=()
  FAILED_PRS=()
  FAIL_KIND=(); FAIL_DETAIL=(); FAIL_SRC=(); FAIL_SUBJECT=(); FAIL_BRANCH=(); FAIL_WHO=(); FAIL_AUTHOR=()
  git fetch -q origin next
  git -C "${WT}" reset -q --hard origin/next
  local pr
  while read -r pr; do
    process_pr "${pr}"
  done < <(jq -c '.[]' "${PRS}")
}

# ---------------------------------------------------------------------------------------
# Issues
# ---------------------------------------------------------------------------------------
issue_exists() { # title
  local titles
  titles=$(gh issue list --repo "${REPO}" --label "${ISSUE_LABEL}" --state all --limit 500 \
             --json title --jq '.[].title' 2>/dev/null || true)
  grep -Fxq -- "$1" <<<"${titles}"
}

issue_body() { # n
  local n=$1
  local url="https://github.com/${REPO}/pull/${n}"
  echo "PR ${url} (branch \`${FAIL_BRANCH[$n]}\`) was merged into \`master\`, but merging it into \`next\` failed."
  if [ -n "${FAIL_AUTHOR[$n]}" ]; then
    echo
    echo "cc @${FAIL_AUTHOR[$n]}"
  fi
  echo
  if [ "${FAIL_KIND[$n]}" = conflict ]; then
    echo "Conflicts in:"
    echo
    echo '```'
    echo "${FAIL_DETAIL[$n]}"
    echo '```'
    echo
    echo "Most likely a PR that is in \`next\` but not yet in \`master\` touches the same files."
  else
    echo "Error:"
    echo
    echo '```'
    echo "${FAIL_DETAIL[$n]}"
    echo '```'
  fi
  echo
  echo "Until this is fixed, \`next\` is missing a PR that is in \`master\`. Please merge the PR"
  echo "branch into \`next\` by hand, resolving the conflicts, with the same message body as the"
  echo "\`master\` merge, and \`into next\` in the first line:"
  echo
  echo '```'
  echo "git fetch origin && git checkout -B next origin/next"
  echo "git merge --no-ff ${FAIL_SRC[$n]}   # first line: ${FAIL_SUBJECT[$n]}"
  echo "git push origin next"
  echo '```'
  if [ -n "${RUN_URL}" ]; then
    echo
    echo "Workflow run: ${RUN_URL}"
  fi
}

NEW_FAILURES=0
report_failures() {
  local n title
  [ ${#FAILED_PRS[@]} -gt 0 ] || return 0
  [ "${DRY_RUN}" = 1 ] || gh label create "${ISSUE_LABEL}" --repo "${REPO}" --color B60205 \
      --description "Automatic merge of a PR into next failed" >/dev/null 2>&1 || true
  for n in "${FAILED_PRS[@]}"; do
    title="Automatic merge of PR #${n} into next failed"
    if [ "${DRY_RUN}" = 1 ]; then
      log "DRY RUN: would open issue '${title}'"
      issue_body "${n}"
      NEW_FAILURES=$((NEW_FAILURES + 1))
      continue
    fi
    if issue_exists "${title}"; then
      log "PR #${n}: issue already exists, not opening another"
      continue
    fi
    issue_body "${n}" > "${TMP}/issue-${n}.md"
    # Fall back to fewer options if the label or the assignee is not accepted
    gh issue create --repo "${REPO}" --title "${title}" --body-file "${TMP}/issue-${n}.md" \
        --label "${ISSUE_LABEL}" ${FAIL_WHO[$n]:+--assignee "${FAIL_WHO[$n]}"} \
    || gh issue create --repo "${REPO}" --title "${title}" --body-file "${TMP}/issue-${n}.md" \
        --label "${ISSUE_LABEL}" \
    || gh issue create --repo "${REPO}" --title "${title}" --body-file "${TMP}/issue-${n}.md"
    NEW_FAILURES=$((NEW_FAILURES + 1))
  done
}

# ---------------------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------------------
git fetch -q origin next master
git worktree add -q -B next "${WT}" origin/next
collect_prs

attempt=1
while :; do
  sweep
  if [ ${#MERGED_PRS[@]} -eq 0 ]; then
    break
  fi
  if [ "${DRY_RUN}" = 1 ]; then
    log "DRY RUN: would push PRs ${MERGED_PRS[*]} to next"
    break
  fi
  if git -C "${WT}" push -q origin HEAD:next; then
    log "Pushed to next: PRs ${MERGED_PRS[*]}"
    break
  fi
  if [ "${attempt}" -ge "${MAX_PUSH_ATTEMPTS}" ]; then
    echo "::error::Could not push to next after ${attempt} attempts"
    exit 1
  fi
  attempt=$((attempt + 1))
  log "Push to next rejected, retrying (attempt ${attempt}/${MAX_PUSH_ATTEMPTS})"
  sleep $((RANDOM % 10 + 5))
done

report_failures

if [ -n "${GITHUB_STEP_SUMMARY:-}" ]; then
  {
    echo "### merge-to-next"
    echo "- merged into next: ${MERGED_PRS[*]:-none}"
    echo "- failed: ${FAILED_PRS[*]:-none}"
  } >> "${GITHUB_STEP_SUMMARY}"
fi

if [ "${NEW_FAILURES}" -gt 0 ]; then
  exit ${EXIT_REPORTED}
fi
