#!/bin/bash
# Run CIME tests inside an E3SM container image (github.com/E3SM-Project/containers) with
# the image's own CIME machine, e3sm-container. Used by .github/workflows/e3sm-gh-testing.yml
# and adapted from the containers repo's ghci/tests/e3sm-case.sh; see
# docs/dev-guide/gh-hosted-testing.md. Run as the image's default user (uid 1000), under a
# login shell, with these mounted and writable by that user:
#   /projects/e3sm/work/E3SM       E3SM checkout, submodules included
#   /projects/e3sm/data/inputdata  inputdata; missing files are downloaded into it
#   /projects/e3sm/data/baselines  baselines; the machine's BASELINE_ROOT is <this>/$COMPILER
#   /home/e3sm/e3sm_scratch        case directories
# e.g.
# docker run --rm -e TESTS=SMS_Ln5.ne4pg2_oQU480.WCYCL2010NS -e CT_FLAGS="-c -b master" \
#   -v $PWD:/projects/e3sm/work/E3SM -v /path/to/inputdata:/projects/e3sm/data/inputdata \
#   -v /path/to/baselines:/projects/e3sm/data/baselines -v /path/to/scratch:/home/e3sm/e3sm_scratch \
#   ghcr.io/e3sm-project/e3sm-image:gnu13-openmpi4 \
#   bash -l /projects/e3sm/work/E3SM/.github/scripts/e3sm-container-test.sh
#
# TESTS: space-separated test names without machine and compiler (the image's default
# compiler is filled in). CT_FLAGS: baseline options for create_test, e.g. "-c -b master" to
# compare, "-o -g -b master" to generate, empty for neither.
# Results end up in ~/e3sm_scratch/<test>.gh/TestStatus.
set -euo pipefail

: "${TESTS:?TESTS must name at least one test}"
: "${CT_FLAGS:=}"
src=/projects/e3sm/work/E3SM
scratch="$HOME/e3sm_scratch"
step() { echo "== $*"; }

step "CIME machine"
# The env images ship the machine without enabling it; opt in the way a user would.
# shellcheck source=/dev/null
. /projects/e3sm/cime/enable.sh
test "${CIME_MACHINE:-}" = e3sm-container || { echo "CIME_MACHINE is '${CIME_MACHINE:-}', not e3sm-container" >&2; exit 1; }
# The image's default compiler: the first one the machine lists
compiler=$(sed -n 's:.*<COMPILERS>\([^,<]*\).*:\1:p' "$HOME/.cime/config_machines.xml")
echo "machine $CIME_MACHINE, compiler $compiler, flags: ${CT_FLAGS:-none}"
echo "E3SM $(git -C "$src" rev-parse HEAD 2>/dev/null || echo unknown)"
cd "$src/cime/scripts"

# CIME first probes the inputdata server with `wget --spider`, and gives up if that takes
# over 60 s; wget's defaults (15 min read timeout, 20 tries) can spend all of it on one
# stalled connection. Short timeouts and a few tries recover in time.
WGETRC=$(mktemp)
export WGETRC
printf '%s\n' 'timeout = 10' 'tries = 4' 'waitretry = 2' 'retry_connrefused = on' > "$WGETRC"

# TEST.GRID.COMPSET[.TESTMOD] -> TEST.GRID.COMPSET.e3sm-container_<compiler>[.TESTMOD]
tests=()
for t in $TESTS; do
    IFS=. read -r -a p <<< "$t"
    case "${#p[@]}" in
        3) tests+=("$t.${CIME_MACHINE}_$compiler") ;;
        4) tests+=("${p[0]}.${p[1]}.${p[2]}.${CIME_MACHINE}_$compiler.${p[3]}") ;;
        *) echo "test '$t' is not TEST.GRID.COMPSET[.TESTMOD]" >&2; exit 1 ;;
    esac
done

# EAM's P3 builds its lookup table's path itself (micro_p3.F90), so nothing asks CIME to
# download it and the run dies on a missing file.
step "inputdata E3SM does not ask for"
inputdata=/projects/e3sm/data/inputdata
# shellcheck disable=SC2043 # a list of one, for now
for f in atm/cam/physprops/p3_lookup_table_1.dat-v4.1.2; do
    [ -f "$inputdata/$f" ] && continue
    mkdir -p "$(dirname "$inputdata/$f")"
    curl -fsSL -o "$inputdata/$f.part" "https://web.lcrc.anl.gov/public/e3sm/inputdata/$f"
    mv "$inputdata/$f.part" "$inputdata/$f"
    echo "fetched $f"
done

mkdir -p "$scratch"
step "create_test ${tests[*]} $CT_FLAGS"
# --test-id keeps case directory names predictable for the logs CI uploads.
# --proc-pool: create_test will not start a run with more ranks than its pool (cores plus
# 25% by default). The machine oversubscribes, so let 8-rank runs start on a 4-core runner.
cores=$(nproc)
rc=0
# shellcheck disable=SC2086 # CT_FLAGS is a list of options
./create_test "${tests[@]}" $CT_FLAGS --wait --test-id gh \
    --parallel-jobs "$cores" --proc-pool "$(( cores > 8 ? cores : 8 ))" || rc=$?
shopt -s nullglob
for status in "$scratch"/*.gh/TestStatus; do
    echo "-- $status"
    cat "$status"
done
shopt -u nullglob
[ "$rc" -eq 0 ] || { echo "create_test failed ($rc)" >&2; exit "$rc"; }
echo "tests passed"
