# Run this script from the top level directory of the repository.

set -e # Fail if any line fails

printf "1. Run unit tests\n"
printf "==============================================\n"
python -m unittest tests/e3sm_diags/*/test_*.py

printf "\n2. Copy over test data and images from the web\n"
printf "==============================================\n"
# Takes about four minutes
python -m tests.integration.download_data

# To use `scp`, run the following two lines instead. Be sure to change <username>. [Takes about two minutes]
#scp -r <username>@blues.lcrc.anl.gov:/lcrc/group/e3sm/public_html/e3sm_diags_test_data/integration/integration_test_data integration_test_data
#scp -r <username>@blues.lcrc.anl.gov:/lcrc/group/e3sm/public_html/e3sm_diags_test_data/integration/expected/integration_test_images integration_test_images

printf "\n3. Run integration tests\n"
printf "==============================================\n"
# Integration tests must be called individually. `e3sm_diags` uses `cdp_parser.get_parameters()`, which incorrectly picks up `python -m unittest` cmd args.
# This results in `cdp_parser.get_parameters()` raising `error: unrecognized arguments <UNITTEST ARGS>`.
# Running tests using these commands do not work:
#  - `python -m unittest tests/integration/test_*.py`
#  - `python -m unittest discover -s tests/integration --pattern "test_*.py"`
# Related lines of code : https://github.com/CDAT/cdp/blob/2917603097112f7db52be0bf7a2e766e14cc2e16/cdp/cdp_parser.py#L498
command='python -m unittest tests/integration/'
declare -a files=("test_all_sets" "test_dataset" "test_diags" "test_run")
for file in "${files[@]}"
    do
        file="${file}.py"
        printf "\n${file}\n"
        printf "~~~~~~~~~~~~~~~~~~~~\n"
        eval "${command}${file}"
done
