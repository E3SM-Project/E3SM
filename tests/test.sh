# Run this script from the top level directory of the repository.

set -e # Fail if any line fails

printf "1. Run unit tests\n"
printf "==============================================\n"
pytest tests/e3sm_diags/

printf "\n2. Copy over test data and images from the web\n"
printf "==============================================\n"
# Takes about four minutes
python -m tests.integration.download_data

# To use `scp`, run the following two lines instead. Be sure to change <username>. [Takes about two minutes]
#scp -r <username>@blues.lcrc.anl.gov:/lcrc/group/e3sm/public_html/e3sm_diags_test_data/integration/integration_test_data integration_test_data
#scp -r <username>@blues.lcrc.anl.gov:/lcrc/group/e3sm/public_html/e3sm_diags_test_data/integration/expected/integration_test_images integration_test_images

printf "\n3. Run integration tests\n"
printf "==============================================\n"
pytest tests/integration/
