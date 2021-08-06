# Run this script from the top level directory of the repository.

set -e # Fail if any line fails

# 1. Run unit tests
python -m unittest tests/e3sm_diags/*/test_*.py

# 2. Enter tests directory
cd tests

# 3. Copy over test data and images from the web. [Takes about four minutes]
cd integration && python download_data.py && cd ..

# To use `scp`, run the following two lines instead. Be sure to change <username>. [Takes about two minutes]
#scp -r <username>@blues.lcrc.anl.gov:/lcrc/group/e3sm/public_html/e3sm_diags_test_data/integration/integration_test_data integration_test_data
#scp -r <username>@blues.lcrc.anl.gov:/lcrc/group/e3sm/public_html/e3sm_diags_test_data/integration/expected/integration_test_images integration_test_images

# 4. Run integration tests
cd integration && python -m unittest && cd ..
