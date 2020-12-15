# Run from top level directory of repository.

# Copy over test data [takes about two minutes]
# TODO: How will an automated CI tool access the E3SM data server? We could set up a public key for the CI server.
# Change <username>
scp -r <username>@blues.lcrc.anl.gov:/lcrc/group/e3sm/public_html/e3sm_diags_test_data/unit_test_data tests/unit_test_data

# Use "&&" so the test will fail if either python call fails.
cd tests/system && python -m unittest && cd .. && python -m unittest && cd ..

# Uncomment if you want to delete the test data copied from blues after running the tests.
#rm -rf tests/unit_test_data
