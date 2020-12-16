# Run from top level directory of repository.

# Uncomment if you want to remove the results directory from the previous test run.
#rm -rf tests/system/all_sets_results_test/

# Copy over test data and images [takes about two minutes]
# TODO: CI tool will access data from the web
# Change <username>
scp -r <username>@blues.lcrc.anl.gov:/lcrc/group/e3sm/public_html/e3sm_diags_test_data/unit_test_data tests/unit_test_data
scp -r <username>@blues.lcrc.anl.gov:/lcrc/group/e3sm/public_html/e3sm_diags_test_data/unit_test_images tests/unit_test_images

# Use "&&" so the test will fail if either python call fails.
cd tests/system && python -m unittest && cd .. && python -m unittest && cd ..

# Uncomment if you want to delete the test data/images copied from Blues after running the tests.
#rm -rf tests/unit_test_data
#rm -rf tests/unit_test_images
