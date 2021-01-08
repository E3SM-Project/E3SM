# Run from top level directory of repository.

# Uncomment if you want to remove the results directory from the previous test run.
#rm -rf tests/system/all_sets_results_test/

cd tests
# Copy over test data and images from the web. [Takes about four minutes]
python download_data.py
# To use `scp`, run the following two lines instead. Be sure to change <username>. [Takes about two minutes]
#scp -r <username>@blues.lcrc.anl.gov:/lcrc/group/e3sm/public_html/e3sm_diags_test_data/unit_test_data unit_test_data
#scp -r <username>@blues.lcrc.anl.gov:/lcrc/group/e3sm/public_html/e3sm_diags_test_data/unit_test_images unit_test_images

# Use "&&" so the test will fail if either python call fails.
cd system && python -m unittest && cd .. && python -m unittest && cd ..

# Uncomment if you want to delete the test data/images copied from Blues after running the tests.
#rm -rf tests/unit_test_data
#rm -rf tests/unit_test_images
