# Running the testing framework

Before running anything, be sure to `cd` to the right directory and download the test data:

* `cd E3SM/components/cam/src/physics/crm/samxx/test/build`
* `./download_data.sh`

The test system configures itself from the test data, so the CMake configure script must be given two parameters: (1) the 2-D file and (2) the 3-D file

First, source the appropriate environment file, e.g.:

* `source summit_gpu.sh`  OR
* `source summit_cpu.sh`

You may wish to change the number of CRMs being run in those files before sourcing them.

Next, run the CMake configure script

* `./cmakescript.sh crmdata_nx32_ny1_nz28_nxrad2_nyrad1.nc crmdata_nx8_ny8_nz28_nxrad2_nyrad2.nc`

Finally, run the tests:

* `./runtest.sh`

## Larger tests

If you want to run with more CRM columns than the existing test files have, you can replicate those columns `N` times with the `ncreplicate.py` script, e.g.:

* `python ncreplicate.py crmdata_nx32_ny1_nz28_nxrad2_nyrad1.nc 80`
* `python ncreplicate.py crmdata_nx8_ny8_nz28_nxrad2_nyrad2.nc 80`
* `./cmakescript.sh crmdata_nx32_ny1_nz28_nxrad2_nyrad1_80x.nc crmdata_nx8_ny8_nz28_nxrad2_nyrad2_80x.nc`
* `./runtests.sh`


