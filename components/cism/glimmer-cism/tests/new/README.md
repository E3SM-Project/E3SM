How to create a new CISM test
=============================

This directory contains a new test template called `runTest.py` which
should be used to create a new test in CISM. This template is based off of the
`tests/higher-order/dome` test, and you should start by playing around with that test
so that you understand how the options work, and how the tests are run. If you
maintain the structure of this test template, your new test will work with the
`tests/regression/build_and_test.py` structure (BATS). 

Creating a new test
----------------------------

Say we want to create a higher-order test for Antarctica, which we will 
call `antarctica` The first thing we do is make the Antarctica test
directory, and then move into it:

```bash
$ mkdir $CISM/tests/higher-order/antarctica
$ cd $CISM/tests/higher-order/antarctica
```

Next, we copy the `runTest.py` python script template into your new test
directory, the `*.config` files, and the netCDF tools module `netCDF.py` from
the `$CISM/tests/new` directory. 

```bash
$ cp $CISM/tests/new/runTest.py runAntarctica.py
$ cp $CISM/tests/new/test.config antarctica.config
$ cp $CISM/tests/new/test-forcing.config antarctica-forcing.config
$ cp $CISM/test/new/netCDF.py ./
```

Next, we open up the `runAntarctica.py` script in our favorite editor 
and create our test. Most of the structure of the test has been created
from the template. To finish the test, we find the `#FIXME:` statements
in the `runAntarctica.py`, read the associated comment text, and create
what we need for the test case. 


_Last updated: August 13, 2015 by Joseph H Kennedy at ORNL_
