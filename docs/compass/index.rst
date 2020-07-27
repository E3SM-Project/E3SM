COMPASS
=======

Overview
--------

The COMPASS (Configuration Of Model for Prediction Across Scales Setups)
infrastructure provides a capability for defining simple test-case workflows.
It is intended to house a small number of files which can describe the steps to
setup and configure a test case.

It provides four utility python scripts:

* ``clean_testcase.py``
* ``list_testcases.py``
* ``setup_testcase.py``
* ``manage_regression_suite.py``

and two configuration file templates:

  * ``general.config.test``
  * ``general.config.ocean``

Each of the python scripts can be run with a ``-h`` argument only to get usage
information.

Additionally, each core has a directory at the top level (e.g. ocean for the
ocean test cases). There is also a templates directory where a core can place
template files that are intended to be available for it's test cases.

An example test case is placed in ``ocean/baroclinic_channel/10km``
An example template is placed in ``templates/ocean/global_stats.xml``

Test cases are described by XML files. Each test case can have an arbitrary
number of XML files that configure the steps for setting up the test case.

The various XML files that can be used with this test case infrastructure are
described in the README files contained in the doc directory.

.. toctree::
   :titlesonly:

   details
   scripts
   ocean
   ocean_testcases/index