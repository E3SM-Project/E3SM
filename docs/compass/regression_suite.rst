.. _compass_regression_suite:

regression\_suite
=================

A ``regression_suite`` file is used to define a regression suite, which
involves a set of tests that should be run. This file contains information
describing a set of tests that should be setup and run as part of a regression
test suite.

Below, you will see text describing the various XML tags available in a
``regression_suite`` file. Each will describe the tag itself, any attributes the
tag can have, and what children can be placed below the tag.

``<regression_suite>`` - This is the overarching parent tag in a ``regression_suite``
file. It defines the suite that will be setup

    - Attributes:
        * ``name``: This attribute defines the name of the regression suite. A
          script will be generated named ``<name>.py`` that will run the entire
          regression suite in the location the regression suite is setup (i.e.
          ``work_dir``)

    - Children:
        * ``<test>``

``<test>`` - This tag defines a test that will be included as part of the
regression suite.

    - Attributes:
        * ``name``: This attribute defines the name of the test as part of this
          regression suite. NOTE: This name is only used in this
          regression suite, so multiple suites can name the same test in
          different ways.

        * ``core``: This attribute defines the core that would be passed to
          ``setup_testcases.py`` to setup this test.

        * configuration: This attribute defines the configuration that would be passed to
          ``setup_testcases.py`` to setup this test.

        * resolution: This attribute defines the resolution that would be passed to
          ``setup_testcases.py`` to setup this test.

    - Children:
        * ``<script>``

``<script>`` - This tag defines a script that will be run to perform the
specified test.

    - Attributes:
        * ``name``: This attribute defines the name of the script that will be run
          to perform the specified test. Typically this is a driver script.
