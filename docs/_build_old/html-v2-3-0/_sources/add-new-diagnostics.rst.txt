How to Add New Diagnostics Runs
-------------------------------


Adding Derived Variables
~~~~~~~~~~~~~~~~~~~~~~~~

We have a set of built-in derived variables for the E3SM model
diagnostics
`here <https://github.com/E3SM-Project/e3sm_diags/blob/master/acme_diags/derivations/acme.py>`__
(search for ``derived_variables``). The diagnostics software looks into
the ``derived_variables`` dictionary for variable keys and operations
needed for deriving new variables (renaming, unit conversions,
calculations, etc).

If users want to, they can add their own derived variables, which is
added to the default list during runtime and overwrites any default
values if there's a collision. Since derived variables require code,
such functionality cannot be added to json/cfg files. You can do the
following in the parameters script, which is a Python script (ex: the
Python script is ``run_e3sm_diags.py`` in ``python run_e3sm_diags.py -d mydiags.cfg`` or
``myparams.py`` in ``e3sm_diags -p myparams.py -d mydiags.cfg``).

Format of the ``derived_variables`` dictionary
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: python

    derived_variables = {
        'user_inputted_var': {
            ('output_var1', 'output_var2'): function_to_call
        }
    }

Above is how a ``derived_variables`` dictionary is formatted.
``'user_inputted_var'`` is the variable that is defined in the
``variables`` part of the json file. ``'output_var1'`` and
``'output_var2'`` are the variables inside the ``test_name`` (model)
file.

**Example of adding derived variables to a parameters script**

.. code:: python

    # in run_e3sm_diags.py or myparams.py

    def albedo_obs(rsdt, rsut):
        """ TOA (top-of-atmosphere) albedo, (solin - fsntoa) / solin, unit is nondimension """
        var = rsut / rsdt
        var.units = "dimensionless"
        var.long_name = "TOA albedo"
        return var

    derived_variables = {
        'New_ALBEDO': {
            ('rsdt', 'rsut'): albedo_obs
        }
    }

The above code will allow it so that if ``variables = "New_ALBEDO"`` in the
diagnostics file (the json or cfg file) and ``rsdt`` and ``rsut`` are
variables in the test (model) file, the ``albedo_obs()`` function is run
on the ``rsdt`` and ``rsut`` variables from the test (model) file.

**Note**: Please do check first if our built-in derived variable list already has included what you would like to
calculate. Also, for more advanced users, you can clone our repo and make direct changes to the source code and maybe
create a pull request if you think the derived variable is used frequently. We really appreciate it!

Example
~~~~~~~

The example below will do one diagnostics run globally with the
``New_ALBEDO`` variable, annually. Below is the json file, call it
``mydiags.cfg``.

::

    [Diags]
    sets = ['lat_lon']
    case_id = "lat_lon_CERES"
    variables = ["New_ALBEDO"]
    ref_name = "edition_4_ceres_ebaf_toa"
    reference_name = "edition_4_ceres_ebaf_toa"
    seasons = ["ANN"]
    regions = ["global"]
    contour_levels = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95]
    diff_levels = [-0.25, -0.2, -0.15, -0.1, -0.07, -0.05, -0.03, 0.0, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.25]

And below is the parameters file, call it ``run_e3sm_diags.py``. This is
to run on aims4. To run on another machine, please edit the
``reference_data_path``, ``test_data_path``, and ``test_name``
accordingly.

.. code:: python

    import os
    from acme_diags.parameter.core_parameter import CoreParameter
    from acme_diags.run import runner

    param = CoreParameter()

    param.reference_data_path = '/space1/test_data/CERES-EBAF/'
    param.test_data_path = '/space/golaz1/ACME_simulations/20160520.A_WCYCL1850.ne30_oEC.edison.alpha6_01/pp/clim_rgr/0070-0099/'

    param.test_name = '20160520.A_WCYCL1850.ne30'

    param.backend = 'vcs'
    param.diff_title = 'Test - Reference'
    param.results_dir = 'myresults'

    def albedo_obs(rsdt, rsut):
        """TOA (top-of-atmosphere) albedo, (solin - fsntoa) / solin, unit is nondimension"""
        var = rsut / rsdt
        var.units = "dimensionless"
        var.long_name = "TOA albedo"
        return var

    param.derived_variables = {
        'New_ALBEDO': {
            ('rsdt', 'rsut'): albedo_obs
        }
    }

    runner.sets_to_run = ['lat_lon']
    runner.run_diags([param])

Run the command like so:
``python run_e3sm_diags.py -d mydiags.cfg``
