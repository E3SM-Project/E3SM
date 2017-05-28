Mass and Heat Budgets
=====================

Mass and heat are conserved in the coupler to several digits over centuries.
Several steps have been taken to ensure this level of conservation, and these are described in other sections of the document. 
In addition, efforts have been made to make sure each component is internally conservative with respect to mass and heat.

The budgets can be turned on and off using the namelist variable ``do_budgets``.
The value of that namelist is set by the ``$CASEROOT/env_run.xml`` variable, ``BUDGETS``.

The driver coupler can diagnose heat and mast budgets at several levels and over different periods.
The periods are *instantenous*, *daily average*, *monthly average*, *annual average*, or since the start of the run. 
The budget output for each of these periods is controlled by the driver namelist variables ``budget_inst``, ``budget_daily``, ``budget_month``, ``budget_ann``, ``budget_ltann``, and ``budget_ltend``. 
``budget_ltann`` and ``budget_ltend`` are used to write the long term budget at either the end of every year or the end of every run.
Other budgets are written at their period interval. 

The namelist input is an integer specifying what to write.
The budget flags are controlled by ``$CASEROOT/env_run.xml`` variables ``BUDGET_INST``, ``BUDGET_DAILY``, ``BUDGET_MONTHLY``, ``BUDGET_ANNUAL``, ``BUDGET_LONGTER_EOY``, and ``BUDGET_LONGTERM_STOP`` respectively. 
Valid values are 0, 1, 2, or 3.
If 0 is set, no budget data is written. 
The value 1 generates a net heat and water budget for each component, 2 adds a detailed heat and water budget for each component, and 3 adds a detailed heat and water budget of the different conmponents on the atmosphere grid.
Normally values of 0 or 1 are specified. 
Values of 2 or 3 are generally used only when debugging problems involving conservation.
