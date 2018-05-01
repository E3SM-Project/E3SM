
Modifying driver namelists
-------------------------------------------

Driver namelist variables belong in two groups:

1. Those that are set directly from ``$CASEROOT`` xml variables.

2. Those that are set by the driver utility **$CIMEROOT/src/drivers/mct/cime_config/buildnml**.

All driver namelist variables are defined in the file **$CIMEROOT/src/drivers/mct/cime_config/namelist_definition_drv.xml**.
The variables that can be changed only by modifying xml variables appear with the *entry* attribute ``modify_via_xml="xml_variable_name"``.

All other variables that appear in the **namelist_definition_drv.xml** file can be modified by adding a keyword value pair at the end of ``user_nl_cpl``.
For example, to change the driver namelist value of ``eps_frac`` to ``1.0e-15``, add the following line to the end of the ``user_nl_cpl``:
::

   eps_frac = 1.0e-15

To see the result of change, call **preview_namelists** and verify that the new value appears in **CaseDocs/drv_in**.
