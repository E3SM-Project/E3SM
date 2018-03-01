.. _esp:

################################
CIME ESP Data and Stub XML Files
################################

External System Processing **ESP** component XML files for data, stub, and dead components. 

.. toctree::
   :maxdepth: 1

***************************************************
CIMEROOT/src/components/data_comps/desp/cime_config
***************************************************

ESP data model, **desp**, XML files and settings.

XML variables and component descriptions specific to desp. 

.. literalinclude:: ../../../src/components/data_comps/desp/cime_config/config_component.xml

XML namelist definitions for desp.

.. literalinclude:: ../../../src/components/data_comps/desp/cime_config/namelist_definition_desp.xml


***************************************************
CIMEROOT/src/components/stub_comps/sesp/cime_config
***************************************************

The ESP stub model, **sesp**, does not output any files in the RUNDIR nor
does it have any namelist settings.

XML variables and component descriptions specific to sesp. 

.. literalinclude:: ../../../src/components/stub_comps/sesp/cime_config/config_component.xml






