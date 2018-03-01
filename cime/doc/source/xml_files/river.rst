.. _river:

#########################################
CIME River Runoff Data and Stub XML Files
#########################################

River runoff component XML files for data, stub, and dead components. 

.. toctree::
   :maxdepth: 1

***************************************************
CIMEROOT/src/components/data_comps/drof/cime_config
***************************************************

River runoff data model, **drof**, XML files and settings.

XML specification for archiving drof output files.

.. literalinclude:: ../../../src/components/data_comps/drof/cime_config/config_archive.xml

XML variables and component descriptions specific to drof. 

.. literalinclude:: ../../../src/components/data_comps/drof/cime_config/config_component.xml

XML variables and component descriptions specific to drof. 

.. literalinclude:: ../../../src/components/data_comps/drof/cime_config/config_component.xml

XML namelist definitions for drof.

.. literalinclude:: ../../../src/components/data_comps/drof/cime_config/namelist_definition_drof.xml


***************************************************
CIMEROOT/src/components/stub_comps/srof/cime_config
***************************************************

The river runoff stub model, **srof**, does not output any files in the RUNDIR nor
does it have any namelist settings.

XML variables and component descriptions specific to srof. 

.. literalinclude:: ../../../src/components/stub_comps/srof/cime_config/config_component.xml


***************************************************
CIMEROOT/src/components/xcpl_comps/xrof/cime_config
***************************************************

The river runoff dead model, **xrof**, does not output any files in the RUNDIR nor
does it have any namelist settings.

XML variables and component descriptions specific to xrof.

.. literalinclude:: ../../../src/components/xcpl_comps/xrof/cime_config/config_component.xml





