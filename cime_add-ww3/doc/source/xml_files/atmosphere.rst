.. _atmosphere:

#######################################
CIME Atmosphere Data and Stub XML Files
#######################################

Atmosphere component XML files for data, stub, and dead components. 

.. toctree::
   :maxdepth: 1

***************************************************
CIMEROOT/src/components/data_comps/datm/cime_config
***************************************************

Atmosphere data model, **datm**, XML files and settings.

XML specification for archiving datm output files.

.. literalinclude:: ../../../src/components/data_comps/datm/cime_config/config_archive.xml

XML variables and component descriptions specific to datm. 

.. literalinclude:: ../../../src/components/data_comps/datm/cime_config/config_component.xml

XML namelist definitions for datm.

.. literalinclude:: ../../../src/components/data_comps/datm/cime_config/namelist_definition_datm.xml


***************************************************
CIMEROOT/src/components/stub_comps/satm/cime_config
***************************************************

The atmosphere stub model, **satm**, does not output any files in the RUNDIR nor
does it have any namelist settings.

XML variables and component descriptions specific to satm. 

.. literalinclude:: ../../../src/components/stub_comps/satm/cime_config/config_component.xml


***************************************************
CIMEROOT/src/components/xcpl_comps/xatm/cime_config
***************************************************

The atmosphere dead model, **xatm**, does not output any files in the RUNDIR nor
does it have any namelist settings.

XML variables and component descriptions specific to xatm. 

.. literalinclude:: ../../../src/components/xcpl_comps/xatm/cime_config/config_component.xml





