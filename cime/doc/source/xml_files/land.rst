.. _land:

#################################
CIME Land Data and Stub XML Files
#################################

Land component XML files for data, stub, and dead components. 

.. toctree::
   :maxdepth: 1

***************************************************
CIMEROOT/src/components/data_comps/dlnd/cime_config
***************************************************

Land data model, **dlnd**, XML files and settings.

XML specification for archiving dlnd output files.

.. literalinclude:: ../../../src/components/data_comps/dlnd/cime_config/config_archive.xml

XML variables and component descriptions specific to dlnd. 

.. literalinclude:: ../../../src/components/data_comps/dlnd/cime_config/config_component.xml

XML namelist definitions for dlnd.

.. literalinclude:: ../../../src/components/data_comps/dlnd/cime_config/namelist_definition_dlnd.xml


***************************************************
CIMEROOT/src/components/stub_comps/slnd/cime_config
***************************************************

The land stub model, **slnd**, does not output any files in the RUNDIR nor
does it have any namelist settings.

XML variables and component descriptions specific to slnd. 

.. literalinclude:: ../../../src/components/stub_comps/slnd/cime_config/config_component.xml


***************************************************
CIMEROOT/src/components/xcpl_comps/xlnd/cime_config
***************************************************

The land dead model, **xlnd**, does not output any files in the RUNDIR nor
does it have any namelist settings.

XML variables and component descriptions specific to xlnd. 

.. literalinclude:: ../../../src/components/xcpl_comps/xlnd/cime_config/config_component.xml





