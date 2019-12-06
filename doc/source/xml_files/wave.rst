.. _wave:

#################################
CIME Wave Data and Stub XML Files
#################################

Wave component XML files for data, stub, and dead components. 

.. toctree::
   :maxdepth: 1

***************************************************
CIMEROOT/src/components/data_comps/dwav/cime_config
***************************************************

Wave data model, **dwav**, XML files and settings.

XML specification for archiving dwav output files.

.. literalinclude:: ../../../src/components/data_comps/dwav/cime_config/config_archive.xml

XML variables and component descriptions specific to dwav. 

.. literalinclude:: ../../../src/components/data_comps/dwav/cime_config/config_component.xml

XML variables and component descriptions specific to dwav. 

.. literalinclude:: ../../../src/components/data_comps/dwav/cime_config/config_component.xml

XML namelist definitions for dwav.

.. literalinclude:: ../../../src/components/data_comps/dwav/cime_config/namelist_definition_dwav.xml


***************************************************
CIMEROOT/src/components/stub_comps/swav/cime_config
***************************************************

The wave stub model, **swav**, does not output any files in the RUNDIR nor
does it have any namelist settings.

XML variables and component descriptions specific to swav. 

.. literalinclude:: ../../../src/components/stub_comps/swav/cime_config/config_component.xml


***************************************************
CIMEROOT/src/components/xcpl_comps/xwav/cime_config
***************************************************

The wave dead model, **xwav**, does not output any files in the RUNDIR nor
does it have any namelist settings.

XML variables and component descriptions specific to xwav. 

.. literalinclude:: ../../../src/components/xcpl_comps/xwav/cime_config/config_component.xml





