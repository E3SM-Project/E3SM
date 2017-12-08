.. _seaice:

####################################
CIME Sea Ice Data and Stub XML Files
####################################

Sea Ice component XML files for data, stub, and dead components. 

.. toctree::
   :maxdepth: 1

***************************************************
CIMEROOT/src/components/data_comps/dice/cime_config
***************************************************

Sea Ice data model, **dice**, XML files and settings.

XML specification for archiving dice output files.

.. literalinclude:: ../../../src/components/data_comps/dice/cime_config/config_archive.xml

XML variables and component descriptions specific to dice. 

.. literalinclude:: ../../../src/components/data_comps/dice/cime_config/config_component.xml

XML variables and component descriptions specific to dice. 

.. literalinclude:: ../../../src/components/data_comps/dice/cime_config/config_component.xml

XML namelist definitions for dice.

.. literalinclude:: ../../../src/components/data_comps/dice/cime_config/namelist_definition_dice.xml


***************************************************
CIMEROOT/src/components/stub_comps/sice/cime_config
***************************************************

The sea ice stub model, **sice**, does not output any files in the RUNDIR nor
does it have any namelist settings.

XML variables and component descriptions specific to sice. 

.. literalinclude:: ../../../src/components/stub_comps/sice/cime_config/config_component.xml


***************************************************
CIMEROOT/src/components/xcpl_comps/xice/cime_config
***************************************************

The sea ice dead model, **xice**, does not output any files in the RUNDIR nor
does it have any namelist settings.

XML variables and component descriptions specific to satm. 

.. literalinclude:: ../../../src/components/xcpl_comps/xice/cime_config/config_component.xml





