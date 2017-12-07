.. _ocean:

##################################
CIME Ocean Data and Stub XML Files
##################################

Ocean component XML files for data, stub, and dead components. 

.. toctree::
   :maxdepth: 1

***************************************************
CIMEROOT/src/components/data_comps/docn/cime_config
***************************************************

Ocean data model, **docn**, XML files and settings.

XML specification for archiving docn output files.

.. literalinclude:: ../../../src/components/data_comps/docn/cime_config/config_archive.xml

XML variables and component descriptions specific to docn. 

.. literalinclude:: ../../../src/components/data_comps/docn/cime_config/config_component.xml

XML namelist definitions for docn.

.. literalinclude:: ../../../src/components/data_comps/docn/cime_config/namelist_definition_docn.xml


***************************************************
CIMEROOT/src/components/stub_comps/socn/cime_config
***************************************************

The ocean stub model, **socn**, does not output any files in the RUNDIR nor
does it have any namelist settings.

XML variables and component descriptions specific to socn. 

.. literalinclude:: ../../../src/components/stub_comps/socn/cime_config/config_component.xml


***************************************************
CIMEROOT/src/components/xcpl_comps/xocn/cime_config
***************************************************

The ocean dead model, **xocn**, does not output any files in the RUNDIR nor
does it have any namelist settings.

XML variables and component descriptions specific to xocn. 

.. literalinclude:: ../../../src/components/xcpl_comps/xocn/cime_config/config_component.xml





