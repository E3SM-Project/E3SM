.. _drivers:

#################
Driver XML Files
#################

Driver/Coupler XML files in CIMEROOT/src/drivers/mct/cime_config

.. toctree::
   :maxdepth: 1

************************************
CIMEROOT/src/drivers/mct/cime_config
************************************

The Model Coupling Toolkit (MCT) based driver/coupler is
treated as a component by CIME with associated XML files
to define behavior. 

XML specification for archiving coupler output files.

.. literalinclude:: ../../../src/drivers/mct/cime_config/config_archive.xml

XML variables and component descriptions specific to the driver/coupler

.. literalinclude:: ../../../src/drivers/mct/cime_config/config_component.xml

XML variables and component descriptions specific to the E3SM driver/coupler

.. literalinclude:: ../../../src/drivers/mct/cime_config/config_component_e3sm.xml

XML variables and component descriptions specific to the E3SM driver/coupler

.. literalinclude:: ../../../src/drivers/mct/cime_config/config_component_cesm.xml

XML settings for driver/coupler defined component set (compset) configurations.

.. literalinclude:: ../../../src/drivers/mct/cime_config/config_compsets.xml

XML settings for driver/coupler defined component set (compset) PE layouts.

.. literalinclude:: ../../../src/drivers/mct/cime_config/config_pes.xml


**************************************************************
CIMEROOT/src/drivers/mct/cime_config/namelist_definition_*.xml 
**************************************************************

XML namelist definitions for the driver/coupler.

.. literalinclude:: ../../../src/drivers/mct/cime_config/namelist_definition_drv.xml

XML namelist definitions for the driver/coupler fields.

.. literalinclude:: ../../../src/drivers/mct/cime_config/namelist_definition_drv_flds.xml

XML namelist definitions for the driver/coupler model input/output settings.

.. literalinclude:: ../../../src/drivers/mct/cime_config/namelist_definition_modelio.xml



