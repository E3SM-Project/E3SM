.. _cesm:

#############################
CESM Coupled Model XML Files
#############################

XML files for CESM in CIMEROOT/config/cesm. 

.. toctree::
   :maxdepth: 1

********************
CIMEROOT/config/cesm
********************

CESM XML settings for short term archiver.

.. literalinclude:: ../../../config/cesm/config_archive.xml

CESM XML settings for defining CASEROOT env_*.xml file entries.

.. literalinclude:: ../../../config/cesm/config_files.xml

CESM XML settings for defining supported grids.

.. literalinclude:: ../../../config/cesm/config_grids.xml


*****************************
CIMEROOT/config/cesm/machines
*****************************

CESM XML settings for supported batch queuing systems.

.. literalinclude:: ../../../config/cesm/machines/config_batch.xml

CESM XML settings for supported compilers.

.. literalinclude:: ../../../config/cesm/machines/config_compilers.xml

CESM XML settings for supported machines.

.. literalinclude:: ../../../config/cesm/machines/config_machines.xml

CESM XML settings for Parallel Input/Output (PIO) library. 

.. literalinclude:: ../../../config/cesm/machines/config_pio.xml


******************************
allactive SRCROOT/cime_config
******************************

The CESM all-active model settings are stored in the `CESM cime_config github repository
<https://github.com/CESM-Development/cime_config>`_. That repository includes
the following XML files.

CESM XML settings for all-active component set (compset) configurations.

.. literalinclude:: ../../../../cime_config/config_compsets.xml

CESM XML settings for all-active test configurations.

.. literalinclude:: ../../../../cime_config/testlist_allactive.xml

CESM XML settings for optimized processor elements (PEs) layout configurations.

.. literalinclude:: ../../../../cime_config/config_pes.xml

