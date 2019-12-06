.. _cime-dir:

******************
Directory content
******************

If you use CIME as part of a climate model or standalone, the content of the **cime** directory is the same.

If you are using it as part of a climate model, **cime** is usually one of the first subdirectories under the main directory.

.. table:: **CIME directory in a climate model**

   ====================== ===================================
   Directory or Filename               Description
   ====================== ===================================
   README, etc.           typical top-level directory content
   components/            source code for active models
   cime/                  All of CIME code
   ====================== ===================================

CIME's content is split into several subdirectories. Users should start in the **scripts/** subdirectory.

.. table::  **CIME directory content**

   ========================== ==================================================================
   Directory or Filename               Description
   ========================== ==================================================================
   CMakeLists.txt	      For building with CMake
   ChangeLog		      Developer-maintained record of changes to CIME
   ChangeLog_template	      Template for an entry in ChangeLog
   LICENSE.TXT		      The CIME license
   README		      Brief intro to CIME
   README.md		      README in markdown language
   README.unit_testing	      Instructions for running unit tests with CIME
   **config/**		      **Shared and model-specific configuration files**
   config/cesm/	              CESM-specific configuration options
   config/e3sm/	              E3SM-specific configuration options
   **scripts/**		      **The CIME user interface**
   scripts/lib/  	      Infrastructure source code for CIME scripts and functions
   scripts/Tools/	      Auxiliary tools; scripts and functions
   **src/**		      **Model source code provided by CIME**
   src/components/	      CIME-provided components including data and stub models
   src/drivers/  	      CIME-provided main driver for a climate model
   src/externals/	      Software provided with CIME for building a climate model
   src/share/    	      Model source code provided by CIME and used by multiple components
   **tests/**		      **Tests**
   **tools/**		      **Standalone climate modeling tools**
   utils/		      Some Perl source code needed by some prognostic components
   ========================== ==================================================================
