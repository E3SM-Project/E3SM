.. _introduction-and-overview:


*************
Introduction
*************

Part 1 of this guide explains the basic commands
that are needed to get a model running.

Prerequisites
=============

Part 1 of this guide assumes that CIME and the necessary input files have been installed on
the computer you are using. If that is not the case, see Installing CIME.

Other prerequisites:

- Familiarity with basic climate modeling concepts.

- Familiarity with UNIX command line terminals and the UNIX development environment.

- A correct version of the Python interpreter.

CIME's commands are Python scripts and require a correct version of
the Python interpreter to be installed. The Python version must be
greater than 2.7 but less than 3.0. Determine which version you have
like this:
::

   > python --version


Key Terms and concepts
======================

The following key terms and concepts are ingrained in CIME and are used frequently in this documentation.
See the glossary for a more complete list of terms.

**active, data and stub models**
   *active models*: Components of a model that solve a complex set of equations to describe the model's behavior are called
   *active* models. Sometimes they are called *prognostic* or *full* models.

   CIME recognizes 7 different active models. They are:

       atmosphere, ocean, sea-ice, land surface, river, glacier, wave

   An external system processing (ESP) stub-only component is also allowed.

   *data models*: For some climate problems, it is necessary to reduce feedbacks within the system by replacing an active model with a
   version that sends and receives the same variables to and from other models, but with the values read from files rather
   than computed from the equations. The values received are ignored. These active-model substitutes are called *data models*.
   CIME provides data models for each of the supported components.

   *stub models*: For some configurations, no data model is needed, so CIME provides *stub* versions that simply occupy the
   required place in the climate execution sequence and do not send or receive any data.

**case**:
    The most important concept in CIME is a *case*. To build and execute a CIME-enabled climate model, you have to
    make choices of compset, model grid, machine and compiler. A collection of these choices, and any additional
    customizations you may make, is called the *case*.

**compiler**:
   CIME controls compiling of your model's source code (Fortran, C and C++) into an executable.
   Some machines support multiple compilers, so you may need to specify which one to use.

**component set** or **compset**:
   CIME allows several sub-models and other tools to be linked together into a climate model. These sub-models and
   tools are called *components* of the climate model. For example, a climate model has an atmosphere component, an
   ocean component, and so on. The resulting set of components is called the *component set* or *compset*.

**grid** or **model grid**:
   Each active model must solve its equations on a numerical grid. CIME allows models within the system to have
   different grids. The resulting set of numerical grids is called the *model grid* or sometimes just the *grid*, where
   *grid* is a unique name that denotes a set of numerical grids. Sometimes the *resolution* also refers to a specific set
   of grids with different resolutions.

**machine**:
   The *machine* is the computer you are using to run CIME and build and run the climate model. It could be a workstation
   or a national supercomputer. The *machine* is typically the UNIX hostname but it could be any string.

**out-of-the-box**:
   Any case or capability available with a basic installation of a CIME-driven model.

CIME and your environment
=========================

Before using any CIME commands, set the ``CIME_MODEL`` environment variable. In bash, use **export** as shown and replace
**<model>** with the appropriate text. Current possibilities are "e3sm" or "cesm."
::

   > export CIME_MODEL=<model>

There are a number of possible ways to set CIME variables.
For variables that can be set in more than one way, the order of precedence is:

- variable appears in a command line argument to a CIME command

- variable is set as an environment variable

- variable is set in ``$HOME/.cime/config``

- variable is set in a ``$CASEROOT`` xml file

Directory content
==================

If you use CIME as part of a climate model or standalone, the content of the **cime** directory is the same.

If you are using it as part of a climate model, **cime** is usually one of the first subdirectories under the main directory:

.. csv-table:: **CIME directory in a climate model**
   :header: "Directory or Filename", "Description"
   :widths: 200, 300

   "README, etc.", "typical top-level directory content"
   "components/", "source code for active models"
   "cime/", "All of CIME code"

CIME's content is split into several subdirectories. Users should start in the **scripts** subdirectory.

.. csv-table:: **CIME directory content**
   :header: "Directory or Filename", "Description"
   :widths: 150, 300

   "CMakeLists.txt", "For building with CMake"
   "ChangeLog", "Developer-maintained record of changes to CIME"
   "ChangeLog_template", "Template for an entry in ChangeLog"
   "LICENSE.TXT", "The CIME license"
   "README", "Brief intro to CIME"
   "README.md", "README in markdown language"
   "README.unit_testing", "Instructions for running unit tests with CIME"
   "config/", "Shared and model-specific configuration files"
   "scripts/", "The CIME user interface"
   "src/", "Model source code provided by CIME"
   "tests/", "tests"
   "tools/", "Standalone climate modeling tools"
   "utils/", "Some Perl source code for CIME scripts; see **scripts/lib** for Python version"

Here are some other key subdirectories, down one level in the
directory structure.

.. csv-table:: **Content of some key CIME subdirectories**
   :header: "Directory or Filename", "Description"
   :widths: 150, 300

   "config/cesm/", "CESM-specific configuration options"
   "config/e3sm/", "E3SM-specific configuration options"
   "src/components/", "CIME-provided components including data and stub models"
   "src/drivers/", "CIME-provided main driver for a climate model"
   "src/externals/", "Software provided with CIME for building a climate model"
   "src/share/", "Model source code provided by CIME and used by multiple components"
   "scripts/lib/", "Infrastructure source code for CIME scripts and functions"
   "scripts/Tools/", "Auxiliary tools; scripts and functions"

Discovering available cases with **query_config**
=================================================


Use the utility **$CIMEROOT/scripts/query_config** to see which out-of-the-box compsets, components, grids and machines are available for a model.

Optional arguments include the following::

  --compsets
  --components
  --grids
  --machines

If CIME is downloaded in standalone mode, only standalone CIME compsets can be queried. If CIME is part of CIME
-driven model, **query_config** will allow you to query all prognostic component compsets.
To see lists of available compsets, components, grids and machines, look at the **help** text::

  > query_config --help

**Usage examples**

To run **query_config** for compset information, follow this example, where **drv** is the component name::

  > query_config --compsets drv

The output will be similar to this::

     --------------------------------------
     Compset Short Name: Compset Long Name
     --------------------------------------
   A                    : 2000_DATM%NYF_SLND_DICE%SSMI_DOCN%DOM_DROF%NYF_SGLC_SWAV
   ADWAV                : 2000_SATM_SLND_SICE_SOCN_SROF_SGLC_DWAV%CLIMO
   S                    : 2000_SATM_SLND_SICE_SOCN_SROF_SGLC_SWAV_SESP
   ADLND                : 2000_SATM_DLND%SCPL_SICE_SOCN_SROF_SGLC_SWAV
   ADESP_TEST           : 2000_DATM%NYF_SLND_DICE%SSMI_DOCN%DOM_DROF%NYF_SGLC_SWAV_DESP%TEST
   X                    : 2000_XATM_XLND_XICE_XOCN_XROF_XGLC_XWAV
   ADESP                : 2000_DATM%NYF_SLND_DICE%SSMI_DOCN%DOM_DROF%NYF_SGLC_SWAV_DESP
   AIAF                 : 2000_DATM%IAF_SLND_DICE%IAF_DOCN%IAF_DROF%IAF_SGLC_SWAV

Each model component specifies its own definitions of what can appear after the ``%`` modifier in the compset longname (for example, ``DOM`` in ``DOCN%DOM``).

To see what supported modifiers are for ``DOCN``, run **query_config** as in this example::

  > query_config --component docn

The output will be similar to this::

     =========================================
     DOCN naming conventions
     =========================================

         _DOCN%AQP1 : docn prescribed aquaplanet sst - option 1
        _DOCN%AQP10 : docn prescribed aquaplanet sst - option 10
         _DOCN%AQP2 : docn prescribed aquaplanet sst - option 2
         _DOCN%AQP3 : docn prescribed aquaplanet sst - option 3
         _DOCN%AQP4 : docn prescribed aquaplanet sst - option 4
         _DOCN%AQP5 : docn prescribed aquaplanet sst - option 5
         _DOCN%AQP6 : docn prescribed aquaplanet sst - option 6
         _DOCN%AQP7 : docn prescribed aquaplanet sst - option 7
         _DOCN%AQP8 : docn prescribed aquaplanet sst - option 8
         _DOCN%AQP9 : docn prescribed aquaplanet sst - option 9
          _DOCN%DOM : docn prescribed ocean mode
          _DOCN%IAF : docn interannual mode
         _DOCN%NULL : docn null mode
          _DOCN%SOM : docn slab ocean mode
       _DOCN%SOMAQP : docn aquaplanet slab ocean mode
    _DOCN%SST_AQUAP : docn aquaplanet mode:


For more details on how CIME determines the output for **query_config**, see :ref:`cime-internals`.

Quick start
==================

To see an example of how a case is created, configured, built and run with CIME, execute the following commands for an example. (This assumes that CIME has been ported to your current machine).
::

   > cd cime/scripts
   > ./create_newcase --case mycase --compset X --res f19_g16
   > cd mycase
   > ./case.setup
   > ./case.build
   > ./case.submit

The output from each command is explained in the following sections.

After you submit the case, you can follow the progress of your run by monitoring the **CaseStatus** file.

::

   > tail CaseStatus

Repeat the command until you see the message ``case.run success``.
