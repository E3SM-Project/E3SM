.. _introduction-and-overview:

.. role:: red


*************
Introduction
*************

Part 1 of this guide explains the basic commands in the CIME Case Control System
that are needed to get a model running.

Prerequisites
=============

Part 1 of this guide assumes that CIME or a CIME-driven model and the necessary input files
have been installed on the computer you are using. If that is not the case, see :ref:`Porting CIME<porting>`.

Other prerequisites:

- Familiarity with basic climate modeling concepts.

- Familiarity with UNIX command line terminals and the UNIX development environment.

- A correct version of the Python interpreter.

CIME's commands are Python scripts and require a correct version of
the Python interpreter to be installed. The Python version must be
greater than 2.7.  Determine which version you have
like this:
::

   > python --version

Consult your local documentation if you need to update your python version.

Key Terms and concepts
======================

The following key terms and concepts are ingrained in CIME and used frequently in this documentation.
See the :ref:`glossary` for a more complete list of terms.

**components**

   In CIME, a coupled earth system model is made up of *components* that interact through a coupler and are all controlled by a driver.

   In the current version of CIME, there are 7 physical components allowed.  They are:

       atmosphere, ocean, sea-ice, land surface, river, ice sheet, ocean waves

   Components are also referred to as "models".  The choice of 7 is partly historical and partly determined by the physics of the
   Earth system: these 7 components
   occupy physically distinct domains in the Earth system and/or require different numerical grids for solving.


**component types**

   For each of the 7 physical components (models), there can be three different implementations in a CIME-driven coupled model.

   *active*: Solve a complex set of equations to describe the model's behavior. Also called *prognostic* or *full* models.
   These can be full General Circulation Models. Multiple active models might be available (for example POP and MPAS-ocean to represent the global ocean) but only one ocean or atmosphere model at a time can be used in a component set.

   *data*: For some climate problems, it is necessary to reduce feedbacks within the system by replacing an active model with a
   version that sends and receives the same variables to and from other models, but with the values read from files rather
   than computed from the equations. The values received are ignored. These active-model substitutes are called *data models*.
   CIME provides data models for each of the possible components.  You could add your own data model implementation of a component
   but as for active models only one at a time can be used.

   *stub*: For some configurations, no data model is needed, so CIME provides *stub* versions that simply occupy the
   required place in the driver and do not send or receive any data.

**component set** or **compset**:   The particular combination of active, data and stub versions of the 7 components is referred to
   as a *component set* or  *compset*.  The Case Control System allows one to define
   several possible compsets and configure and run them on supported platforms. See :ref:`Component Sets<compsets>` for more information.

**grid** or **model grid**:
   Each active model must solve its equations on a numerical grid. CIME allows models within the system to have
   different grids. The resulting set of all numerical grids is called the *model grid* or sometimes just the *grid*, where
   *grid* is a unique name that denotes a set of numerical grids. Sometimes the *resolution* also refers to a specific set
   of grids.

**machine and compilers**:
   The *machine* is the computer you are using to run CIME and build and run the climate model. It could be a workstation
   or a national supercomputer. The exact name of  *machine* is typically the UNIX hostname but it could be any string.  A machine
   may have one more more versions of Fortran, C and C++ *compilers* that are needed to compile the model's source code and CIME

**case**:
    To build and execute a CIME-enabled climate model, you have to make choices of compset, model grid,
    machine and compiler. The collection of these choices, and any additional
    customizations you may make, is called the *case*.

**out-of-the-box**:
   Any case that can be defined by the coupled model's CIME configuration files and built with only basic commands in the
   CIME Case Control System is an "out-of-the-box" case.  Since CIME and its configuration files are kept with
   the model source code and version-controlled together, its possible to match supported out-of-the-box cases with specific
   versions of the model source code, promoting reproducibility and provenance.  An out-of-the-box case is also called a *base case*

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

- variable is set in ``$HOME/.cime/config`` as explained further :ref:`here<customizing-cime>`.

- variable is set in a ``$CASEROOT`` xml file

Quick start
==================

To see an example of how a case is created, configured, built and run with CIME, execute the following commands. (This assumes that CIME has been ported to your current machine).
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


Discovering available cases with **query_config**
=================================================

Your CIME-driven model has many more possible cases besides the simple one in the above Quick Start.

Use the utility `query_config <../Tools_user/query_config.html>`_  to see which out-of-the-box compsets, components, grids and machines are available for your model.

If CIME is downloaded in standalone mode, only standalone CIME compsets can be queried.

If CIME is part of a CIME-driven model, `query_config <../Tools_user/query_config.html>`_ will allow you to query all prognostic component compsets.

To see lists of available compsets, components, grids and machines, look at the **help** text::

  > query_config --help

To see all available component sets, try::

  > query_config --compsets all

**Usage examples**

To run `query_config <../Tools_user/query_config.html>`_ for compset information, follow this example, where **drv** is the component name::

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

Each model component specifies its own definitions of what can appear after the **%**  modifier in the compset longname (for example, **DOM** in **DOCN%DOM**).

To see what supported modifiers are for **DOCN**, run `query_config <../Tools_user/query_config.html>`_ as in this example::

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

For more details on how CIME determines the output for query_config, see :ref:`Component Sets<compsets>`.
