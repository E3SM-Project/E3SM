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

Setting defaults
=================

Before using any CIME commands, set the ``CIME_MODEL`` environment variable. In csh, use **setenv** as shown and replace 
**<model>** with the appropriate text. Current possibilities are "acme" or "cesm."
::

   > setenv CIME_MODEL <model>


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
   "config/acme/", "ACME-specific configuration options"
   "src/components/", "CIME-provided components including data and stub models"
   "src/drivers/", "CIME-provided main driver for a climate model"
   "src/externals/", "Software provided with CIME for building a climate model"
   "src/share/", "Model source code provided by CIME and used by multiple components"
   "scripts/lib/", "Infrastructure source code for CIME scripts and functions"
   "scripts/Tools/", "Auxiliary tools; scripts and functions"

Discovering available cases
==============================

To identify which compsets, grids and machines your CIME-enabled model supports, use the **query_config** command found in **cime/scripts**.  See the **help** text for more information.

::

   > ./query_config --help

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
