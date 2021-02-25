.. _creating-a-case:

*********************************
Creating a Case
*********************************

This and following sections provide more detail about the basic commands of the CIME Case Control System: **create_newcase**,
**case.setup**, **case.build** and **case.submit**. On a supported system, you can configure, build and run many complex
climate model configurations with only these 4 commands.

To see if your machine is supported try::

  > query_config --machines

If you are not on an out-of-the box CIME-supported platform, you will need to :ref:`port <porting>` CIME to your system before proceeding.

===================================
Calling **create_newcase**
===================================

The first step in creating a CIME-based experiment is to use `create_newcase  <../Tools_user/create_newcase.html>`_.

See the options for `create_newcase  <../Tools_user/create_newcase.html>`_ in the  **help** text.::

  > create_newcase --help

The only required arguments to `create_newcase  <../Tools_user/create_newcase.html>`_ are::

  > create_newcase --case CASENAME --compset COMPSET --res GRID

Creating a CIME experiment or *case* requires, at a minimum, specifying a compset and a model grid and a case directory.
CIME supports out-of-the-box *component sets*, *model grids* and *hardware platforms* (machines).

.. warning::
   The ``--case`` argument must be a string and may not contain any of the following special characters
   ::
      > + * ? < > { } [ ] ~ ` @ :

The ``--case`` argument is used to define the name of your case, a very important piece of
metadata that will be used in filenames, internal metadata and directory paths. The
``CASEROOT`` is a directory create_newcase will create with the same name as the
``CASENAME``. If ``CASENAME`` is simply a name (not a path), ``CASEROOT`` is created in
the directory where you execute create_newcase. If ``CASENAME`` is a relative or absolute
path, ``CASEROOT`` is created there, and the name of the case will be the last component
of the path.

======================================
Results of calling **create_newcase**
======================================

Suppose **create_newcase** was called as follows.
Here, $CIMEROOT is the full pathname of the root directory of the CIME distribution::

  > cd $CIMEROOT/scripts
  > create_newcase --case ~/cime/example1 --compset A --res f09_g16_rx1

In the example, the command creates a ``$CASEROOT`` directory: ``~/cime/example1``.
If that directory already exists, a warning is printed and the command aborts.

In the argument to ``--case``, the case name is taken from the string after the last slash
--- so here the case name is ``example1``.

The output from create_newcase includes information such as.

- The compset longname is ``2000_DATM%NYF_SLND_DICE%SSMI_DOCN%DOM_DROF%NYF_SGLC_SWAV``
- The model resolution is ``a%0.9x1.25_l%0.9x1.25_oi%gx1v6_r%r05_m%gx1v6_g%null_w%null``

`create_newcase  <../Tools_user/create_newcase.html>`_ installs files in ``$CASEROOT`` that will build and run the model and to optionally archive the case on the target platform.

Running `create_newcase  <../Tools_user/create_newcase.html>`_ creates the following scripts, files and directories in ``$CASEROOT``:

**User Scripts**

- `case.build  <../Tools_user/case.build.html>`_
     Script to build component and utility libraries and model executable.

- `case.setup  <../Tools_user/case.setup.html>`_
    Script used to set up the case (create the case.run script, Macros file and user_nl_xxx files).

- `case.st_archive <../Tools_user/case.st_archive.html>`_
     Script to perform short term archiving to disk for your case output. Note that this script is run automatically by the normal CIME workflow.

- `case.submit <../Tools_user/case.submit.html>`_
     Script to submit the case to run using the machine's batch queueing system.

- `case.cmpgen_namelist <../Tools_user/case.submit.html>`_
     Script to perform namelist baseline operations (compare, generate, or both)."

- `xmlchange <../Tools_user/xmlchange.html>`_
     Script to modify values in the xml files.

- `xmlquery <../Tools_user/xmlquery.html>`_
     Script to query values in the xml files.

- `preview_namelists <../Tools_user/preview_namelists.html>`_
     Script for users to see their component namelists in ``$CASEROOT/CaseDocs`` before running the model.

- `preview_run <../Tools_user/preview_run.html>`_
     Script for users to see batch submit and mpirun command."

- `check_input_data <../Tools_user/check_input_data.html>`_
     Script for checking for various input data sets and moving them into place.

- `check_case <../Tools_user/check_case.html>`_
     Script to verify case is set up correctly.

- `pelayout <../Tools_user/pelayout.html>`_
     Script to query and modify the NTASKS, ROOTPE, and NTHRDS for each component model.
     This a convenience script that can be used in place of `xmlchange <../Tools_user/xmlchange.html>`_ and `xmlquery <../Tools_user/xmlquery.html>`_.

**XML Files**

- env_archive.xml
   Defines patterns of files to be sent to the short-term archive.
   You can edit this file at any time. You **CANNOT** use `xmlchange <../Tools_user/xmlchange.html>`_  to modify variables in this file."

- env_mach_specific.xml
   Sets a number of machine-specific environment variables for building and/or running.
   You **CANNOT** use `xmlchange <../Tools_user/xmlchange.html>`_  to modify variables in this file.

- env_build.xml
   Sets model build settings. This includes component resolutions and component compile-time configuration options.
   You must run the case.build command after changing this file.

- env_run.xml
   Sets runtime settings such as length of run, frequency of restarts, output of coupler diagnostics, and short-term and long-term archiving.
   This file can be edited at any time before a job starts.

- env_mach_pes.xml
   Sets component machine-specific processor layout (see changing pe layout ).
   The settings in this are critical to a well-load-balanced simulation (see :ref:`load balancing <optimizing-processor-layout>`).

- env_batch.xml
   Sets batch system settings such as wallclock time and queue name."

**User Source Mods Directory**

- SourceMods
   Top-level directory containing subdirectories for each compset component where you can place modified source code for that component.
   You may also place modified buildnml and buildlib scripts here."

**Provenance**

- README.case
   File detailing `create_newcase  <../Tools_user/create_newcase.html>`_ usage.
   This is a good place to keep track of runtime problems and changes."

- CaseStatus
   File containing a list of operations done in the current case.


**Non-modifiable work directories**

- Buildconf,
   Work directory containing scripts to generate component namelists and component and utility libraries (PIO or MCT, for example). You should never have to edit the contents of this directory.

- LockedFiles/
   Work directory that holds copies of files that should not be changed. Certain xml files are *locked* after their variables have been used by should no longer be changed (see below).

- Tools/
   Work directory containing support utility scripts. You should never need to edit the contents of this directory."

===================================
Locked files in your case directory
===================================

The ``$CASEROOT`` xml files are organized so that variables can be
locked at certain points after they have been resolved (used) in other
parts of the scripts system.

CIME does this by *locking* a file in ``$CASEROOT/LockedFiles`` and
not permitting you to modify that file unless, depending on the file,
you call `case.setup --clean <../Tools_user/case.setup.html>`_ or
`case.build --clean <../Tools_user/case.build.html>`_ .

CIME locks your ``$CASEROOT`` files according to the following rules:

- Locks variables in **env_case.xml** after `create_newcase  <../Tools_user/create_newcase.html>`_.
   The **env_case.xml** file can never be unlocked.

- Locks variables in **env_mach_pes.xml** after `case.setup  <../Tools_user/case.setup.html>`_.
   To unlock **env_mach_pes.xml**, run `case.setup --clean <../Tools_user/case.setup.html>`_.

- Locks variables in **env_build.xml** after completion of `case.build  <../Tools_user/case.build.html>`_.
   To unlock **env_build.xml**, run `case.build --clean  <../Tools_user/case.build.html>`_

- Variables in **env_run.xml**, **env_batch.xml** and **env_archive.xml** are never locked, and most can be changed at any time.

- There are some exceptions in the **env_batch.xml** file.

===================================
Adding a --user-mods-dir argument to **create_newcase**
===================================

A user may want to customize a target case with a combination of
``user_nl_xxx`` file modifications and/or ``SourceMods`` for some
components and/or **xmlchange** commands. As an example, the user
might want to carry out a series of experiments based on a common set
of changes to the namelists, source code and/or case xml settings.
Rather than make these changes each time a new experimental
``CASEROOT`` is generated, the user can create a directory on local
disk with a set of changes that will be applied to each case.

As an example, the directory could contain the following files: ::

  > user_nl_cpl
  > shell_commands  (this would contain ./xmlchange commands)
  > SourceMods/src.cam/dyncomp.F90

When the user calls **create_newcase** with the ``--user-mods-dir`` pointing to the
full pathname of the directory containing these changes, then the ``CASEROOT`` will be
created with these changes applied.
