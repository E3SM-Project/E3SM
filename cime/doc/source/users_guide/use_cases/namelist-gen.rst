.. _namelist-gen:

Changing namelist values
=========================

All CIME-compliant components generate their namelist settings using a **buildnml** file located in the component's **cime_config** directory
For example, the CIME data atmosphere model (DATM) generates namelists using the script **$CIMEROOT/components/data_comps/datm/cime_config/buildnml**.

User-specific component namelist changes should be made only by:

- editing the **$CASEROOT/user_nl_xxx** files.

- using :ref:`xmlchange<modifying-an-xml-file>` to modify xml variables in **env_run.xml**, **env_build.xml** or **env_mach_pes.xml**.

You can preview the component namelists by running **preview_namelists** from ``$CASEROOT``.
This results in the creation of component namelists (for example, atm_in, lnd_in, and so on) in **$CASEROOT/CaseDocs/**. The namelist files are there only for user reference and SHOULD NOT BE EDITED since they are overwritten every time **preview_namelists**  and  **case.submit** are called.

Here are two examples of how to invoke **xmlchange**:

::

   xmlchange <entry id>=<value>
   -- OR --
   xmlchange -id <entry id> -val <name> -file <filename>

The ``-id`` argument identifies the variable to be changed, and ``-val`` is the intended value of that variable. See the **help** text for more usage information:
::

   xmlchange --help
