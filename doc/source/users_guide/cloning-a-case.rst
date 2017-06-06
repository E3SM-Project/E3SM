.. _cloning-a-case:

**************************
Cloning a Case
**************************

If you have access to a run that you want to clone, the **create_clone** command will create a new case while preserving local modifications to the case.

Here is a simple example:
::

   > cd $CIMEROOT/scripts
   > create_clone --case $CASEROOT --clone $CLONEROOT

The **create_clone** script preserves any local namelist modifications made in the **user_nl_xxxx** files as well as any source code modifications in the **SourceMods** tree. Otherwise, **$CASEROOT/** directory will appear as if **create_newcase** had just been run.

**Important**: Do not change anything in the **env_case.xml** file.

See the **help** text for more usage information.
::

   > create_clone --help

Another approach to duplicating a case is to use the information in the case's **README.case** and **CaseStatus** files to create a new case and duplicate the relevant **xlmchange** commands that were issued in the original case. This alternative will *not* preserve any local modifications that were made to the original case, such as source-code or build-script revisions; you will need to import those changes manually.
