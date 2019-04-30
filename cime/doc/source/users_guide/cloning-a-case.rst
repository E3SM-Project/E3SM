.. _cloning-a-case:

**************************
Cloning a Case
**************************

If you have access to a run that you want to clone, the
`create_clone  <../Tools_user/create_clone.html>`_ command will create a new case and run `case.setup  <../Tools_user/case.setup.html>`_
while preserving local modifications to the case.

Here is a simple example:
::

   > cd $CIMEROOT/scripts
   > create_clone --case $CASEROOT --clone $CLONEROOT
   > cd $CASEROOT
   > case.build
   > case.submit

The `create_clone <../Tools_user/create_clone.html>`_ script preserves any local namelist modifications
made in the **user_nl_xxxx** files as well as any source code
modifications in the **SourceMods/** directory tree. Otherwise, your **$CASEROOT** directory
directory will appear as if `create_newcase  <../Tools_user/create_newcase.html>`_ had just been run.

**Important**: Do not change anything in the **env_case.xml** file.

See the **help** text for more usage information.

::

   > create_clone --help

`create_clone  <../Tools_user/create_clone.html>`_ has several useful optional arguments. To point to
the executable of the original case you are cloning from.

::

   > create_clone --case $CASEROOT --clone $CLONEROOT --keepexe
   > cd $CASEROOT
   > case.submit

If the ``--keepexe`` optional argument is used, then no SourceMods
will be permitted in the cloned directory.  A link will be made when
the cloned case is created pointing the cloned SourceMods/ directory
to the original case SourceMods directory.

.. warning:: No changes should be made to ``env_build.xml`` or ``env_mach_pes.xml`` in the cloned directory.

`create_clone <../Tools_user/create_clone.html>`_ also permits you to invoke the ``shell_commands``
 and ``user_nl_xxx`` files in a user_mods directory by calling:

::

   > create_clone --case $CASEROOT --clone $CLONEROOT --user-mods-dir USER_MODS_DIR [--keepexe]

Note that an optional ``--keepexe`` flag can also be used in this case.

.. warning:: If there is a ``shell_commands`` file, it should not have any changes to xml variables in either ``env_build.xml`` or ``env_mach_pes.xml``.

Another approach to duplicating a case is to use the information in
the case's **README.case** and **CaseStatus** files to create a new
case and duplicate the relevant `xmlchange <../Tools_user/xmlchange.html>`_ commands that were
issued in the original case. This alternative will *not* preserve any
local modifications that were made to the original case, such as
source-code or build-script revisions; you will need to import those
changes manually.
