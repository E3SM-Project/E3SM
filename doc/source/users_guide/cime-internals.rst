.. _cime-internals:

========================
Main Configuration File
========================

The file **$CIMEROOT/config/[cesm,e3sm]/config_files.xml** contains all model-specific information that CIME uses to determine compsets, compset component settings, model grids, machines, batch queue settings, and compiler settings. It contains the following xml nodes, which are discussed below or in subsequent sections of this guide.
::

   compset definitions:
      <entry id="COMPSETS_SPEC_FILE">

   component specific compset settings:
      <entry id="CONFIG_CPL_FILE">
      <entry id="CONFIG_CPL_FILE_MODEL_SPECIFIC">
      <entry id="CONFIG_ATM_FILE">
      <entry id="CONFIG_LND_FILE">
      <entry id="CONFIG_ROF_FILE">
      <entry id="CONFIG_ICE_FILE">
      <entry id="CONFIG_OCN_FILE">
      <entry id="CONFIG_GLC_FILE">
      <entry id="CONFIG_WAV_FILE">
      <entry id="CONFIG_ESP_FILE">

   pe-settings:
      <entry id="PES_SPEC_FILE">

   grid definitions:
      <entry id="GRIDS_SPEC_FILE">

   machine specific definitions:
      <entry id="MACHINES_SPEC_FILE">
      <entry id="BATCH_SPEC_FILE">
      <entry id="COMPILERS_SPEC_FILE">
      <entry id="PIO_SPEC_FILE">

   testing:
      <entry id="CONFIG_TESTS_FILE">
      <entry id="TESTS_SPEC_FILE">
      <entry id="TESTS_MODS_DIR">
      <entry id="SYSTEM_TESTS_DIR">

   archiving:
      <entry id="LTARCHIVE_SPEC_FILE">

   CIME components namelists definitions:
      <entry id="NAMELIST_DEFINITION_FILE">

   user-mods directories:
      <entry id="USER_MODS_DIR">

.. _customizing-cime:

Customizing CIME For Your Needs
-------------------------------

CIME recognizes a user-created custom configuration directory, ``$HOME/.cime``. The contents of this directory may include any one of the following list of files:

* ``config``

   This file must have a format which follows the python config format. See `Python Config Parser Examples <https://wiki.python.org/moin/ConfigParserExamples>`_

   In the [main] block you can set the following variables:

   * ``CIME_MODEL=[e3sm, cesm]``

   * ``PROJECT=<account number>``

     This is your project account code for batch submission and/or directory priveleges

   * ``CHARGE_ACCOUNT=<account number>``

     An alternative to PROJECT for batch charging>

   * ``MAIL_USER=<email address>``

     Used request a non-default email for batch summary output

   * ``MAIL_TYPE=[never,all,begin,fail,end]``

    Any **or** all the above valid values can be set to list the batch events that emails will be sent for.

   * **create_test** input arguments

     Any argument to the **create_test** script can have its default changed by listing it here with the new default.

* ``config_machines.xml``

  This file must the same format as ``$CIMEROOT/config/$model/machines/config_machines.xml`` with the appropriate definitions for your machine.

  If you have a customized version of this file in ``$HOME/.cime``, it will **append** to the file in ``$CIMEROOT/config/$model/machines/config_machines.xml``.

* ``config_compilers.xml``

  .. todo:: Add content for config_compilers.xml

* ``config_batch.xml``

  .. todo:: Add content for config_batch.xml
