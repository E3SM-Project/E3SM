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

