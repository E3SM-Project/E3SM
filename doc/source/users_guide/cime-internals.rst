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

.. _defining-pes:

Where are pe-settings defined for the target compset and model grid?
--------------------------------------------------------------------

CIME looks at the xml element ``PES_SPEC_FILE`` in the **config_files.xml** file to determine where
to find the supported out-of-the-box model grids for the target component.

Each component that sets compsets has an associated **config_pes.xml** file that specifies an out-of-the-box pe-layout for those compsets.
The pe-layout might also have dependencies on the model grid and the target machine.
Finally, there might be more than one out-of-the-box pe-layout that could be used for a compset/grid/machine combination: one for a low processor setting and one for a high processor setting.

A typical entry in **config_pes.xml** looks like this:

::

  <grid name="a%T62">
    <mach name="cheyenne">
      <pes pesize="any" compset="DATM%IAF">
      .......
      </pes>
    </mach>
  </grid>

Given the various dependencies, CIME uses an order of precedence to determine the optimal match. This order is as follows:

1. grid match

   CIME first searches the grid nodes for a grid match in **config_grids.xml**.
   The search is based on a regular expression match for the grid longname.
   All nodes that have a grid match are used in the subsequent search. If there is no grid match, all nodes that have ``<grid name="any">`` are used in the subsequent search.


2. machine match

   CIME next uses the list of nodes obtained in the grid match to search for the machine name using the ``<mach>`` nodes. If there is no machine match, then all nodes with ``<machine name="any">`` are used in the subsequent search.


3. pesize and compset match

   CIME next uses the list of nodes obtained in the machine match to search for pesize and compset using the ``<pes>`` nodes. If there is no match, the node with ``<pes pesize="any" compset="any">`` is used.

The **create_newcase** script outputs the matches that are found in determining the best out-of-the-box pe-layout.

.. _defining-machines:

Where are machines defined?
---------------------------

CIME looks at the xml node ``MACHINE_SPEC_FILE`` in the **config_files.xml** file to identify supported out-of-the-box machines for the target model. The node has the following contents:
::

   <entry id="MACHINES_SPEC_FILE">
     <type>char</type>
     <default_value>$CIMEROOT/cime_config/$MODEL/machines/config_machines.xml</default_value>
     <group>case_last</group>
     <file>env_case.xml</file>
     <desc>file containing machine specifications for target model primary component (for documentation only - DO NOT EDIT)</desc>
     <schema>$CIMEROOT/cime_config/xml_schemas/config_machines.xsd</schema>
   </entry>

When porting, you will need to :ref:`customize the config_machines.xml file <customizing-machine-file>`.

.. _defining-the-batch-system:

Where are batch system settings defined?
----------------------------------------

CIME looks at the xml node ``BATCH_SPEC_FILE`` in the **config_files.xml** file to identify supported out-of-the-box batch system details for the target model. The node has the following contents:
::

   <entry id="BATCH_SPEC_FILE">
     <type>char</type>
     <default_value>$CIMEROOT/cime_config/$MODEL/machines/config_batch.xml</default_value>
     <group>case_last</group>
     <file>env_case.xml</file>
     <desc>file containing batch system details for target system  (for documentation only - DO NOT EDIT)</desc>
     <schema>$CIMEROOT/cime_config/xml_schemas/config_batch.xsd</schema>
   </entry>

When porting, you will need to :ref:`customize the config_batch.xml file <customizing-batch-file>`.

.. _defining-compiler-settings:

Where are compiler settings defined?
------------------------------------

CIME looks at the xml element ``COMPILERS_SPEC_FILE`` in the **config_files.xml** file to identify supported out-of-the-box compiler details for the target model. The node has the following contents:
::

  <entry id="COMPILERS_SPEC_FILE">
    <type>char</type>
    <default_value>$CIMEROOT/cime_config/$MODEL/machines/config_compilers.xml</default_value>
    <group>case_last</group>
    <file>env_case.xml</file>
    <desc>file containing compiler specifications for target model primary component (for documentation only - DO NOT EDIT)</desc>
    <schema>$CIMEROOT/cime_config/xml_schemas/config_compilers_v2.xsd</schema>
  </entry>

When porting, you will need to :ref:`customize the config_compilers.xml file <customizing-compiler-file>`.

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
