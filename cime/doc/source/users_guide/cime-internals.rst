.. _cime-internals:

==============
CIME internals
==============

The file **$CIMEROOT/config/[cesm,acme]/config_files.xml** contains all model-specific information that CIME uses to determine compsets, compset component settings, model grids, machines, batch queue settings, and compiler settings. It contains the following xml nodes, which are discussed below or in subsequent sections of this guide.
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

.. _defining-compsets:

Where are compsets defined?
---------------------------

CIME looks at the xml element ``COMPSETS_SPEC_FILE`` in the **config_files.xml** file to determine which xml file to read in order to satisfy the configuration specifications defined in the compset longname.

In the case of CESM, this xml element has the contents shown here, where ``$SRCROOT`` is the root of your CESM sandbox and contains ``$CIMEROOT`` as a subdirectory:

::

     <entry id="COMPSETS_SPEC_FILE">
       <type>char</type>
       <default_value>unset</default_value>
       <values>
         <value component="allactive">$SRCROOT/cime_config/config_compsets.xml</value>
         <value component="drv"      >$CIMEROOT/src/drivers/mct/cime_config/config_compsets.xml</value>
         <value component="cam"      >$SRCROOT/components/cam/cime_config/config_compsets.xml</value>
         <value component="cism"     >$SRCROOT/components/cism/cime_config/config_compsets.xml</value>
         <value component="clm"      >$SRCROOT/components/clm/cime_config/config_compsets.xml</value>
         <value component="cice"     >$SRCROOT/components/cice/cime_config/config_compsets.xml</value>
         <value component="pop"      >$SRCROOT/components/pop/cime_config/config_compsets.xml</value>
       </values>
       <group>case_last</group>
       <file>env_case.xml</file>
       <desc>file containing specification of all compsets for primary component (for documentation only - DO NOT EDIT)</desc>
       <schema>$CIMEROOT/config/xml_schemas/config_compsets.xsd</schema>
     </entry>

.. _defining-component-specific-compset-settings:

Where are component-specific settings defined for the target compset?
---------------------------------------------------------------------

Every model component contains a **config_component.xml** file that has two functions:

1. Specifying the component-specific definitions of what can appear after the ``%`` in the compset longname, (for example, ``DOM`` in ``DOCN%DOM``).

2. Specifying the compset-specific ``$CASEROOT`` xml variables.

CIME first parses the following nodes to identify appropriate **config_component.xml** files for the driver. There are two such files; one is model-independent and the other is model-specific.
::

   <entry id="CONFIG_CPL_FILE">
      ...
      <default_value>$CIMEROOT/driver_cpl/cime_config/config_component.xml</default_value>
      ..
      </entry>

     <entry id="CONFIG_CPL_FILE_MODEL_SPECIFIC">
        <default_value>$CIMEROOT/driver_cpl/cime_config/config_component_$MODEL.xml</default_value>
     </entry>

CIME then parses each of the nodes listed below, using using the value of the *component* attribute to determine which xml files to use for the requested compset longname.
::

     <entry id="CONFIG_ATM_FILE">
     <entry id="CONFIG_ESP_FILE">
     <entry id="CONFIG_ICE_FILE">
     <entry id="CONFIG_GLC_FILE">
     <entry id="CONFIG_LND_FILE">
     <entry id="CONFIG_OCN_FILE">
     <entry id="CONFIG_ROF_FILE">
     <entry id="CONFIG_WAV_FILE">

As an example, the possible atmosphere components for CESM have the following associated xml files.
::

     <entry id="CONFIG_ATM_FILE">
       <type>char</type>
       <default_value>unset</default_value>
       <values>
         <value component="cam" >$SRCROOT/components/cam/cime_config/config_component.xml</value>
         <value component="datm">$CIMEROOT/components/data_comps/datm/cime_config/config_component.xml</value>
         <value component="satm">$CIMEROOT/components/stub_comps/satm/cime_config/config_component.xml</value>
         <value component="xatm">$CIMEROOT/components/xcpl_comps/xatm/cime_config/config_component.xml</value>
       </values>
       <group>case_last</group>
       <file>env_case.xml</file>
       <desc>file containing specification of component specific definitions and values(for documentation only - DO NOT EDIT)</desc>
       <schema>$CIMEROOT/cime_config/xml_schemas/entry_id.xsd</schema>
     </entry>

If the compset's atm component attribute is ``datm``, the file ``$CIMEROOT/components/data_comps/datm/cime_config/config_component.xml`` specifies all possible component settings for ``DATM``.

The schema for every **config_component.xml** file has a ``<description>`` node that specifies all possible values that can follow the ``%`` character in the compset name. To list the possible values, use the **query_case** command with ``--components`` as shown in this example for CAM:
::

  query_case --components cam

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

.. _defining-model-grids:

Where are model grids defined?
------------------------------

CIME looks at the xml node ``GRIDS_SPEC_FILE`` in the **config_files.xml** file to identify supported out-of-the-box model grids for the target model. The node has the following contents:
::

   <entry id="GRIDS_SPEC_FILE">
     <type>char</type>
     <default_value>$CIMEROOT/cime_config/$MODEL/config_grids.xml</default_value>
     <group>case_last</group>
     <file>env_case.xml</file>
     <desc>file containing specification of all supported model grids, domains and mapping files (for documentation only - DO NOT EDIT)</desc>
     <schema>$CIMEROOT/cime_config/xml_schemas/config_grids_v2.xsd</schema>
   </entry>

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


