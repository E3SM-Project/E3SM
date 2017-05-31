.. _cime-internals:

==============
CIME internals
==============

The file ``$CIMEROOT/config/[cesm,acme]/config_files.xml`` contains all model specfic information that CIME utilizes to determine
compsets, compset component settings, model grids, machines, batch queue settings, and compiler settings. ``config_files.xml`` contains
the following xml nodes that are discussed below.

::

   compset definitions:
      <entry id="COMPSETS_SPEC_FILE">

   component specfic compset settings: 
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
   
   machine specfic definitions:
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

CIME looks at the xml element ``COMPSETS_SPEC_FILE`` in the file ``$CIMEROOT/config/[cesm,acme]/config_files.xml`` to determine where
to find the list of supported compsets lists all model components that set compsets. 

In the case of CESM, this xml element has the following contents

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

where ``$SRCROOT`` is the root of your CESM sandbox and contains ``$CIMEROOT`` as a sub-directory. 

.. note:: CIME searches each of the above ``<config_compsets.xml>`` files in the ``<value>`` nodes to determine if there is a match for the compset longname. If a match is found then the value of the ``component`` attribute is set for all other searches that appear below. 

.. _defining-component-specfic-compset-settings:

Where are component-specific settings defined for the target compset?
---------------------------------------------------------------------

Every model component contains a ``config_component.xml`` file that has two functions:

1. Specifying the component-specfici definitions of what can appear after the ``%`` in the compset longname, (e.g. ``DOM`` in ``DOCN%DOM``). 

2. Specifying the compset-specific ``$CASEROOT`` xml variables.  


CIME first parses the following nodes to determine the ``config_component.xml`` files for the driver. There are two such files, one is model independent
::
   
   <entry id="CONFIG_CPL_FILE">
      ... 
      <default_value>$CIMEROOT/driver_cpl/cime_config/config_component.xml</default_value>
      ..
   </entry>

and the other is model specific
::

     <entry id="CONFIG_CPL_FILE_MODEL_SPECIFIC">
        <default_value>$CIMEROOT/driver_cpl/cime_config/config_component_$MODEL.xml</default_value>
     </entry>

CIME the parses each of the following nodes 
::
     
     <entry id="CONFIG_ATM_FILE">
     <entry id="CONFIG_ESP_FILE">
     <entry id="CONFIG_ICE_FILE">
     <entry id="CONFIG_GLC_FILE">
     <entry id="CONFIG_LND_FILE">
     <entry id="CONFIG_OCN_FILE">
     <entry id="CONFIG_ROF_FILE">
     <entry id="CONFIG_WAV_FILE">

using the value of the ``component`` attribute determined 
to determine the list of ``config_component.xml`` files to use for the requested compset longname.

As an example, the possible atmosphere components for CESM have the following associated ``config_component.xml`` files.  
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

If the compset's atm ``component`` attribute is ``datm`` then the file ``$CIMEROOT/components/data_comps/datm/cime_config/config_component.xml`` will specify all possible component settings for ``DATM``.

The schema for every ``config_component.xml`` file has a ``<description>`` node that specfies all possible values that can follow the ``%`` character in the compset name.

As an example, use the command **manage_case --query-component-name cam** to query the possibly ``%`` modifiers for CESM cam.

.. _defining-pes:

Where are pe-settings defined for the target compset and model grid?
--------------------------------------------------------------------

CIME looks at the xml element ``PES_SPEC_FILE`` in the file ``$CIMEROOT/config/[cesm,acme]/config_files.xml`` to determine where
to find the supported out-of-the-box model grids for the target ``component``. Currently, this node has the following contents (where ``MODEL`` can be ``[acme, cesm]``.

Each ``component`` that sets compsets, must also have an associated ``config_pes.xml`` file that specifies an out-of-the-box pe-layout for those compsets. 
In addition, this out-of-the-box pe-layout might also have dependencies on the model grid and the target machine. 
Finally, there might be more than one-of-the-box pe-layout that could be used for a compset/grid/machine combination: one for a low processor setting and one for a high processor setting.

A typical entry in ``config_pes.xml`` looks like the following:

::

  <grid name="a%T62">
    <mach name="yellowstone|pronghorn">
      <pes pesize="any" compset="DATM%IAF">
      .......
      </pes>
    </mach>
  </grid>

Given these various dependencies, an order of precedence has been established to determine the optimal match. This order is as follows:

1. grid match

   CIME first searches the grid nodes and tries to find a grid match in ``config_grids.xml``. 
   The search is based on a regular expression match for the grid longname.
   All nodes that have a grid match are then used in the subsequent search. If there is no grid match, than all nodes that have 
   ::
   
      <grid name="any">

   are used in the subsequent search.


2. machine match

   CIME next uses the list of nodes obtained in (1.) to match on the machine name using the ``<mach>`` nodes. If there is no machine match, then all nodes in (1.)
   that have 
   ::

      <machine name="any"> 

   are used in the subsequent search.


3. pesize and compset match

   CIME next uses the list of nodes obtained in (2.) to match on pesize and compset using the ``<pes>`` nodes. If there is no match, then the node in (2.) that has
   ::

      <pes pesize="any" compset="any">

   will be the one used.


**create_newcase** outputs the matches that are found in determining the best out-of-the-box peylayout. 

.. _defining-model-grids:

Where are model grids defined?
------------------------------

CIME looks at the xml node ``GRIDS_SPEC_FILE`` in the file ``$CIMEROOT/config/[cesm,acme]/config_files.xml`` to determine where
to find the supported out-of-the-box model grids for the target model. Currently, this node has the following contents (where ``MODEL`` can be ``[acme, cesm]``.
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

CIME looks at the xml node ``MACHINE_SPEC_FILE`` in the file ``$CIMEROOT/config/[cesm,acme]/config_files.xml`` to determine where
to find the supported out-of-the-box machines for the target model. Currently, this node has the following contents (where ``MODEL`` can be ``[acme, cesm]``.
::

   <entry id="MACHINES_SPEC_FILE">
     <type>char</type>
     <default_value>$CIMEROOT/cime_config/$MODEL/machines/config_machines.xml</default_value>
     <group>case_last</group>
     <file>env_case.xml</file>
     <desc>file containing machine specifications for target model primary component (for documentation only - DO NOT EDIT)</desc>
     <schema>$CIMEROOT/cime_config/xml_schemas/config_machines.xsd</schema>
   </entry>

As part of porting, you will need to :ref:`customize the config_machines.xml file <customizing-machine-file>`. 

.. _defining-the-batch-system:

Where are batch system settings defined?
----------------------------------------

CIME looks at the xml node ``BATCH_SPEC_FILE`` in the file ``$CIMEROOT/config/[cesm,acme]/config_files.xml`` to determine where
to find the supported out-of-the-box  batch system details for the target model. Currently, this node has the following contents (where ``MODEL`` can be ``[acme, cesm]``.
::

   <entry id="BATCH_SPEC_FILE">
     <type>char</type>
     <default_value>$CIMEROOT/cime_config/$MODEL/machines/config_batch.xml</default_value>
     <group>case_last</group>
     <file>env_case.xml</file>
     <desc>file containing batch system details for target system  (for documentation only - DO NOT EDIT)</desc>
     <schema>$CIMEROOT/cime_config/xml_schemas/config_batch.xsd</schema>
   </entry>

As part of porting, you will need to :ref:`customize the config_batch.xml file <customizing-batch-file>`. 

.. _defining-compiler-settings:

Where are compiler settings defined?
------------------------------------

CIME looks at the xml element ``COMPILERS_SPEC_FILE`` in the file ``$CIMEROOT/config/[cesm,acme]/config_files.xml`` to determine where
to find the supported out-of-the-box  batch system details for the target model. Currently, this node has the following contents (where ``MODEL`` can be ``[acme, cesm]``.
::

  <entry id="COMPILERS_SPEC_FILE">
    <type>char</type>
    <default_value>$CIMEROOT/cime_config/$MODEL/machines/config_compilers.xml</default_value>
    <group>case_last</group>
    <file>env_case.xml</file>
    <desc>file containing compiler specifications for target model primary component (for documentation only - DO NOT EDIT)</desc>
    <schema>$CIMEROOT/cime_config/xml_schemas/config_compilers_v2.xsd</schema>
  </entry>

As part of porting, you will need to :ref:`customize the config_batch.xml file <customizing-compiler-file>`. 

