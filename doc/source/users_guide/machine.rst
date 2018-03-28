.. _machine:

========================
Defining the machine
========================

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

Batch system definition
-----------------------

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

Compiler settings
-----------------

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

