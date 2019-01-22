.. _compsets:

===============
Component sets
===============

In CIME, multiple components can define compsets that are targeted to their model development needs.

Each component supports a set of compset longnames that are used in testing and supported in out of the box configurations.

To determine if the compset name to `create_newcase  <../Tools_user/create_newcase.html>`_ is a supported component, CIME looks in the **config_files.xml** file and parses the
the xml element ``COMPSETS_SPEC_FILE`` in order to determine which component is defining the compset.

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


Every file listed in COMPSETS_SPEC_FILE will be searched for the compset specified in the call to create_newcase.

CIME will note which component's config_compsets.xml had the matching compset name and that component will be treated as
the **primary component** As an example, the primary component for a compset that has a prognostic atmosphere,
land and cice (in prescribed mode) and a data ocean is the atmosphere component (for cesm this is CAM) because the compset
is defined, using the above example, in ``$SRCROOT/components/cam/cime_config/config_compsets.xml``
In a compset where all components are prognostic, the primary component will be **allactive**.

.. _defining-compsets:

Compset longname
-------------------

Each config_compsets.xml file has a list of allowed component sets in the form of a longname and an alias.

A compset longname has this form::

  TIME_ATM[%phys]_LND[%phys]_ICE[%phys]_OCN[%phys]_ROF[%phys]_GLC[%phys]_WAV[%phys]_ESP[_BGC%phys]

Supported values for each element of the longname::

  TIME = model time period (e.g. 1850, 2000, 20TR, RCP8...)

  CIME supports the following values for ATM,LND,ICE,OCN,ROF,GLC,WAV and ESP.
  ATM  = [DATM, SATM, XATM]
  LND  = [DLND, SLND, XLND]
  ICE  = [DICE, SICE, SICE]
  OCN  = [DOCN, SOCN, XOCN]
  ROF  = [DROF, SROF, XROF]
  GLC  = [SGLC, XGLC]
  WAV  = [SWAV, XWAV]
  ESP  = [SESP]

A CIME-driven model may have other options available.  Use `query_config  <../Tools_user/query_config.html>`_ to determine the available options.

The OPTIONAL %phys attributes specify sub-modes of the given system.
For example, DOCN%DOM is the DOCN data ocean (rather than slab-ocean) mode.
**All** the possible %phys choices for each component are listed by calling `query_config --compsets <../Tools_user/query_config.html>`_.
**All** data models have a %phys option that corresponds to the data model mode.

.. _defining-component-specific-compset-settings:

Component specific settings in a compset
-----------------------------------------

Every model component also contains a **config_component.xml** file that has two functions:

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

The schema for every **config_component.xml** file has a ``<description>`` node that specifies all possible values that can follow the ``%`` character in the compset name.

To list the possible values, use the `query_config --component datm <../Tools_user/query_config.html>`_ command.

.. _creating-new-compsets:

Creating New Compsets
-----------------------

A description of how CIME interprets a compset name is given in the section :ref:`defining-compsets` .

To create a new compset, you will at a minimum have to:

1. edit the approprite ``config_components.xml`` file(s) to add your new requirements
2. edit associate ``namelist_definitions_xxx.xml`` in the associated ``cime_config`` directories.
   (e.g. if a change is made to the the ``config_components.xml`` for ``DOCN`` then ``namelist_definitions_docn.xml`` file will also need to be modified).

It is important to point out, that you will need expertise in the target component(s) you are trying to modify in order to add new compset functionality for that particular component.
We provide a few examples below that outline this process for a few simple cases.


Say you want to add a new mode, ``FOO``,  to the data ocean model, ``DOCN``. Lets call this mode, ``FOO``.
This would imply when parsing the compset longname, CIME would need to be able to recognize the string ``_DOCN%FOO_``.
To enable this, you will need to do the following:

1. edit ``$CIMEROOT/src/components/data_comps/docn/cime_config/config_component.xml`` (see the ``FOO`` additions below).

   * add an entry to the ``<description modifier block="1">`` block as shown below ::

       <description modifier_mode="1">
          <desc ocn="DOCN...[%FOO]">DOCN </desc>
          ...
          <desc option="FOO"> new  mode</desc>
          ....
       </description>

   * add an entry to the ``<entry id="DOCN_MODE">`` block as shown below::

       <entry id="DOCN_MODE">
          ....
          <values match="last">
          ....
          <value compset="_DOCN%FOO_" >prescribed</value>
          ...
       </entry>

   * modify any of the other xml entries that need a new dependence on ``FOO``

2. edit ``$CIMEROOT/src/components/data_comps/docn/cime_config/namelist_definition_docn.xml`` (see the ``FOO`` additions below).

   * add an entry to the ``datamode`` block as shown below. ::

       <entry id="datamode">
          ....
          <valid_values>...FOO</valid_values>
          ...
       </entry>

   * add additional changes to ``namelist_definition_docn.xml`` for the new mode


.. todo:: Add additional examples for creating a case
