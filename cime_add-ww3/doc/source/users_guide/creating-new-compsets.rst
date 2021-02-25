.. _creating-new-compsets:

=====================
Creating New Compsets
=====================

A description of how CIME interprets a compset name is given in the section :ref:`defining-compsets` .

To create a new compset, you will at a minimum have to:

1. edit the approprite ``config_components.xml`` file(s) to add your new requirements
2. edit associate ``namelist_definitions_xxx.xml`` in the associated ``cime_config`` directories.
   (e.g. if a change is made to the the ``config_components.xml`` for ``DOCN`` then ``namelist_definitions_docn.xml`` file will also need to be modified).

It is important to point out, that you will need expertise in the target component(s) you are trying to modify in order to add new compset functionality for that particular component.
We provide a few examples below that outline this process for a few simple cases.

Example 1:
----------

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
