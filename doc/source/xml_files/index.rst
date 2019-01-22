.. _xml_files:

#########
XML Files
#########

CIME includes a number of XML files used for configuration of a case.
General Users should never need to modify the XML files in the
CIMEROOT.  Modifcations to XML settings is case specific and the tools
**:ref:`xmlquery`** and **:ref:`xmlchange`** in every CASEROOT can be used to query
and modify these settings while ensuring the continued schema
integrity of the XML.

For advanced CIME developers, there are XML schema definition files 
in the CIMEROOT/config/xml_schemas directory that can be used with
**xmllint** to verify the XML. 

.. toctree::
   :maxdepth: 2

   e3sm.rst
   cesm.rst
   common.rst
   components.rst
   drivers.rst




