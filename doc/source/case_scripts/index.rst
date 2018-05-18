.. _templates:

#########
Templates
#########

CIME includes a script files in the CASEROOT that are created 
by **create_newcase**, **create_clone**, **create_test** and **case.setup**.
The **.case.run** and **.case.test**
scripts are hidden files in the CASEROOT because these scripts should
only be run via **case.submit**.

.. toctree::
   :maxdepth: 1

   case.run
   case.test
   case.st_archive
