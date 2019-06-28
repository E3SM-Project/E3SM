.. _customizing-cime:

===============================
CIME user config directory
===============================

CIME recognizes a user-created custom configuration directory, ``$HOME/.cime``. The contents of this directory may include any of the following files:

* ``config``

  This file must have a format which follows the python config format. See `Python Config Parser Examples <https://wiki.python.org/moin/ConfigParserExamples>`_

  In the [main] block you can set the following variables:

  * ``CIME_MODEL=[e3sm, cesm]``

  * ``PROJECT=<account number>``

    Used to specify a project id for compute accounting and directory permissions when on a batch system.

  * ``CHARGE_ACCOUNT=<account number>``

    Used to override the accounting (only) aspect of PROJECT

  * ``MAIL_USER=<email address>``

    Used to request a non-default email for batch summary output

  * ``MAIL_TYPE=[never,all,begin,fail,end]``

    Any **or** all the above valid values can be set to list the batch events that emails will be sent for.

  * **create_test** input arguments

    Any argument to the **create_test** script can have its default changed by listing it here with the new default.

  * The following is an example ``config`` file:

    ::

       [main]
       CIME_MODEL=cesm
       SRCROOT=$CIMEROOT/..
       MAIL_TYPE=end
       [create_test]
       MAIL_TYPE=fail

* ``config_machines.xml``

  This file must the same format as **$CIMEROOT/config/$model/machines/config_machines.xml** with the appropriate definitions for your machine.

  If you have a customized version of this file in the directory ``$HOME/.cime``, it will **append** to the file in ``$CIMEROOT/config/$model/machines/config_machines.xml``.

  For an example of a **config_machines.xml** file for a linux cluster, look at **$CIMEROOT/config/xml_schemas/config_machines_template.xml**.

* ``config_compilers.xml``

  This file permits you to customize compiler settings for your machine and is appended to the file **$CIMEROOT/config/$model/machines/config_compilers.xml**.

  The following is an example of what would be needed for customized a ibm compiler flags on a BlueGeneQ machine.

  ::

     <?xml version="1.0" encoding="UTF-8"?>
     <config_compilers version="2.0">
        <compiler COMPILER="ibm" OS="BGQ">
           <FFLAGS> -g -qfullpath -qmaxmem=-1 -qspillsize=2500 -qextname=flush </FFLAGS>
	   <ADD_FFLAGS DEBUG="FALSE"> -O3 -qstrict -qinline=auto </ADD_FFLAGS>
	   <ADD_FFLAGS DEBUG="FALSE" compile_threaded="TRUE"> -qsmp=omp </ADD_FFLAGS>
	   <ADD_FFLAGS DEBUG="TRUE" compile_threaded="TRUE"> -qsmp=omp:noopt </ADD_FFLAGS>
	   <ADD_CPPDEFS> -DLINUX  </ADD_CPPDEFS>
	   <CONFIG_ARGS> --build=powerpc-bgp-linux --host=powerpc64-suse-linux </CONFIG_ARGS>
	   <LDFLAGS>  -Wl,--relax -Wl,--allow-multiple-definition </LDFLAGS>
        </compiler>
     </config_compilers>

* ``config_batch.xml``

  This file permits you to customize batch settings for you machine and is appended to the file **$CIMEROOT/config/$model/machines/config_batch.xml**.

  The following is an example of what would be needed to add batch settings for pbs on the machine brutus.

  ::

     <?xml version="1.0"?>
     <config_batch version="2.0">
        <batch_system type="pbs" MACH="brutus" >
	   <directives>
	      <directive default="/bin/bash" > -S {{ shell }}  </directive>
	   </directives>
           <queues>
	     <queue walltimemax="00:59:00" nodemin="1" nodemax="624" default="true">batch</queue>
	   </queues>
	</batch_system>
     </config_batch>
