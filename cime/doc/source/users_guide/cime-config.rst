.. _customizing-cime:

===============================
CIME user config directory
===============================

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

