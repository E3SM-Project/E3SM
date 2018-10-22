.. _changing-data-model-namelists:

Customizing data model namelists and stream files
---------------------------------------------------

Data Atmosphere (DATM)
~~~~~~~~~~~~~~~~~~~~~~

DATM is discussed in detail in :ref:`data atmosphere overview <data-atm>`.
DATM can be user-customized by changing either its  *namelist input files* or its *stream files*.
The namelist file for DATM is **datm_in** (or **datm_in_NNN** for multiple instances).

- To modify **datm_in** or **datm_in_NNN**, add the appropriate keyword/value pair(s) for the namelist changes that you want at the end of the **user_nl_datm** file or the **user_nl_datm_NNN** file in ``$CASEROOT``.

- To modify the contents of a DATM stream file, first run **preview_namelists** to list the *streams.txt* files in the **CaseDocs/** directory. Then, in the same directory:

  1. Make a *copy* of the file with the string *"user_"* prepended.
        ``> cp datm.streams.txt.[extension] user_datm.streams.txt[extension.``
  2. **Change the permissions of the file to be writeable.** (chmod 644)
        ``chmod 644 user_datm.streams.txt[extension``
  3. Edit the **user_datm.streams.txt.*** file.

**Example**

If the stream txt file is **datm.streams.txt.CORE2_NYF.GISS**, the modified copy should be **user_datm.streams.txt.CORE2_NYF.GISS**.
After calling **preview_namelists** again, your edits should appear in **CaseDocs/datm.streams.txt.CORE2_NYF.GISS**.

Data Ocean (DOCN)
~~~~~~~~~~~~~~~~~~~~~~

DOCN is discussed in detail in :ref:`data ocean overview <data-ocean>`.
DOCN can be user-customized by changing either its namelist input or its stream files.
The namelist file for DOCN is **docn_in** (or **docn_in_NNN** for multiple instances).

- To modify **docn_in** or **docn_in_NNN**, add the appropriate keyword/value pair(s) for the namelist changes that you want at the end of the file in ``$CASEROOT``.

- To modify the contents of a DOCN stream file, first run **preview_namelists** to list the *streams.txt* files in the **CaseDocs/** directory. Then, in the same directory:

  1. Make a *copy* of the file with the string *"user_"* prepended.
        ``> cp docn.streams.txt.[extension] user_docn.streams.txt[extension.``
  2. **Change the permissions of the file to be writeable.** (chmod 644)
        ``chmod 644 user_docn.streams.txt[extension``
  3. Edit the **user_docn.streams.txt.*** file.

**Example**

As an example, if the stream text file is **docn.stream.txt.prescribed**, the modified copy should be **user_docn.streams.txt.prescribed**.
After changing this file and calling **preview_namelists** again, your edits should appear in **CaseDocs/docn.streams.txt.prescribed**.

Data Sea-ice (DICE)
~~~~~~~~~~~~~~~~~~~~~~

DICE is discussed in detail in :ref:`data sea-ice overview <data-seaice>`.
DICE can be user-customized by changing either its namelist input or its stream files.
The namelist file for DICE is ``dice_in`` (or ``dice_in_NNN`` for multiple instances) and its values can be changed by editing the ``$CASEROOT`` file ``user_nl_dice`` (or ``user_nl_dice_NNN`` for multiple instances).

- To modify **dice_in** or **dice_in_NNN**, add the appropriate keyword/value pair(s) for the namelist changes that you want at the end of the file in ``$CASEROOT``.

- To modify the contents of a DICE stream file, first run **preview_namelists** to list the *streams.txt* files in the **CaseDocs/** directory. Then, in the same directory:

  1. Make a *copy* of the file with the string *"user_"* prepended.
        ``> cp dice.streams.txt.[extension] user_dice.streams.txt[extension.``
  2. **Change the permissions of the file to be writeable.** (chmod 644)
        ``chmod 644 user_dice.streams.txt[extension``
  3. Edit the **user_dice.streams.txt.*** file.

Data Land (DLND)
~~~~~~~~~~~~~~~~~~~~~~

DLND is discussed in detail in :ref:`data land overview <data-lnd>`.
DLND can be user-customized by changing either its namelist input or its stream files.
The namelist file for DLND is ``dlnd_in`` (or ``dlnd_in_NNN`` for multiple instances) and its values can be changed by editing the ``$CASEROOT`` file ``user_nl_dlnd`` (or ``user_nl_dlnd_NNN`` for multiple instances).

- To modify **dlnd_in** or **dlnd_in_NNN**, add the appropriate keyword/value pair(s) for the namelist changes that you want at the end of the file in ``$CASEROOT``.

- To modify the contents of a DLND stream file, first run **preview_namelists** to list the *streams.txt* files in the **CaseDocs/** directory. Then, in the same directory:

  1. Make a *copy* of the file with the string *"user_"* prepended.
        ``> cp dlnd.streams.txt.[extension] user_dlnd.streams.txt[extension.``
  2. **Change the permissions of the file to be writeable.** (chmod 644)
        ``chmod 644 user_dlnd.streams.txt[extension``
  3. Edit the **user_dlnd.streams.txt.*** file.

Data River (DROF)
~~~~~~~~~~~~~~~~~~~~~~

DROF is discussed in detail in :ref:`data river overview <data-river>`.
DROF can be user-customized by changing either its namelist input or its stream files.
The namelist file for DROF is ``drof_in`` (or ``drof_in_NNN`` for multiple instances) and its values can be changed by editing the ``$CASEROOT`` file ``user_nl_drof`` (or ``user_nl_drof_NNN`` for multiple instances).

- To modify **drof_in** or **drof_in_NNN**, add the appropriate keyword/value pair(s) for the namelist changes that you want at the end of the file in ``$CASEROOT``.

- To modify the contents of a DROF stream file, first run **preview_namelists** to list the *streams.txt* files in the **CaseDocs/** directory. Then, in the same directory:

  1. Make a *copy* of the file with the string *"user_"* prepended.
        ``> cp drof.streams.txt.[extension] user_drof.streams.txt[extension.``
  2. **Change the permissions of the file to be writeable.** (chmod 644)
        ``chmod 644 user_drof.streams.txt[extension``
  3. Edit the **user_drof.streams.txt.*** file.

