.. _pesthreads:

==================================
Controlling processors and threads
==================================

Once a compset and resolution for a case has been defined, CIME provides ways to define how many processors and
threads the case will use.


.. _defining-pes:

pe-settings for a case
-------------------------

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

Threading control
-------------------------

.. todo:: Add threading control info

SMT
-------------------------

.. todo:: Add SMT info
