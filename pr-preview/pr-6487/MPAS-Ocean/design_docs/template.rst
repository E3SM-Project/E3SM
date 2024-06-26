
Some descriptive title
======================

date: YYYY/MM/DD

Contributors:



Summary
-------

The purpose of this section is to summarize what capability is to be added to
MPAS-Model through this design process. It should be clear what new code will do
that the current code does not. Summarizing the primary challenges with respect
to software design and implementation is also appropriate for this section.
Finally, this statement should contain a general statement with regard to what
is "success."


Requirements
------------

Requirement: name-of-topic-here
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: YYYY/MM/DD

Contributors: (add your name to this list if it does not appear)


Each requirement is to be listed under a "section" heading, as there will be a
one-to-one correspondence between requirements, design, proposed implementation
and testing. Requirements should not discuss technical software issues, but
rather focus on model capability. To the extent possible, requirements should
be relatively independent of each other, thus allowing a clean algorithm design,
implementation and testing plan.


Algorithm Design (optional)
---------------------------

Algorithm Design: name-of-topic-here (same as Requirement)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: YYYY/MM/DD

Contributors: (add your name to this list if it does not appear)

For each requirement, you can provide an algorithm design that is intended to
meet that requirement. Algorithm can include detailed technical discussions of
PDEs, algorithms, solvers and similar, as well as technical discussion of
performance issues. In general, this section should steer away from a detailed
discussion of low-level software issues such as variable declarations,
interfaces and sequencing.

If you want to add math, the syntax is almost identical to Latex:

.. math::

   (a + b)^2  &=  (a + b)(a + b) \\
              &=  a^2 + 2ab + b^2

If you want to add an image, keep in mind that these should be quite small
(jpegs are preferred but small pngs are okay) to keep the size of the repo
small.  Here is an example:

.. image:: small_images/ocean.jpg
   :width: 500 px
   :align: center


Implementation
--------------

Implementation: name-of-topic-here (same as Requirement)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: YYYY/MM/DD

Contributors: (add your name to this list if it does not appear)

This section should detail the plan for implementing the algorithm design for
this requirement. In general, this section is software-centric with a focus on
software implementation. Pseudo code is appropriate in this section. Links to
actual source code are appropriate. Project management items, such as git
branches, timelines and staffing are also appropriate. Pseudo code can be
included via blocks like

.. code-block:: python

   def example_function(foo):
       return foo**2.0


Testing
-------

Testing and Validation: name-of-topic-here (same as Requirement)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: YYYY/MM/DD

Contributors: (add your name to this list if it does not appear)

How will the implementation of this requirement be tested, showing that we have
met the requirement? Which tests from the regression suites are appropriate?
How would they need to be configured or modified to test that the new software
is working properly?  What additions or modifications to the nightly (or
another) regression suite might be made to ensure that the new capability
continues to work as expected?
