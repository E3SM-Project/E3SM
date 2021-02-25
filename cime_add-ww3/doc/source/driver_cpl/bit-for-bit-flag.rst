Bit-for-bit flag
============================

The driver namelist variable ``bfbflag`` provides the option of preserving bit-for-bit results on different coupler processor counts.
This flag has no impact on other components and their ability to generate bit-for-bit results on different pe counts. 
When this flag is set, all mappings become "X" types where the source data is rearranged to the destination processor and then local mapping is carried out.
The order of operations of this mapping is independent of the pe count or decomposition of the grids. 

The other feature that is changed by the ``bfbflag`` is the global sum diagnostics.

- When  ``bfbflag`` is set to *.false.*, a partial sum is done on each processors and those partial sums are added together to form a global sum.  This is generally not order of operations independent for different pe counts or decompositions.
- When ``bfbflag`` is set to *.true.*, the global sums are computed by gathering the global field on the root processor and doing an ordered sum there.
