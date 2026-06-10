# EAMxx new campaign instructions.

## Ultimate goal
The ultimate goal of a series of campaigns that I will run here 
is to implement water isotope and tracer capability in EAMxx, 
following conventions used in existing Earth system model implementation.
This campaign spec will configure some of the base-level infrastructure
needed after discussing with EAMxx developers.

## Campaign outline
1. A new process or processes need to be defined and added to the EAMxx DAG
that represent where all of the functions related to water isotopes/tracers 
will be called. I'd prefer one to be a subset of the other (e.g., water 
isotopes are a special case of water tracer that have undergo water 
isotope fractionation physics.) Add the process to the dag, and add a folder to 
src/physics to hold the implementation files.
2. After discussion with development team, we decided to implement specific
water tracer arrays using the field infrastructure in EAMxx. These arrays should be
something like qv_trace, qi_trace, (and for isotopes, qv_iso, qi_iso, ...) etc., 
they should have dimensions (col, cmp, lev), where the size of the cmp dimension can be specified at runtime
with an input yaml file. In this task, the spec design should fan out subagents to
to understand how to achieve different tasks along the way (e.g., how to specify yaml,
how to set up new fields, etc.)
3. The surface water fluxes into EAMxx will also need to be extended in a similar 
way to the atmospheric water arrays. 