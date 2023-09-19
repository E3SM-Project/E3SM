# OMEGA src directory

This directory contains most of the source code for the
Ocean Model for E3SM Global Applications (OMEGA). These
are sorted into several subdirectories, including:

- **base:**
  The lowest level infrastructure for decomposition,
  message passing, basic data types, etc.

- **infra:**
  Other infrastructure shared by many routines, including
  time managers, configuration, metadata, etc.

- **ocn:**
   Source code for actual ocean tendencies, operators and
   physical parameterizations.

- **drivers:**
   Top-level drivers for both standalone and coupled execution

- **analysis:**
   Source code for ocean analyses
