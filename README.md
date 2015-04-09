MPAS-v2.1
====

The Model for Prediction Across Scales (MPAS) is a collaborative project for
developing atmosphere, ocean, and other earth-system simulation components for
use in climate, regional climate, and weather studies. The primary development
partners are the climate modeling group at Los Alamos National Laboratory
(COSIM) and the National Center for Atmospheric Research. Both primary
partners are responsible for the MPAS framework, operators, and tools common to
the applications; LANL has primary responsibility for the ocean model, and NCAR
has primary responsibility for the atmospheric model.

The MPAS framework facilitates the rapid development and prototyping of models
by providing infrastructure typically required by model developers, including
high-level data types, communication routines, and I/O routines. By using MPAS,
developers can leverage pre-existing code and focus more on development of
their model.

BUILDING
========

This README is provided as a brief introduction to the MPAS framework. It does
not provide details about each specific model, nor does it provide building
instructions.

For information about building and running each core, please refer to each
core's user's guide, which can be found at the following web sites:

[MPAS-Atmosphere](http://mpas-dev.github.io/atmosphere/atmosphere_download.html)

[MPAS-Land Ice](http://mpas-dev.github.io/land_ice/download.html)

[MPAS-Ocean](http://mpas-dev.github.io/ocean/releases.html)

Code Layout
----------

Within the MPAS repository, code is laid out as follows. Sub-directories are
only described below the src directory.

	MPAS
	└── src
	    ├── registry -- Code for building Registry.xml parser (Shared)
	    ├── driver -- Main driver for MPAS in stand-alone mode (Shared)
	    ├── external -- External software for MPAS (Shared)
	    ├── framework -- MPAS Framework (Includes DDT Descriptions, and shared routines. Shared)
	    ├── operators -- MPAS Opeartors (Includes Operators for MPAS meshes. Shared)
	    ├── inc -- Empty directory for include files that Registry generates (Shared)
	    └── core_* -- Individual model cores.

Model cores are typically developed independently. For information about
building and running a particular core, please refer to that core's user's
guide.
