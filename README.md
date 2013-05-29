MPAS
====

MPAS is a collaborative project for the rapid development and prototyping of
dynamical cores. A shared framework provides infrastructure typically required
by model developers, including communication routines, and I/O routines. By
using MPAS, developers can leverage pre-existing code and focus more on
development of their model.

Code Layout
----------

Within the MPAS repository code is laid out as follows. Sub-directories are only described below the src directory.

	MPAS
	├── graphics
	│   ├── dx -- Graphics for OpenDX
	│   ├── matlab -- Graphicx for MATLAB
	│   └── ncl -- Graphics for NCAR Command Language
	└── src
	    ├── registry -- Code for building Registry.xml parser (Shared)
	    ├── driver -- Main driver for MPAS in stand-alone mode (Shared)
	    ├── external -- External software for MPAS (Shared)
	    ├── framework -- MPAS Framework (Includes DDT Descriptions, and shared routines. Shared)
	    ├── operators -- MPAS Opeartors (Includes Operators for MPAS meshes. Shared)
	    ├── inc -- Empty directory for include files that Registry generates (Shared)
	    └── core_* -- Individual dynamical cores. (Private)


Dynamical cores are private and typically developed independently. Each core is stored in it's own directory under src, with an abbreviated name. For example:

src/core_sw houses the shallow water core.

For information about building and running each core, please refer to each core's user's guide.
