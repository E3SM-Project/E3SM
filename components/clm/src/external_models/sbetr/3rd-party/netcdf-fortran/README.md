Unidata NetCDF Fortran Library
========================================

[![Build Status](https://travis-ci.org/Unidata/netcdf-fortran.svg?branch=master)](https://travis-ci.org/Unidata/netcdf-fortran)

Overview 
--------


The Unidata network Common Data Form (netCDF) is an interface for scientific data access and a set of freely-distributed software libraries that provide an implementation of the interface.  The netCDF library also defines a machine-independent format for representing scientific data.  Together, the interface, library, and format support the creation, access, and sharing of scientific data.  This software package provides Fortran application interfaces for accessing netCDF data.  It depends on the netCDF C library, which must be installed first.  Other separate software distributions available from Unidata provide C Java, C++, and Python.  They have been tested on various common platforms.

NetCDF files are self-describing, network-transparent, directly accessible, and extendible.  `Self-describing` means that a netCDF file includes information about the data it contains.  `Network-transparent` means that a netCDF file is represented in a form that can be accessed by computers with different ways of storing integers, characters, and floating-point numbers.  `Direct-access` means that a small subset of a large dataset may be accessed efficiently, without first reading through all the preceding data.  `Extendible` means that data can be appended to a netCDF dataset without copying it or redefining its structure.

NetCDF is useful for supporting access to diverse kinds of scientific data in heterogeneous networking environments and for writing application software that does not depend on application-specific formats.  For information about a variety of analysis and display packages that have been developed to analyze and display data in netCDF form, see 

* http://www.unidata.ucar.edu/netcdf/software.html

For more information about netCDF, see the netCDF Web page at

* http://www.unidata.ucar.edu/netcdf/

Getting NetCDF
--------------

You can obtain a copy of the latest released versions of netCDF software from

* http://github.com/Unidata/netcdf-c 
* http://github.com/Unidata/netcdf-fortran
* http://github.com/Unidata/netcdf-cxx4
* http://github.com/Unidata/netcdf4-python

Copyright and licensing information can be found here, as well as in the COPYRIGHT file accompanying the software

* http://www.unidata.ucar.edu/software/netcdf/copyright.html

To install this package, please see the file INSTALL in the distribution, or the (usually more up-to-date) document:

* http://www.unidata.ucar.edu/netcdf/docs/netcdf-fortran-install.html

Additional Documentation
------------------------

The netCDF-3 C and FORTRAN-77 interfaces are documented in man(1) pages at

* http://www.unidata.ucar.edu/netcdf/docs/netcdf-man-3.html
* http://www.unidata.ucar.edu/netcdf/docs/netcdf-man-3f.html 

User's Guides are also available in several forms from the same location.

Mailing List
------------

A mailing list, netcdfgroup@unidata.ucar.edu, exists for discussion of the netCDF interface and announcements about netCDF bugs, fixes, and enhancements.  For information about how to subscribe, see the URL

* http://www.unidata.ucar.edu/netcdf/mailing-lists.html

We appreciate feedback from users of this package.  Please send comments, suggestions, and bug reports to <support-netcdf@unidata.ucar.edu>.  Please identify the version of the package.
