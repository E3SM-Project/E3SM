#ifndef _DEFINITIONS_H_
#define _DEFINITIONS_H_
#if defined(_MSC_VER)
#pragma once
#endif

/*
 * LEGAL NOTICE
 * This computer software was prepared by Battelle Memorial Institute,
 * hereinafter the Contractor, under Contract No. DE-AC05-76RL0 1830
 * with the Department of Energy (DOE). NEITHER THE GOVERNMENT NOR THE
 * CONTRACTOR MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
 * LIABILITY FOR THE USE OF THIS SOFTWARE. This notice including this
 * sentence must appear on any copies of this computer software.
 * 
 * EXPORT CONTROL
 * User agrees that the Software will not be shipped, transferred or
 * exported into any country or used in any manner prohibited by the
 * United States Export Administration Act or any other applicable
 * export laws, restrictions or regulations (collectively the "Export Laws").
 * Export of the Software may require some form of license or other
 * authority from the U.S. Government, and failure to obtain such
 * export control license may result in criminal liability under
 * U.S. laws. In addition, if the Software is identified as export controlled
 * items under the Export Laws, User represents and warrants that User
 * is not a citizen, or otherwise located within, an embargoed nation
 * (including without limitation Iran, Syria, Sudan, Cuba, and North Korea)
 *     and that User is not otherwise prohibited
 * under the Export Laws from receiving the Software.
 * 
 * Copyright 2011 Battelle Memorial Institute.  All Rights Reserved.
 * Distributed as open-source under the terms of the Educational Community 
 * License version 2.0 (ECL 2.0). http://www.opensource.org/licenses/ecl2.php
 * 
 * For further details, see: http://www.globalchange.umd.edu/models/gcam/
 * 
 */


/*!
* \file definitions.h
* \ingroup Objects
* \brief A set of standard definitions which should be included in all files.
* \details This is a set of definitions, used mainly to work around platform and
*          compiler specific bugs. The intention is to hide many of the hacks
*          used to avoid compiler bugs.
* \author Josh Lurz
*/

#include <vector>
#include <map>
#include <list>
#include <string>
#include <iostream>

// Configuration constants.

//! A flag which tells whether to attempt linking of Fortran portions.
#define __HAVE_FORTRAN__ 0

#define __ACCESS_DB_OVERRIDE__ 0

//! A flag which turns on or off compilation of database code. Database
// cannot be compiled on VC8.
#if( defined(_MSC_VER) && _MSC_VER < 1400 && !__ACCESS_DB_OVERRIDE__)
#define __HAVE_DB__ 1
#else
#define __HAVE_DB__ 0
#endif

//! A flag which turns on or off the compilation of the XML database code.
#define __USE_XML_DB__ 1

// This allows for memory leak debugging.
#if defined(_MSC_VER)
#   ifdef _DEBUG
// usually the following two lines are defined in Microsoft's generated stdafx.h
#       define VC_EXTRALEAN // do not include rarely used parts
// extra definition for check whether all needed headers are included
#       undef SEARCH_MEMORY_LEAKS_ENABLED
#       define SEARCH_MEMORY_LEAKS_ENABLED
#   endif   // _DEBUG
#endif  // _MSC_VER

// Remove the _stdcall needed for WIN32 from externs
#if !defined(WIN32) && !defined(_stdcall)
#define _stdcall
#endif

#endif // _DEFINITIONS_H_
