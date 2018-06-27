#ifndef _IROUND_TRIPPABLE_H_
#define _IROUND_TRIPPABLE_H_
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
* \file iround_trippable.h
* \ingroup Objects
* \brief IRoundTrippable class header file.
* \author Josh Lurz
*/

class Tabs;
#include <iosfwd>

/*! \brief The IRoundTrippable interface allows the data necessary to create an object 
*          to be written to an XML file so that it can be read back in by a later scenario.
* \details Interface which specifies that the object must implement an interface
*          which writes out the data necessary to create this object so it can be read back 
*          in and used in a subsequent scenario. The interface specifies a single method, 
*          toInputXML which must output the object to the given stream. 
*/

class IRoundTrippable {
public:
	//! Virtual destructor so that instances of the interface may be deleted
    //! correctly through a pointer to the interface.
    inline virtual ~IRoundTrippable();

	/*! \brief Serialize the object to an output stream in an XML format.
	* \details Function which writes out all data members of an object which are
    *          necessary to duplicate a model run. This should not include
    *          internal state variables, only variables that were read-in or
    *          changed by calibration.
	* \param aOut Stream into which to write.
	* \param aTabs Object which controls formatting of the file.
    */
	virtual void toInputXML( std::ostream& aOut, Tabs* aTabs ) const = 0;
};

IRoundTrippable::~IRoundTrippable(){
}

#endif // _IROUND_TRIPPABLE_H_
