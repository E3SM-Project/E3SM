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
* \file male.cpp
* \ingroup Objects-SGM
* \brief Male class source file.
* \author Sonny Kim
*/

#include "util/base/include/definitions.h"
#include <xercesc/dom/DOMNode.hpp>
#include "demographics/include/male.h"
#include "util/base/include/xml_helper.h"
#include "util/base/include/ivisitor.h"

using namespace std;
using namespace xercesc;
// static initialize
const string Male::XML_NAME = "male";

//! Parse xml file for data
bool Male::XMLDerivedClassParse( const string& nodeName, const DOMNode* curr ) {
    return false;
}

//! For derived classes to output XML data
void Male::toInputXMLDerived( ostream& out, Tabs* tabs ) const {
    // nothing
}

//! Output debug info for derived class
void Male::toDebugXMLDerived(ostream& out, Tabs* tabs ) const {
    // nothing
}

/*! \brief Get the XML node name for output to XML.
*
* This public function accesses the private constant string, XML_NAME.
* This way the tag is always consistent for both read-in and output and can be easily changed.
* This function may be virtual to be overridden by derived class pointers.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME.
*/
const std::string& Male::getXMLName() const {
    return XML_NAME;
}

/*! \brief Get the XML node name in static form for comparison when parsing XML.
*
* This public function accesses the private constant string, XML_NAME.
* This way the tag is always consistent for both read-in and output and can be easily changed.
* The "==" operator that is used when parsing, required this second function to return static.
* \note A function cannot be static and virtual.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME as a static.
*/
const std::string& Male::getXMLNameStatic(){
    return XML_NAME;
}

/*! \brief Update a visitor with information about a Male.
* \param aVisitor Visitor to update.
* \param aPeriod Period for which to update.
*/
void Male::accept( IVisitor* aVisitor, const int aPeriod ) const {
    aVisitor->startVisitMale( this, aPeriod );
    // Update the parent class.
    Gender::accept( aVisitor, aPeriod );
    aVisitor->endVisitMale( this, aPeriod );
}
