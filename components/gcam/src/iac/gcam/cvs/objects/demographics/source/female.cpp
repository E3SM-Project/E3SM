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
* \file female.cpp
* \ingroup Objects-SGM
* \brief Female class source file.
* \author Sonny Kim
*/

#include "util/base/include/definitions.h"
#include <xercesc/dom/DOMNode.hpp>
#include "demographics/include/female.h"
#include "util/base/include/xml_helper.h"
#include "util/base/include/ivisitor.h"

using namespace std;
using namespace xercesc;
// static initialize
const string Female::XML_NAME = "female";

//! Default constructor
Female::Female() : Gender(){
    mFertilityRate = 0;
    mMaleBirth = 0;
    mFemaleBirth = 0;
    mMaleBirthFrac = 0;
}

//! Parse xml file for data
bool Female::XMLDerivedClassParse( const string &nodeName, const DOMNode* curr ) {
    if ( nodeName == "fertilityRate" ) {
        mFertilityRate = XMLHelper<double>::getValue( curr );
    }
    else if ( nodeName == "maleBirthFrac" ) {
        mMaleBirthFrac = XMLHelper<double>::getValue( curr );
    }
    else {
        return false;
    }
    return true;
}

//! For derived classes to output XML data
void Female::toInputXMLDerived( ostream& out, Tabs* tabs ) const {
    XMLWriteElementCheckDefault( mFertilityRate, "fertilityRate", out, tabs );
    XMLWriteElementCheckDefault( mMaleBirthFrac, "maleBirthFrac", out, tabs );
    XMLWriteElementCheckDefault( mFemaleBirth, "femaleBirth", out, tabs );
    XMLWriteElementCheckDefault( mMaleBirth, "maleBirth", out, tabs );
}

//! Output debug info for derived class
void Female::toDebugXMLDerived( ostream& out, Tabs* tabs ) const {
    XMLWriteElement( mFertilityRate, "fertilityRate", out, tabs );
    XMLWriteElement( mMaleBirthFrac, "maleBirthFrac", out, tabs );
    XMLWriteElement( mFemaleBirth, "femaleBirth", out, tabs );
    XMLWriteElement( mMaleBirth, "maleBirth", out, tabs );
}

// calculate male births
double Female::calcMaleBirth() {
    assert( mPopulation != -1 );
    mMaleBirth = mPopulation * mFertilityRate * mMaleBirthFrac;
    return mMaleBirth;
}

// calculate female births
double Female::calcFemaleBirth() {
    assert( mPopulation != -1 );
    mFemaleBirth = mPopulation * mFertilityRate * ( 1 - mMaleBirthFrac );
    return mFemaleBirth;
}

/*! \brief Get the XML node name for output to XML.
*
* This public function accesses the private constant string, XML_NAME.
* This way the tag is always consistent for both read-in and output and can be easily changed.
* This function may be virtual to be overridden by derived class pointers.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME.
*/
const std::string& Female::getXMLName() const {
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
const std::string& Female::getXMLNameStatic(){
    return XML_NAME;
}

/*! \brief Update a visitor with information about a Female.
* \param aVisitor Visitor to update.
* \param aPeriod Period for which to update.
*/
void Female::accept( IVisitor* aVisitor, const int aPeriod ) const {
    aVisitor->startVisitFemale( this, aPeriod );
    // Update the parent class.
    Gender::accept( aVisitor, aPeriod );
    aVisitor->endVisitFemale( this, aPeriod );
}
