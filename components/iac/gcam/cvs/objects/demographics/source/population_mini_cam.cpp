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
* \file population_mini_cam.cpp
* \ingroup Objects
* \brief PopulationMiniCAM class source file.
* \author Sonny Kim
* \author Katherine Chung
*/

#include "util/base/include/definitions.h"
#include <string>
#include <map>
#include <cassert>
#include <vector>
#include <cmath>
#include <xercesc/dom/DOMNode.hpp>

#include "containers/include/scenario.h"
#include "demographics/include/population_mini_cam.h"
#include "util/base/include/model_time.h"
#include "util/base/include/xml_helper.h"
#include "util/base/include/configuration.h"
#include "util/base/include/ivisitor.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;
// static initialize.
const string PopulationMiniCAM::XML_NAME = "populationMiniCAM";

//! Default constructor.
PopulationMiniCAM::PopulationMiniCAM():
mFractionWorking( 1.0 ) {
}

Population* PopulationMiniCAM::cloneAndInterpolate( const int aNewYear, const Population* aNextPopulation ) const {
    /*!
     * \pre aNewYear is between the year of this population and aNextPopulation.
     */
    assert( this->mYear < aNewYear < aNextPopulation->getYear() );
    
    PopulationMiniCAM* newInterpolatedPop = new PopulationMiniCAM();
    newInterpolatedPop->mYear = aNewYear;
    newInterpolatedPop->mTotalPop = util::linearInterpolateY( aNewYear, this->mYear, aNextPopulation->getYear(),
                                                              this->mTotalPop, aNextPopulation->getTotal() );
    newInterpolatedPop->mIsParsed = false;
    return newInterpolatedPop;
}

//! parses the rest of PopuationMiniCAM xml object
bool PopulationMiniCAM::XMLDerivedClassParse( const string &nodeName, const xercesc::DOMNode* curr ){
    if ( nodeName == "totalPop" ) {
        mTotalPop = XMLHelper<double>::getValue( curr );
    }
    else if ( nodeName == "fraction-working" ) {
        mFractionWorking = XMLHelper<double>::getValue( curr );
    }
    else {
        return false;
    }
    return true;
}

//! returns total working age population (ages 15-65)
double PopulationMiniCAM::getWorkingAgePop() const { // ages 15-65
    // for population MiniCAM we just take the total and multiply
    // by an exogenous working age fraction which defaults to 1
    return mTotalPop * mFractionWorking;
}

//! Write out data members to XML output stream.
void PopulationMiniCAM::toInputXMLDerived( ostream& out, Tabs* tabs ) const {
    // note that mTotalPop is actually written out by Population
    XMLWriteElementCheckDefault( mFractionWorking, "fraction-working", out, tabs, 1.0 );
}

//! Write out XML for debugging purposes.
void PopulationMiniCAM::toDebugXMLDerived( ostream& out, Tabs* tabs ) const {
    XMLWriteElement( mFractionWorking, "fraction-working", out, tabs );
}

//! Complete the initialization.
void PopulationMiniCAM::completeInit( const vector<double>& femalePopFromPrev, const vector<double>& malePopFromPrev ){

}

//! initialize anything that won't change during the calcuation
void PopulationMiniCAM::initCalc(){

}

/*! \brief Get the XML node name for output to XML.
*
* This public function accesses the private constant string, XML_NAME.
* This way the tag is always consistent for both read-in and output and can be easily changed.
* This function may be virtual to be overriden by derived class pointers.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME.
*/
const string& PopulationMiniCAM::getXMLName() const {
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
const string& PopulationMiniCAM::getXMLNameStatic() {
    return XML_NAME;
}

void PopulationMiniCAM::csvSGMOutputFile( ostream& aFile, const int period ) const {
    Population::csvSGMOutputFile( aFile, period );
}

/*! \brief Update a visitor with information about a MiniCAM population.
* \param aVisitor Visitor to update.
* \param aPeriod Period for which to update.
*/
void PopulationMiniCAM::accept( IVisitor* aVisitor, const int aPeriod ) const {
    aVisitor->startVisitPopulationMiniCAM( this, aPeriod );
    // Call the parent class visit.
    Population::accept( aVisitor, aPeriod );
    aVisitor->endVisitPopulationMiniCAM( this, aPeriod );
}
