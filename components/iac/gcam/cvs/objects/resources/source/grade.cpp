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
* \file grade.cpp
* \ingroup Objects
* \brief Grade class source file.
* \author Sonny Kim
*/

#include "util/base/include/definitions.h"
#include <string>
#include <iostream>
#include <cassert>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>

#include "containers/include/scenario.h"
#include "resources/include/grade.h"
#include "util/base/include/model_time.h"
#include "util/base/include/xml_helper.h"
#include "containers/include/info_factory.h"
#include "containers/include/iinfo.h"
#include "util/base/include/ivisitor.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;


//! Default constructor
Grade::Grade():
mAvailable( 0.0 ),
mExtractCost( 0.0 ),
mTotalCost( scenario->getModeltime()->getmaxper(), 0.0 )
{
}

//! Initialize data members from XML.
void Grade::XMLParse( const DOMNode* tempNode ) {
    /*! \pre assume we are passed a valid node. */
    assert( tempNode );

    // get the name attribute.
    mName = XMLHelper<string>::getAttr( tempNode, "name" );
    DOMNodeList* tempNodeLst = tempNode->getChildNodes();

    for( unsigned int i = 0; i < tempNodeLst->getLength(); ++i ) {
        DOMNode* tNode = tempNodeLst->item( i );
        string tNodeName = XMLHelper<string>::safeTranscode( tNode->getNodeName() );

        if( tNodeName == "#text" ) {
            continue;
        }

        else if( tNodeName == "available" ){
            mAvailable = XMLHelper<double>::getValue( tNode );
        }
        else if( tNodeName == "extractioncost" ){
            mExtractCost = XMLHelper<double>::getValue( tNode );
        }
        else {
            cout << "Unrecognized text string: " << tNodeName << " found while parsing grade." << endl;
        }
    }
}

/*! \brief Complete the initialization
*
* This routine is only called once per model run
*
* \author Josh Lurz, Sonny Kim
* \warning markets are not necessarily set when completeInit is called
*/
void Grade::completeInit( const IInfo* aSubresourceInfo ) {
    mGradeInfo.reset( InfoFactory::constructInfo( aSubresourceInfo, mName ) ); 
}


//! Write data members to data stream in XML format for replicating input file.
void Grade::toInputXML( ostream& out, Tabs* tabs ) const {
    XMLWriteOpeningTag( getXMLName(), out, tabs, mName );
    XMLWriteElementCheckDefault( mAvailable, "available", out, tabs, 0.0 );
    XMLWriteElementCheckDefault( mExtractCost, "extractioncost", out, tabs, 0.0 );
    XMLWriteClosingTag( getXMLName(), out, tabs );
}

//! Write data members to debugging data stream in XML format.
void Grade::toDebugXML( const int period, ostream& out, Tabs* tabs ) const {
    
    XMLWriteOpeningTag( getXMLName(), out, tabs, mName );
    
    XMLWriteElement( mAvailable, "available", out, tabs );
    XMLWriteElement( mExtractCost, "extractioncost", out, tabs );
    XMLWriteElement( mTotalCost[period], "totalcost", out, tabs );
    
    XMLWriteClosingTag( getXMLName(), out, tabs );
}

/*! \brief Get the XML node name for output to XML.
*
* This public function accesses the private constant string, XML_NAME.
* This way the tag is always consistent for both read-in and output and can be easily changed.
* This function may be virtual to be overridden by derived class pointers.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME.
*/
const std::string& Grade::getXMLName() const {
    return getXMLNameStatic();
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
const std::string& Grade::getXMLNameStatic() {
    const static string XML_NAME = "grade";
    return XML_NAME;
}

//! Total cost of each grade.
void Grade::calcCost( const double aTax, const double aCumTechChange, const double aEnvironCost, const int aPeriod ) {
    mTotalCost[ aPeriod ] = ( mExtractCost + aEnvironCost ) / aCumTechChange + aTax;
}

//! Return available amount in each Grade.
double Grade::getAvail() const {
    return mAvailable;
}

//! Return the total cost.
double Grade::getCost( const int aPeriod ) const {
    return mTotalCost[ aPeriod ];
}

//! Return the extraction cost.
double Grade::getExtCost() const {
    return mExtractCost;
}

//! Get the name.
const string& Grade::getName() const {
    return mName;
}

/*! \brief Perform any initializations needed before each period.
* \details Any initializations or calcuations that only need to be done once per
*          period(instead of every iteration) should be placed in this function.
* \author Sonny Kim
* \param aRegionName Region name.
* \param aResourceName Resource name.
* \param aPeriod Model period
*/
void Grade::initCalc( const string& aRegionName, const string& aResourceName, const int aPeriod ) {
    // do nothing
}

/*! \brief Perform any initializations needed after each period.
* \details Any initializations or calcuations that only need to be done once per
*          period after the model is solved should be placed in this function.
* \author Sonny Kim
* \param aRegionName Region name.
* \param aResourceName Resource name.
* \param aPeriod Model period
*/
void Grade::postCalc( const string& aRegionName, const string& aResourceName, const int aPeriod ) {
    // do nothing
}

/*! \brief Update a visitor for a SubResource.
* \param aVisitor Visitor to update.
* \param aPeriod Period to update.
*/
void Grade::accept( IVisitor* aVisitor, const int aPeriod ) const {
    aVisitor->startVisitGrade( this, aPeriod );
    aVisitor->endVisitGrade( this, aPeriod );
}
