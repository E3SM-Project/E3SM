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
* \file subresource.cpp
* \ingroup Objects
* \brief SubResource class source file.
* \author Sonny Kim
*/

#include "util/base/include/definitions.h"
#include <vector>
#include <string>
#include <iostream>
#include <cassert>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>

#include "containers/include/scenario.h"
#include "util/base/include/model_time.h"
#include "resources/include/subresource.h"
#include "resources/include/grade.h"
#include "util/base/include/xml_helper.h"
#include "containers/include/info_factory.h"
#include "containers/include/iinfo.h"
#include "util/base/include/ivisitor.h"
#include "sectors/include/sector_utils.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;
// static initialize.
const string SubResource::XML_NAME = "subresource";
const double GDP_EXPANS_DEFAULT = 1;
const Value VALUE_DEFAULT = 0.0;

//! Default constructor.
SubResource::SubResource():
//mPriceAdder( scenario->getModeltime()->getmaxper() , 0.0 ),
mEffectivePrice( scenario->getModeltime()->getmaxper() , -1.0 ),
mCalProduction( scenario->getModeltime()->getmaxper() , -1.0 ),
mCumulativeTechChange( scenario->getModeltime()->getmaxper() , 1.0 ),
mAnnualProd( scenario->getModeltime()->getmaxper() , 0.0 ),
mAvailable( scenario->getModeltime()->getmaxper() , 0.0 ),
mCumulProd( scenario->getModeltime()->getmaxper() , 0.0 )
{
}

//! Destructor.
SubResource::~SubResource() {
    for ( vector<Grade*>::iterator outerIter = mGrade.begin(); outerIter != mGrade.end(); outerIter++ ) {
        delete *outerIter;
    }
}

//! Initialize member variables from xml data
void SubResource::XMLParse( const DOMNode* node ){  
    // make sure we were passed a valid node.
    assert( node );

    // get the name attribute.
    mName = XMLHelper<string>::getAttr( node, "name" );
    
    // get all child nodes.
    DOMNodeList* nodeList = node->getChildNodes();
    const Modeltime* modeltime = scenario->getModeltime();

    // loop through the child nodes.
    for( unsigned int i = 0; i < nodeList->getLength(); ++i ){
        DOMNode* curr = nodeList->item( i );
        string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );

        if( nodeName == "#text" ) {
            continue;
        }
        else if( nodeName == Grade::getXMLNameStatic() ){
            parseContainerNode( curr, mGrade, mGradeNameMap, new Grade );
        }
        else if( nodeName == "annualprod" ){
            XMLHelper<double>::insertValueIntoVector( curr, mAnnualProd, modeltime );
        }
        else if( nodeName == "techChange" ){
            XMLHelper<Value>::insertValueIntoVector( curr, mTechChange, modeltime );
        }
        else if( nodeName == "environCost" ){
            XMLHelper<Value>::insertValueIntoVector( curr, mEnvironCost, modeltime );
        }
        else if( nodeName == "severanceTax" ){
            XMLHelper<Value>::insertValueIntoVector( curr, mSeveranceTax, modeltime );
        }
        else if( nodeName == "cal-production" ){
            XMLHelper<double>::insertValueIntoVector( curr, mCalProduction, modeltime );
        }
        else if( nodeName == "price-adder" ){
            XMLHelper<Value>::insertValueIntoVector( curr, mPriceAdder, modeltime );
        }
        else if( !XMLDerivedClassParse( nodeName, curr ) ){
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Unrecognized text string: " << nodeName 
                << " found while parsing " << getXMLName() << "." << endl;
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
void SubResource::completeInit( const IInfo* aResourceInfo ) {
    mSubresourceInfo.reset( InfoFactory::constructInfo( aResourceInfo, mName ) ); 
    // update the available resource for period 0
    // this function must be called after all the grades have been parsed and nograde set
    updateAvailable( 0 );

    // call completeInit for grades
    for( vector<Grade*>::iterator gradeIter = mGrade.begin(); gradeIter != mGrade.end(); gradeIter++ ) {
        ( *gradeIter )->completeInit( mSubresourceInfo.get() );
    }

    // Interpolate and initialized exogenous resource variables created dynamically.
    const Modeltime* modeltime = scenario->getModeltime();

    // If unitialize, initialize following variables to null for final calibration period.
    const int FinalCalPer = modeltime->getFinalCalibrationPeriod();
    if( !mTechChange[ FinalCalPer ].isInited() ) {
        mTechChange[ FinalCalPer ] = VALUE_DEFAULT;
    }
    if( !mEnvironCost[ FinalCalPer ].isInited() ) {
        mEnvironCost[ FinalCalPer ] = VALUE_DEFAULT;
    }
    if( !mSeveranceTax[ FinalCalPer ].isInited() ) {
        mSeveranceTax[ FinalCalPer ] = VALUE_DEFAULT;
    }

    // decrement from terminal period to copy backward the technical change for missing periods
    for( int per = modeltime->getmaxper() - 1; per > modeltime->getFinalCalibrationPeriod(); --per ) {
        if( !mTechChange[ per ].isInited() ){
            // Copy backwards.
            if( per < modeltime->getmaxper() - 1 ){
                mTechChange[ per ] = mTechChange[ per + 1 ];
            }
            // For unitialized terminal period.
            // This may occur if running GCAM time periods beyond read-in dataset and
            // fillout is not used.
            else{
                mTechChange[ per ] = VALUE_DEFAULT;
            }
        }
    }

    SectorUtils::fillMissingPeriodVectorInterpolated( mEnvironCost );
    SectorUtils::fillMissingPeriodVectorInterpolated( mSeveranceTax );
}

/*! \brief Perform any initializations needed for each period.
* \details Any initializations or calculations that only need to be done once per
*          period(instead of every iteration) should be placed in this function.
* \author Sonny Kim
* \param aRegionName Region name.
* \param aResourceName Resource name.
* \param aPeriod Model aPeriod
*/
void SubResource::initCalc( const string& aRegionName, const string& aResourceName,
                           const int aPeriod )
{
    // call grade initializations
    for (unsigned int i = 0; i < mGrade.size(); i++) {
        mGrade[i]->initCalc( aRegionName, aResourceName, aPeriod );
    }
    // calculate total extraction cost for each grade
    for ( unsigned int gr=0; gr< mGrade.size(); gr++) {
        if ( aPeriod > 0) {
            const Modeltime* modeltime = scenario->getModeltime();
            mCumulativeTechChange[ aPeriod ] = mCumulativeTechChange[ aPeriod - 1 ] * 
                pow( ( 1.0 + mTechChange[ aPeriod ] ), modeltime->gettimestep( aPeriod ) );
        }
        // Determine cost
        mGrade[gr]->calcCost( mSeveranceTax[ aPeriod ], mCumulativeTechChange[ aPeriod ],
            mEnvironCost[ aPeriod ], aPeriod );
    }

    // Fill price added after it is calibrated.  This will interpolate to any
    // price adders read in the future or just copy forward if there is nothing
    // to interpolate to.
    if( aPeriod == scenario->getModeltime()->getFinalCalibrationPeriod() + 1 ) {
        SectorUtils::fillMissingPeriodVectorInterpolated( mPriceAdder );
    }
    
}

/*! \brief Perform any initializations needed after each period.
* \details Any initializations or calculations that only need to be done once
*          after each period(instead of every iteration) should be placed in
*          this function.
* \author Sonny Kim
* \param aRegionName Region name.
* \param aResourceName Resource name.
* \param period Model aPeriod
*/
void SubResource::postCalc( const string& aRegionName, const string& aResourceName, const int aPeriod ) {

    // Available is the total resource (stock) initialized in initCalc and
    // is the initial amount at the beginning of the period.
    // It does not subtract the amount used in that period.
    updateAvailable( aPeriod ); // reinitialize available amount
    if( aPeriod > 0 ) {
        mAvailable[ aPeriod ] -= mCumulProd[ aPeriod - 1 ];
        mAvailable[ aPeriod ] = max( mAvailable[ aPeriod ], 0.0 );
    }

    // call grade post calculations.
    for( unsigned int i = 0; i < mGrade.size(); i++ ) {
        mGrade[i]->postCalc( aRegionName, aResourceName, aPeriod );
    }
}

//! Blank definition so that don't have to define in derived classes if there is nothing to write out
void SubResource::toXMLforDerivedClass( ostream& out, Tabs* tabs ) const {   
}   

//! Write data members to data stream in XML format for replicating input file.
void SubResource::toInputXML( ostream& out, Tabs* tabs ) const {

    const Modeltime* modeltime = scenario->getModeltime();

    XMLWriteOpeningTag( getXMLName(), out, tabs, mName );

    // write the xml for the class members.
    const Value VALUE_DEFAULT = 0.0; // enables template function to recognize Value Class
    XMLWriteVector( mEnvironCost, "environCost", out, tabs, modeltime, VALUE_DEFAULT );
    XMLWriteVector( mSeveranceTax, "severanceTax", out, tabs, modeltime, VALUE_DEFAULT );
    XMLWriteVector( mTechChange, "techChange", out, tabs, modeltime, VALUE_DEFAULT );
    
    // for base year only
    XMLWriteElementCheckDefault(mAnnualProd[0],"annualprod",out, tabs, 0.0 , modeltime->getper_to_yr(0)); 

    XMLWriteVector( mCalProduction, "cal-production", out, tabs, modeltime, -1.0 );
    XMLWriteVector( mPriceAdder, "price-adder", out, tabs, modeltime, VALUE_DEFAULT  );
    // finished writing xml for the class members.

    // write out anything specific to the derived classes
    toXMLforDerivedClass( out, tabs );

    // write out the grade objects.
    for( vector<Grade*>::const_iterator i = mGrade.begin(); i != mGrade.end(); i++ ){ 
        ( *i )->toInputXML( out, tabs );
    }

    XMLWriteClosingTag( getXMLName(), out, tabs );
}

void SubResource::toDebugXML( const int period, ostream& out, Tabs* tabs ) const {

    XMLWriteOpeningTag( getXMLName(), out, tabs, mName );

    // Write out data for the period we are in from the vectors.
    XMLWriteElement( mAvailable[ period ], "available", out, tabs );
    XMLWriteElement( mAnnualProd[ period ], "annualprod", out, tabs );
    XMLWriteElement( mCumulProd[ period ], "cumulprod", out, tabs );
    XMLWriteElement( mTechChange[ period ], "techChange", out, tabs );
    XMLWriteElement( mEnvironCost[ period ], "environCost", out, tabs );
    XMLWriteElement( mSeveranceTax[ period ], "severanceTax", out, tabs );
    XMLWriteElement( mCalProduction[ period ], "cal-production", out, tabs );
    XMLWriteElement( mEffectivePrice[ period ], "effective-price", out, tabs );
    XMLWriteElement( mPriceAdder[ period ], "price-adder", out, tabs );

    // write out the grade objects.
    for( int i = 0; i < static_cast<int>( mGrade.size() ); i++ ){    
        mGrade[ i ]->toDebugXML( period, out, tabs );
    }

    // finished writing xml for the class members.

    // write out anything specific to the derived classes
    toXMLforDerivedClass( out, tabs );

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
const std::string& SubResource::getXMLName() const {
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
const std::string& SubResource::getXMLNameStatic() {
    return XML_NAME;
}

//! return SubResource name
string SubResource::getName() const {
    return mName;
}

void SubResource::cumulsupply( double aPrice, int aPeriod )
{  
    if ( aPeriod == 0 ) {
        mCumulProd[ aPeriod ] = 0.0;
        mAnnualProd[ aPeriod ] = mCalProduction[ aPeriod ];
        mPriceAdder[ aPeriod ] = 0.0;
    }

    // Always calculate the effective price
    mEffectivePrice[ aPeriod ] = aPrice + mPriceAdder[ aPeriod ];

    if ( aPeriod > 0 ) {
        // Case 1
        // if market price is less than cost of first grade, then zero cumulative 
        // production
        if ( mEffectivePrice[ aPeriod ] <= mGrade[0]->getCost( aPeriod )) {
            mCumulProd[ aPeriod ] = mCumulProd[ aPeriod - 1 ];
        }

        // Case 2
        // if market price is in between cost of first and last grade, then calculate 
        // cumulative production in between those grades
        if ( mEffectivePrice[ aPeriod ] > mGrade[0]->getCost( aPeriod ) && mEffectivePrice[ aPeriod ] <= mGrade[ mGrade.size() - 1 ]->getCost( aPeriod )) {
            mCumulProd[ aPeriod ] = 0;
            int i = 0;
            int iL = 0;
            int iU = 0;
            while ( mGrade[ i ]->getCost( aPeriod ) < mEffectivePrice[ aPeriod ] ) {
                iL=i; i++; iU=i;
            }
            // add subrsrcs up to the lower grade
            for ( i = 0; i <= iL; i++ ) {
                mCumulProd[ aPeriod ] += mGrade[i]->getAvail();
            }
            // price must reach upper grade cost to produce all of lower grade
            double slope = mGrade[iL]->getAvail()
                / ( mGrade[iU]->getCost( aPeriod ) - mGrade[iL]->getCost( aPeriod ) );
            mCumulProd[ aPeriod ] -= slope * ( mGrade[iU]->getCost( aPeriod ) - mEffectivePrice[ aPeriod ] );
        }

        // Case 3
        // if market price greater than the cost of the last grade, then
        // cumulative production is the amount in all grades
        if ( mEffectivePrice[ aPeriod ] > mGrade[ mGrade.size() - 1 ]->getCost( aPeriod ) ) {
            mCumulProd[ aPeriod ] = 0;
            for ( unsigned int i = 0; i < mGrade.size(); i++ ) {
                mCumulProd[ aPeriod ] += mGrade[i]->getAvail();
            }
        }
    }
}

double SubResource::getCumulProd( const int aPeriod ) const {
    return mCumulProd[ aPeriod ];
}

/*! Update the sub-resource availability for a period
* Resource depletion by grade is not calculated.
* This function only returns the maximum amount of resource
* available by grade.  
*
*/
void SubResource::updateAvailable( const int aPeriod ){
    mAvailable[ aPeriod ] = 0;
    for ( unsigned int i = 0; i < mGrade.size(); ++i ) {
        mAvailable[ aPeriod ] += mGrade[ i ]->getAvail();
    }
}

//! calculate annual supply
/*! Takes into account short-term capacity limits.
Note that cumulsupply() must be called before calling this function. */
void SubResource::annualsupply( int aPeriod, const GDP* aGdp, double aPrice, double aPrev_price ) {
    const Modeltime* modeltime = scenario->getModeltime();
    // For period 0 use initial annual supply. Cumulative production is 0 for period 0
    if ( aPeriod >= 1) {
        // Calculate the annual production given that the cumulative production
        // for the period is known. Cumulative production for the current period
        // is equal to the cumulative production of the previous period plus the
        // current annual production times the timestep.  This assumes that
        // the production in the current period is the average production for
        // the whole timestep.
        mAnnualProd[ aPeriod ] = ( mCumulProd[ aPeriod ] - mCumulProd[ aPeriod - 1 ] ) 
                                / modeltime->gettimestep( aPeriod );


        if( mAnnualProd[ aPeriod ] <= 0) {
            mCumulProd[ aPeriod ] = mCumulProd[ aPeriod - 1 ];
            mAnnualProd[ aPeriod ] = 0.0;
        } 

        // mAvailable is the total resource (stock) remaining
        mAvailable[ aPeriod ] = mAvailable[ aPeriod - 1 ] - ( mAnnualProd[ aPeriod ] * modeltime->gettimestep( aPeriod ) );
        mAvailable[ aPeriod ] = max( mAvailable[ aPeriod ], 0.0 );
    }
}

//! return annual production for period
double SubResource::getAnnualProd( int aPeriod ) const {
    return mAnnualProd[ aPeriod ];
}

/*! \brief Update an output container for a SubResource.
* \param aVisitor Output container to update.
* \param aPeriod Period to update.
*/
void SubResource::accept( IVisitor* aVisitor, const int aPeriod ) const {
    aVisitor->startVisitSubResource( this, aPeriod );

    // Update the output container for the subresources.
    for( unsigned int i = 0; i < mGrade.size(); ++i ){
        mGrade[ i ]->accept( aVisitor, aPeriod );
    }
    aVisitor->endVisitSubResource( this, aPeriod );
}

//! return available resource for period
double SubResource::getAvailable(int per) const {
    return mAvailable[per];
}

//! write SubResource output to database
void SubResource::dbOutput( const string &regname, const string& secname ){
    // function protocol
    void dboutput4(string var1name,string var2name,string var3name,string var4name,
        string uname,vector<double> dout);

    int m=0;
    const Modeltime* modeltime = scenario->getModeltime();
    const int maxper = modeltime->getmaxper();
    const string outputUnit = mSubresourceInfo->getString( "output-unit", true );
    const string priceUnit = mSubresourceInfo->getString( "price-unit", true );
    vector<double> temp(maxper);
    string tssname = mName; // tempory subsector name

    // function arguments are variable name, double array, db name, table name
    // the function writes all years
    // total subsector output
    //dboutput4(regname,"Primary Energy", "Production for " + secname,name,outputUnit,annualprod);
    //    dboutput4(regname,"Resource",secname,str,"EJ",available);
    dboutput4(regname,"Resource","Available "+secname,mName,outputUnit,mAvailable);
    dboutput4(regname,"Resource","CummProd "+secname,mName,outputUnit,mCumulProd);

    // do for all grades in the sector
    for ( unsigned int i=0;i< mGrade.size();i++) {
        string str = tssname + "_" + mGrade[i]->getName();
        // grade cost
        for (m=0;m<maxper;m++) {
            temp[m] = mGrade[i]->getCost(m);
        }
        dboutput4(regname,"Price",secname,str,priceUnit,temp);
        // grade extraction cost
        for (m=0;m<maxper;m++) {
            temp[m] = mGrade[i]->getExtCost();
        }
        dboutput4(regname,"Price ExtCost",secname,str,priceUnit,temp);
        // available resource for each grade
        for (m=0;m<maxper;m++) {
            temp[m] = mGrade[i]->getAvail();
        }
        dboutput4(regname,"Resource",secname,str,outputUnit,temp);
    }
}

//! write SubResource output to file
void SubResource::csvOutputFile( const string &regname, const string& sname) {
    const Modeltime* modeltime = scenario->getModeltime();
    // function protocol
    void fileoutput3( string var1name,string var2name,string var3name,
        string var4name,string var5name,string uname,vector<double> dout);

    const int maxper = modeltime->getmaxper();
    const string outputUnit = mSubresourceInfo->getString( "output-unit", true );
    vector<double> temp(maxper);

    // function arguments are variable name, double array, db name, table name
    // the function writes all years
    // total subsector output
    fileoutput3( regname,sname,mName," ","production",outputUnit,mAnnualProd);
    fileoutput3( regname,sname,mName," ","resource",outputUnit,mAvailable);

}


// ************************************************************
// Definitions for two of the derived classes below.
// Since these are very small changes, keep in same file for simplicity
// ************************************************************

//! Parses any input variables specific to this derived class
bool SubDepletableResource::XMLDerivedClassParse( const string& nodeName, const DOMNode* node ) {
    return false;
}

//! Parses any input variables specific to this derived class
bool SubFixedResource::XMLDerivedClassParse( const string& nodeName, const DOMNode* node ) {
    return false;
}
//! get variance
/*! do nothing here.  Applies to derived subrenewableresource
* \author Marshall Wise
*/
double SubResource::getVariance() const {
    return 0.0;
}

//! get resource capacity factor
/*! do nothing here.  Applies to derived subrenewableresource
* \author Marshall Wise
*/
double SubResource::getAverageCapacityFactor() const {
    return 0.0;
}
