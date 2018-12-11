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
 * \file land_allocator.cpp
 * \ingroup Objects
 * \brief LandAllocator class source file.
 * \author James Blackwood, Kate Calvin
 */

#include "util/base/include/definitions.h"
#include "util/base/include/xml_helper.h"

#include "land_allocator/include/land_allocator.h"
#include "containers/include/scenario.h"
#include "containers/include/iinfo.h"
#include "util/base/include/model_time.h"
#include "ccarbon_model/include/carbon_model_utils.h"
#include "util/base/include/configuration.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

/*!
 * \brief Constructor.
 * \author James Blackwood
 */
LandAllocator::LandAllocator()
: LandNode( 0 ),
  mCarbonPriceIncreaseRate( 0.0 ),
  mSoilTimeScale( CarbonModelUtils::getSoilTimeScale() )
{
}

//! Destructor
LandAllocator::~LandAllocator() {
}

const string& LandAllocator::getXMLName() const {
    return getXMLNameStatic();
}

/*! \brief Get the XML node name in static form for comparison when parsing XML.
*
* This public function accesses the private constant string, XML_NAME.
* This way the tag is always consistent for both read-in and output and can be easily changed.
* The "==" operator that is used when parsing, required this second function to return static.
* \note A function cannot be static and virtual.
* \author James Blackwood
* \return The constant XML_NAME as a static.
*/
const string& LandAllocator::getXMLNameStatic() {
    const static string XML_NAME = "LandAllocatorRoot";    // original XML tag
    return XML_NAME;
}

bool LandAllocator::XMLParse( const DOMNode* aNode ){
    // Call the XML parse.
    return LandNode::XMLParse( aNode );
}

bool LandAllocator::XMLDerivedClassParse( const string& aNodeName, const DOMNode* aCurr ){
    if( aNodeName == "landAllocation" ){
        XMLHelper<double>::insertValueIntoVector( aCurr, mLandAllocation, scenario->getModeltime() );
    }
    else if( aNodeName == "carbonPriceIncreaseRate" ){
        XMLHelper<double>::insertValueIntoVector( aCurr, mCarbonPriceIncreaseRate, scenario->getModeltime() );
    }
    else if( aNodeName == "soilTimeScale" ){
        mSoilTimeScale = XMLHelper<int>::getValue( aCurr );
    }
    else {
        return false;
    }
    return true;
}

void LandAllocator::toDebugXML( const int aPeriod, std::ostream& aOut, Tabs* aTabs ) const {
    // Call the node toDebugXML
    ALandAllocatorItem::toDebugXML( aPeriod, aOut, aTabs );  
}

void LandAllocator::toInputXML( std::ostream& aOut, Tabs* aTabs ) const {
    // Call the node toInputXML
    LandNode::toInputXML( aOut, aTabs );
}

void LandAllocator::toInputXMLDerived( ostream& aOut, Tabs* aTabs ) const {
    XMLWriteVector( mLandAllocation, "landAllocation", aOut, aTabs, scenario->getModeltime() );
    XMLWriteVector( mCarbonPriceIncreaseRate, "carbonPriceIncreaseRate", aOut, aTabs, scenario->getModeltime() );
    XMLWriteElement( mSoilTimeScale, "soilTimeScale", aOut, aTabs );
}

void LandAllocator::initCalc( const string& aRegionName, const int aPeriod )
{
    // In calibration periods, check land area and set calibration values
    if ( aPeriod <= scenario->getModeltime()->getFinalCalibrationPeriod() ) {
        checkLandArea( aRegionName, aPeriod);
        calibrateLandAllocator( aRegionName, aPeriod );
    } 
    else {
        /* In non-calibration periods, adjust the profit scalers to account for new technologies. 
           If a new technology is added to a node, it will increase the node's land share because 
           it represents another "throw of the dice".  The adjustment procedure modifies the 
           profit scalers to correct for this problem.  The adjustment ensures that if a new 
           technology that is identical to an existing technology is added to a node, then the node
           will be allocated the same amount of land as when the new technology didn't exist.
        */
        adjustProfitScalers( aRegionName, aPeriod );
    }

    // Call land node's initCalc
    LandNode::initCalc( aRegionName, aPeriod );
    
    // Ensure that carbon price increase rate is positive
    if ( mCarbonPriceIncreaseRate[ aPeriod ] < 0 ) {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "Carbon price increase rate is negative in region "
                << aRegionName << endl;
        exit( -1 );
    }
    setCarbonPriceIncreaseRate( mCarbonPriceIncreaseRate[ aPeriod ] , aPeriod );
}

void LandAllocator::completeInit( const string& aRegionName, 
                                      const IInfo* aRegionInfo )
{
    // Call generic node method (since LandAllocator is just a specialized node)
    LandNode::completeInit( aRegionName, aRegionInfo );

    // Ensure that soil time scale is positive
    if ( mSoilTimeScale < 0 ) {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "Soil time scale is negative in region "
                << aRegionName << endl;
        exit( -1 );
    }

    // Set the soil time scale
    setSoilTimeScale( mSoilTimeScale );
}


/*!
 * \brief Ensures consistency in land area
 * \details Checks if the sum of land leafs equals the total land allocation
 * \param aRegionName Region name.
 * \param aPeriod model period.
 */
void LandAllocator::checkLandArea( const string& aRegionName, const int aPeriod ) {

    // Check to make sure that the sum of the area of the leafs is equal
    // to the total land area read in.
    double excessLand = getCalLandAllocation( eManaged, aPeriod )
        + getCalLandAllocation( eUnmanaged, aPeriod )
        - mLandAllocation[ aPeriod ];

    // If the difference between the total and the sum of the leafs
    // is too large, print an error message and exit. 
    double fractionDiff = excessLand / mLandAllocation[ aPeriod ];
    if ( abs( fractionDiff ) > 0.01 ) {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "The sum of land areas in region " << aRegionName 
            << " exceeds the land allocation by " 
            << excessLand << " (" << fractionDiff * 100 << "%) "
            << " in period " << aPeriod
            << endl;
		if ( abs( fractionDiff ) > 0.05 ) {
			exit( -1 );
		}
    }
}


/*!
 * \brief Calibrate the land allocator.
 * \details Sets initial land shares, sets unmanaged land profit rates,
 *          calculate calibration prfot rates, then calculates share profit scalers
 *          that will be used for scaling the profit rates for future year sharing
 * \param aRegionName Region name.
 * \param aPeriod model period.
 * \author Marshall Wise and Katherine V Calvin
 */
void LandAllocator::calibrateLandAllocator( const string& aRegionName, const int aPeriod ) {
    
/*  Step 1. Calculate and set initial land shares based on read in data for a 
    calibration period. */
    
    setInitShares( aRegionName,
                   0, // No land allocation above this node.
                   aPeriod );


/* Step 2. Set the profit rate of unmanaged land leafs equal to the read in land price
   (which is also the marginal profit rate) for that region or land node subregion. This
   is the only way the unmanaged land will have a profit rate at this point. It is implied
   by the read in price of land.  */

    setUnmanagedLandProfitRate( aRegionName, mUnManagedLandValue, aPeriod );

/* For these steps to work, the profit rates of managed land leaves will have been computed before
   this method call (i.e., calibrateLandAllocator) in the initCalc() of the agTechnology Class
   and set in the leafs based on read in calibration prices, yields, and non-land variable 
   costs. These might also be called "observed profits" since they are based on this simple 
   calculation. Also, the node profits do not need to be calculated here as they are not a 
   part of the calibration.  All the info needed is in the leaves. */


/* Step 3. Calculate the profit rate implied by the shares in the calibration data. These rates 
   are what the profit rates would have to be based on the actual shares, the logit exponent, and
   the average profit of the containing node. These are equivalent to what was called "intrinsic
   rates" in the 2008 version of the code based on Sands and Leimbech. */

    calculateCalibrationProfitRate( aRegionName, mUnManagedLandValue, 
                                    0, // logit exponent is zero at this level
                                    aPeriod );

/* Step 4. Calculate profit scalers. Because the calibration profit rate computed in Step 4
   will most likely differ from the profit rate computed using the yield times price - cost, a
   scaler or multiplier is solved for which makes the profit equal to the calibration profit. In
   future periods, this scaler is then applied to the computed profit for use in the sharing
   and node profit equations.  It is analagous to the share weight calibration approach. Also, it
   will work the same for unmanaged land leafs with the land price being used as the profit.
   
   All of the calibration is captured in the leaves, so the share profit scalers for nodes are
   set equal to 1.  */

    calculateProfitScalers( aRegionName, aPeriod );
}

/*!
 * \brief Calculates share profit scalers
 * \param aRegionName Region name.
 * \param aPeriod model period.
 */
void LandAllocator::calculateProfitScalers( const string& aRegionName, 
                                          const int aPeriod ) {
    LandNode::calculateProfitScalers( aRegionName, aPeriod );

    // Land allocator gets a shareweight of 1
    mProfitScaler[ aPeriod ] = 1.0;
}

/*!
 * \brief Adjusts share profit scalers
 * \param aRegionName Region name.
 * \param aPeriod model period.
 */
void LandAllocator::adjustProfitScalers( const string& aRegionName, 
                                          const int aPeriod ) {
    LandNode::adjustProfitScalers( aRegionName, aPeriod );
}

double LandAllocator::getLandAllocation( const string& aProductName,
                                         const int aPeriod ) const
{
    const ALandAllocatorItem* node = findChild( aProductName, eLeaf );

    if( node ){
        return node->getLandAllocation( aProductName, aPeriod );
    }
    return 0;
}

void LandAllocator::setCarbonPriceIncreaseRate( const double aCarbonPriceIncreaseRate, const int aPeriod )
{
    for ( unsigned int i = 0; i < mChildren.size(); i++ ) {
        mChildren[ i ]->setCarbonPriceIncreaseRate( aCarbonPriceIncreaseRate, aPeriod );
    }
}

/*!
* \brief Set the number of years needed to for soil carbons emissions/uptake
* \details This method sets the soil time scale into the carbon calculator
*          for each land leaf.
* \param aTimeScale soil time scale (in years)
* \author Kate Calvin
*/
void LandAllocator::setSoilTimeScale( const int aTimeScale ) {

    for ( unsigned int i = 0; i < mChildren.size(); i++ ) {
        mChildren[ i ]->setSoilTimeScale( aTimeScale );
    }

}

void LandAllocator::setProfitRate( const string& aRegionName,
                                   const string& aProductName,
                                   const double aProfitRate,
                                   const int aPeriod )
{
    ALandAllocatorItem* node = findChild( aProductName, eLeaf );
    if( node ){
        node->setProfitRate( aRegionName, aProductName,
                                aProfitRate, aPeriod );
    } 
}

void LandAllocator::setInitShares( const string& aRegionName,
                                       const double aLandAllocationAbove,
                                       const int aPeriod )
{
    // Calculating the shares
    for ( unsigned int i = 0; i < mChildren.size(); i++ ) {
        mChildren[ i ]->setInitShares( aRegionName,
                                       mLandAllocation[ aPeriod ],
                                       aPeriod );
    }

    // This is the root node so its share is 100%.
    mShare[ aPeriod ] = 1;
}

double LandAllocator::calcLandShares( const string& aRegionName,
                                          const double aLogitExpAbove,
                                          const int aPeriod ){

    // First set value of unmanaged land leaves
    setUnmanagedLandProfitRate( aRegionName, mUnManagedLandValue, aPeriod );

    LandNode::calcLandShares( aRegionName, aLogitExpAbove, aPeriod );
 
    // This is the root node so its share is 100%.
    mShare[ aPeriod ] = 1;
    return 1;
}

void LandAllocator::calcLandAllocation( const string& aRegionName,
                                            const double aLandAllocationAbove,
                                            const int aPeriod ){
    for ( unsigned int i = 0; i < mChildren.size(); ++i ){
        mChildren[ i ]->calcLandAllocation( aRegionName, mLandAllocation[ aPeriod ], aPeriod );
    }
}

void LandAllocator::calcLUCEmissions( const string& aRegionName, const int aPeriod,
                                      const int aEndYear )
{
    // Calculate emissions for all years in this model period.  Note that in
    // period 0 historical emissions are also calculated.
    for ( unsigned int i = 0; i < mChildren.size(); ++i ){
        mChildren[ i ]->calcLUCEmissions( aRegionName, aPeriod, aEndYear );
    } 
}


void LandAllocator::calcFinalLandAllocation( const string& aRegionName,
                                                 const int aPeriod ){
    // Calculate land shares
    calcLandShares( aRegionName,
                    0, // No logit exponent above the root.
                    aPeriod );

    // Calculate land allocation
    calcLandAllocation( aRegionName,
                        0, // No land allocation above the root.
                        aPeriod );

    // Calculate land-use change emissions but only to the end of this model
    // period for performance reasons.
    calcLUCEmissions( aRegionName,
                      aPeriod,
                      scenario->getModeltime()->getper_to_yr( aPeriod ) );
}

void LandAllocator::postCalc( const string& aRegionName, const int aPeriod ) {
    // Calculate land-use change emissions for the entire model time horizon.
    calcLUCEmissions( aRegionName, aPeriod, CarbonModelUtils::getEndYear() );
}

ALandAllocatorItem* LandAllocator::findProductLeaf( const string& aProductName ) {
    return findChild( aProductName, eLeaf );
}

void LandAllocator::accept( IVisitor* aVisitor, const int aPeriod ) const {
    LandNode::accept( aVisitor, aPeriod );
}
