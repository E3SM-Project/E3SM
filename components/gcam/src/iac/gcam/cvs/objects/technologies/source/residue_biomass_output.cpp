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

#include "util/base/include/definitions.h"
#include "containers/include/dependency_finder.h"
#include "containers/include/iinfo.h"
#include "containers/include/scenario.h"
#include "marketplace/include/marketplace.h"
#include "technologies/include/residue_biomass_output.h"
#include "util/base/include/ivisitor.h"
#include "util/base/include/xml_helper.h"
#include "util/base/include/TValidatorInfo.h"
#include "util/curves/include/point_set_curve.h"
#include "util/curves/include/curve.h"
#include "util/curves/include/explicit_point_set.h"
#include "util/curves/include/xy_data_point.h"
#include "land_allocator/include/aland_allocator_item.h"

#include <xercesc/dom/DOMNodeList.hpp>

#include <cstdio>

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

// static initialize.
const string ResidueBiomassOutput::XML_REPORTING_NAME = "output-residue-biomass";

ResidueBiomassOutput::ResidueBiomassOutput( const std::string& sectorName ) : parent(),
        mPhysicalOutputs( scenario->getModeltime()->getmaxper(), 0 ),
        mName( sectorName ),
        mCachedCO2Coef( ),
        mProductLeaf( 0 ),
        mHarvestIndex( 0 ),
        mErosCtrl( 0 ),
        mMassConversion( 0 ),
        mMassToEnergy( 0 ),
        mCostCurve( ),
        mWaterContent( 0 )

{
}

ResidueBiomassOutput::ResidueBiomassOutput( const ResidueBiomassOutput& other ) : parent(),
        mPhysicalOutputs( scenario->getModeltime()->getmaxper(), 0 ), // Do not copy results
        mName( other.mName ),
        mCachedCO2Coef( other.mCachedCO2Coef ),
        mProductLeaf( other.mProductLeaf ),
        mHarvestIndex( other.mHarvestIndex ),
        mErosCtrl( other.mErosCtrl ),
        mMassConversion( other.mMassConversion ),
        mMassToEnergy( other.mMassToEnergy ),
        mWaterContent( other.mWaterContent )
{
    mCostCurve.reset( other.mCostCurve->clone() );
}


ResidueBiomassOutput::~ResidueBiomassOutput( ) 
{
}

ResidueBiomassOutput& ResidueBiomassOutput::operator =( const ResidueBiomassOutput& other ){
    if ( &other != this ) {
        // Note mPhysicalOutputs is not copied as results should never be copied.
        mName = other.mName;
        mCachedCO2Coef = other.mCachedCO2Coef;
        mProductLeaf = other.mProductLeaf;
        mHarvestIndex = other.mHarvestIndex;
        mErosCtrl = other.mErosCtrl;
        mMassConversion = other.mMassConversion;
        mMassToEnergy = other.mMassToEnergy;
        mCostCurve.reset( other.mCostCurve->clone() );
        mWaterContent = other.mWaterContent;
    }
    return *this;
}


/*! \brief Get the XML name for reporting to XML file.
*
* This public function accesses the private constant string, XML_NAME. This way
* the tag is always consistent for reporting outputs and can be easily
* changed.
* \author Sonny Kim
* \return The constant XML_NAME.
*/
const string& ResidueBiomassOutput::getXMLReportingName() const{
    return XML_REPORTING_NAME;
}

void ResidueBiomassOutput::accept( IVisitor* aVisitor, const int aPeriod ) const
{
   if ( aVisitor )
   {
      aVisitor->startVisitOutput( this, aPeriod );
      aVisitor->endVisitOutput( this, aPeriod );
   }
}

/*!
 * \brief Calculates the crop or forest residue bioenergy production
 * \details Calculates the potential bioenergy production which
 *          is equal to the energy content of the primary output of 
 *          the crop or forest less the biomass needed to prevent
 *          erosion. The fraction of the potential bioenergy that is
 *          actually produced is determined by the bioenergy price.
 *          Returns the residue bioenergy produced at the current
 *          market price and primary output quantity
 * \param aPrimaryOutput Primary Production for the Parent Technology
 * \param aRegionName Region
 * \param aCaptureComponent Emissions Capture Component for the Parent Technology
 * \param aPeriod Model Period
 * \return Biomass Production
 */
IOutput::OutputList ResidueBiomassOutput::calcPhysicalOutput( const double aPrimaryOutput,
                                                              const std::string& aRegionName,
                                                              const ICaptureComponent* aCaptureComponent,
                                                              const int aPeriod ) const 
{
    // create the list to return and default it to zero
    OutputList outputList;
    outputList.push_back( make_pair( mName, 0 ) );

    // If primary output is less than or equal to zero, then residue biomass production
    // is zero.
    if ( aPrimaryOutput <= 0 ) {
        return outputList;
    }

    const Marketplace* marketplace = scenario->getMarketplace();
    double price = marketplace->getPrice( getName(), aRegionName, aPeriod, true );

    // If there is no market price, return
    if ( price == Marketplace::NO_MARKET_PRICE ) {
        return outputList;
    }

    // If there is an erosion control parameter, then we need
    // the land area to calculate the biomass needed to prevent
    // erosion. ( mErosCtrl is measured in tonnes/kHa ).
    double landArea = 0;
    if ( mErosCtrl > 0 ) {
        landArea = mProductLeaf->getLandAllocation( "", aPeriod );
    }
   
    // Compute the amount of crop produced in tonnes
    // Primary Output is in billion m^3 or in ECal
    // mMassConversion is in tonnes/billion m^3 or tonnes/ECal
    mCropMass = aPrimaryOutput * mMassConversion;

    // Compute the above-ground biomass that is not used for food
    // production. Harvest index is the fraction of total biomass produced for food
    // Thus, the amount of biomass available for bioenergy production
    // is food production * [ ( 1 / harvest index ) - 1 ]
    if ( mHarvestIndex > 0.0 ) {
        mResMass = mCropMass * ( std::pow( mHarvestIndex, double( -1 ) ) - double( 1 ) );
    }
    else {
        mResMass = 0.0;
    }

    // Compute the mass of residue that are required to sustain
    // agriculture in terms of soil nutrients and erosion 
    // mErosCtrl is in tonnes/kHa, landArea is in kHa
    // So, mMeanErosCtrl is in tonnes
    mMeanErosCtrl = landArea * mErosCtrl;

    // Compute the amount of residue biomass available to the energy sector
    // The amount available is the total biomass not produced for food
    // less the amount of biomass needed for erosion control
    // multiply by moisture content to get DRY residue available; residue from model is wet
    mResAvail = (1 - mWaterContent) * ( mResMass - mMeanErosCtrl);

    // If there is no biomass available for production after
    // erosion control demands are met, then return. Residue
    // biomass production is zero.
    if ( mResAvail <= 0 ) {
        return outputList;
    }

    // Compute energy content of the available biomass 
    // mResAvail is measured in tonnes
    // mMassToEnergy is measured in EJ/tonne
    // mMaxBioEnergySupply is measured in EJ
    mMaxBioEnergySupply = mResAvail * mMassToEnergy;

    // Compute the fraction of the total possible supply that is
    // produced at the current biomass price 
    mFractProduced = mCostCurve->getY( price );

    // Compute the quantity of a crop residue biomass produced
    double resEnergy = mMaxBioEnergySupply * mFractProduced;

    outputList.front().second = resEnergy;
    return outputList;
}

void ResidueBiomassOutput::completeInit( const std::string& aSectorName, DependencyFinder* aDependencyFinder,
                                         const IInfo* aTechInfo, const bool aIsTechOperating )
{
    #if !defined( _MSC_VER )
        using std::exit;
    #endif

    // If erosion control is positive, but a land allocator
    // has not been assigned, then print an error message
    // Erosion control equations need land allocations.
    if ( !mProductLeaf  && mErosCtrl > 0 ) {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "Invalid land allocator for " << getXMLNameStatic()
                << " in sector " << aSectorName << std::endl;
        exit( -1 );
    }

    // Validate input parameters
    typedef ObjECTS::TValidatorInfo<> validator_type;
    validator_type validator[] = {
                    validator_type( mCostCurve->getMaxY(), "max-harvest-fraction",
                                    mCostCurve->getMaxY() <= 1 ),
                    validator_type( mErosCtrl, "eros-ctrl", mErosCtrl >= 0 ),
                    validator_type( mHarvestIndex, "harvest-index", mHarvestIndex >= 0 ),
                    validator_type( mMassConversion, "mass-conversion", mMassConversion >= 0 ),
                    validator_type( mMassToEnergy, "mass-to-energy", mMassToEnergy >= 0 ) };
    unsigned short numParams = sizeof( validator ) / sizeof( validator[0] );
    std::string msg = ObjECTS::getInvalidNames( &validator[0], &validator[numParams] );

    // If input parameters are invalid, print an error message and exit
    if ( msg.length() ) {
        // Invalid input parameter
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "Invalid input parameter(s) to " << getXMLNameStatic()
                << " in sector " << aSectorName << ": " << msg << std::endl;
        exit( -1 );
    }
}

double ResidueBiomassOutput::getEmissionsPerOutput( const std::string& aGHGName, const int aPeriod ) const
{
    // Currently other GHGs do not use output emissions coefficients.
    assert( aGHGName == "CO2" );
    return mCachedCO2Coef;
}

double ResidueBiomassOutput::getPhysicalOutput( const int aPeriod ) const
{
    /*if ( !mPhysicalOutputs.size() ) {
        const Modeltime* modeltime = scenario->getModeltime();
        mPhysicalOutputs.resize( modeltime->getmaxper() );
    }*/

    return mPhysicalOutputs[ aPeriod ];
}

double ResidueBiomassOutput::getValue( const std::string& aRegionName,
                                       const ICaptureComponent* aCaptureComponent,
                                       const int aPeriod ) const
{
    // TODO: need to change code to make it residue per unit of crop,
    // rather than total residue
    return 0;
}

const std::string& ResidueBiomassOutput::getXMLNameStatic( void )
{
    static const std::string XMLName = "residue-biomass-production";
    return XMLName;
}

void ResidueBiomassOutput::initCalc( const std::string& aRegionName, const std::string& aSectorName,
                                     const int aPeriod )
{
    assert( scenario != 0 );
    const Marketplace* marketplace = scenario->getMarketplace();
    const IInfo* productInfo = marketplace->getMarketInfo( getName(), aRegionName, aPeriod, false );

    mCachedCO2Coef.set( productInfo ? productInfo->getDouble( "CO2Coef", false ) : 0 );
}

void ResidueBiomassOutput::postCalc( const std::string& aRegionName, const int aPeriod )
{
   // Not used
}

void ResidueBiomassOutput::scaleCoefficient( const double aScaler )
{
   // Not used
}

void ResidueBiomassOutput::sendLandAllocator( const ILandAllocator* aLandAllocator,
                                              const std::string& aName )
{
    // TODO: maybe the technology should just pass the product leaf
    mProductLeaf = const_cast<ILandAllocator*>( aLandAllocator )->findProductLeaf( aName );
}

void ResidueBiomassOutput::setPhysicalOutput( const double aPrimaryOutput, const std::string& aRegionName,
                                              ICaptureComponent* aCaptureComponent, const int aPeriod )
{
    /*
    if ( !mPhysicalOutputs.size() ) {
        const Modeltime* modeltime = scenario->getModeltime();
        mPhysicalOutputs.resize( modeltime->getmaxper() );
    }*/

    // Set the physical output for the specified period
    // calcPhysicalOutput only returns a single value.
    OutputList outputList = calcPhysicalOutput( aPrimaryOutput, aRegionName, aCaptureComponent, aPeriod );
    mPhysicalOutputs[ aPeriod ].set( outputList.front().second );

    // Add output to the supply
    Marketplace* marketplace = scenario->getMarketplace();
    marketplace->addToSupply( getName(), aRegionName, mPhysicalOutputs[ aPeriod ], aPeriod, true );
}

void ResidueBiomassOutput::toDebugXML( const int aPeriod, std::ostream& aOut, Tabs* aTabs ) const
{
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs, getName() );
    if ( value_vector_type::size_type( aPeriod ) < mPhysicalOutputs.size() ) {
        XMLWriteElement( mPhysicalOutputs[ aPeriod ], "output", aOut, aTabs );
    }
    XMLWriteElement( mErosCtrl, "eros-ctrl", aOut, aTabs );
    XMLWriteElement( mHarvestIndex, "harvest-index", aOut, aTabs );
    XMLWriteElement( mMassConversion, "mass-conversion", aOut, aTabs );
    XMLWriteElement( mMassToEnergy, "mass-to-energy", aOut, aTabs );
    XMLWriteElement( mWaterContent, "water-content", aOut, aTabs );
    XMLWriteElement( mResMass, "Residue-Mass", aOut, aTabs );
    XMLWriteElement( mCropMass, "Crop-Mass", aOut, aTabs );
    XMLWriteElement( mResAvail, "resAvail", aOut, aTabs );
    XMLWriteElement( mMeanErosCtrl, "meanErosCtrl", aOut, aTabs );
    XMLWriteElement( mMaxBioEnergySupply, "maxBioEnergySupply", aOut, aTabs );
    XMLWriteElement( mFractProduced, "fraction-produced", aOut, aTabs );
    
    const vector<pair<double,double> > pairs = mCostCurve->getSortedPairs();
    typedef vector<pair<double, double> >::const_iterator PairIterator;
    map<string, double> attrs;
    for( PairIterator currPair = pairs.begin(); currPair != pairs.end(); ++currPair ) {
        attrs[ "price" ] = currPair->first;
        XMLWriteElementWithAttributes( currPair->second, "fract-harvested", aOut, aTabs, attrs );
    }
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

void ResidueBiomassOutput::toInputXML( std::ostream& aOut, Tabs* aTabs ) const
{
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs, getName() );
    XMLWriteElementCheckDefault( mErosCtrl, "eros-ctrl", aOut, aTabs, 0.0 );
    XMLWriteElement( mHarvestIndex, "harvest-index", aOut, aTabs );
    XMLWriteElement( mMassConversion, "mass-conversion", aOut, aTabs );
    XMLWriteElement( mMassToEnergy, "mass-to-energy", aOut, aTabs );
    XMLWriteElement( mWaterContent, "water-content", aOut, aTabs );

    
    const vector<pair<double,double> > pairs = mCostCurve->getSortedPairs();
    typedef vector<pair<double, double> >::const_iterator PairIterator;
    map<string, double> attrs;
    for( PairIterator currPair = pairs.begin(); currPair != pairs.end(); ++currPair ) {
        attrs[ "price" ] = currPair->first;
        XMLWriteElementWithAttributes( currPair->second, "fract-harvested", aOut, aTabs, attrs );
    }

    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

bool ResidueBiomassOutput::XMLParse( const xercesc::DOMNode* aNode )
{
    // assume we are passed a valid node.
    if ( !aNode ) {
        return false;
    }

    // get all the children.
    xercesc::DOMNodeList* nodeList = aNode->getChildNodes();
    if ( !nodeList ) {
        return false;
    }

    ExplicitPointSet* currPoints = new ExplicitPointSet();

    // get the sector name attribute.
    std::string sectorName = XMLHelper<std::string>::getAttr( aNode, "name" );
    if ( !sectorName.length() ) {
        return false;
    }
    setName( sectorName );

    XMLSize_t n = nodeList->getLength();
    for ( XMLSize_t i = 0; i != n; ++i ) {
        const xercesc::DOMNode* curr = nodeList->item( i );
        if ( !curr ) {
            return false;
        }

        const std::string nodeName = XMLHelper<std::string>::safeTranscode( curr->getNodeName() );
        if( nodeName == "#text" ) {
            continue;
        }
        else if( nodeName == "eros-ctrl" ) {
            mErosCtrl = XMLHelper<double>::getValue( curr );
        }
        else if( nodeName == "harvest-index" ) {
            mHarvestIndex = XMLHelper<double>::getValue( curr );
        }
        else if( nodeName == "mass-conversion" ) {
            mMassConversion = XMLHelper<double>::getValue( curr );
        }
        else if( nodeName == "mass-to-energy" ) {
            mMassToEnergy = XMLHelper<double>::getValue( curr );
        }
        else if( nodeName == "water-content" ) {
            mWaterContent = XMLHelper<double>::getValue( curr );
        }
        else if ( nodeName == "fract-harvested" ){
            double price = XMLHelper<double>::getAttr( curr, "price" );  
            double fractionHarvested = XMLHelper<double>::getValue( curr );
            XYDataPoint* currPoint = new XYDataPoint( price, fractionHarvested );
            currPoints->addPoint( currPoint );
        }
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Unrecognized text string: " << nodeName << " found while parsing "
                    << getXMLNameStatic() << "." << std::endl;
            return false;
        }
    }
    
    mCostCurve.reset( new PointSetCurve( currPoints ) );

    return true;
}

void ResidueBiomassOutput::doInterpolations( const int aYear, const int aPreviousYear,
                                             const int aNextYear, const IOutput* aPreviousOutput,
                                             const IOutput* aNextOutput )
{
    // TODO: what to interpolate?
}

// end of residue_biomass_output.cpp ***********************************


