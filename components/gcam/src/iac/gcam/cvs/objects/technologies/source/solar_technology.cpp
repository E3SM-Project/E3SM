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
 * solar_technology.cpp
 * Created: 03/06/2007
 * Version: 04/24/2007
 *
 * This software, which is provided in confidence, was prepared by employees
 * of Pacific Northwest National Laboratory operated by Battelle Memorial
 * Institute. Battelle has certain unperfected rights in the software
 * which should not be copied or otherwise disseminated outside your
 * organization without the express written authorization from Battelle.
 * All rights to the software are reserved by Battelle.   Battelle makes no
 * warranty, express or implied, and assumes no liability or responsibility
 * for the use of this software.
 */

// include files ***********************************************************

#include "util/base/include/definitions.h"
#include "technologies/include/solar_technology.h"
#include "technologies/include/marginal_profit_calculator.h"
#include "technologies/include/iproduction_state.h"
#include "marketplace/include/marketplace.h"
#include "containers/include/iinfo.h"
#include "sectors/include/sector_utils.h"
#include "util/base/include/TValidatorInfo.h"
#include "util/base/include/util.h"
#include "util/base/include/xml_helper.h"
#include "functions/include/non_energy_input.h"

#include <algorithm>
using namespace std;

// namespaces **************************************************************

// constants ***************************************************************

const double SolarTechnology::kWhrtoGJ = 0.0036;
const std::string SolarTechnology::ELECTRIC_SECTOR_NAME_KEY = "electricSectorName";
const std::string SolarTechnology::NO_SUN_DAYS_KEY = "no-sun-days";
//const std::string SolarTechnology::TOTAL_ANNUAL_IRRADIANCE_KEY = "total-annual-irradiance";
const double DEFAULT_TOTAL_ANNUAL_IRRADIANCE = 1000.0;

// Constructors: SolarTechnology *********************************************

/*! Default constructor
 *  \param aName the name of the technology
 *  \param aYear the year
 */
SolarTechnology::SolarTechnology(
   const std::string& aName,
   const int          aYear )
   : parent( aName, aYear ),
     mCapitalCost( 3486 ),
     mCSPCapacityFactor( 1 ),
     mElectricSectorName(),
     mFCR( 0.0856 ),
     mGridConnectionCost( 1500 ),
     mOM( 47.87 ),
     mRegionName(),
     mSectorName(),
     mSolarFieldFraction( 0.3 ),
     mSolarFieldArea( 6.9 ),
     mMaxLoss( 0.55 ),
     mEfficiencyLossExponent( 3.0 ),
     mMaxSectorLoadServed( 0.15 ),
     mPlantAvailability( 0.98 ),
     mScheduledMaintenance( 0.15 ),
     mRandomMaintenanceFraction( 0.50 ),
     mCSPEfficiency( 0.001 ) // non-sensical value -- this should always be read-in
{
}

/*! Copy constructor
 *  \param other the instance to copy
 */
SolarTechnology::SolarTechnology( const SolarTechnology& other )
   : parent( other ),
     mCapitalCost( other.mCapitalCost ),
     mCSPCapacityFactor( other.mCSPCapacityFactor ),
     mElectricSectorName( other.mElectricSectorName ),
     mFCR( other.mFCR ),
     mGridConnectionCost( other.mGridConnectionCost ),
     mOM( other.mOM ),
     mRegionName( other.mRegionName ),
     mSectorName( other.mSectorName ),
     mSolarFieldFraction( other.mSolarFieldFraction ),
     mSolarFieldArea( other.mSolarFieldArea )
{
}

// Destructor: SolarTechnology ***********************************************

SolarTechnology::~SolarTechnology(void)
{
}

// SolarTechnology::calcCost *************************************************

// Documentation is inherited
void SolarTechnology::calcCost(
   const std::string& aRegionName,
   const std::string& aSectorName,
   const int          aPeriod )
{
    if( !mProductionState[ aPeriod ]->isOperating() ){
        return;
    }

   // Var to hold capital and operating costs
   double totalTechCapOMCost = 0;

   static const double denom = 24 * 365;

   // Get marketplace and calculate costs
   Marketplace*       pMarketplace = scenario->getMarketplace();
   const IInfo*       pInfo        = pMarketplace->getMarketInfo( ( *mResourceInput )->getName(), aRegionName, aPeriod, true );

   double totalAnnualIrradiance = pInfo->hasValue( mTotalAnnualIrradianceKey ) ? pInfo->getDouble( mTotalAnnualIrradianceKey, true ) : DEFAULT_TOTAL_ANNUAL_IRRADIANCE;
   double dConnect              = pMarketplace->getPrice( ( *mResourceInput )->getName(), aRegionName, aPeriod );
   double CSPEfficiency         = getSolarEfficiency( aPeriod );

   // Outage Fraction -- amount of time plant is not operating either on solar or backup mode
   // This is equal to the amount of time required for scheduled maintenance plus unscheduled outages
   // All other times plant is operating in solar or hybrid mode, so all costs are amortized over the remaining time
   double nonRandomMaintenance = 1 - mRandomMaintenanceFraction;
   double outageFraction         = mScheduledMaintenance +
                                   mPlantAvailability * ( 1 - mScheduledMaintenance * nonRandomMaintenance ); 
                                   
   if( !mProductionState[ aPeriod ]->isNewInvestment() ||
       mFixedOutput != IProductionState::fixedOutputDefault() )
   {
      // If an existing stock, only have O&M costs.
      totalTechCapOMCost = mOM / ( totalAnnualIrradiance * mSolarFieldArea * CSPEfficiency * kWhrtoGJ * outageFraction );
   }
   else {
       // Calculate total cost, equal to generation plus connection and backup cost. 
       // Equation 1.
       mCGeneration = ( mFCR * mCapitalCost + mOM ) / ( totalAnnualIrradiance * mSolarFieldArea * CSPEfficiency * kWhrtoGJ * outageFraction );
       // Equation 3.
       mCSPCapacityFactor = ( totalAnnualIrradiance * mSolarFieldArea * CSPEfficiency * outageFraction ) / denom;
       // Equation 2.
       mCConnect = mFCR * dConnect * mGridConnectionCost / ( mCSPCapacityFactor * 1000 * kWhrtoGJ * 24.0 * 365.0 );

       // Total cost is the sum of the generation and connection cost
       totalTechCapOMCost = std::max( mCGeneration + mCConnect, util::getSmallNumber() );
   }

   // Set tech capital and operating costs to input object
   ( *mTechCostInput )->setPrice( aRegionName, totalTechCapOMCost, aPeriod );
    
   // Call parent function to set costs and coefficients
   parent::calcCost( aRegionName, aSectorName, aPeriod );
   
}

// SolarTechnology::calcResourceArea *****************************************

/*! Calculate the resource area in km^2
 *  \param aRegionName the region name
 *  \param aSectorName the sector name
 *  \param aVariableDemand the variable demand
 *  \param aPeriod the period
 *  \return the resource area
 */
double SolarTechnology::calcResourceArea(
   const std::string& aRegionName,
   const std::string& aSectorName,
   double             aVariableDemand,
   const int          aPeriod )
{
/*
 * Equation 4 from "First Phase Implementation of CSP" is as follows:
 *
 * CSPGeneration / (kWhrtoGJ • 10^-9) = totalAnnualIrradiance • (ResourceArea •1000^2)• solarFieldFraction • CSPEfficiency
 *
 * Solving for CSPGeneration, we get:
 *
 * CSPGeneration  = totalAnnualIrradiance • ResourceArea • solarFieldFraction • CSPEfficiency • kWhrtoGJ • 10^-9 •1000^2
 *
 * And finally, solving for ResourceArea, we get:
 *
 * ResourceArea = CSPGeneration / ( totalAnnualIrradiance  • solarFieldFraction • CSPEfficiency • kWhrtoGJ • 10^-9 •1000^2 )
 */
   static const double conversionFact = kWhrtoGJ * 1e-3;

   // Set market demand for km^2
   Marketplace*       pMarketplace = scenario->getMarketplace();
   const IInfo*       pInfo        = pMarketplace->getMarketInfo( ( *mResourceInput )->getName(), aRegionName, aPeriod, true );

   double totalAnnualIrradiance = pInfo->hasValue( mTotalAnnualIrradianceKey ) ? pInfo->getDouble( mTotalAnnualIrradianceKey, true ) : DEFAULT_TOTAL_ANNUAL_IRRADIANCE;
   double CSPGeneration         = aVariableDemand;
   double CSPEfficiency         = getSolarEfficiency( aPeriod );

   // Efficiency should be positive
   assert( CSPEfficiency > 0 );

   // Equation 4.
   double resourceArea          = CSPGeneration / ( totalAnnualIrradiance * mSolarFieldFraction * CSPEfficiency * conversionFact );

   return resourceArea;
}

// SolarTechnology::calcShare ************************************************

// SolarTechnology::clone ****************************************************

// Documentation is inherited
SolarTechnology* SolarTechnology::clone( void ) const
{
   return new SolarTechnology( *this );
}

// SolarTechnology::completeInit *********************************************

// Documentation is inherited
void SolarTechnology::completeInit(
   const std::string&              aRegionName,
   const std::string&              aSectorName,
   const std::string&              aSubsectorName,
   DependencyFinder*               aDepFinder,
   const IInfo*                    aSubsectorIInfo,
   ILandAllocator*                 aLandAllocator )
{
   // Initialize a non-energy input to hold technology costs
   if( util::searchForValue( mInputs, getTechCostName() ) == mInputs.end() ){
        mInputs.push_back( new NonEnergyInput( getTechCostName() ) );
    }

   parent::completeInit( aRegionName, aSectorName, aSubsectorName, aDepFinder, aSubsectorIInfo, aLandAllocator );

   // Validate input parameters
   typedef ObjECTS::TValidatorInfo<> validator_type;
   validator_type   validator[] =
   {
      validator_type( mCapitalCost, "capital-cost", mCapitalCost > 0 ),
      validator_type( mFCR, "fcr", mFCR > 0 ),
      validator_type( mSolarFieldFraction, "solar-field-fraction", mSolarFieldFraction > 0 ),
      validator_type( mSolarFieldArea, "solar-field-area", mSolarFieldArea > 0 )
   };

   unsigned short numParams = sizeof( validator ) / sizeof( validator[0] );
   std::string    msg       = ObjECTS::getInvalidNames(
      &validator[0],
      &validator[numParams] );

   if ( msg.length() )
   // Invalid input parameter
   {
      ILogger& mainLog = ILogger::getLogger( "main_log" );
      mainLog.setLevel( ILogger::ERROR );
      mainLog << "Invalid input parameter(s) to "
         << getXMLNameStatic()
         << " in sector " << aSectorName
         << ": " << msg << std::endl;
      exit( -1 );
   }
}

/*! \brief Return penetration level of this technology
* \author Steve Smith
* \param aPeriod Model period.
*/
double SolarTechnology::getSolarPenetration( const int aPeriod ) const
{
   // Compute the CSP penetration level
   double SolarPenetration = 0;
   if ( aPeriod > 0 && mMaxSectorLoadServed > 0 ){
      SolarPenetration = SectorUtils::getTrialSupply( mRegionName, mTrialMarketName, aPeriod )
                         / mMaxSectorLoadServed;
   }

   return SolarPenetration;
}

// SolarTechnology::getSolarEfficiency ********************************************

/*! \brief Return efficiency of the CSP plant
*  The efficiency is the unitless net solar energy conversion efficiency of the plant
*  taking into account any solar losses due to system interactions
* \author Kevin Walker and Steve Smith
* \param aPeriod Model period.
*/
double SolarTechnology::getSolarEfficiency( const int aPeriod ) const
{
   // Compute the amount of solar generation lost
   double CSPLoss = std::min( mMaxLoss * std::pow( getSolarPenetration( aPeriod ), mEfficiencyLossExponent ), 0.99 );

   // Compute the net efficiency
   double result = mCSPEfficiency * ( 1.0 - CSPLoss );

   return result;
}

// SolarTechnology::getXMLName1D *********************************************

// Documentation is inherited
const std::string& SolarTechnology::getXMLName1D( void ) const
{
   return getXMLNameStatic();
}

// SolarTechnology::getXMLNameStatic1D ***************************************

const std::string& SolarTechnology::getXMLNameStatic1D( void )
{
   static const std::string XML_NAME1D = "solar-technology";
   return XML_NAME1D;
}

/*! \brief Return name to be used for input object containing technology costs.
*
* This input object will contain technology capital, operation, and any other costs
* exclusive of backup or fuel costs. Setting to blank indicates that this object does
* not use this cost.
*
* \author Steve Smith
* \return The constant XML_NAME as a static.
*/
const std::string& SolarTechnology::getTechCostName( ) const {
   const static std::string TECH_COST_NAME = "turbine-and-connect-costs";
   return TECH_COST_NAME;
}

// SolarTechnology::initCalc *************************************************

// Documentation is inherited
void SolarTechnology::initCalc(
   const std::string& aRegionName,
   const std::string& aSectorName,
   const IInfo*       aSubsectorIInfo,
   const Demographic* aDemographics,
   PreviousPeriodInfo& aPrevPeriodInfo,
   const int          aPeriod )
{
   // Need resource location so call this here.
   initializeInputLocations( aRegionName, aSectorName, aPeriod );

   // Get marketplace and make sure we have the total annual irradiance
   Marketplace*       pMarketplace = scenario->getMarketplace();
   const IInfo*       pInfo        = pMarketplace->getMarketInfo( ( *mResourceInput )->getName(), aRegionName, aPeriod, true );

   if ( !pInfo || !pInfo->hasValue( mTotalAnnualIrradianceKey ) )
   // Invalid input parameter
   {
      ILogger& mainLog = ILogger::getLogger( "main_log" );
      mainLog.setLevel( ILogger::ERROR );
      mainLog << "Invalid input parameter(s) to "
         << getXMLNameStatic()
         << " in sector " << aSectorName
         << ": " << mTotalAnnualIrradianceKey << std::endl;
   }

   // Get number of no sun days
   if ( pInfo || !pInfo->hasValue( "no-sun-days" ) )
   // Invalid input parameter
   {
      mNoSunDays = pInfo->getDouble( "no-sun-days", false );
   }
   else {
      mNoSunDays = 0;
   }

   // Put variables into an info object so can be passed into backup calculator
   mIntermittTechInfo->setDouble( "max-sector-load-served", mMaxSectorLoadServed );
   mIntermittTechInfo->setDouble( "no-sun-days", mNoSunDays );
   mIntermittTechInfo->setDouble( "random-maintenance-fraction", mRandomMaintenanceFraction );
   mIntermittTechInfo->setDouble( "scheduled-maintenance", mScheduledMaintenance );
   
   parent::initCalc( aRegionName, aSectorName, aSubsectorIInfo, aDemographics, 
       aPrevPeriodInfo, aPeriod );

   // Cache the region name
   mRegionName = aRegionName;

   // Cache the sector name
   mSectorName = aSectorName;

   // Cache the electric sector name
   if ( aSubsectorIInfo->hasValue( ELECTRIC_SECTOR_NAME_KEY ) )
   {
      mElectricSectorName = aSubsectorIInfo->getString( ELECTRIC_SECTOR_NAME_KEY, true );
   }
   else
   {
      mElectricSectorName = "electricity";
   }

}

/*! \brief Return amount of resource needed per unit of energy output.
*  This method should be used when a technology uses a resource that is
*  not in energy units. 
* \author Steve Smith
* \param aPeriod Model period.
*/
double SolarTechnology::getResourceToEnergyRatio( const std::string& aRegionName,
                                                  const std::string& aSectorName,
                                                  const int aPeriod )
{
    // Default assumpion is that resource is in energy units
   // Need to return the amount of resource (in km^2) needed per EJ of energy output
   return calcResourceArea( aRegionName, aSectorName, 1.0, aPeriod );
}

// SolarTechnology::toDebugXMLDerived ****************************************

// Documentation is inherited
void SolarTechnology::toDebugXMLDerived(
   const int     period,
   std::ostream& out,
   Tabs*         tabs ) const
{
   parent::toDebugXMLDerived( period, out,  tabs );
   XMLWriteElement( mCapitalCost, "capital-cost", out, tabs );
   XMLWriteElement( mCConnect, "c-connect", out, tabs );
   XMLWriteElement( mCGeneration, "c-generation", out, tabs );
   XMLWriteElement( mCSPCapacityFactor, "csp-capacity-factor", out, tabs );
   XMLWriteElement( mPlantAvailability, "plant-availability-fraction", out, tabs );
   XMLWriteElement( mScheduledMaintenance, "plant-scheduled-maintenance-fraction", out, tabs );
   XMLWriteElement( mRandomMaintenanceFraction, "random-maintence-fraction", out, tabs );
   XMLWriteElement( mFCR, "fcr", out, tabs );
   XMLWriteElement( mGridConnectionCost, "grid-connection-cost", out, tabs );
   XMLWriteElement( mOM, "om", out, tabs );
   XMLWriteElement( mSolarFieldFraction, "solar-field-fraction", out, tabs );
   XMLWriteElement( mSolarFieldArea, "solar-field-area", out, tabs );
   XMLWriteElement( mMaxLoss, "max-solar-loss", out, tabs );
   XMLWriteElement( mEfficiencyLossExponent, "loss-exponent", out, tabs );
   XMLWriteElement( mCSPEfficiency, "net-solar-conversion-efficiency", out, tabs );
   XMLWriteElement( mMaxSectorLoadServed, "max-sector-load-served", out, tabs );
   XMLWriteElement( mTotalAnnualIrradianceKey, "irradiance-tagname", out, tabs );
   XMLWriteElement( getSolarEfficiency( period ), "net-efficiency", out, tabs );
   XMLWriteElement( getSolarPenetration( period ), "estimated-sector-penetration", out, tabs );
}

// SolarTechnology::toInputXMLDerived ****************************************

// Documentation is inherited
void SolarTechnology::toInputXMLDerived(
   std::ostream& out,
   Tabs*         tabs ) const
{
   parent::toInputXMLDerived( out,  tabs );
   XMLWriteElementCheckDefault( mCapitalCost, "capital-cost", out, tabs, double( 3486 ) );
   XMLWriteElementCheckDefault( mFCR, "fcr", out, tabs, double( 0.0856 ) );
   XMLWriteElementCheckDefault( mGridConnectionCost, "grid-connection-cost", out, tabs, double( 1500 ) );
   XMLWriteElementCheckDefault( mOM, "om", out, tabs, double( 47.87 ) );
   XMLWriteElementCheckDefault( mSolarFieldFraction, "solar-field-fraction", out, tabs, double( 0.3 ) );
   XMLWriteElementCheckDefault( mSolarFieldArea, "solar-field-area", out, tabs, double( 6.9 ) );
   XMLWriteElementCheckDefault( mMaxLoss, "max-solar-loss", out, tabs, double( 0.55 ) );
   XMLWriteElementCheckDefault( mEfficiencyLossExponent, "loss-exponent", out, tabs, double( 3.0 ) );
   XMLWriteElementCheckDefault( mPlantAvailability, "plant-availability-fraction", out, tabs, double( 0.98 ) );
   XMLWriteElementCheckDefault( mScheduledMaintenance, "plant-scheduled-maintenance-fraction", out, tabs, double( 0.15 ) );
   XMLWriteElementCheckDefault( mRandomMaintenanceFraction, "random-maintence-fraction", out, tabs, double( 0.50 ) );
   XMLWriteElement( mCSPEfficiency, "net-solar-conversion-efficiency", out, tabs );
   XMLWriteElement( mMaxSectorLoadServed, "max-sector-load-served", out, tabs );
   XMLWriteElementCheckDefault<std::string>( mTotalAnnualIrradianceKey, "irradiance-tagname", out, tabs, "total-annual-irradiance" ); 
}

// SolarTechnology::XMLDerivedClassParse *************************************

// Documentation is inherited
bool SolarTechnology::XMLDerivedClassParse(
   const std::string&      nodeName,
   const xercesc::DOMNode* curr )
{
   if ( nodeName == "capital-cost" )
   {
      mCapitalCost = XMLHelper<double>::getValue( curr );
   }
   else if ( nodeName == "plant-availability-fraction" )
   {
      mPlantAvailability = XMLHelper<double>::getValue( curr );
   }
   else if ( nodeName == "plant-scheduled-maintenance-fraction" )
   {
      mScheduledMaintenance = XMLHelper<double>::getValue( curr );
   }
   else if ( nodeName == "random-maintence-fraction" )
   {
      mRandomMaintenanceFraction = XMLHelper<double>::getValue( curr );
   }
   else if ( nodeName == "fcr" )
   {
      mFCR = XMLHelper<double>::getValue( curr );
   }
   else if ( nodeName == "grid-connection-cost" )
   {
      mGridConnectionCost = XMLHelper<double>::getValue( curr );
   }
   else if ( nodeName == "om" )
   {
      mOM = XMLHelper<double>::getValue( curr );
   }
   else if ( nodeName == "solar-field-fraction" )
   {
      mSolarFieldFraction = XMLHelper<double>::getValue( curr );
   }
   else if ( nodeName == "solar-field-area" )
   {
      mSolarFieldArea = XMLHelper<double>::getValue( curr );
   }
   else if ( nodeName == "max-solar-loss" ){
       mMaxLoss = XMLHelper<double>::getValue( curr );
   }
   else if ( nodeName == "loss-exponent" ){
       mEfficiencyLossExponent = XMLHelper<double>::getValue( curr );
   }
   else if ( nodeName == "net-solar-conversion-efficiency" ){
       mCSPEfficiency = XMLHelper<double>::getValue( curr );
   }
   else if ( nodeName == "irradiance-tagname" ){
       mTotalAnnualIrradianceKey = XMLHelper<std::string>::getValue( curr );
   }
   else if ( nodeName == "max-sector-load-served" ){
       mMaxSectorLoadServed = XMLHelper<double>::getValue( curr );
   }
   else if ( IntermittentTechnology::XMLDerivedClassParse( nodeName, curr ) ) {
   }
   else
   {
      return false;
   }

   return true;
}

// end of solar_technology.cpp ***********************************************
