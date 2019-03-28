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
* \file interm_subsector.cpp
* \ingroup Objects
* \brief IntermittentSubsector class source file.
* \author Marshall Wise, Sonny Kim
*/

#include "util/base/include/definitions.h"
#include <string>
#include <cassert>
#include <vector>
#include <xercesc/dom/DOMNode.hpp>

#include "technologies/include/intermittent_technology.h"
#include "containers/include/scenario.h"
#include "containers/include/iinfo.h"
#include "containers/include/info_factory.h"
#include "util/base/include/model_time.h"
#include "util/base/include/xml_helper.h"
#include "marketplace/include/marketplace.h"
#include "sectors/include/ibackup_calculator.h"
#include "sectors/include/backup_calculator_factory.h"
#include "sectors/include/sector_utils.h"
#include "functions/include/iinput.h"
#include "functions/include/non_energy_input.h"
#include "technologies/include/iproduction_state.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

/*!
 * \brief Constructor.
 * \author Marshall Wise, Sonny Kim
 */

IntermittentTechnology::IntermittentTechnology( const string& aName, const int aYear ) 
:Technology( aName, aYear ),
mElectricSectorName("electricity"),
mBackupCapacityFactor(0.05),
mBackupCapitalCost( 0.0 ),
mElecReserveMargin( 0.15 ),
mAveGridCapacityFactor( 0.60 )
{
}

/*! 
 * \brief Copy constructor.
 * \param aOther Technology from which to copy data.
 */
IntermittentTechnology::IntermittentTechnology( const IntermittentTechnology& aOther )
: Technology( aOther ),mElectricSectorName( aOther.mElectricSectorName ),
mTrialMarketNameParsed( aOther.mTrialMarketNameParsed ),
mBackupCapacityFactor( aOther.mBackupCapacityFactor ),
mBackupCapitalCost( aOther.mBackupCapitalCost )
// Only copy member variables that are read-in. The rest will be filled in by
// initialization methods.
{
    if( aOther.mBackupCalculator.get() ){
        mBackupCalculator.reset( aOther.mBackupCalculator->clone() );
    }
}
    
IntermittentTechnology* IntermittentTechnology::clone() const {
    return new IntermittentTechnology( *this );
}

const string& IntermittentTechnology::getXMLName() const {
    return getXMLNameStatic();
}

/*! \brief Get the XML node name in static form for comparison when parsing XML.
*
* This public function accesses the private constant string, XML_NAME. This way
* the tag is always consistent for both read-in and output and can be easily
* changed. The "==" operator that is used when parsing, required this second
* function to return static.
* \note A function cannot be static and virtual.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME as a static.
*/
const string& IntermittentTechnology::getXMLNameStatic() {
    const static string XML_NAME = "intermittent-technology";
    return XML_NAME;
}

/*! \brief Return name to be used for input object containing backup capital costs.
*
* \author Steve Smith
* \return The constant XML_NAME as a static.
*/
const string& IntermittentTechnology::getBackupCapCostName( ) const {
   const static string BACKUP_CAPCOST_NAME = "backup-cap-cost";
   return BACKUP_CAPCOST_NAME;
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
const string& IntermittentTechnology::getTechCostName( ) const {
   const static string TECH_COST_NAME = ""; // no tech cost implmented in this class
   return TECH_COST_NAME;
}

bool IntermittentTechnology::XMLDerivedClassParse( const string& aNodeName,
                                                   const DOMNode* aCurr )
{
    if( aNodeName == "electric-sector-name" ){
        mElectricSectorName = XMLHelper<string>::getValue( aCurr );
    }
    // Dependent trial market name to allow multiple and different intermittent
    // technologies to contribute to same trial market.
    else if( aNodeName == "trial-market-name" ){
        mTrialMarketNameParsed = XMLHelper<string>::getValue( aCurr );
    }
    // Backup capacity factor.
    else if( aNodeName == "backup-capacity-factor" ){
        mBackupCapacityFactor = XMLHelper<double>::getValue( aCurr );
    }
    // Backup capital cost.
    else if( aNodeName == "backup-capital-cost" ){
        mBackupCapitalCost = XMLHelper<double>::getValue( aCurr );
    }
    else if( BackupCalculatorFactory::isOfType( aNodeName ) ){
        // Check if a new backup calculator needs to be created because
        // there is not currently one or the current type does not match the
        // new type.
        if( !mBackupCalculator.get() || !mBackupCalculator->isSameType( aNodeName ) ){
            mBackupCalculator = BackupCalculatorFactory::create( aNodeName );
        }
        mBackupCalculator->XMLParse( aCurr );
    }
    else if( aNodeName == "trial-market-price" ){
        mTrialMarketPrice.set( XMLHelper<double>::getValue(aCurr) );
    }
    else {
        return false;
    }
    return true;
}

void IntermittentTechnology::toInputXMLDerived( ostream& aOut, Tabs* aTabs ) const {
    XMLWriteElement( mElectricSectorName, "electric-sector-name", aOut, aTabs);
    XMLWriteElementCheckDefault( mTrialMarketNameParsed, "trial-market-name", aOut, aTabs, string("") );
    if( mTrialMarketPrice.isInited() ){
        XMLWriteElementCheckDefault( mTrialMarketPrice.get(), "trial-market-price", aOut, aTabs, 0.001 );
    }
    if( mBackupCapacityFactor.isInited() ){
        XMLWriteElement( mBackupCapacityFactor.get(), "backup-capacity-factor", aOut, aTabs);
    }
    if( mBackupCapitalCost.isInited() ){
        XMLWriteElement( mBackupCapitalCost.get(), "backup-capital-cost", aOut, aTabs);
    }
    if( mBackupCalculator.get() ){
        mBackupCalculator->toInputXML( aOut, aTabs );
    }
}   

void IntermittentTechnology::toDebugXMLDerived( const int period, ostream& aOut, Tabs* aTabs ) const {
    XMLWriteElement( mElectricSectorName, "electric-sector-name", aOut, aTabs);
    XMLWriteElementCheckDefault( mTrialMarketNameParsed, "trial-market-name", aOut, aTabs, string("") );
    if( mTrialMarketPrice.isInited() ){
        XMLWriteElementCheckDefault( mTrialMarketPrice.get(), "trial-market-price", aOut, aTabs, 0.001 );
    }
    if( mBackupCapacityFactor.isInited() ){
        XMLWriteElement( mBackupCapacityFactor.get(), "backup-capacity-factor", aOut, aTabs );
    }
    if( mBackupCapitalCost.isInited() ){
        XMLWriteElement( mBackupCapitalCost.get(), "backup-capital-cost", aOut, aTabs );
    }
    XMLWriteElement( calcEnergyFromBackup(), "energy-to-backup", aOut, aTabs );
    if( mBackupCalculator.get() ){
        mBackupCalculator->toDebugXML( period, aOut, aTabs );
    }
}

/*! \brief Create a trial market for the technology and complete the initialization.
*
* This routine is only called once per model run
* \param aRegionName Region name.
* \param aSectorName Sector name.
* \param aSubsectorName Subsector name.
* \param aSectorInfo Sector information object.
* \param aDependencyFinder Regional dependency finder.
* \param aLandAllocator Regional land allocator.
* \author Marshall Wise, Sonny Kim
* \detail A trial market for the intermittent technology is created here,  
*         as well as the completetion of technology initialization. 
*         A trial market is created for each intermittent technology if common
*         trial market name is not given.
*         Use techInfo to pass infomation to SectorUtils for setting
*         market parameters, such as units.
* \todo Member constant values (Backup Capacity Factor, Backup Cost,
*       Electricity Reserve Margin, and Ave grid Capacity Factor) could
*       be dynamically calculated and utilized by intermittent technology.
*/
void IntermittentTechnology::completeInit( const string& aRegionName,
                                           const string& aSectorName,
                                           const string& aSubsectorName,
                                           DependencyFinder* aDepFinder,
                                           const IInfo* aSubsectorInfo,
                                           ILandAllocator* aLandAllocator )
{
	// The parent method must be called first due to sequence issues
    Technology::completeInit( aRegionName, aSectorName, aSubsectorName, aDepFinder, aSubsectorInfo,
							 aLandAllocator );	
	
	// Initialize electric reserve margin and average grid capacity factor from the Sector.
    mElecReserveMargin = aSubsectorInfo->getDouble( "electricity-reserve-margin", true );
    mAveGridCapacityFactor = aSubsectorInfo->getDouble( "average-grid-capacity-factor", true );

    // Initialize a non-energy input to hold backup capacity charges
    // This needs be be done before calling other methods so that 1) input vector size is fixed 
    // (so references won't change) and 2) so that initCalc() methods for new objects can be called.
    if( util::searchForValue( mInputs, getBackupCapCostName() ) == mInputs.end() ){
        mInputs.push_back( new NonEnergyInput( getBackupCapCostName() ) );
    }

    // Inititalize info object
    mIntermittTechInfo.reset( InfoFactory::constructInfo( 0, getName() ) );
    // Output unit for intermittent technology to be for market.
    mIntermittTechInfo->setString( "output-unit", "share" );

    // If trial market name has not been readin use technology name as default trial market name.
    if( mTrialMarketNameParsed.empty() ){
        mTrialMarketName = getName(); // technology name
    }
    else {
        mTrialMarketName = mTrialMarketNameParsed;
    }

    // Create trial market for intermettent technology if backup exists and needs to be
    // calculated.
    if( mBackupCalculator.get() ){
        SectorUtils::createTrialSupplyMarket( aRegionName, mTrialMarketName, mIntermittTechInfo.get() );
    }

    // Warn if a backup calculator was not read-in.
    if( !mBackupCalculator.get() ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::NOTICE );
        mainLog << "Intermittent technology " << mName << " in sector " << aSectorName
                << " in region " << aRegionName
                << " did not read in a backup calculator. Backup costs will default to zero. " << endl;
    }
    // Convert technology year to period
    int period = scenario->getModeltime()->getyr_to_per( year );
    // Initialize trial market price here instead of in initCalc to avoid
    // price override.
    if( mTrialMarketPrice.isInited() ){
        scenario->getMarketplace()->setPrice( 
            SectorUtils::getTrialMarketName(mTrialMarketName),
            aRegionName, mTrialMarketPrice.get(), period, true );
    }
}

void IntermittentTechnology::initCalc( const string& aRegionName,
                                       const string& aSectorName,
                                       const IInfo* aSubsectorInfo,
                                       const Demographic* aDemographics,
                                       PreviousPeriodInfo& aPrevPeriodInfo,
                                       const int aPeriod )
{
    // Note: initCalc is called for all past, current and future technologies.
    Technology::initCalc( aRegionName, aSectorName, aSubsectorInfo,
        aDemographics, aPrevPeriodInfo, aPeriod );
    if ( mBackupCalculator.get() ) {
        mBackupCalculator->initCalc( mIntermittTechInfo.get() );
    }
    initializeInputLocations( aRegionName, aSectorName, aPeriod );
}

void IntermittentTechnology::postCalc( const string& aRegionName,
                                      const int aPeriod )
{
    // Use the base class postCalc.
    Technology::postCalc( aRegionName, aPeriod );
    // Store trial market solution price 
    if( mProductionState[ aPeriod ]->isNewInvestment() ) { // do only once
        mTrialMarketPrice.set( scenario->getMarketplace()->getPrice(
            SectorUtils::getTrialMarketName(mTrialMarketName),
            aRegionName, aPeriod, true) );
    }
    // TODO: hack since interpolated intermittent technologies are written out
    if( !mParsedShareWeight.isInited() ) {
        mParsedShareWeight = mShareWeight;
    }
}

/*! \brief Set intermittent technology outputs and inputs and shares for trial market.
* \author Sonny Kim
* \param aRegionName name of the region.
* \param aSectorName name of the sector.
* \param aVariableDemand current demand for the technology.
* \param aFixedOutputScaleFactor
* \param aGDP GDP
* \param aPeriod Model period
*/
void IntermittentTechnology::production( const string& aRegionName,
                                         const string& aSectorName, 
                                         double aVariableDemand,
                                         double aFixedOutputScaleFactor,
                                         const GDP* aGDP,
                                         const int aPeriod )
{
    // Use the base class production to set outputs and inputs.
    Technology::production( aRegionName, aSectorName, aVariableDemand,
                            aFixedOutputScaleFactor, aGDP, aPeriod );
    
    // For the trial intermittent technology market, set the trial supply amount to
    // the ratio of intermittent-technology output to the trial electricity output.
    // Trial electricity market must exist.
    double dependentSectorOutput = scenario->getMarketplace()->getPrice( 
        SectorUtils::getTrialMarketName(mElectricSectorName), aRegionName, aPeriod, true );

    double currentTechRatio = 0;
    if ( dependentSectorOutput > 0 ){
        currentTechRatio = getOutput( aPeriod ) / dependentSectorOutput;
    }

    // Multiple vintaged intermittent technology ratios are additive. This gives one 
    // share for backup calculation and proper behavior for vintaging intermittent technologies.
    SectorUtils::addToTrialDemand( aRegionName, mTrialMarketName, currentTechRatio, aPeriod );
}

double IntermittentTechnology::calcShare( const string& aRegionName,
                                          const string& aSectorName, 
                                          const GDP* aGDP,
                                          const double aLogitExp,
                                          const int aPeriod ) const
{
    // Use the base class calcShare.
    return Technology::calcShare( aRegionName, aSectorName, aGDP, aLogitExp, aPeriod );
}

/*! \brief Set tech shares based on backup energy needs for an intermittent
*          resource.
* \author Marshall Wise
* \param aPeriod Model period.
*/
void IntermittentTechnology::setCoefficients( const string& aRegionName,
                                              const string& aSectorName,
                                              const int aPeriod )
{
    // Convert backup capacity per unit of resource energy to energy required
    // (in EJ) per unit of resource energy (in EJ) using backup capacity factor.
    // Based on average backup capacity as this is multiplied by sector output
    // to get total backup electricity.
    double backupEnergyFraction = getAverageBackupCapacity( aRegionName, aSectorName,
                                  aPeriod ) * calcEnergyFromBackup();

    /*! \invariant Backup energy fraction must be positive. */
    assert( util::isValidNumber( backupEnergyFraction ) &&
              backupEnergyFraction >= 0 );

    // The inputs will only be invalid if the dataset was invalid.
    if( mResourceInput != mInputs.end() && mBackupInput != mInputs.end() ){
        // Normalize coefficients so that the energy output of the sector is apportioned appropriately to the two energy inputs.
        //
        // NOTE: the version currently in the multi-inputs branch only changes the coefficient for 
        // backup, leaving the coefficient for the resource 1. This, essentially, assumes that as backup comes in 
        // some of the resource energy is lost. This assumption makes relatively little difference for wind
        // (since the amount of energy from backup is small), but is not correct for CSP, where sector output is split between
        // CSP and backup mode. 
        double newCoefficient = getResourceToEnergyRatio(aRegionName, aSectorName, aPeriod ) /
                                ( 1.0 + backupEnergyFraction );
        ( *mResourceInput )->setCoefficient( newCoefficient, aPeriod );
        ( *mBackupInput )->setCoefficient( backupEnergyFraction / ( 1.0 + backupEnergyFraction ), aPeriod );
    }
}

/*! \brief Return amount of resource needed per unit of energy output.
*  This method should be used when a technology uses a resource that is
*  not in energy units. 
* \author Steve Smith
* \param aPeriod Model period.
*/
double IntermittentTechnology::getResourceToEnergyRatio( const string& aRegionName,
                                                         const string& aSectorName,
                                                         const int aPeriod )
{
    // Default assumpion is that resource is in energy units
   return 1.0;
}

/*! \brief Computes weighted cost of all technologies in Subsector plus backup
*          costs.
* \details Computes a total cost of the subsector by adding the weighted
*          technology costs and adding the additional backup cost based on the
*          backup capacity required.
* \author Marshall Wise, Sonny Kim
* \param aGDP Regional GDP container.
* \param aPeriod Model period.
*/
void IntermittentTechnology::calcCost( const string& aRegionName,
                                       const string& aSectorName,
                                       const int aPeriod )
{
    // Set marginal cost for backup to the input object set asside for this
    ( *mBackupCapCostInput )->setPrice( aRegionName, 
                              getMarginalBackupCapCost( aRegionName, mTrialMarketName, aPeriod ), 
                              aPeriod );
   
    // Set the coefficients for energy and backup in the production function.
    // Must call this after costs for backup capital and technology have been set.
    setCoefficients( aRegionName, mTrialMarketName, aPeriod );

    // Calculate the base technology cost. This will use the standard leontief
    // production function with updated coefficients for the fuel and the
    // backup.
    Technology::calcCost( aRegionName, aSectorName, aPeriod );
}

/*! \brief Returns marginal cost for backup capacity
* \author Marshall Wise, Steve Smith
* \param aSectorName Sector name.
* \param aRegionName Region name.
* \param aPeriod Model period.
*/
double IntermittentTechnology::getMarginalBackupCapCost( const string& aRegionName,
                                                         const string& aSectorName,
                                                         const int aPeriod ) const
{
    // Add per unit cost of backup capacity to subsector price backup capacity
    // is in GW/EJ, so have to convert to kW/GJ (multiply numerator by 1E6 and
    // denominator by 1E9 to get * 1/1000) to make consistent with market price
    // which is in $/GJ. BackupCost is in $/kw/yr.
    double backupCost = getMarginalBackupCapacity( aRegionName, aSectorName, aPeriod )
                        / 1000 * mBackupCapitalCost;   
   return backupCost;
}

/*!
 * \brief Get the marginal backup capacity required per unit of energy output.
 * \details Uses the internal backup calculator to determine the marginal backup
 *          capacity per unit output. If a backup calculator was not read-in,
 *          this is assumed to be zero. 
 * \author Marshall Wise, Steve Smith, Sonny Kim
 * \param aPeriod Model period.
 * \return Marginal backup capacity per unit of energy output.
 */
double IntermittentTechnology::getMarginalBackupCapacity( const string& aRegionName,
                                                          const string& aSectorName,
                                                          const int aPeriod ) const {
    double backupCapacity = 0;
    if( mBackupCalculator.get() && mResourceInput != mInputs.end() ){
        const string& resourceName = ( *mResourceInput )->getName();
        backupCapacity = mBackupCalculator->getMarginalBackupCapacity( aSectorName,
                         mElectricSectorName, resourceName, aRegionName,
                         mElecReserveMargin, mAveGridCapacityFactor, aPeriod );
    }

    /*! \post Backup capacity is a valid number and positive. */
    assert( backupCapacity >= 0 && util::isValidNumber( backupCapacity ) );
    return backupCapacity;
}

/*!
 * \brief Get the average backup capacity per unit output for the intermittent
 *        subsector.
 * \details Uses the internal backup calculator to determine the average backup
 *          capacity per unit output. If a backup calculator was not read-in,
 *          this is assumed to be zero.
 * \author Marshall Wise, Steve Smith, Sonny Kim
 * \param aPeriod Model period.
 * \return Average backup capacity per unit output.
 */
double IntermittentTechnology::getAverageBackupCapacity( const string& aRegionName,
                                                         const string& aSectorName,
                                                         const int aPeriod ) const
{
    double backupCapacity = 0;
    if( mBackupCalculator.get() ){
        const string& resourceName = ( *mResourceInput )->getName();
        backupCapacity = mBackupCalculator->getAverageBackupCapacity( aSectorName,
                         mElectricSectorName, resourceName, aRegionName,
                         mElecReserveMargin, mAveGridCapacityFactor, aPeriod );
    }

    /*! \post Backup capacity is a valid number and positive. */
    assert( backupCapacity >= 0 && util::isValidNumber( backupCapacity ) );
    return backupCapacity;
}

/*! 
 * \brief Initialize the cached locations of the resource and backup inputs.
 * \details Determines and caches the locations of the resource and backup
 *          inputs. The resource input is assumed to be the input with a
 *          variance. The backup input is assumed to be the remaining energy
 *          input.
 * \param aRegionName Name of the containing region.
 * \param aSectorName Name of the containing sector.
 * \param aPeriod Period.
 * \warning If the input vector changes size, these positions will not be valid.
 */
void IntermittentTechnology::initializeInputLocations( const string& aRegionName,
                                                       const string& aSectorName,
                                                       const int aPeriod )
{
    // Set the inputs to the error value.
    mBackupCapCostInput = mTechCostInput = mResourceInput = mBackupInput = mInputs.end();

    const Marketplace* marketplace = scenario->getMarketplace();
    for( InputIterator i = mInputs.begin(); i != mInputs.end(); ++i ){
        // Parse location for non-energy inputs.
        if( !( *i )->hasTypeFlag( IInput::ENERGY ) ){
            if ( ( *i )->getName() == getBackupCapCostName() && ( getBackupCapCostName() != "" ) ) {
               mBackupCapCostInput = i;
               continue;
            } 
            else if ( ( *i )->getName() == getTechCostName() && ( getTechCostName() != "" ) ) {
               mTechCostInput = i;
               continue;
            }
        }

        // Otherwise, if this is any other non-energy input then skip this input
        if( !( *i )->hasTypeFlag( IInput::ENERGY ) ){
            continue;
        }
        
        // Determine the location of the resource and backup input. Resource input
        // is known to be the input with a variance. Backup input is assumed to be
        // the other energy input.
        
        // Use period 0 marketInfo object to determine which input has a resource variance and is therefore a resource.
        const IInfo* info = marketplace->getMarketInfo( ( *i )->getName(), aRegionName, 0, true );
        if( !info ){
            continue;
        }

        // If the good has a variance it must be the resource input.
        if( info->hasValue( "resourceVariance" ) ){
            if( mResourceInput != mInputs.end() ){
                // There already was a resource input.
                ILogger& mainLog = ILogger::getLogger( "main_log" );
                mainLog.setLevel( ILogger::NOTICE );
                mainLog << "Intermittent technology " << mName << " in sector " << aSectorName
                        << " in region " << aRegionName << " has more than one variable resource input." << endl;
            }
            else {
                mResourceInput = i;
            }
        }
        // If it is an energy input that is not the resource it must be the backup.
        else {
            if( mBackupInput != mInputs.end() ){
              // There already was a resource input.
              ILogger& mainLog = ILogger::getLogger( "main_log" );
                mainLog.setLevel( ILogger::NOTICE );
                mainLog << "Intermittent technology " << mName << " in sector " << aSectorName
                        << " in region " << aRegionName
                        << " has more than one energy input that is not the resource." << endl;
            }
            else {
                mBackupInput = i;
            }
        }
    }

    // Check that both the resource and backup input were set.
    if( mResourceInput == mInputs.end() || mBackupInput == mInputs.end() ){
        // There already was a resource input.
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::NOTICE );
        mainLog << "Intermittent technology " << mName << " in sector " << aSectorName
                << " in region " << aRegionName << " does not have the required resource and backup inputs."
                << endl;
    }
}

/*!
 * \brief Determine the amount of energy produced by the backup per unit of
 *        capacity.
 * \details Determines the amount of energy produced per unit of backup capacity
 *          by adjusting for the backup capacity factor and converting the
 *          resulting operating backup into energy.
 * \return Capacity to energy production conversion factor.
 * \todo Move all units conversion to utilities.
 */
double IntermittentTechnology::calcEnergyFromBackup() const {
    // Conversion: 1 gigaWattHour of electricity = 3.6E-6 ExaJoules
    const double EJ_PER_GWH = 0.0000036;
    // Number of hours in a year.
    const int HOURS_PER_YEAR = 8760;
    return HOURS_PER_YEAR * mBackupCapacityFactor * EJ_PER_GWH;
}
