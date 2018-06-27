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
* \file sgm_gen_table.cpp
* \ingroup Objects
* \brief The SGMGenTable class source file.
*
* \author Sonny Kim
*/

#include "util/base/include/definitions.h"
#include <iostream>
#include <string>
#include <map>

#include "reporting/include/sgm_gen_table.h"

#include "containers/include/scenario.h"
#include "marketplace/include/marketplace.h"
#include "containers/include/region.h"
#include "containers/include/region_cge.h"
#include "demographics/include/demographic.h"
#include "sectors/include/sector.h"
#include "sectors/include/production_sector.h"
#include "sectors/include/subsector.h"
#include "emissions/include/aghg.h"
#include "containers/include/national_account.h"
#include "technologies/include/expenditure.h"
#include "technologies/include/base_technology.h"
#include "consumers/include/consumer.h"
#include "consumers/include/household_consumer.h"
#include "consumers/include/govt_consumer.h"
#include "consumers/include/trade_consumer.h"
#include "consumers/include/invest_consumer.h"
#include "technologies/include/production_technology.h"
#include "functions/include/iinput.h"
#include "sectors/include/factor_supply.h"
#include "util/base/include/model_time.h"
#include "functions/include/function_utils.h"
#include "technologies/include/ioutput.h"
#include "containers/include/iinfo.h"

using namespace std;
extern Scenario* scenario;

//! Default Constructor
SGMGenTable::SGMGenTable( const string& aName, const string& aHeader, const Modeltime* aModeltime ):
mName( aName ), mHeader( aHeader ), mModeltime( aModeltime ), mFile( 0 ) {
}

/*! \brief Set the file the table will print to.
* \todo Remove this.
* \param aOutputFile File to which to output.
*/
void SGMGenTable::setOutputFile( ostream& aOutputFile ){
    mFile = &aOutputFile;
}

//! Add to the value for the DCT specified by the account type key.
void SGMGenTable::addToType( const int aTypeRow, const string aTypeCol, const double value ){
    // add to column and row totals here
    // Do not add to totals anywhere else
    mTable[ aTypeRow ][ aTypeCol ] += value;
    // add to total for this type of table only
    mTable[ aTypeRow ][ "zTotal" ] += value;
}

//! set the value for the DCT specified by the account type key.
void SGMGenTable::setType( const int aTypeRow, const string aTypeCol, const double value ){
    mTable[ aTypeRow ][ aTypeCol ] = value;
}

//! Get the value for the DCT specified by the account type key.
double SGMGenTable::getValue( const int aTypeRow, const string aTypeCol ) const {
    return util::searchForValue( util::searchForValue( mTable, aTypeRow ), aTypeCol );
}

/*! \brief For outputting SGM data to a flat csv File
*
*/
void SGMGenTable::finish() const {
    /*! \pre The output file has been set. */
    assert( mFile );
    if ( !mTable.empty() ){ // don't print if empty

        *mFile << mHeader << endl;
        // Note: This is structurally different from SAM. This goes through the rows and prints
        // out each of the category values.
        // write column labels
        *mFile << "Year" << ',';
        for( map< string, double>::const_iterator strIter = (*mTable.begin()).second.begin(); strIter != (*mTable.begin()).second.end(); ++strIter ) {
            *mFile << (*strIter).first << ',';
        }
        *mFile << endl;

        for ( map<int, map<string, double> >::const_iterator mapIter = mTable.begin(); mapIter != mTable.end(); ++mapIter ) {
            *mFile << (*mapIter).first;
            for( map< string, double >::const_iterator strIter = ((*mapIter).second).begin(); strIter != ((*mapIter).second).end(); strIter++ ) {
                *mFile << ',' << (*strIter).second;
            }
            *mFile << endl;
        }
        *mFile << endl;
    }
}

void SGMGenTable::startVisitRegionCGE( const RegionCGE* regionCGE, const int aPeriod ) {
    // Store the current region name.
    mCurrentRegionName = regionCGE->getName();
    // if table output is by sector, then add sector names to map
    if( mName == "PECbySector" ) {
        for( std::vector<Sector*>::const_iterator secNameIter = regionCGE->supplySector.begin(); secNameIter != regionCGE->supplySector.end(); ++secNameIter ) {
            for ( int per = 0; per < mModeltime->getmaxper(); per++ ) {
                int year = mModeltime->getper_to_yr( per );
                addToType( year, (*secNameIter)->getName(), 0 );
            }
        }
    }
}

void SGMGenTable::endVisitRegionCGE( const RegionCGE* aRegionCGE, const int aPeriod ){
    // Clear the stored region name.
    mCurrentRegionName.clear();
}

void SGMGenTable::startVisitSector( const Sector* aSector, const int aPeriod ){
    // Store the current sector name.
    mCurrentSectorName = aSector->getName();
}

void SGMGenTable::endVisitSector( const Sector* aSector, const int aPeriod ){
    // Clear the stored sector name.
    mCurrentSectorName.clear();
}

// Write sector market prices to the SGM gen file
void SGMGenTable::startVisitProductionSector( const ProductionSector* aProductionSector, const int aPeriod ) {
    if( mName == "PRICE" ) {
        addToType( mModeltime->getper_to_yr( aPeriod ), aProductionSector->getName(),
            scenario->getMarketplace()->getPrice( aProductionSector->getName(), aProductionSector->regionName, aPeriod ) );
    }
}

// Write factor supply market prices to the SGM gen file
void SGMGenTable::startVisitFactorSupply( const FactorSupply* aFactorSupply, const int aPeriod ) {
    if( mName == "PRICE" ) {
        addToType( mModeltime->getper_to_yr( aPeriod ), aFactorSupply->getName(),
            scenario->getMarketplace()->getPrice( aFactorSupply->getName(), aFactorSupply->marketName, aPeriod ) );
    }
}

// Write National Account information to the SGM gen file
void SGMGenTable::startVisitNationalAccount( const NationalAccount* aNationalAccount, const int aPeriod ) {
    if( mName == "GNPREAL" ) {
        addToType( mModeltime->getper_to_yr( aPeriod ), "Consumption",
            aNationalAccount->getAccountValue( NationalAccount::CONSUMPTION_REAL ) );
        addToType( mModeltime->getper_to_yr( aPeriod ), "Investment",
            aNationalAccount->getAccountValue( NationalAccount::INVESTMENT_REAL ) );
        addToType( mModeltime->getper_to_yr( aPeriod ), "Government",
            aNationalAccount->getAccountValue( NationalAccount::GOVERNMENT_REAL ) );
        addToType( mModeltime->getper_to_yr( aPeriod ), "Trade Balance",
            aNationalAccount->getAccountValue( NationalAccount::NET_EXPORT_REAL ) );
        addToType( mModeltime->getper_to_yr( aPeriod ), "GNP",
            aNationalAccount->getAccountValue( NationalAccount::GNP_REAL ) );
    }
    else if( mName == "GNPNOM" ) {
        addToType( mModeltime->getper_to_yr( aPeriod ), "Consumption",
            aNationalAccount->getAccountValue( NationalAccount::CONSUMPTION_NOMINAL ) );
        addToType( mModeltime->getper_to_yr( aPeriod ), "Investment",
            aNationalAccount->getAccountValue( NationalAccount::INVESTMENT_NOMINAL ) );
        addToType( mModeltime->getper_to_yr( aPeriod ), "Government",
            aNationalAccount->getAccountValue( NationalAccount::GOVERNMENT_NOMINAL ) );
        addToType( mModeltime->getper_to_yr( aPeriod ), "Trade Balance",
            aNationalAccount->getAccountValue( NationalAccount::NET_EXPORT_NOMINAL ) );
        addToType( mModeltime->getper_to_yr( aPeriod ), "GNP",
            aNationalAccount->getAccountValue( NationalAccount::GNP_NOMINAL ) );
    }
}

// Write demographics results to the SGM gen file
void SGMGenTable::startVisitDemographic( const Demographic* aDemographic, const int aPeriod ) {
    if( mName == "DEM" ) {
        addToType( mModeltime->getper_to_yr( aPeriod ), "Tot Pop", aDemographic->getTotal( aPeriod ) );
        addToType( mModeltime->getper_to_yr( aPeriod ), "Working Age Male",
            aDemographic->getWorkingAgePopulationMales( aPeriod ) );
        addToType( mModeltime->getper_to_yr( aPeriod ), "Working Age Female",
            aDemographic->getWorkingAgePopulationFemales( aPeriod ) );
        addToType( mModeltime->getper_to_yr( aPeriod ), "Working Age Tot",
            aDemographic->getWorkingAgePopulation( aPeriod ) );
    }
}

// Write to SGM general table.
// This routine assumes that only and all operating technology vintages are passed in as
// an argument.
void SGMGenTable::startVisitConsumer( const Consumer* aConsumer, const int aPeriod )
{
    // Only update the current consumer.
    if( aConsumer->getYear() != mModeltime->getper_to_yr( aPeriod ) ){
        return;
    }

    // add output of each technology for each period
    if( ( mName == "CO2" ) || ( mName == "CO2bySec" ) || ( mName == "CO2byTech" ) ){
        unsigned int CO2index = util::searchForValue( aConsumer->mGhgNameMap, string( "CO2" ) );
        double CO2Emiss = aConsumer->mGhgs[ CO2index ]->getEmission( aPeriod );

        if( mName == "CO2" ) {
            addToType( mModeltime->getper_to_yr( aPeriod ), "CO2",  CO2Emiss );
        }
        else if( mName == "CO2bySec" ) {
            addToType( mModeltime->getper_to_yr( aPeriod ), aConsumer->getXMLName(), CO2Emiss );
        }
        else if( mName == "CO2byTech" ) {
            addToType( mModeltime->getper_to_yr( aPeriod ), aConsumer->getName(), CO2Emiss );
        }
    }
}

void SGMGenTable::startVisitHouseholdConsumer( const HouseholdConsumer* householdConsumer, const int aPeriod ) {
    if( mName == "CAP" ) {
        // add only current year consumer
        if( householdConsumer->year == mModeltime->getper_to_yr( aPeriod ) ) {
            addToType( mModeltime->getper_to_yr( aPeriod ), "Savings",
                householdConsumer->expenditures[ aPeriod ].getValue( Expenditure::SAVINGS ) );
        }
    }
    else if( mName == "DEM" ) {
        if( householdConsumer->year == mModeltime->getper_to_yr( aPeriod ) ) {
            addToType( mModeltime->getper_to_yr( aPeriod ), "Labor Supply",
                householdConsumer->getLaborSupply() );
        }
    }
}

void SGMGenTable::startVisitGovtConsumer( const GovtConsumer* govtConsumer, const int aPeriod ) {
    if( mName == "CAP" ) {
        // add only current year consumer
        if( govtConsumer->year == mModeltime->getper_to_yr( aPeriod ) ) {
            // note the negative sign to get
            addToType( mModeltime->getper_to_yr( aPeriod ), "Govt Deficit",
                govtConsumer->expenditures[ aPeriod ].getValue( Expenditure::SAVINGS ) );
        }
    }
}

void SGMGenTable::startVisitTradeConsumer( const TradeConsumer* tradeConsumer,
                                          const int aPeriod )
{
    // net energy trade
    if( mName == "ETRADE" ) {
        // add only current year consumer
        if( tradeConsumer->getYear() == mModeltime->getper_to_yr( aPeriod ) ) {
            // get energy inputs only
            for( unsigned int i=0; i<tradeConsumer->mLeafInputs.size(); i++ ){
                if( tradeConsumer->mLeafInputs[ i ]->hasTypeFlag( IInput::ENERGY ) ){
                    addToType( mModeltime->getper_to_yr( aPeriod ), tradeConsumer->mLeafInputs[ i ]->getName(),
                        tradeConsumer->mLeafInputs[ i ]->getPhysicalDemand( aPeriod ) );
                }
            }
        }
    }
    else if( mName == "EmissBySource" ){
        // add or remove emissions only for the current consumer.
        if( tradeConsumer->getYear() != mModeltime->getper_to_yr( aPeriod ) ){
            return;
        }

        // Loop through the inputs and find primary goods.
        for( unsigned int i = 0; i < tradeConsumer->mLeafInputs.size(); ++i ){
            // Skip non-primary inputs
            if( !tradeConsumer->mLeafInputs[ i ]->hasTypeFlag( IInput::PRIMARY ) ){
                continue;
            }

            // Calculate the amount of emissions that are being traded.
            const double tradedEmissions = tradeConsumer->mLeafInputs[ i ]->getPhysicalDemand( aPeriod )
                * tradeConsumer->mLeafInputs[ i ]->getCO2EmissionsCoefficient( "CO2", aPeriod );
            // Add or remove the emissions to the column for the sector and year. Check that the sign is right.
            addToType( mModeltime->getper_to_yr( aPeriod ),
                tradeConsumer->mLeafInputs[ i ]->getName(), -1 * tradedEmissions );
        }
    }
}

void SGMGenTable::startVisitInvestConsumer( const InvestConsumer* investConsumer, const int aPeriod ) {
    // add only current year consumer
    if( investConsumer->year == mModeltime->getper_to_yr( aPeriod ) ) {
    }
}

// Write to SGM general table.
// This routine assumes that only and all operating technology vintages are passed in as
// an argument.
void SGMGenTable::startVisitProductionTechnology( const ProductionTechnology* prodTech,
                                                 const int aPeriod )
{
    if( aPeriod == -1 || ( prodTech->isAvailable( aPeriod ) && !prodTech->isRetired( aPeriod ) ) ) {
        // add output of each technology for each period
        if( ( mName == "CO2" ) || ( mName == "CO2bySec" ) || ( mName == "CO2byTech" ) ){
            unsigned int CO2index = util::searchForValue( prodTech->mGhgNameMap, string( "CO2" ) );
            double CO2Emiss = prodTech->mGhgs[CO2index]->getEmission( aPeriod );
            if( mName == "CO2" ) {
                addToType( mModeltime->getper_to_yr( aPeriod ), "CO2", CO2Emiss );
            }
            else if( mName == "CO2bySec" ) {
                addToType( mModeltime->getper_to_yr( aPeriod ), mCurrentSectorName, CO2Emiss );
            }
            else if( mName == "CO2byTech" ) {
                addToType( mModeltime->getper_to_yr( aPeriod ), prodTech->getName(), CO2Emiss );
            }
        }
        else if( mName == "EmissBySource" ){
            unsigned int CO2index = util::searchForValue( prodTech->mGhgNameMap, string( "CO2" ) );
            const double emissFuel = prodTech->mGhgs[ CO2index ]->getEmissFuel( aPeriod );
            // Only primary energy sectors will have getEmissFuel > 0.
            if( emissFuel > 0 ){
                addToType( mModeltime->getper_to_yr( aPeriod ), mCurrentSectorName, emissFuel );
            }
        }
        else if( mName == "PEC" ) {
            for( unsigned int i=0; i<prodTech->mLeafInputs.size(); i++ ){
                // get primary energy input only
                if( prodTech->mLeafInputs[ i ]->hasTypeFlag( IInput::PRIMARY ) ){
                    // isn't this the same as getting physical demand?
                    addToType( mModeltime->getper_to_yr( aPeriod ), prodTech->mLeafInputs[ i ]->getName(),
                        prodTech->mLeafInputs[ i ]->getPhysicalDemand( aPeriod ) );
                }
            }
            // special code to add renewable, nuclear and hydro electricity to primary energy consumption
            if( mCurrentSectorName == "ElectricityGeneration" ) {
                // TODO: use average fossil efficiency instead of hard-coded 0.3
                double fossilEfficiency = 0.3;
                if( prodTech->categoryName == "Renewable"){
                    addToType( mModeltime->getper_to_yr( aPeriod ), prodTech->name,
                        prodTech->getOutput( aPeriod ) / fossilEfficiency );
                }
            }
        }
        // primary energy production
        else if( mName == "PEP" ) {
            // get primary energy input only
            if( isPrimaryEnergyGood( mCurrentRegionName, mCurrentSectorName ) ){
                addToType( mModeltime->getper_to_yr( aPeriod ), mCurrentSectorName,
                    prodTech->getOutput( aPeriod ) );
            }
            // special code to add renewable, nuclear and hydro electricity to primary energy production
            if( mCurrentSectorName == "ElectricityGeneration" ) {
                // TODO: use average fossil efficiency instead of hard-coded 0.3
                double fossilEfficiency = 0.3;
                if( prodTech->categoryName == "Renewable"){
                    addToType( mModeltime->getper_to_yr( aPeriod ), prodTech->name,
                        prodTech->getOutput( aPeriod )
                        / fossilEfficiency );
                }
            }
        }
        // secondary energy production
        else if( mName == "SEP" ) {
            // get secondary energy goods only
            if( isSecondaryEnergyGood( mCurrentRegionName, mCurrentSectorName  ) ){
                addToType( mModeltime->getper_to_yr( aPeriod ), mCurrentSectorName,
                    prodTech->getOutput( aPeriod ) );
            }
        }
        // non-energy sector output
        else if( mName == "NEP" ) {
            // get non-energy goods only
            if( !isEnergyGood( mCurrentRegionName, mCurrentSectorName ) ){
                addToType( mModeltime->getper_to_yr( aPeriod ), mCurrentSectorName,
                    prodTech->getOutput( aPeriod ) );
            }
        }
        // electricity generation by technology
        else if( mName == "ELEC" ) {
            if( mCurrentSectorName == "ElectricityGeneration" ) {
                addToType( mModeltime->getper_to_yr( aPeriod ), prodTech->getName(),
                    prodTech->getOutput(aPeriod) );
            }
        }
        // fuel consumption for electricity generation
        else if( mName == "ElecFuel" ) {
            if( mCurrentSectorName == "ElectricityGeneration" ) {
                for( unsigned int i=0; i<prodTech->mLeafInputs.size(); i++ ){
                    // get energy input only
                    if( prodTech->mLeafInputs[ i ]->hasTypeFlag( IInput::ENERGY ) ){
                        addToType( mModeltime->getper_to_yr( aPeriod ),
                            prodTech->mLeafInputs[ i ]->getName(),
                            prodTech->mLeafInputs[ i ]->getPhysicalDemand( aPeriod ) );
                    }
                }
            }
        }
        // capital stock and other related output
        else if( mName == "CAP" ) {
            addToType( mModeltime->getper_to_yr( aPeriod ), "CapitalStock", prodTech->getCapitalStock() );
            /*
            addToType( mModeltime->getper_to_yr( aPeriod ), "CapitalStock/1000 Worker", prodTech->getCapitalStock() /
                scenario->getMarketplace()->getSupply( "Labor", mCurrentRegionName, aPeriod ) );
                */
            addToType( mModeltime->getper_to_yr( aPeriod ), "Profits", prodTech->mProfits[ aPeriod ] );
            addToType( mModeltime->getper_to_yr( aPeriod ), "Retained Earnings",
                prodTech->expenditures[ aPeriod ].getValue( Expenditure::RETAINED_EARNINGS ) );
        }
        // energy investments annual
        else if( mName == "EINV" ) {
            // get energy technologies only
            if( isEnergyGood( mCurrentRegionName, mCurrentSectorName ) ){
                addToType( mModeltime->getper_to_yr( aPeriod ), mCurrentSectorName,
                    prodTech->getAnnualInvestment( aPeriod ) );
            }
        }
        // non-energy investments annual
        else if( mName == "NEINV" ) {
            // get non-energy technologies only
            if( !isEnergyGood( mCurrentRegionName, mCurrentSectorName ) ){
                addToType( mModeltime->getper_to_yr( aPeriod ), mCurrentSectorName,
                    prodTech->getAnnualInvestment( aPeriod ) );
            }
        }
        // write out all passenger transportation sector results
        else if( (mName == "PASSTRAN") || (mName == "PASSTRANFC") || (mName == "PASSTRANFCM") ||
            (mName == "PASSTRANFCT")  || (mName == "PASSTRANMPG") || (mName == "PASSTRANCOST") ) {
                // get passenger transport technologies only that have non zero production
                if( (prodTech->categoryName == "PassTransport") && (prodTech->mOutputs[ 0 ]->getCurrencyOutput( aPeriod ) != 0) ){
                    double conversionFactor = FunctionUtils::getMarketConversionFactor( mCurrentRegionName, mCurrentSectorName );
                    if( mName == "PASSTRAN" ) {
                        addToType( mModeltime->getper_to_yr( aPeriod ), mCurrentSectorName, prodTech->mOutputs[ 0 ]->getPhysicalOutput( aPeriod ) );
                    }
                    else if ( mName == "PASSTRANCOST" ) {
                        addToType( mModeltime->getper_to_yr( aPeriod ), mCurrentSectorName, scenario->getMarketplace()->getPrice(mCurrentSectorName, mCurrentRegionName, aPeriod) * conversionFactor );
                    }
                    // for all other tables that require inputs
                    else {
                        for( unsigned int i=0; i<prodTech->mLeafInputs.size(); i++ ){
                            // get secondary energy input only
                            string inputName = prodTech->mLeafInputs[ i ]->getName();
                            // *** skip if not refined oil input ***
                            // this is problematic for other vehicles that do not use refined oil
                            if( inputName != "RefinedOil" ) {
                                continue;
                            }
                            double fuelConsumption = prodTech->mLeafInputs[ i ]->getPhysicalDemand( aPeriod );
                            // passenger transportation technology fuel consumption by fuel
                            if( mName == "PASSTRANFC" ) {
                                addToType( mModeltime->getper_to_yr( aPeriod ), inputName, fuelConsumption );
                            }
                            // passenger transportation technology fuel consumption by mode
                            else if( mName == "PASSTRANFCM" ) {
                                addToType( mModeltime->getper_to_yr( aPeriod ), mCurrentSectorName, fuelConsumption );
                            }
                            // passenger transportation technology fuel consumption by technology
                            else if( mName == "PASSTRANFCT" ) {
                                addToType( mModeltime->getper_to_yr( aPeriod ), prodTech->name, fuelConsumption );
                            }
                            // passenger transportation technology fuel economy
                            else if( mName == "PASSTRANMPG" ) {
                                double mpg = prodTech->mOutputs[ 0 ]->getPhysicalOutput( aPeriod ) / fuelConsumption;
                                addToType( mModeltime->getper_to_yr( aPeriod ), prodTech->name, mpg );
                            }
                        }
                    }
                }
            }
    }
}

/*! \brief Return whether a good is an energy good.
* \param aGoodName Good name.
* \return Whether the good is an energy price good.
*/
bool SGMGenTable::isEnergyGood( const string& aRegionName, const string& aGoodName ){
    const IInfo* marketInfo = scenario->getMarketplace()->getMarketInfo( aGoodName, aRegionName,
                                                                         0, false );
    return marketInfo && marketInfo->getBoolean( "IsEnergyGood", false );
}

/*! \brief Return whether a good is a primary energy good.
* \param aRegionName Region name.
* \param aGoodName Good name.
* \return Whether the good is a primary energy price good.
*/
bool SGMGenTable::isPrimaryEnergyGood( const string& aRegionName, const string& aGoodName ){
    const IInfo* marketInfo = scenario->getMarketplace()->getMarketInfo( aGoodName, aRegionName,
                                                                         0, false );
    return marketInfo && marketInfo->getBoolean( "IsPrimaryEnergyGood", false );
}

/*! \brief Return whether a good is a secondary energy good.
* \param aRegionName Region name.
* \return Whether the good is a secondary energy price good.
*/
bool SGMGenTable::isSecondaryEnergyGood( const string& aRegionName, const string& aGoodName ){
    const IInfo* marketInfo = scenario->getMarketplace()->getMarketInfo( aGoodName, aRegionName,
                                                                         0, false );
    return marketInfo && marketInfo->getBoolean( "IsSecondaryEnergyGood", false );
}
