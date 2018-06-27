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
* \file social_accounting_matrix.cpp
* \ingroup Objects
* \brief The SocialAccountingMatrix class source file.
*
* \author Praelit Patel
*/

#include "util/base/include/definitions.h"
#include <iostream>
#include <string>

#include "reporting/include/social_accounting_matrix.h"

#include "containers/include/scenario.h"
#include "util/base/include/model_time.h"
#include "marketplace/include/marketplace.h"
#include "containers/include/region.h"
#include "sectors/include/sector.h"
#include "sectors/include/subsector.h"
#include "technologies/include/base_technology.h"
#include "consumers/include/consumer.h"
#include "consumers/include/household_consumer.h"
#include "consumers/include/govt_consumer.h"
#include "consumers/include/trade_consumer.h"
#include "consumers/include/invest_consumer.h"
#include "technologies/include/production_technology.h"
#include "functions/include/iinput.h"
#include "sectors/include/factor_supply.h"
#include "reporting/include/storage_table.h"
#include "functions/include/function_utils.h"

using namespace std;
extern Scenario* scenario;

//! Default Constructor
SocialAccountingMatrix::SocialAccountingMatrix( const string& aRegionName, ostream& aFile ):
mFile( aFile ), mRegionName( aRegionName ), mTable( new StorageTable ){
}

/*! \brief For outputting SGM data to a flat csv File
 * 
 * \author Pralit Patel
 */
void SocialAccountingMatrix::finish() const {
    mFile << "Social Accounting Matrix" << endl
          << "Receipts" << ',' << "Activites" << ',' << "Commoditites" << ',' << "Land" << ',' 
          << "Labor" << ',' << "Capital" << ',' << "Households" << ',' << "Enterprises" << ',' 
          << "Government" << ','<< "Carbon Tax" << ',' << "Capital Account" << ','
          << "Rest Of World" << ',' << "Totals" << endl;

    for ( int categoryIndex = ACTIVITIES; categoryIndex <= TOTALS; categoryIndex++ ) {
        mTable->addColumn( getString( CategoryType( categoryIndex ) ) );
        mFile << getString( CategoryType( categoryIndex ) );
        for( int sectorIndex = ACTIVITIES; sectorIndex <= TOTALS; sectorIndex++ ) {
            if( !( sectorIndex == TOTALS && categoryIndex == TOTALS ) ) {
                mFile << ',' << mTable->getValue( getString( CategoryType ( sectorIndex ) ),
                    getString( CategoryType( categoryIndex ) ) );
            }
        }
        mFile << endl;
    }
    mFile << endl;
}

void SocialAccountingMatrix::startVisitHouseholdConsumer( const HouseholdConsumer* householdConsumer,
                                                      const int aPeriod ) 
{
    if( householdConsumer->year == scenario->getModeltime()->getper_to_yr( aPeriod ) ){
        // for household column of SAM
        addToType( HOUSEHOLDS, COMMODITIES,
            householdConsumer->expenditures[ aPeriod ].getValue( Expenditure::CONSUMPTION ) );
        addToType( HOUSEHOLDS, GOVERNMENT,
            householdConsumer->expenditures[ aPeriod ].getValue( Expenditure::SOCIAL_SECURITY_TAX ) );
        addToType( HOUSEHOLDS, GOVERNMENT,
            householdConsumer->expenditures[ aPeriod ].getValue( Expenditure::DIRECT_TAXES ) );
        addToType( HOUSEHOLDS, CAPITAL_ACCOUNT,
            householdConsumer->expenditures[ aPeriod ].getValue( Expenditure::SAVINGS ) );
    }
}

void SocialAccountingMatrix::startVisitGovtConsumer( const GovtConsumer* govtConsumer,
                                                 const int aPeriod )
{
    if( govtConsumer->year == scenario->getModeltime()->getper_to_yr( aPeriod ) ){
        // for government column of SAM
        addToType( GOVERNMENT, ACTIVITIES,
            govtConsumer->expenditures[ aPeriod ].getValue( Expenditure::SUBSIDY ) );
        addToType( GOVERNMENT, COMMODITIES,
            govtConsumer->expenditures[ aPeriod ].getValue( Expenditure::CONSUMPTION ) );
        // we are spliting out carbon tax into its own category
        addToType( GOVERNMENT, HOUSEHOLDS,
            govtConsumer->expenditures[ aPeriod ].getValue( Expenditure::TRANSFERS )
            - govtConsumer->expenditures[ aPeriod ].getValue( Expenditure::CARBON_TAX ) );
        addToType( CARBON_TAX, HOUSEHOLDS,
            govtConsumer->expenditures[ aPeriod ].getValue( Expenditure::CARBON_TAX ) );
        addToType( GOVERNMENT, CAPITAL_ACCOUNT,
            govtConsumer->expenditures[ aPeriod ].getValue( Expenditure::SAVINGS ) );   
    }
}

void SocialAccountingMatrix::startVisitTradeConsumer( const TradeConsumer* tradeConsumer, const int aPeriod )
{
    if ( tradeConsumer->year == scenario->getModeltime()->getper_to_yr( aPeriod ) ) {
        // for net export or rest of world column of SAM
        addToType( REST_OF_WORLD, ACTIVITIES,
            tradeConsumer->expenditures[ aPeriod ].getValue( Expenditure::TOTAL_IMPORTS ) );

        // we have to get the leaves again because they would have been cleared by now
        // maybe we should store this value in the expenditure somewhere
        IInput* capInput = FunctionUtils::getCapitalInput( FunctionUtils::getLeafInputs( tradeConsumer->mNestedInputRoot ) );
        // add the capital flow
        addToType( REST_OF_WORLD, CAPITAL_ACCOUNT,
            capInput->getCurrencyDemand( aPeriod ) );
    }
}

void SocialAccountingMatrix::startVisitInvestConsumer( const InvestConsumer* investConsumer, 
                                                   const int aPeriod )
{
    if( investConsumer->year == scenario->getModeltime()->getper_to_yr( aPeriod ) ){
        // for capital account column of SAM
        addToType( CAPITAL_ACCOUNT, COMMODITIES,
            investConsumer->expenditures[ aPeriod ].getValue( Expenditure::INVESTMENT ) );
    }
}

void SocialAccountingMatrix::startVisitProductionTechnology( const ProductionTechnology* prodTech,
                                                         const int aPeriod ) 
{
    if( aPeriod == -1 || ( prodTech->isAvailable( aPeriod ) && !prodTech->isRetired( aPeriod ) ) ) {
        // for activities column of SAM
        // split out carbon tax
        addToType( ACTIVITIES, COMMODITIES,
            prodTech->expenditures[ aPeriod ].getValue( Expenditure::INTERMEDIATE_INPUTS )
            - prodTech->expenditures[ aPeriod ].getValue( Expenditure::CARBON_TAX ) );
        addToType( ACTIVITIES, CARBON_TAX, prodTech->expenditures[ aPeriod ].getValue( Expenditure::CARBON_TAX ) );
        addToType( ACTIVITIES, FACTORS_LABOR, prodTech->expenditures[ aPeriod ].getValue( Expenditure::WAGES ) );
        addToType( ACTIVITIES, FACTORS_LAND, prodTech->expenditures[ aPeriod ].getValue( Expenditure::LAND_RENTS ) );
        addToType( ACTIVITIES, FACTORS_CAPITAL, prodTech->expenditures[ aPeriod ].getValue( Expenditure::RENTALS ) );
        addToType( FACTORS_CAPITAL, ENTERPRISES, prodTech->expenditures[ aPeriod ].getValue( Expenditure::RENTALS ) );
        addToType( ACTIVITIES, GOVERNMENT, prodTech->expenditures[ aPeriod ].getValue( Expenditure::INDIRECT_TAXES ) );

        // for commodities column of SAM
        addToType( COMMODITIES, ACTIVITIES, prodTech->expenditures[ aPeriod ].getValue( Expenditure::SALES ) );
        addToType( COMMODITIES, GOVERNMENT, prodTech->expenditures[ aPeriod ].getValue( Expenditure::TARIFFS ) );
        addToType( COMMODITIES, REST_OF_WORLD, prodTech->expenditures[ aPeriod ].getValue( Expenditure::IMPORTS ) );

        // for enterprise column of SAM
        addToType( ENTERPRISES, HOUSEHOLDS, prodTech->expenditures[ aPeriod ].getValue( Expenditure::DIVIDENDS ) );
        addToType( ENTERPRISES, GOVERNMENT, prodTech->expenditures[ aPeriod ].getValue( Expenditure::DIRECT_TAXES ) );
        addToType( ENTERPRISES, CAPITAL_ACCOUNT,
            prodTech->expenditures[ aPeriod ].getValue( Expenditure::RETAINED_EARNINGS ) );
    }
}

void SocialAccountingMatrix::startVisitFactorSupply( const FactorSupply* factorSupply,
                                                 const int period )
{   
    const Marketplace* marketplace = scenario->getMarketplace();
    double tempSupply = factorSupply->getSupply( mRegionName, period )
        * marketplace->getPrice( factorSupply->getName(), mRegionName, period );

    if( factorSupply->getName() == "Land" ) {
        // for land column of SAM
        addToType( FACTORS_LAND, HOUSEHOLDS, tempSupply );
    }
    else if( factorSupply->getName() == "Labor" ) {
        // for labor column of SAM
        addToType( FACTORS_LABOR, HOUSEHOLDS, tempSupply );
    }
    // for factors_capital column see startVisitProductionTechnology
    // Is this right?
}

//! Helper function which converts enums to strings and writes to the internal table.
void SocialAccountingMatrix::addToType( CategoryType aRowCat, CategoryType aColCat, double aValue ){
    // TODO: I think rows a columns are switched somewhere but everything turns
    // out right the way it is now
    assert( mTable.get() );
    mTable->addToType( getString( aRowCat ), getString( aColCat ), aValue );

    // need to sum rows as well
    const string total = getString( TOTALS );
    mTable->addToType( total, getString( aColCat ), aValue );
}

const string& SocialAccountingMatrix::getString( const CategoryType aType ) {
    // Create a static array of labels in order of their enum.
    const static string labels[] = { "Activities", "Commodities", "Factors: Land", "Factors: Labor",
                                     "Factors: Capital", "Households", "Enterprises", "Government",
                                     "Carbon Tax", "Capital Account", "Rest of world", "Total" };
    // Return the correct label based on the enum value.
    return labels[ aType ];
}
