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
* \file demand_components_table.cpp
* \ingroup Objects
* \brief The DemandComponentsTable source file.
*
*  Detailed description.
*
* \author Pralit Patel
* \author Sonny Kim
* \author Josh Lurz
*/
#include "util/base/include/definitions.h"
#include "reporting/include/demand_components_table.h"
#include "containers/include/region_cge.h"
#include "sectors/include/factor_supply.h"
#include "consumers/include/household_consumer.h"
#include "consumers/include/govt_consumer.h"
#include "consumers/include/trade_consumer.h"
#include "consumers/include/invest_consumer.h"
#include "technologies/include/production_technology.h"
#include "functions/include/iinput.h"
#include "util/base/include/model_time.h"
#include "containers/include/scenario.h"
#include "reporting/include/storage_table.h"
#include "sectors/include/sector.h"

using namespace std;

extern Scenario* scenario;

//! Default Constructor
DemandComponentsTable::DemandComponentsTable( ostream& aFile ):
mFile( aFile ),
mTable( new StorageTable ){
}

/*! \brief For outputing SGM data to a flat csv File
 */
void DemandComponentsTable::finish() const {
    
    mFile << "Demand Components Table" << endl << "Industry" << ',';
    // Now output the table column labels.
    const vector<string> cols = mTable->getColLabels();
    for( vector<string>::const_iterator col = cols.begin(); col != cols.end(); ++col ){
        mFile << *col <<','; // output the label.
    }
    mFile << endl;

    // Note: This is structurally different from SAM. This goes through the rows and prints
    // out each of the category values.
    mFile.precision(0);

    // Get the row labels.
    const vector<string> rows = mTable->getRowLabels();

    // Loop through each row and print all the columns.
    for( vector<string>::const_iterator row = rows.begin(); row != rows.end(); ++row ){
        mFile << *row; // output the row label.
        for( vector<string>::const_iterator col = cols.begin(); col != cols.end(); ++col ){
            mFile << ',' << mTable->getValue( *row, *col ); // output the value.
        }
        mFile << endl;
    }
    mFile << endl;
    // reset format to default
    mFile.precision(3);
}

void DemandComponentsTable::startVisitRegionCGE( const RegionCGE* regionCGE, const int aPeriod ) {
    // Add columns to the table.
    for( int i = TOTAL; i <= TRADE; ++i ){
        mTable->addColumn( getLabel( CategoryType( i ) ) );
    }
    // Add the rows in the right order.
    // This has to be done explicitly because the rows are added later by input, so the
    // ordering would be wrong.
    for( unsigned int i = 0; i < regionCGE->supplySector.size(); ++i ){
        mTable->addToType( regionCGE->supplySector[ i ]->getName(), "Trade", 0 ); // Column doesn't matter.
    }
    // Add the factor supplies.
    for( unsigned int i = 0; i < regionCGE->factorSupply.size(); ++i ){
        mTable->addToType( regionCGE->factorSupply[ i ]->getName(), "Trade", 0 );
    }
}

void DemandComponentsTable::startVisitHouseholdConsumer( const HouseholdConsumer* householdConsumer,
                                                     const int aPeriod )
{
    // add only current year consumer
    if( householdConsumer->year == scenario->getModeltime()->getper_to_yr( aPeriod ) ) {
        for( unsigned int i = 0; i < householdConsumer->mLeafInputs.size(); ++i ){
            mTable->addToType( householdConsumer->mLeafInputs[ i ]->getName(), getLabel( CONSUMPTION ),
                               householdConsumer->mLeafInputs[ i ]->getCurrencyDemand( aPeriod ) );
        }
    }
}

void DemandComponentsTable::startVisitGovtConsumer( const GovtConsumer* govtConsumer,
                                                const int aPeriod )
{
    // add only current year consumer
    if( govtConsumer->year == scenario->getModeltime()->getper_to_yr( aPeriod ) ) {
        for( unsigned int i = 0; i < govtConsumer->mLeafInputs.size(); ++i ){
            mTable->addToType( govtConsumer->mLeafInputs[ i ]->getName(), getLabel( GOVERNMENT ),
                               govtConsumer->mLeafInputs[ i ]->getCurrencyDemand( aPeriod ) );
        }
    }
}

void DemandComponentsTable::startVisitTradeConsumer( const TradeConsumer* tradeConsumer, const int aPeriod )
{
    // add only current year consumer
    if( tradeConsumer->year == scenario->getModeltime()->getper_to_yr( aPeriod ) ) {
        for( unsigned int i = 0; i< tradeConsumer->mLeafInputs.size(); i++ ){
            mTable->addToType( tradeConsumer->mLeafInputs[ i ]->getName(), getLabel( TRADE ),
                               tradeConsumer->mLeafInputs[ i ]->getCurrencyDemand( aPeriod ) );
        }
    }
}

void DemandComponentsTable::startVisitInvestConsumer( const InvestConsumer* investConsumer,
                                                  const int aPeriod )
{
    // add only current year consumer
    if( investConsumer->year == scenario->getModeltime()->getper_to_yr( aPeriod ) ) {
        for( unsigned int i = 0; i < investConsumer->mLeafInputs.size(); ++i ){
            mTable->addToType( investConsumer->mLeafInputs[ i ]->getName(), getLabel( INVESTMENT ),
                               investConsumer->mLeafInputs[ i ]->getCurrencyDemand( aPeriod ) );
        }
    }
}

void DemandComponentsTable::startVisitProductionTechnology( const ProductionTechnology* prodTech,
                                                        const int aPeriod ) 
{
    if( aPeriod == -1 || ( prodTech->isAvailable( aPeriod ) &&
            !prodTech->isRetired( aPeriod ) ) ) {
        for( unsigned int i=0; i<prodTech->mLeafInputs.size(); i++ ){
            // if capital, add only current vintage demand
            if( prodTech->mLeafInputs[i]->hasTypeFlag( IInput::CAPITAL ) ) {
                if( prodTech->getYear() == scenario->getModeltime()->getper_to_yr( aPeriod ) ) {
                    mTable->addToType( prodTech->mLeafInputs[ i ]->getName(), getLabel( INTERMED ),
                                       prodTech->mLeafInputs[ i ]->getCurrencyDemand( aPeriod ) );
                }
            }
            // for all other inputs, add all vintages
            else {
                mTable->addToType( prodTech->mLeafInputs[ i ]->getName(), getLabel( INTERMED ),
                                   prodTech->mLeafInputs[ i ]->getCurrencyDemand( aPeriod ) );
            }
        }
    }
}

/*! \brief Determine the category label based on the enum type value.
* \param aType The type of category.
* \return The label for the category. 
*/
const string& DemandComponentsTable::getLabel( const CategoryType aType ) const {
    // Setup a static array with the labels in the correct positions.
    // This will only be done on the first call to the function.
    const static string labels[] = 
    { "Total", "Intermediate Production", "Consumption", "Investment", "Government", "Trade" };

    // Return the correct label.
    return labels[ aType ];
}
