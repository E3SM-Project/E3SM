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
* \file input_output_table.cpp
* \ingroup Objects
* \brief The InputOutputTable class source file.
*
* \author Josh Lurz
*/

#include "util/base/include/definitions.h"

#include <string>
#include <vector>

#include "reporting/include/input_output_table.h"
#include "functions/include/iinput.h"
#include "functions/include/demand_input.h"
#include "functions/include/production_input.h"
#include "technologies/include/production_technology.h"
#include "technologies/include/ioutput.h"
#include "reporting/include/storage_table.h"
#include "sectors/include/sector.h"
#include "sectors/include/factor_supply.h"
#include "containers/include/region_cge.h"
#include "containers/include/scenario.h" // only for modeltime
#include "util/base/include/model_time.h"
#include "consumers/include/consumer.h"

extern Scenario* scenario;

using namespace std;

//! Default Constructor
InputOutputTable::InputOutputTable( const string& aRegionName, ostream& aFile ):
mFile( aFile ),
mRegionName( aRegionName ),
mInternalTable( new StorageTable ),
mParsingConsumer( false ),
mUseInput( false ){
}

/*! \brief Output the IOTable to a CSV
 * 
 * \author Josh Lurz
 */
void InputOutputTable::finish() const {
    mFile << "-----------------------------" << endl;
    mFile << "Regional Input Output Table" << endl;
    mFile << "-----------------------------" << endl << endl;

    // Get the column labels which are sector names.
    const vector<string> colNames = mInternalTable->getColLabels();
    for( vector<string>::const_iterator col = colNames.begin(); col != colNames.end(); ++col ){
        mFile << ',' << *col;
    }   
    mFile << endl;

    // Write out each row of the IOTable.
    const vector<string> rowNames = mInternalTable->getRowLabels();
    for ( vector<string>::const_iterator row = rowNames.begin(); row != rowNames.end(); ++row ){
        mFile << *row;
        // Iterate through the columns.
        for( vector<string>::const_iterator col = colNames.begin(); col != colNames.end(); ++col ) {
            mFile << ',' << mInternalTable->getValue( *row, *col );
        }
        mFile << endl;
    }
    mFile << endl;
}

//! Update the region. 
void InputOutputTable::startVisitRegionCGE( const RegionCGE* aRegion, const int aPeriod ){
    // Loop through the sectors and add a blank item for each just to set the ordering.
    for( unsigned int i = 0; i < aRegion->supplySector.size(); ++i ){
        mInternalTable->setType( aRegion->supplySector[ i ]->getName(), aRegion->supplySector[ i ]->getName(), 0 );
    }
    // Add factor supplies at the end. right?
    for( unsigned int i = 0; i < aRegion->factorSupply.size(); ++i ){
        mInternalTable->setType( aRegion->factorSupply[ i ]->getName(), aRegion->factorSupply[ i ]->getName(), 0 );
    }
}

void InputOutputTable::startVisitSector( const Sector* sector, const int aPeriod ) {
    // Add a column for ourselves.
    mInternalTable->addColumn( sector->getName() );

    // Store the sector name as we'll need it at the technology level.
    mCurrSectorName = sector->getName();
}

void InputOutputTable::startVisitProductionTechnology( const ProductionTechnology* prodTechnology,
                                                   const int aPeriod )
{
    if( aPeriod == -1 || ( prodTechnology->isAvailable( aPeriod ) && 
        !prodTechnology->isRetired( aPeriod ) ) ) {
            mUseInput = true;
            // Add the technologies output as a negative demand on the diagonal.
            mInternalTable->addToType( mCurrSectorName, mCurrSectorName, -1 * prodTechnology->mOutputs[ 0 ]->getCurrencyOutput( aPeriod ) );
            if( prodTechnology->isNewInvestment( aPeriod ) ) {
                mInternalTable->addToType( "Capital", mCurrSectorName, prodTechnology->mAnnualInvestment );
            }

            // Everything else will be updated at the input level.
            mParsingConsumer = false; // set that we aren't currently parsing a consumer.
        }
    else {
        mUseInput = false;
    }
}

void InputOutputTable::startVisitFactorSupply( const FactorSupply* factorSupply, const int aPeriod ){
    // Add factor supplies to the household column.
    mInternalTable->addToType( factorSupply->getName(), "Household", 
        -1 * factorSupply->getSupply( mRegionName, aPeriod ) );
}

//! Update the inputs contribution to the IOTable.
void InputOutputTable::startVisitProductionInput( const ProductionInput* aProdInput, 
                                                 const int aPeriod )
{
    if( mUseInput ) {
        // Add the currency demand.
        // The capital row is not truly capital but other value added.
        // The row is only the OVA row in ProductionSectors, it behaves as capital in Consumers.
        if( aProdInput->hasTypeFlag( IInput::CAPITAL ) && !mParsingConsumer ){
            mInternalTable->addToType( "OVA", mCurrSectorName, 
                aProdInput->getCurrencyDemand( aPeriod ) );
        }
        else {
            mInternalTable->addToType( aProdInput->getName(), mCurrSectorName, 
                aProdInput->getCurrencyDemand( aPeriod ) );
        }
    }
}

//! Update the inputs contribution to the IOTable.
void InputOutputTable::startVisitDemandInput( const DemandInput* aDemandInput, const int aPeriod ){
    if( mUseInput ) {
        mInternalTable->addToType( aDemandInput->getName(), mCurrSectorName, 
            aDemandInput->getCurrencyDemand( aPeriod ) );
    }
}

//! Update the consumer. Set that the current state is parsing a consumer, not a production tech.
void InputOutputTable::startVisitConsumer( const Consumer* aConsumer, const int aPeriod ){
    if( aConsumer->getYear() == scenario->getModeltime()->getper_to_yr( aPeriod ) ) {
        mParsingConsumer = true;
        mUseInput = true;
    }
    else {
        mUseInput = false;
    }
}
