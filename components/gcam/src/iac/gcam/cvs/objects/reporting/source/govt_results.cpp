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
* \file govt_results.cpp
* \ingroup Objects
* \brief The GovtResults class source file.
*
* \author Josh Lurz
*/

#include "util/base/include/definitions.h"

#include <string>
#include <vector>
#include "reporting/include/govt_results.h"
#include "technologies/include/production_technology.h"
#include "reporting/include/storage_table.h"
#include "sectors/include/production_sector.h"
#include "containers/include/scenario.h"
#include "util/base/include/model_time.h"
#include "marketplace/include/marketplace.h"
#include "consumers/include/govt_consumer.h"
#include "consumers/include/household_consumer.h"

using namespace std;

extern Scenario* scenario; // for marketplace. Remove if unneeded and includes.

//! Default Constructor
GovtResults::GovtResults( const string& aRegionName, ostream& aFile ):
mFile( aFile ),
mRegionName( aRegionName ),
mTaxReceipts( new StorageTable ),
mSubsidies( new StorageTable ),
mGovtExpenditures( new StorageTable ),
mGovtTransfers( 0 ),
mParsingGovt( false ){
}

/*! \brief Output the Government Sector Results table to a CSV
 * 
 * \author Josh Lurz
 */
void GovtResults::finish() const {
    mFile << "-----------------------------" << endl;
    mFile << "Government Sector Results Table" << endl;
    mFile << "-----------------------------" << endl << endl;
    
    // Block out the code which writes out tax receipts.
    {
        // Write out the tax receipts table.
        mFile << "Tax receipts by Sector" << endl;

        // Get the column labels for the taxes table which are types of taxes.
        const vector<string> colNames = mTaxReceipts->getColLabels();
        for( vector<string>::const_iterator col = colNames.begin(); col != colNames.end(); ++col ){
            mFile << ',' << *col;
        }   
        mFile << endl;

        // Write out each row of the table, representing one sector
        const vector<string> rowNames = mTaxReceipts->getRowLabels();
        for ( vector<string>::const_iterator row = rowNames.begin(); row != rowNames.end(); ++row ){
            mFile << *row;
            // Iterate through the columns.
            for( vector<string>::const_iterator col = colNames.begin(); col != colNames.end(); ++col ) {
                mFile << ',' << mTaxReceipts->getValue( *row, *col );
            }
            mFile << endl;
        }
        mFile << endl;
    }
    // Block out the code which writes the subsidies table.
    {
        // Write out the subsidies table.
        mFile << "Subsidies by Sector" << endl;
        
        mFile << "Sector";
        // Get the column labels for the subsidies table.
        const vector<string> colNames = mSubsidies->getColLabels();
        for( vector<string>::const_iterator col = colNames.begin(); col != colNames.end(); ++col ){
            mFile << ',' << *col;
        }   
        mFile << endl;

        // Write out each row of the table, representing one sector
        const vector<string> rowNames = mSubsidies->getRowLabels();
        for ( vector<string>::const_iterator row = rowNames.begin(); row != rowNames.end(); ++row ){
            mFile << *row;
            // Iterate through the columns.
            for( vector<string>::const_iterator col = colNames.begin(); col != colNames.end(); ++col ) {
                mFile << ',' << mSubsidies->getValue( *row, *col );
            }
            mFile << endl;
        }
        mFile << endl;
    }
    
    // Write out total transfers.
    mFile << "Government Transfers, " << mGovtTransfers << endl << endl;

    // Write out expenditures.
}

//! Update the CGE Region
void GovtResults::startVisitRegionCGE( const RegionCGE* aRegionCGE, const int aPeriod ){
    // Add the column labels to all the tables.
    mTaxReceipts->addColumn( "Proportional Tax" ); // what is the AND in Legacy? Additive taxes?
    mTaxReceipts->addColumn( "Social Security Tax" );
    mTaxReceipts->addColumn( "Corporate Tax" ); // whats the "And" in Legacy?
    mTaxReceipts->addColumn( "IBT" );
    mTaxReceipts->addColumn( "Carbon" );
    mTaxReceipts->addColumn( "ITC" );

    mSubsidies->addColumn( "Amount" );

    mGovtExpenditures->addColumn( "Price" );
    mGovtExpenditures->addColumn( "Physical Quantity" );
    mGovtExpenditures->addColumn( "Amount" );
}

void GovtResults::startVisitSector( const Sector* aSector, const int aPeriod ){    
    // Store the sector name as we'll need it at the technology level.
    // This has to be done here instead of updateProductionSector because that is called
    // after all the technologies are updated.
    mCurrSectorName = aSector->getName();
}

void GovtResults::startVisitProductionTechnology( const ProductionTechnology* aProdTechnology, 
                                                 const int aPeriod )
{
    if( aPeriod == -1 || ( aProdTechnology->isAvailable( aPeriod ) && 
        !aProdTechnology->isRetired( aPeriod ) ) ) {
            mParsingGovt = false;
            // Fill up the tax table.
            // This isn't complete yet, its zero anyway right now.
            mTaxReceipts->addToType( mCurrSectorName, "Proportional Tax", 0 );

            mTaxReceipts->addToType( mCurrSectorName, "Social Security Tax",
                aProdTechnology->expenditures[ aPeriod ].getValue( Expenditure::SOCIAL_SECURITY_TAX ) );

            // this would be wrong if more was added to DIRECT_TAXES. 
            mTaxReceipts->addToType( mCurrSectorName, "Corporate Tax",
                aProdTechnology->expenditures[ aPeriod ].getValue( Expenditure::DIRECT_TAXES ) );

            // this would be wrong if more was added to indirect taxes.
            mTaxReceipts->addToType( mCurrSectorName, "IBT",
                aProdTechnology->expenditures[ aPeriod ].getValue( Expenditure::INDIRECT_TAXES ) );

            // these two aren't finished yet.
            mTaxReceipts->addToType( mCurrSectorName, "Carbon", 0 );
            mTaxReceipts->addToType( mCurrSectorName, "ITC", 0 );

            // Fill up the subsidy table. This isn't done per sector in the model yet.
            mSubsidies->addToType( mCurrSectorName, "Amount", 0 );
        }
}    

void GovtResults::startVisitGovtConsumer( const GovtConsumer* aGovtConsumer, const int aPeriod ){
    if( aGovtConsumer->getYear() == scenario->getModeltime()->getper_to_yr( aPeriod ) ) {
        // Set the total transfers.
        mGovtTransfers = aGovtConsumer->expenditures[ aPeriod ].getValue( Expenditure::TRANSFERS );

        // Need to setup the expenditure table.
        mParsingGovt = true;
    }
}

void GovtResults::startVisitHouseholdConsumer( const HouseholdConsumer* aHouseholdConsumer, const int aPeriod ){
    if( aHouseholdConsumer->getYear() == scenario->getModeltime()->getper_to_yr( aPeriod ) ) {
        // Need to get social security tax added.
        mTaxReceipts->addToType( mCurrSectorName, "Social Security Tax",
            aHouseholdConsumer->expenditures[ aPeriod ].getValue( Expenditure::SOCIAL_SECURITY_TAX ) );

        // Also corporate income tax?
        mTaxReceipts->addToType( mCurrSectorName, "Corporate Tax",
            aHouseholdConsumer->expenditures[ aPeriod ].getValue( Expenditure::DIRECT_TAXES ) ); // this might be wrong.
    }
}
