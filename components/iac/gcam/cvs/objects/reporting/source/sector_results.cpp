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
* \file sector_results.cpp
* \ingroup Objects
* \brief The SectorResults class source file.
*
* \author Josh Lurz
*/

#include "util/base/include/definitions.h"

#include <string>
#include <vector>
#include "reporting/include/sector_results.h"
#include "technologies/include/production_technology.h"
#include "reporting/include/storage_table.h"
#include "sectors/include/production_sector.h"
#include "containers/include/scenario.h"
#include "marketplace/include/marketplace.h"
#include "functions/include/function_utils.h"
#include "technologies/include/ioutput.h"

using namespace std;

extern Scenario* scenario; // for marketplace

//! Default Constructor
SectorResults::SectorResults( const string& aRegionName, ostream& aFile ):
mCurrentRegionName( aRegionName ),
mFile( aFile ),
mInternalTable( new StorageTable ){
}

/*! \brief Output the Sector Results table to a CSV
 * 
 * \author Josh Lurz
 */
void SectorResults::finish() const {
    mFile << "-----------------------------" << endl;
    mFile << "Regional Sector Results Table" << endl;
    mFile << "-----------------------------" << endl << endl;
    
    // The table is stored inverted so that totals are easily calculated.
    // Use row labels to print the columns.
    // Get the column labels which are the types of outputs.
    mFile << "Sector";
    const vector<string> colNames = mInternalTable->getRowLabels();
    for( vector<string>::const_iterator col = colNames.begin(); col != colNames.end(); ++col ){
        mFile << ',' << *col;
    }   
    mFile << endl;

    //aFile.precision(0);
    // Write out each row of the table, representing one sector
    const vector<string> rowNames = mInternalTable->getColLabels();
    for ( vector<string>::const_iterator row = rowNames.begin(); row != rowNames.end(); ++row ){
        mFile << *row;
        // Iterate through the columns.
        for( vector<string>::const_iterator col = colNames.begin(); col != colNames.end(); ++col ) {
            mFile << ',' << mInternalTable->getValue( *col, *row );
        }
        mFile << endl;
    }
    mFile << endl;
}
void SectorResults::startVisitSector( const Sector* aSector, const int aPeriod ){    
    // Store the sector name as we'll need it at the technology level.
    // This has to be done here instead of updateProductionSector because that is called
    // after all the technologies are updated.
    mCurrentSectorName = aSector->getName();
}

void SectorResults::startVisitProductionSector( const ProductionSector* aProdSector, const int aPeriod ) {
    // Add a column for ourselves.
    mInternalTable->addColumn( aProdSector->getName() );
}

void SectorResults::startVisitProductionTechnology( const ProductionTechnology* aProdTechnology,
                                                const int aPeriod ) 
{
    if( aPeriod == -1 || ( aProdTechnology->isAvailable( aPeriod ) &&
            !aProdTechnology->isRetired( aPeriod ) ) ) {
        // Add to the physical output column.
        mInternalTable->addToType( "Phys. Out.", mCurrentSectorName, aProdTechnology->mOutputs[ 0 ]->getPhysicalOutput( aPeriod ) );

        // Add to the sales column
        mInternalTable->addToType( "Sales", mCurrentSectorName, aProdTechnology->expenditures[ aPeriod ].getValue( Expenditure::SALES ) );

        double priceReceived = FunctionUtils::getPriceReceived( mCurrentRegionName, mCurrentSectorName, aPeriod ); 
        // Add to the revenue column.
        mInternalTable->addToType( "Revenue", mCurrentSectorName,
                                    aProdTechnology->expenditures[ aPeriod ].getValue( Expenditure::SALES ) * priceReceived );

        // Add to the profits column.
        mInternalTable->addToType( "Profits", mCurrentSectorName, aProdTechnology->mProfits[ aPeriod ] );

        // Add to the retained earnings column.
        mInternalTable->addToType( "Retained Earnings", mCurrentSectorName,
                                   aProdTechnology->expenditures[ aPeriod ].getValue( Expenditure::RETAINED_EARNINGS ) );

        // Add to the costs column.
        mInternalTable->addToType( "Costs", mCurrentSectorName,
                                   aProdTechnology->mCostsReporting[ aPeriod ] );
    }
}
