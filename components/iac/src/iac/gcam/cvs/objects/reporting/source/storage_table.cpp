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
* \file storage_table.cpp
* \ingroup Objects
* \brief The StorageTable class source file.
*
* \author Josh Lurz
*/

#include "util/base/include/definitions.h"
#include <string>
#include <vector>
#include <algorithm>
#include <iostream> // until we get logger.
#include "reporting/include/storage_table.h"
#include "util/base/include/util.h"

using namespace std;

//! Default Constructor
StorageTable::StorageTable(){
}

//! Clear the data in the table. 
// This does not clear column labels.
void StorageTable::clear() {
    mInternalTable.rows.clear();
}

//! Return if the table has any rows.
bool StorageTable::isEmpty() const {
    return mInternalTable.rows.empty();
}

/*! \brief Add a column label to the list of columns.
* \note This will only add to the column if it is unique.
*/
void StorageTable::addColumn( const string& aColLabel ){
    // Check if the label is unique.
    if( find( mColLabels.begin(), mColLabels.end(), aColLabel ) == mColLabels.end() ){
        // Add the label.
        mColLabels.push_back( aColLabel );
    }
    else {
        // When logging is added this could be at DEBUG level.
        // Currently it is too common to print.
        // cout << aColLabel << " is not unique. Not adding column." << endl;
    }
}

//! Add to a type which has a year as a row.
void StorageTable::addToType( const int aRow, const string& aCol, const double aValue ){
    // Convert the integer row to a string and call the standard addToType.
    addToType( util::toString( aRow ), aCol, aValue );
}

//! Add to the value for the table specified by the account type key. 
void StorageTable::addToType( const string& aRow, const string& aCol, const double aValue ){
    // Find the correct column.
    int rowIndex = getRowIndex( aRow );
    
    // If the row does not exist insert one on the end.
    if( rowIndex == NO_ITEM_FOUND ){
        mInternalTable.rows.push_back( Row( aRow ) );
        // Set the row index to the new row so we can add the value.
        rowIndex = static_cast<int>( mInternalTable.rows.size() ) - 1;
    }

    // Find the column index.
    int colIndex = mInternalTable.rows[ rowIndex ].getColIndex( aCol );

    // If the column does not exist add it onto the end of the row.
    if( colIndex == NO_ITEM_FOUND ){
        mInternalTable.rows[ rowIndex ].data.push_back( Item( aCol ) );
        // Set the col index to the new col so we can add the value.
        colIndex = static_cast<int>( mInternalTable.rows[ rowIndex ].data.size() ) - 1;
    }
    // Add the value to the correct position.
    mInternalTable.rows[ rowIndex ].data[ colIndex ].value += aValue;
    // Add to the total.
    mInternalTable.rows[ rowIndex ].total += aValue;
}

//! set the value for the table specified by the account type key.
// This function should be replaced with a clear row and then use addToType.
// Right now it is massive copy-paste.
void StorageTable::setType( const string& aRow, const string& aCol,
                              const double aValue )
{
 	// Find the correct column.
    int rowIndex = getRowIndex( aRow );
    
    // If the row does not exist insert one on the end.
    if( rowIndex == NO_ITEM_FOUND ){
        mInternalTable.rows.push_back( Row( aRow ) );
        // Set the row index to the new row so we can add the value.
        rowIndex = static_cast<int>( mInternalTable.rows.size() ) - 1;
    }

    // Find the column index.
    int colIndex = mInternalTable.rows[ rowIndex ].getColIndex( aCol );

    // If the column does not exist add it onto the end of the row.
    if( colIndex == NO_ITEM_FOUND ){
        mInternalTable.rows[ rowIndex ].data.push_back( Item( aCol ) );
        // Set the col index to the new col so we can add the value.
        colIndex = static_cast<int>( mInternalTable.rows[ rowIndex ].data.size() ) - 1;
    }
    // Add the value to the correct position.
    mInternalTable.rows[ rowIndex ].data[ colIndex ].value = aValue;
}

//! Get the value for the table specified by the account type key. 
double StorageTable::getValue( const string& aRow, const string& aCol ) const {
    const int rowIndex = getRowIndex( aRow );
    if( rowIndex != NO_ITEM_FOUND ){
        // Special case total here.
        if( aCol == "Total" ){
            return mInternalTable.rows[ rowIndex ].total;
        }
        // Find the correct column.
        const int colIndex = mInternalTable.rows[ rowIndex ].getColIndex( aCol );
        if( colIndex != NO_ITEM_FOUND ){
            return mInternalTable.rows[ rowIndex ].data[ colIndex ].value;
        }
    }
    // Is this an error?
    return 0;
}

//! Get the list of all row labels in order.
const vector<string> StorageTable::getRowLabels() const {
    vector<string> rowLabels;

    // Loop through the rows and add the label for each to the vector.
    for( unsigned int row = 0; row < mInternalTable.rows.size(); ++row ){
        rowLabels.push_back( mInternalTable.rows[ row ].label );
    }

    // Return the list of row labels.
    return rowLabels;
}

//! Get the list of all column labels.
const vector<string> StorageTable::getColLabels() const {
    // Create a copy of the internal labels and tack total onto it if it does not exist.
    // Maybe users should have to request total explicitly?
    vector<string> colLabels( mColLabels );
    if( find( colLabels.begin(), colLabels.end(), "Total" ) == colLabels.end() ){
        colLabels.push_back( "Total" );
    }
    return colLabels;
}

/*! \brief Get the row index for a given string which represents a row label.
* \param aTypeRow The row label string.
* \return The index of the row, NO_ITEM_FOUND if it is not found.
* \author Josh Lurz
*/
int StorageTable::getRowIndex( const string& aRow ) const {
    // Search the vector.
    for( unsigned int row = 0; row < mInternalTable.rows.size(); ++row ){
        if( mInternalTable.rows[ row ].label == aRow ){ // The row label matches the search label.
            return row;
        }
    }
    // The index does not exist.
    return NO_ITEM_FOUND;
}

//! Constructor for the Item.
StorageTable::Item::Item( const string& aLabel ):
label( aLabel ),
value( 0 ){
}

//! Constructor for a Row
StorageTable::Row::Row( const string& aLabel ):
label( aLabel ),
total( 0 ){
}

//! Get the index of a column within a row.
int StorageTable::Row::getColIndex( const string& aCol ) const {
    // Search the vector.
    for( unsigned int col = 0; col < data.size(); ++col ){
        if( data[ col ].label == aCol ){ // The column label matches the search label.
            return col;
        }
    }
    // The index does not exist.
    return NO_ITEM_FOUND;
}
