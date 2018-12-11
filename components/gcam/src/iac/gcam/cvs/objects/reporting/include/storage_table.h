#ifndef _STORAGE_TABLE_H_
#define _STORAGE_TABLE_H_
#if defined(_MSC_VER)
#pragma once
#endif

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
* \file storage_table.h
* \ingroup Objects
* \brief StorageTable class header file.
*
*  Detailed description.
*
* \author Josh Lurz
*/

#include <string>
#include <vector>

/*! 
* \ingroup Objects
* \brief A datastructure which stores in a row column format which is referenced by the column and row names.
* \details This is a sparse table as all columns do not have to contain the same rows. SWAP. For this reason the columns
* can be iterated by the rows cannot.
* \author Josh Lurz
* \todo Fix handling of total row so that a consumer of this class must request it.
*/

class StorageTable {
public:
	StorageTable();
	void clear();
    bool isEmpty() const;
    void addColumn( const std::string& aCol );
    void addToType( const int aRow, const std::string& aCol, const double aValue );
    void addToType( const std::string& aRow, const std::string& aCol, const double aValue );
    void setType( const std::string& aRow, const std::string& aCol, const double aValue ); 
    double getValue( const std::string& aRow, const std::string& aCol ) const;
    const std::vector<std::string> getRowLabels() const;
    const std::vector<std::string> getColLabels() const;
private:
    int getRowIndex( const std::string& aRow ) const;
    const static int NO_ITEM_FOUND = -1;
    std::vector<std::string> mColLabels;
    //! Structure for each Item
    struct Item {
        explicit Item( const std::string& aLabel );
        std::string label;
        double value;
    };

    //! Structure for each Column.
    struct Row {
        explicit Row( const std::string& aLabel );
        int getColIndex( const std::string& aCol ) const;
        std::string label;
        std::vector<Item> data;
        double total;
    };
    //! Structure for the internal storage.
    struct InternalTable {
        std::string label;
        std::vector<Row> rows;
    };
	InternalTable mInternalTable; //!< The internal storage.
};

#endif // _STORAGE_TABLE_H_

