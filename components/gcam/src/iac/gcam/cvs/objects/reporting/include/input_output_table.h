#ifndef _INPUT_OUTPUT_TABLE_H_
#define _INPUT_OUTPUT_TABLE_H_
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
* \file input_output_table.h
* \ingroup Objects
* \brief InputOutputTable class header file.
*
*  Detailed description.
*
* \author Josh Lurz
*/

#include <string>
#include <iosfwd>
#include <memory>

#include "util/base/include/default_visitor.h"


class Sector;
class ProductionTechnology;
class FactorSupply;
class StorageTable;
class RegionCGE;
class ProductionInput;
class DemandInput;
class Consumer;

/*! 
* \ingroup Objects
* \brief A visitor that represents the IO table for a single Region.
* \details The InputOutputTable is instantiated at the Region level where it is passed to all sectors
* so that each may add its additions to the InputOutputTable. FactorSupplies and FinalDemands are also
* added to the table. The table is in currency values.
* \author Josh Lurz
*/
class InputOutputTable : public DefaultVisitor {
public:
    InputOutputTable( const std::string& aRegionName, std::ostream& aFile );
    void finish() const;
    void startVisitRegionCGE( const RegionCGE* aRegion, const int aPeriod );
    void startVisitSector( const Sector* sector, const int aPeriod );
    void startVisitProductionTechnology( const ProductionTechnology* prodTech, const int aPeriod );
    void startVisitProductionInput( const ProductionInput* aProdInput, const int aPeriod );
    void startVisitDemandInput( const DemandInput* aDemandInput, const int aPeriod );
    void startVisitFactorSupply( const FactorSupply* factorSupply, const int aPeriod );
    void startVisitConsumer( const Consumer* aConsumer, const int aPeriod );
private:
    //! The name of the file to which to write.
    std::ostream& mFile;
    const std::string mRegionName;
    std::auto_ptr<StorageTable> mInternalTable; //!< The internal storage structure in row-column order.
    std::string mCurrSectorName; //!< The cached name of the current sector we are adding values for.
    bool mParsingConsumer; //!< Whether a consumer is currently being parsed.
    
    //! Whether an input should be visited since they should only be visited when the technology
    //! is active.
    bool mUseInput;
};

#endif // _INPUT_OUTPUT_TABLE_H_


