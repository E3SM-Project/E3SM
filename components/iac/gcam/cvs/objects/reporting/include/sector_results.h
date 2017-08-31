#ifndef _SECTOR_RESULTS_H_
#define _SECTOR_RESULTS_H_
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
* \file sector_results.h
* \ingroup Objects
* \brief SectorResults class header file.
*
*  Detailed description.
*
* \author Josh Lurz
*/

#include <string>
#include <iosfwd>
#include <memory>

#include "util/base/include/default_visitor.h"

/*! 
* \ingroup Objects
* \brief A visitor that represents the total results for all sectors within a region.
* \details ADD HERE
* \author Josh Lurz
*/
class Sector;
class ProductionTechnology;
class ProductionSector;
class StorageTable;

class SectorResults : public DefaultVisitor {
public:
    SectorResults( const std::string& aRegionName, std::ostream& aFile );
    void finish() const;
    void startVisitSector( const Sector* aSector, const int aPeriod );
    void startVisitProductionSector( const ProductionSector* aProdSector, const int aPeriod );
    void startVisitProductionTechnology( const ProductionTechnology* prodTech, const int period );
private:
    //! The current region name.
    const std::string mCurrentRegionName;

    //! The current sector name name.
    std::string mCurrentSectorName;
    
    //! The file to which to write.
    std::ostream& mFile;

    std::auto_ptr<StorageTable> mInternalTable; //!< The internal storage structure in row-column order.
};

#endif // _SECTOR_RESULTS_H_


