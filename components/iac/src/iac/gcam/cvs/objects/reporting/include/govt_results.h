#ifndef _GOVT_RESULTS_H_
#define _GOVT_RESULTS_H_
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
* \file govt_results.h
* \ingroup Objects
* \brief GovtResults class header file.
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
class ProductionSector;
class RegionCGE;
class GovtConsumer;
class StorageTable;

/*! 
* \ingroup Objects
* \brief A visitor which collects the total results for the government final
*        demand sector.
* \details ADD HERE
* \author Josh Lurz
*/
class GovtResults : public DefaultVisitor {
public:
    GovtResults( const std::string& aRegionName, std::ostream& aFile );
    void finish() const;
    void startVisitRegionCGE( const RegionCGE* aRegionCGE, const int aPeriod );
    void startVisitSector( const Sector* aSector, const int aPeriod );
    void startVisitProductionTechnology( const ProductionTechnology* prodTech, const int period );
    void startVisitGovtConsumer( const GovtConsumer* aGovtConsumer, const int aPeriod );
    void startVisitHouseholdConsumer( const HouseholdConsumer* aHouseholdConsumer, const int aPeriod );
private:
    //! The file to which to write.
    std::ostream& mFile;
    const std::string mRegionName; //!< The name of the region this container is reporting for.
    std::auto_ptr<StorageTable> mTaxReceipts; //!< The tax receipts storage table.
    std::auto_ptr<StorageTable> mSubsidies; //!< The subsidies storage table.
    std::auto_ptr<StorageTable> mGovtExpenditures; //!< The government expenditures storage table.
    std::string mCurrSectorName; //!< The cached name of the current sector we are adding values for.
    double mGovtTransfers; //!< Total government transfers.
    bool mParsingGovt; //!< Whether we are currently in the govt consumer.
};

#endif // _GOVT_RESULTS_H_


