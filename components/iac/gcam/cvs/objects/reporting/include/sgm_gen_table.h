#ifndef _SGM_GEN_TABLE_H_
#define _SGM_GEN_TABLE_H_
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
* \file sgm_gen_table.h
* \ingroup Objects
* \brief SGMGenTable class header file.
*
*  Detailed description.
*
* \author Katherine Chung, Sonny Kim
*/

#include <string>
#include <map>
#include <iosfwd>
#include "util/base/include/default_visitor.h"

class Region;
class NationalAccount;
class Demographic;
class Sector;
class ProductionSector;
class Subsector;
class BaseTechnology;
class Consumer;
class HouseholdConsumer;
class GovtConsumer;
class InvestConsumer;
class TradeConsumer;
class ProductionTechnology;
class FactorSupply;
class Modeltime;

/*!
 * \brief A visitor which constructs the general SGM output tables.
 */
class SGMGenTable : public DefaultVisitor {
public:
    SGMGenTable( const std::string& aName, const std::string& aHeader,
                 const Modeltime* aModeltime );
    void finish() const;
    void startVisitRegionCGE( const RegionCGE* regionCGE, const int aPeriod );
    void endVisitRegionCGE( const RegionCGE* aRegionCGE, const int aPeriod );
    void startVisitSector( const Sector* aSector, const int aPeriod );
    void endVisitSector( const Sector* aSector, const int aPeriod );
    void startVisitProductionSector( const ProductionSector* aProductionSector, const int aPeriod );
    void startVisitDemographic( const Demographic* aDemographic, const int aPeriod );
    void startVisitConsumer( const Consumer* aConsumer, const int aPeriod );
    void startVisitHouseholdConsumer( const HouseholdConsumer* householdConsumer, const int aPeriod );
    void startVisitGovtConsumer( const GovtConsumer* govtConsumer, const int aPeriod );
    void startVisitTradeConsumer( const TradeConsumer* tradeConsumer, const int aPeriod );
    void startVisitInvestConsumer( const InvestConsumer* investConsumer, const int aPeriod );
    void startVisitProductionTechnology( const ProductionTechnology* prodTech, const int aPeriod );
    void startVisitFactorSupply( const FactorSupply* aFactorSupply, const int aPeriod );
    void startVisitNationalAccount( const NationalAccount* aNationalAccount, const int aPeriod );

    // Non-interface function.
    void setOutputFile( std::ostream& aOutputFile );
private:

    void addToType( const int aTypeRow, const std::string aTypeCol, const double value );
    void setType( const int aTypeRow, const std::string aTypeCol, const double value );
    double getValue( const int aTypeRow, const std::string aTypeCol) const;
    
    static bool isEnergyGood( const std::string& aRegionName,
                              const std::string& aGoodName );

    static bool isPrimaryEnergyGood( const std::string& aRegionName,
                                     const std::string& aGoodName );

    static bool isSecondaryEnergyGood( const std::string& aRegionName,
                                       const std::string& aGoodName );

    const std::string mName;
    const std::string mHeader;

    //! The current region name.
    std::string mCurrentRegionName;
    
    //! The current sector name.
    std::string mCurrentSectorName;
    
    //! The file to which to write.
    std::ostream* mFile;

    std::map<int, std::map< std::string, double> > mTable;
    const Modeltime* mModeltime;
};

#endif // _SGM_GEN_TABLE_H_


