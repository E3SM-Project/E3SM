#ifndef _SOCIAL_ACCOUNTING_MATRIX_H_
#define _SOCIAL_ACCOUNTING_MATRIX_H_
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
* \file social_accounting_matrix.h
* \ingroup Objects
* \brief SocialAccountingMatrix class header file.
* \author Pralit Patel
* \author Katherine Chung
*/

#include <string>
#include <memory>
#include <iosfwd>
#include "util/base/include/default_visitor.h"



class Region;
class Sector;
class Subsector;
class BaseTechnology;
class HouseholdConsumer;
class GovtConsumer;
class InvestConsumer;
class TradeConsumer;
class ProductionTechnology;
class FactorSupply;
class Input;
class StorageTable;

/*! 
* \ingroup Objects
* \brief An object which outputs a social accounting matrix.
* \details TODO
* \author Pralit Patel, Katherine Chung
*/
class SocialAccountingMatrix : public DefaultVisitor {
public:
    SocialAccountingMatrix( const std::string& aRegionName, std::ostream& aFile );
    void finish() const;
    void startVisitHouseholdConsumer( const HouseholdConsumer* householdConsumer, const int aPeriod );
    void startVisitGovtConsumer( const GovtConsumer* govtConsumer, const int aPeriod );
    void startVisitTradeConsumer( const TradeConsumer* tradeConsumer, const int aPeriod );
    void startVisitInvestConsumer( const InvestConsumer* investConsumer, const int aPeriod );
    void startVisitProductionTechnology( const ProductionTechnology* prodTech, const int period );
    void startVisitFactorSupply( const FactorSupply* factorSupply, const int period );
private:
    enum CategoryType {
        ACTIVITIES,
        COMMODITIES,
        FACTORS_LAND,
        FACTORS_LABOR,
        FACTORS_CAPITAL,
        HOUSEHOLDS,
        ENTERPRISES,
        GOVERNMENT,
        CARBON_TAX,
        CAPITAL_ACCOUNT,
        REST_OF_WORLD,
        TOTALS
    };
    void addToType( CategoryType aRowCat, CategoryType aColCat, double aValue );
    static const std::string& getString( const CategoryType aType );
    
    std::auto_ptr<StorageTable> mTable;
    const std::string mRegionName;

    //! File to which to write.
    std::ostream& mFile;
};

#endif // _SOCIAL_ACCOUNTING_MATRIX_H_
