#ifndef _EXPENDITURE_H_
#define _EXPENDITURE_H_
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
* \file expenditure.h
* \ingroup Objects
* \brief Expenditure class header file.
* \author Pralit Patel
* \author Sonny Kim
*/

#include <string>
#include <vector>

#include "util/base/include/ivisitable.h"


class Tabs;
class IVisitor;

/*! 
* \ingroup Objects
* \brief A class to track various types of expenditures by Consumers or
*        ProductionTechnologies.  This mostly keeps track of values for
*        reporting.
* \author Pralit Patel, Sonny Kim
*/
class Expenditure : IVisitable
{
public:
    //! An enumeration of all possible types of expenditure.
    enum ExpenditureType {
        SOCIAL_SECURITY_TAX,
        SAVINGS,
        TAXABLE_INCOME,
        DIRECT_TAXES,
        TRANSFERS,
        DISPOSABLE_INCOME,
        CONSUMPTION,
        INCOME,
        BUDGET,
        SUBSIDY,
        INVESTMENT,
        TOTAL_IMPORTS,
        CARBON_TAX,
        // for production sectors
        DIVIDENDS,
        RETAINED_EARNINGS,
        INDIRECT_TAXES,
        INTERMEDIATE_INPUTS,
        WAGES,
        LAND_RENTS,
        RENTALS,
        TARIFFS,
        IMPORTS,
        SALES,
        COSTS,
        // Insert new values before END marker.
        END
    };

    Expenditure();
    void reset();
    void setType( const ExpenditureType aType, const double aValue );
    void addToType( const ExpenditureType aType, const double aValue );
    double getValue( const ExpenditureType aType ) const;
    void toDebugXML( const int period, std::ostream& out, Tabs* tabs ) const;
    void csvSGMOutputFile( std::ostream& aFile, const int period ) const;
    const std::string& enumToName( const ExpenditureType aType ) const;
    const std::string& enumToXMLName( const ExpenditureType aType ) const;
    virtual void accept( IVisitor* aVisitor, const int aPeriod ) const;
private:
    //! Vector of expenditures indexed by type.
    std::vector<double> mExpenditures;
};

#endif // _EXPENDITURE_H_
