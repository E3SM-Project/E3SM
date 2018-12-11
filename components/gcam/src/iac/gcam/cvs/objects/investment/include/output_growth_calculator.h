#ifndef _OUTPUT_GROWTH_CALCULATOR_H_
#define _OUTPUT_GROWTH_CALCULATOR_H_
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
 * \file output_growth_calculator.h
 * \ingroup Objects
 * \brief The OutputGrowthCalculator class header file.
 * \author Josh Lurz
 */
#include <string>
#include <vector>
#include "investment/include/igrowth_calculator.h"

class Tabs;
class Demographic;
class IInvestable;
class NationalAccount;

/*! 
* \ingroup Objects
* \brief This object contains the standard methodology and variables to
*        calculate the growth of investment based on an output growth rate.
* \details This class is used to determine a scalar which when applied to the
*          previous period's investment will grow a sector through investment to
*          reach a proscribed level of output. This growth rate is read-in to
*          the model and may vary by period. Note that since the trial capital
*          required to grow output is only calculated once per period, the
*          equilibrium level of output for this sector will not exactly equal
*          the previous period's output times the growth rate.
* \note This object encapsulates the DEMK functionality of Legacy-SGM.
* \author Josh Lurz
*/
class OutputGrowthCalculator: public IGrowthCalculator
{
public:
    OutputGrowthCalculator();
    static const std::string& getXMLNameStatic();
    void XMLParse( const xercesc::DOMNode* aCurr );
    void toDebugXML( const int period, std::ostream& aOut, Tabs* aTabs ) const;
    void toInputXML( std::ostream& aOut, Tabs* aTabs ) const;
    double calcInvestmentDependencyScalar( const std::vector<IInvestable*>& aInvestables,
                                           const Demographic* aDemographic,
                                           const NationalAccount& aNationalAccount,
                                           const std::string& aGoodName,
                                           const std::string& aRegionName,
                                           const double aPrevInvestment,
                                           const double aInvestmentLogitExp,
                                           const int aPeriod );
protected:
    //! The read-in per period output growth rate.
    std::vector<double> mOutputGrowthRate;

    //! The cached level of trial capital per period so that it is not
    //! recalculated in each iteration.
    mutable std::vector<double> mTrialCapital;

    //! Fraction of total regional investment to use for new technologies.
    double mAggregateInvestmentFraction;

    double calcTrialCapital( const std::vector<IInvestable*>& aInvestables,
                             const NationalAccount& aNationalAccount,
                             const std::string& aGoodName,
                             const std::string& aRegionName,
                             const double aInvestmentLogitExp,
                             const int aPeriod );

    double calcOutputGap( const std::vector<IInvestable*>& aInvestables,
                          const NationalAccount& aNationalAccount,
                          const std::string& aGoodName,
                          const std::string& aRegionName,
                          const double aInvestmentLogitExp,
                          const int aPeriod ) const;
};

#endif // _OUTPUT_GROWTH_CALCULATOR_H_
