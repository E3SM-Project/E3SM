#ifndef _ACCELERATOR_H_
#define _ACCELERATOR_H_
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
 * \file accelerator.h
 * \ingroup Objects
 * \brief The Accelerator class header file.
 * \author Josh Lurz
 */
#include <string>
#include <vector>
#include <memory>

#include "investment/include/iinvestor.h"

class Tabs;
class Demographic;
class Subsector;
class NationalAccount;
class IGrowthCalculator;
class IExpectedProfitRateCalculator;

/*! 
 * \ingroup Objects
 * \brief This class determines the investment level for a ProductionSector by
 *        accelerating the investment level from the previous period.
 * \details The Accelerator determines the investment level for a
 *          ProductionSector by taking the investment or output from the
 *          previous period and increasing it depending on the expected profit
 *          rate. If there was no investment in the previous period the
 *          investment level is determined by the total investment level in the
 *          region.
 * \author Josh Lurz
 */
class Accelerator: public IInvestor
{
public:
    Accelerator();
    ~Accelerator();
    void XMLParse( const xercesc::DOMNode* node ); 
    void toDebugXML( const int period, std::ostream& out, Tabs* tabs ) const;
    void toInputXML( std::ostream& out, Tabs* tabs ) const;
    static const std::string& getXMLNameStatic();
    void completeInit( const std::string& aRegionName, const std::string& aSectorName );
    void initCalc( std::vector<IInvestable*>& aInvestables,
                   NationalAccount& aNationalAccount, 
                   const Demographic* aDemographic,
                   const int aPeriod );
    double calcAndDistributeInvestment( std::vector<IInvestable*>& aInvestables,
                                        NationalAccount& aNationalAccount, 
                                        const Demographic* aDemographic,
                                        const int aPeriod );
    void setEfficiencyConditions( std::vector<IInvestable*>& aInvestables,
                                  NationalAccount& aNationalAccount, 
                                  const Demographic* aDemographic,
                                  const int aPeriod ) const;
private:

    //! Investment by period.
    std::vector<double> mInvestments;

    //! Fixed(exogenously specified) investment by period.
    std::vector<double> mFixedInvestments;

    //! Region name of the sector for which investment is being calculated.
    std::string mRegionName;

    //! Name of the sector for which investment is being calculated.
    std::string mSectorName;

    //! The type of the current growth calculation object.
    std::string mGrowthCalculatorType;

    //! The type of the expected profit rate calculator.
    std::string mProfitRateCalculatorType;

    //! Object responsible for calculating economic growth scalar.
    std::auto_ptr<IGrowthCalculator> mGrowthCalculator;

    //! Object responsible for calculating the expected profit rates to
    //! distribute investment.
    std::auto_ptr<IExpectedProfitRateCalculator> mProfitRateCalculator;

    //!  The investment logit exponential(RHOINV).
    double mInvestmentLogitExp;

    //! The expected profit rate function exponential(RINV).
    double mProfitElasExp;

    double calcNewInvestment( std::vector<IInvestable*>& aInvestables,
                              NationalAccount& aNationalAccount, 
                              const double aCapDependencyScalar,
                              const int aPeriod ) const;
};

#endif // _ACCELERATOR_H_
