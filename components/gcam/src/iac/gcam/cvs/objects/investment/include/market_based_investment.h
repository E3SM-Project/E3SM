#ifndef _MARKET_BASED_INVESTMENT_H_
#define _MARKET_BASED_INVESTMENT_H_
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
 * \file market_based_investment.h
 * \todo File name and class name don't match.
 * \ingroup Objects
 * \brief The MarketBasedInvestor class header file.
 * \author Josh Lurz
 */
#include <string>
#include <vector>
#include <xercesc/dom/DOMNode.hpp>
#include "investment/include/iinvestor.h"

class Tabs;
class Demographic;
class IInvestable;
class NationalAccount;
class IVisitor;
/*! 
 * \ingroup Objects
 * \brief The MarketBasedInvestor class calculates and distributes new
 *        investment to technologies based on a zero profit condition.
 * \details The MarketBasedInvestor class creates a trial market for each sector
 *          within a region. The markets left hand side is the unit cost, and
 *          the right hand side is the price received of the good. The output of
 *          new investment is the trial value which is adjusted until price
 *          received equals unit cost. The object distributes this output to the
 *          subsectors based on the unit cost. Increasing or decreasing output
 *          of new technologies will then increase or decrease the price and
 *          unit cost. This market will clear when the cost of increasing
 *          production by a single unit is equal to the price of the single
 *          unit.
 * \todo Check over these comments and make sure economics are right.
 * \note Since these markets are solved as a system, there might be problems if
 *       some markets in a region used this method
 * and some did not.
 * \author Josh Lurz
 */
class MarketBasedInvestor: public IInvestor
{
public:
    MarketBasedInvestor(); // should be protected.
    void XMLParse( const xercesc::DOMNode* node ); 
    void toDebugXML( const int period, std::ostream& out, Tabs* tabs ) const;
    void toInputXML( std::ostream& out, Tabs* tabs ) const;
    static const std::string& getXMLNameStatic();
    void completeInit( const std::string& aRegionName, const std::string& aGoodName );
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
protected:
    //! Investment by period.
    std::vector<double> mInvestments;
 
    //! Fixed(Exogenously specified) investment by period.
    std::vector<double> mFixedInvestments;

    //! Region name of the parent sector.
    std::string mRegionName;

    //! Sector name for which the MarketBasedInvestor is determining investment.
    std::string mSectorName;

    //! Name of the market solved for investment.
    std::string mMarketName;

    //! The investment logit exponential(RHOINV).
    double mInvestmentLogitExp;

private:
    void visitInvestables( std::vector<IInvestable*>& aInvestables,
                           IVisitor* aVisitor,
                           const int aPeriod ) const;

};

#endif // _MARKET_BASED_INVESTMENT_
