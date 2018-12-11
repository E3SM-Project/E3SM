#ifndef _INVERSE_CALIBRATION_MARKET_H_
#define _INVERSE_CALIBRATION_MARKET_H_
#if defined(_MSC_VER_)
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
* \file inverse_calibration_market.h
* \ingroup Objects
* \brief The InverseCalibrationMarket class header file.
* \author Josh Lurz
*/

#include "marketplace/include/market.h"

/*!
* \ingroup Objects
* \brief A class which defines a market for solving calibration markets which
*        have an inverse response to the solved term.
* \details Inverse calibration markets allow the calibration of values where as
*          the solved or price term increases, the trial quantity decreases.
*          This is the inverse behavior of the CalibrationMarket. To use a
*          InverseCalibrationMarket initially set for each period supply to
*          equal the known or target value of the variable. Use the price term
*          of the market to calculate the iteration value, which should be
*          stored in the demand side of the market.
* \note The implementation is very similar to a CalibrationMarket. The functions
*       which differ are nullSupply, nullDemand and
*       meetsSpecialSolutionCriteria. These differ due to the constraint being
*       stored in the supply variable instead of the demand variable.
* \author Josh Lurz
*/

class InverseCalibrationMarket: public Market {
public:
    InverseCalibrationMarket( const std::string& goodNameIn, const std::string& regionNameIn, const int periodIn );
    ~InverseCalibrationMarket();
    virtual IMarketType::Type getType() const;

    virtual void initPrice();
    virtual void setPrice( const double priceIn );
    virtual void set_price_to_last_if_default( const double lastPriceIn );
    virtual void set_price_to_last( const double lastPrice );
    virtual double getPrice() const;

    virtual void addToDemand( const double demandIn );
    virtual double getDemand() const;

    virtual void nullSupply();
    virtual void nullDemand();
    virtual double getSupply() const;
    virtual void addToSupply( const double supplyIn );
    virtual bool meetsSpecialSolutionCriteria() const;

    virtual bool shouldSolve() const;
    virtual bool shouldSolveNR() const;
protected:
    virtual void toDebugXMLDerived( std::ostream& out, Tabs* tabs ) const;
};

#endif // _INVERSE_CALIBRATION_MARKET_H_
