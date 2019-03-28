#ifndef _MARKET_NAME_SOLUTION_INFO_FILTER_H_
#define _MARKET_NAME_SOLUTION_INFO_FILTER_H_
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
 * \file market_name_solution_info_filter.h  
 * \ingroup Objects
 * \brief Header file for the MarketNameSolutionInfoFilter class.
 * \author Pralit Patel
 */
#include <xercesc/dom/DOMNode.hpp>
#include <string>

#include "solution/util/include/isolution_info_filter.h"

class SolutionInfo;

/*!
 * \ingroup Objects
 * \brief A solution info filter which will accept any solution info which
 *        match the read in market name.
 * \details Note that it is generally recommended that the user avoid using
 *          this filter in the general case since it would be possible that
 *          market names could get out of sync with the names used in the
 *          dataset.
 *          <b>XML specification for MarketNameSolutionInfoFilter</b>
 *          - XML name: \c market-name-solution-info-filter
 *          - Contained by:
 *          - Parsing inherited from class: None.
 *          - Elements:
 *              - \c market-name string MarketNameSolutionInfoFilter::mAcceptMarketName
 *                      The name that will be directly compared against the market name.
 *                      note that this name should include the market region as well.
 *
 * \author Pralit Patel
 */
class MarketNameSolutionInfoFilter : public ISolutionInfoFilter {
public:
    MarketNameSolutionInfoFilter();
    ~MarketNameSolutionInfoFilter();
    
    static const std::string& getXMLNameStatic();
    
    // ISolutionInfoFilter methods
    virtual bool acceptSolutionInfo( const SolutionInfo& aSolutionInfo ) const;
    
    // IParsable methods
    virtual bool XMLParse( const xercesc::DOMNode* aNode );
    
private:
    //! The name of the market which will be accepted
    std::string mAcceptMarketName;
};

#endif // _MARKET_NAME_SOLUTION_INFO_FILTER_H_
