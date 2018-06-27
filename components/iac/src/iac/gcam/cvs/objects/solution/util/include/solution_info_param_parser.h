#ifndef _SOLUTION_INFO_PARAM_PARSER_H_
#define _SOLUTION_INFO_PARAM_PARSER_H_
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
 * \file solution_info_param_parser.h  
 * \ingroup Objects
 * \brief Header file for the SolutionInfoParamParser class.
 * \author Pralit Patel
 */
#include <xercesc/dom/DOMNode.hpp>
#include <string>

#include "util/base/include/iparsable.h"
#include "util/base/include/time_vector.h"

class Marketplace;

/*!
 * \ingroup Objects
 * \brief A place holder to parse values from xml which can later be set
 *        into the appropriate solution info objects.
 * \details We must store these values while we are parsing since solution info
 *          objects are not created until it is time to solve a given period.  The
 *          values will be transfered from here to the appropriate solution info
 *          during the initialization of the solution info objects by the solver
 *          for the current period.  Note that parameters can be specified in the
 *          following ways and when a solution info matches multiple then parameters
 *          are merged with the first having lowest priority and last has highest:
 *              - Market type for all regions (no region attribute parsed).
 *              - Market type for a single region.
 *              - Fully qualified market by giving a good name and the region.
 *
 *          The solution info parameters which will be read are:
 *              - Solution Tolerance
 *              - Solution Floor
 *              - Bracket Interval
 *              - Max Newton Raphson Price Jump
 *              - Delta Price
 *
 *          <b>XML specification for SolutionInfoParamParser</b>
 *          - XML name: \c solution-info-param-parser
 *          - Contained by: Scenario
 *          - Parsing inherited from class: None.
 *          - Elements:
 *              - \c solution-tolerance double SolutionInfoParamParser::mSolutionInfoParams
 *                  Attributes:
 *                      - \c good The market good name to set for.
 *                      - \c region The market region name to set for.
 *                      - \c period The solution period to set for.
 *                      - \c market-type The market type to set for (only if good is not
 *                                       specified, also when using market type not
 *                                       specifying a region implies global).
 *                  The solution tolerance to set into a solution info.
 *              - \c solution-floor double SolutionInfoParamParser::mSolutionInfoParams
 *                  Attributes:
 *                      - \c good The market good name to set for.
 *                      - \c region The market region name to set for.
 *                      - \c period The solution period to set for.
 *                      - \c market-type The market type to set for (only if good is not
 *                                       specified, also when using market type not
 *                                       specifying a region implies global).
 *                  The solution floor to set into a solution info.
 *              - \c bracket-interval double SolutionInfoParamParser::mSolutionInfoParams
 *                  Attributes:
 *                      - \c good The market good name to set for.
 *                      - \c region The market region name to set for.
 *                      - \c period The solution period to set for.
 *                      - \c market-type The market type to set for (only if good is not
 *                                       specified, also when using market type not
 *                                       specifying a region implies global).
 *                  The bracket interval to set into a solution info.
 *              - \c max-price-change double SolutionInfoParamParser::mSolutionInfoParams
 *                  Attributes:
 *                      - \c good The market good name to set for.
 *                      - \c region The market region name to set for.
 *                      - \c period The solution period to set for.
 *                      - \c market-type The market type to set for (only if good is not
 *                                       specified, also when using market type not
 *                                       specifying a region implies global).
 *                  The max newton raphson price jump to set into a solution info.
 *              - \c delta-price double SolutionInfoParamParser::mDeltaPrice
 *                  Attributes:
 *                      - \c good The market good name to set for.
 *                      - \c region The market region name to set for.
 *                      - \c period The solution period to set for.
 *                      - \c market-type The market type to set for (only if good is not
 *                                       specified, also when using market type not
 *                                       specifying a region implies global).
 *                  The newton raphson delta price to set into a solution info.
 *
 * \author Pralit Patel
 * \todo Do we want to write this back out in toInputXML?
 */
class SolutionInfoParamParser : public IParsable {
public:
    SolutionInfoParamParser();
    ~SolutionInfoParamParser();
    
    static const std::string& getXMLNameStatic();
    
    //! The struct to hold each possible parsable solution info parameter.
    struct SolutionInfoValues {
        SolutionInfoValues();
        
        void mergeValues( const SolutionInfoValues& aMergeValues );
        
        //! Solution tolerance
        double mSolutionTolerance;
        
        //! Solution floor
        double mSolutionFloor;
        
        //! Bracket interval
        double mBracketInterval;
        
        //! Max newton raphson price jump
        double mMaxNRPriceJump;
        
        //! Delta price
        double mDeltaPrice;
    };
    
    SolutionInfoValues getSolutionInfoValuesForMarket( const std::string& aGoodName, const std::string& aRegionName,
                                                       const std::string& aMarketType, const int aPeriod ) const;
    
    // IParsable methods
    virtual bool XMLParse( const xercesc::DOMNode* aNode );
    
private:
    //! A data structure that will go from periods to the pair of good name / market
    //! type and region name which maps to the struct SolutionInfoValues which will
    //! contain each of the possible values for a solution info.
    objects::PeriodVector<std::map<std::pair<std::string, std::string>, SolutionInfoValues> > mSolutionInfoParams;
    
    std::vector<SolutionInfoValues*> getSolutionInfoValuesFromAttrs( const xercesc::DOMNode* aNode );
};

#endif // _SOLUTION_INFO_PARAM_PARSER_H_
