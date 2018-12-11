#ifndef _SOLUTION_INFO_FILTER_FACTORY_H_
#define _SOLUTION_INFO_FILTER_FACTORY_H_
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
 * \file solution_info_filter_factory.h  
 * \ingroup Objects
 * \brief Header file for the SolutionInfoFilterFactory class.
 * \author Pralit Patel
 */
#include <xercesc/dom/DOMNode.hpp>
#include <string>

class ISolutionInfoFilter;

/*!
 * \ingroup Objects
 * \brief A factory which can be used to create instances of solution info filters.
 * \details There are three static methods, one to determine if this factory can create
 *          a filter with the given xml name, one to create and parse a filter with the
 *          given xml name and xml tree, and lastly one that can create them from a syntax
 *          string.  See documentation for createSolutionInfoFilterFromString for more
 *          details on the syntax string format.
 *
 * \author Pralit Patel
 */
class SolutionInfoFilterFactory {
public:
    static bool hasSolutionInfoFilter( const std::string& aXMLName );
    
    static ISolutionInfoFilter* createAndParseSolutionInfoFilter( const std::string& aXMLName,
                                                              const xercesc::DOMNode* aNode );
    
    static ISolutionInfoFilter* createSolutionInfoFilterFromString( const std::string& aFilterString );
    
private:
    static xercesc::DOMNode* buildDOMFromFilterString( const std::string& aFilterString,
                                                       xercesc::DOMDocument* aDocNode );
    
    static std::string::const_iterator findCloseParentheses( const std::string& aFilterString,
                                                             const std::string::const_iterator& aOpenParen );

};

#endif // _SOLUTION_INFO_FILTER_FACTORY_H_
