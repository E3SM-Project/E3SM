#ifndef _TOTAL_SECTOR_EMISSIONS_H_
#define _TOTAL_SECTOR_EMISSIONS_H_
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
* \file total_sector_emissions.h
* \ingroup Objects
* \brief TotalSectorEmissions header file.
* \author James Blackwood
*/

#include <string>
#include <vector>
#include <xercesc/dom/DOMNode.hpp>
#include "util/base/include/value.h"

class Sector;
class Tabs;

/*! \brief Sector level object responsible for calculating an aggregate
*          emissions factor for a set of sectors.
* \details This object allows a total emissions level to be read in for a
*          specified set of sectors. This is then used to calculate an average
*          emissions factor for all of the sectors based on their output in a
*          given year. This object implements the visitor interface so that it
*          can collect calibration information from the supply sectors.
*/
class TotalSectorEmissions
{
public:
    TotalSectorEmissions();

    void XMLParse( const xercesc::DOMNode* aNode );
    
    void toInputXML( std::ostream& aOut,
                     Tabs* aTabs ) const;
    
    void toDebugXML( const int aPeriod,
                     std::ostream& aOut,
                     Tabs* aTabs ) const;

    const std::string& getName() const;

    void setAggregateEmissionFactor( const std::string& aRegionName,
                                     const std::vector<Sector*>& aSectors,
                                     IInfo* aRegionInfo,
                                     const int aPeriod ) const;
    
    double getEmissionFactor() const;

    static const std::string& getXMLNameStatic();
    static const std::string& aggrEmissionsPrefix();
    static const std::string& aggrEmissionsYearPrefix();
private:
    //! Type of sector that this GHG will be emitted from
    std::string mType;

    //! Name of ghg
    std::string mName;

    //! Aggregate emissions of the ghg
    Value mAggregateEmissions;
    
    //! Year in which to set aggregate emissions
    int mApplicableYear;
};

#endif // _TOTAL_SECTOR_EMISSIONS_H_

