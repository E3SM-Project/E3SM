#ifndef _TECHNOLOGY_TYPE_H_
#define _TECHNOLOGY_TYPE_H_
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
* \file technology_type.h
* \ingroup Objects
* \brief The TechnologyType class header file.
* \author Josh Lurz
*/

#include <map>
#include <climits>
// Forward declaration
class BaseTechnology;
class MoreSectorInfo;
class Demographic;
class NationalAccount;

/*! 
* \ingroup Objects
* \brief This class groups vintages of the same type for simplification of subsector level routines.
* \author Josh Lurz
*/

class TechnologyType
{
public:
    TechnologyType();
    bool addVintage( BaseTechnology* aTech );
    double getTotalCapitalStock( const int aUpToYear = INT_MAX ) const;
    void initializeTechsFromBase( const int aTechYear, const MoreSectorInfo* aMoreSectorInfo, const std::string& aRegionName,
                                  const std::string& aSectorName, NationalAccount& aNationalAccount,
                                  const Demographic* aDemographics, const double aCapitalStock );
    BaseTechnology* initOrCreateTech( const int aNewTechYear, const int aCurrTechYear );
    double setTotalInvestment( const std::string& aRegionName, const int aPrevYear, const int aCurrentYear,
                               const double aAnnualInvestment, const int aPeriod );
private:
    std::map<int, BaseTechnology*> mVintages;
    typedef std::map<int, BaseTechnology*>::const_iterator CVintageIterator;
    typedef std::map<int, BaseTechnology*>::iterator VintageIterator;
};

#endif // _TECHNOLOGY_TYPE_H_
