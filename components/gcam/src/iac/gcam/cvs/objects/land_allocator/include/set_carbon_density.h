#ifndef _SET_CARBON_DENSITY_H_
#define _SET_CARBON_DENSITY_H_
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
* \file set_carbon_density.h
* \ingroup Objects
* \brief SetCarbonDensity class header file.
* \author Pralit Patel
* \author Sonny Kim
*/

#include <string>
#include <map>

#include "util/base/include/default_visitor.h"
#include "util/base/include/iparsable.h"
#include "util/base/include/iround_trippable.h"

/*! 
* \ingroup Objects
* \brief A visitor that allows users to scale carbon densitys in land leaves.
* \details A user will have to set values then call the accept to actually
*          push the values into the leaves.
* \author Pralit Patel
* \author Sonny Kim
*/
class SetCarbonDensity : public DefaultVisitor, public IParsable, public IRoundTrippable {
public:
    SetCarbonDensity();

    void doInterpolations();

    void setCarbonDensityToPush( const std::string& aRegionName, 
                                 const int aAEZ,
                                 const std::string& aLandCategory,
                                 const double aCarbonDensityAboveAdj,
                                 const double aCarbonDensityBelowAdj,
				                 const int aYear );

    // IParsable methods
    virtual bool XMLParse( const xercesc::DOMNode* aNode );

    // IRoundTrippable methods
    virtual void toInputXML( std::ostream& aOut, Tabs* aTabs ) const;

    // DefaultVisitor methods
    void startVisitRegion( const Region* aRegion, const int aPeriod );
    void startVisitLandLeaf( const LandLeaf* aLandLeaf, const int aPeriod );
    void startVisitCarbonCalc( const ICarbonCalc* aCarbonCalc, const int aPeriod );
    void startVisitAgProductionTechnology( const AgProductionTechnology* aAgTech, const int aPeriod );
private:
    //! Rename GLM regions (key) to GCAM regions (value)
    std::map<std::string, std::string> mRegionMap;

    //! Rename GLM land type (key) to all GCAM crop types which it should apply (value)
    std::map<std::string, std::vector<std::string> > mLandMap;

    typedef std::map<int, std::pair<double, double> > YearCarbonMap;
    typedef std::map<std::string, YearCarbonMap > LandTypeCarbonMap;
    typedef std::map<std::string, LandTypeCarbonMap > RegionCarbonMap;
    //! The current region being visited.
    RegionCarbonMap::const_iterator mCurrentRegionIter;

    //! The current land type being visited.
    LandTypeCarbonMap::const_iterator mCurrentLandTypeIter;

    //! Carbon density adjustments by GLM region and GLM category.
    RegionCarbonMap mCarbonDensityAdjToSet;

    static std::string convertToGCAMAEZName( const std::string& aLandType, const int aAEZ );
};

#endif // _SET_CARBON_DENSITY_H_

