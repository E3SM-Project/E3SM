#ifndef _GET_GLM_DATA_H__
#define _GET_GLM_DATA_H__
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
* \file get_glm_data.h
* \ingroup Objects
* \brief GetGLMData class header file.
* \author Pralit Patel
* \author Sonny Kim
*/

#include <string>
#include <map>

#include "util/base/include/default_visitor.h"
#include "util/base/include/time_vector.h"

/*! 
* \ingroup Objects
* \brief A visitor to collect and aggregate GCAM data to pass to GLM.
* \details This class will aggregate GCAM land cover data and organize
*          it into GLM land categories.
* \author Pralit Patel
* \author Sonny Kim
*/
class GetGLMData : public DefaultVisitor {
public:
    GetGLMData();

    double getLandCover( const std::string& aGLMRegionName, 
                         const int aAEZ,
                         const std::string& aGLMLandCategory,
                         const int aYear ) const;

    double getProductionInCarbon( const std::string& aGLMRegionName, 
                                  const int aAEZ,
                                  const std::string& aGLMCrop,
                                  const int aYear ) const;

    std::pair<double, double> getCarbonDensity( const std::string& aGLMRegionName, 
                                                const int aAEZ,
                                                const std::string& aGLMLandCategory,
                                                const int aYear ) const;

    // DefaultVisitor methods
    void startVisitRegion( const Region* aRegion, const int aPeriod );
    void startVisitLandLeaf( const LandLeaf* aLandLeaf, const int aPeriod );
    void startVisitAgProductionTechnology( const AgProductionTechnology* aAgTech, const int aPeriod );
private:
    //! Rename GCAM regions (key) to GLM regions (value)
    std::map<std::string, std::string> mRegionMap;

    //! Rename GCAM land types (key) to GLM categories (value)
    std::map<std::string, std::string> mLandMap;

    //! Rename GCAM crop names (key) to GLM categories (value)
    std::map<std::string, std::string> mCropMap;

    //! Rename GCAM land types (key) to GLM categories (value) for carbon densities
    std::map<std::string, std::string> mLandMapForDensities;

    //! The current GLM region name.
    std::string mCurrentGLMRegionName;

    // Typedef the data structure for storing by region / aez / land name / data by period
    typedef std::map<std::string,
              std::map<int,
                std::map<std::string,
                  objects::PeriodVector<double> > > > RegAEZLandMap;

    //! Land cover by GLM region and GLM category.
    RegAEZLandMap mLandCoverData;

    //! Crop production in carbon by GLM region and GLM category.
    RegAEZLandMap mProductionData;

    // Typedef the data structure for storing by region / aez / land name 
    //     / data by year / above and below ground densities
    typedef std::map<std::string,
              std::map<int,
                std::map<std::string,
                  std::map<int,
                    std::pair<double, double> > > > > CarbonDensityMap;

    //! Above and below ground carbon densities by GLM region and category.
    CarbonDensityMap mCarbonDensityData;

    static const std::string& doMapLookup( const std::map<std::string, std::string>& aLookupMap,
                                           const std::string& aKey );

    static std::pair<std::string, int> parseTypeAndAEZFromName( const std::string& aGCAMName );
};

#endif // _GET_GLM_DATA_H__

