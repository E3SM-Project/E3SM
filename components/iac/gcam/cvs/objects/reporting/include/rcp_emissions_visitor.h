#ifndef _RCP_EMISSIONS_VISITOR_H_
#define _RCP_EMISSIONS_VISITOR_H_
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
 * All rights to use the Software are granted on condition that such
 * rights are forfeited if User fails to comply with the terms of
 * this Agreement.
 * 
 * User agrees to identify, defend and hold harmless BATTELLE,
 * its officers, agents and employees from all liability involving
 * the violation of such Export Laws, either directly or indirectly,
 * by User.
 */

/*! 
 * \file rcp_emissions_visitor.h
 * \ingroup Objects
 * \brief RCPEmissionsVisitor class header file.
 * \author Pralit Patel
 */

#include <string>
#include <map>
#include <stack>

#include "util/base/include/default_visitor.h"
#include "util/base/include/time_vector.h"
#include "util/base/include/value.h"

/*! 
 * \ingroup Objects
 * \brief A visitor that retrieves Non-CO2 emissions and attempts to categorize them
 *        by RCP categories.
 * \details Category mappings are currently hard coded in the constructor.  The mappings
 *          are defined such that you could provide as many of sector, subsector, technology, or
 *          GHG names that are mapped to an RCP category.  If any of those names are not
 *          provided they are interpreted as *any*.  The mapping that most closely matches
 *          will be used, emissions that do not match any mappings are ignored.
 * \author Pralit Patel
 */
class RCPEmissionsVisitor : public DefaultVisitor {
public:
    RCPEmissionsVisitor();
    
    double getEmissions( const std::string& aRegionName, 
                         const std::string& aSectorCategory,
                         const std::string& aGasName,
                         const int aYear ) const;
    
    // DefaultVisitor methods
    virtual void startVisitRegion( const Region* aRegion, const int aPeriod );
    
    virtual void startVisitResource( const AResource* aResource, const int aPeriod );
    virtual void endVisitResource( const AResource* aResource, const int aPeriod );
    
    virtual void startVisitSector( const Sector* aSector, const int aPeriod );
    virtual void endVisitSector( const Sector* aSector, const int aPeriod );
    
    virtual void startVisitSubsector( const Subsector* aSubsector, const int aPeriod );
    virtual void endVisitSubsector( const Subsector* aSubsector, const int aPeriod );
    
    virtual void startVisitTechnology( const Technology* aTechnology, const int aPeriod );
    virtual void endVisitTechnology( const Technology* aTechnology, const int aPeriod );
    
    virtual void startVisitGHG( const AGHG* aGHG, const int aPeriod );
private:
    /*!
     * \brief A structure to define filtering of sectors, subsectors, and technologies
     *        to match RCP categories.
     * \details The field mRCPName must be defined.  Should the other fields be left undefined/empty
     *          they are interpreted as *any*.  Also note that the GHG name used
     *          here is the gas name according to GCAM.
     */
    struct RCPMapping {
        /*!
         * \brief Constructor to ensure the mapping contains the RCP category name.
         * \param aRCPName The RCP category mapped to.
         */
        RCPMapping( const std::string& aRCPName ):
        mRCPName( aRCPName )
        {
        }
        
        //! The RCP category name to be mapped to.
        std::string mRCPName;
        
        //! (optional) The GCAM sector name to be mapped from.
        std::string mSectorName;
        
        //! (optional) The GCAM subsector name to be mapped from.
        std::string mSubsectorName;
        
        //! (optional) The GCAM technology name to be mapped from.
        std::string mTechnologyName;
        
        //! (optional) The GCAM GHG name to be mapped from.
        std::string mGHGName;
    private:
        //! Private undefined default constructor to avoid making an invalid mapping.
        RCPMapping();
    };
    
    //! The list of all RCP mapping filters.
    std::vector<RCPMapping> mRCPMappingsList;
    
    //! Rename GCAM GHGs (key) to RCP GHGs (value)
    std::map<std::string, std::string> mGHGMap;
    
    //! The region currently being visited
    std::string mCurrentRegion;
    
    //! The most applicable mapping which is filtered as it moves down the GCAM
    //! heirarchy.
    std::stack<std::vector<const RCPMapping*> > mCurrentMappings;
    
    //! The RCP emissions by the RCP categories.  This is a map from Region ->
    //! RCP category -> GHG (by RCP name) -> emissions by period.
    std::map<std::string, std::map< std::string,
        std::map<std::string, objects::PeriodVector<Value> > > > mRCPEmissions;
    
    const std::string& getRCPGHGName( const std::string& aGasName ) const;
};

#endif // _RCP_EMISSIONS_VISITOR_H_

