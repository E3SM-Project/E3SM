#ifndef _TECHNOLOGY_CONTAINER_H_
#define _TECHNOLOGY_CONTAINER_H_
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
 * \file technology_container.h
 * \ingroup Objects
 * \brief The TechnologyContainer class header file.
 * \author Pralit Patel
 */

#include <map>
#include <xercesc/dom/DOMNode.hpp>
#include "technologies/include/itechnology_container.h"

// Forward declarations
class InterpolationRule;

/*! 
 * \ingroup Objects
 * \brief This class contains technologies of the same type (name) and manages parsing and
 *        access to the contained technologies.
 * \details The container will allow technology vintages to be defined in arbitrary
 *          years.  It will take care of interpolating technologies as necessary and
 *          can enforce that a technology can not exist before  or after a given year.
 *          This technology container also acts as a factory for creating technology
 *          objects.
 *
 *          <b>XML specification for TechnologyContainer</b>
 *          - XML name: -c TechnologyContainer::hasTechnologyType
 *          - Contained by: Subsector, GlobalTechnologyDatabase
 *          - Parsing inherited from class: None
 *          - Attributes: name The name of the technologies contained.
 *          - Elements:
 *              - \c initial-available-year TechnologyContainer::mInitialAvailableYear
 *                   An optional paramater which specifies the first year a
 *                   vintage is allowed.  This does not have to be a model year.
 *              - \c final-available-year TechnologyContainer::mFinalAvailableYear
 *                   An optional paramater which specifies the last year a
 *                   vintage is allowed to be invested in.  This does not have to
 *                   be a model year.
 *              - \c InterpolationRule::getXMLNameStatic() TechnologyContainer::mShareWeightInterpRules
 *                   - Attributes: apply-to
 *                      What the interpolation rule applies to, currently only
 *                      share-weight is recognized.  Note that the delete=1 attribute
 *                      is recognized and in this case will be interpreted as delete
 *                      any rules which have already been parsed.
 *                   Interpolation rules which can be used to control the behavior of
 *                   technology share-weights.  Note that parsing multiple rules will
 *                   stack instead of replacing to allow building complex rules.
 *              - \c Technology::getXMLVintageName() TechnologyContainer::mVintages
 *                   - Attributes: year
 *                      The year of the technology vintage which does not have to be
 *                      a model year.
 *                   The actual technology vintage which will be created through
 *                   TechnologyContainer::createAndParseVintage
 *
 * \author Pralit Patel
 */
class TechnologyContainer : public ITechnologyContainer {
public:
    TechnologyContainer();
    virtual ~TechnologyContainer();
    
    static bool hasTechnologyType( const std::string& aTechNodeName );
    
    // ITechnologyContainer methods
    virtual void completeInit( const std::string& aRegionName, const std::string& aSectorName,
                               const std::string& aSubsectorName, DependencyFinder* aDependencyFinder,
                               const IInfo* aSubsecInfo, ILandAllocator* aLandAllocator );
    
    virtual void initCalc( const std::string& aRegionName, const std::string& aSectorName,
                           const IInfo* aSubsecInfo, const Demographic* aDemographics, const int aPeriod );
    
    virtual void postCalc( const std::string& aRegionName, const int aPeriod );
    
    virtual void toDebugXML( const int aPeriod, std::ostream& aOut, Tabs* aTabs ) const;
    
    virtual ITechnology* getNewVintageTechnology( const int aPeriod );
    
    virtual const ITechnology* getNewVintageTechnology( const int aPeriod ) const;
    
    virtual TechRangeIterator getVintageBegin( const int aPeriod );
    
    virtual CTechRangeIterator getVintageBegin( const int aPeirod ) const;
    
    virtual TechRangeIterator getVintageEnd();
    
    virtual CTechRangeIterator getVintageEnd() const;
    
    // INamed methods
    virtual const std::string& getName() const;
    
    // IParsable methods
    virtual bool XMLParse( const xercesc::DOMNode* aNode );
    
    // IRoundTrippable methods
    virtual void toInputXML( std::ostream& aOut, Tabs* aTabs ) const;
    
    // IVisitable methods
    virtual void accept( IVisitor* aVisitor, const int aPeriod ) const;
    
protected:
    virtual ITechnologyContainer* clone() const;
    
    virtual void interpolateAndParse( const xercesc::DOMNode* aNode );
    
private:
    //! The name of the technology which will be duplicated in each contained
    //! technology.
    std::string mName;
    
    //! The map that will be the primary data structure to contain technology vintages
    //! which do not have to align to model periods
    std::map<int, ITechnology*> mVintages;
    
    // Typedef iterators to help keep code readable
    typedef std::map<int, ITechnology*>::const_iterator CVintageIterator;
    typedef std::map<int, ITechnology*>::iterator VintageIterator;
    
    //! Period vector to organize technologies by model periods.  This will optimize
    //! lookups by period.  Note that all of the technology vintages in the mVintages
    //! map will not necessarily be addressed by this vector.
    objects::PeriodVector<ITechnology*> mVintagesByPeriod;
    
    //! Interpolation rules for technology share weight values.
    std::vector<InterpolationRule*> mShareWeightInterpRules;
    
    // Some typedefs to make using interpolation rules more readable.
    typedef std::vector<InterpolationRule*>::const_iterator CInterpRuleIterator;
    
    //! Optional parameter for the first year in which a vintage should exist.
    int mInitialAvailableYear;
    
    //! Optional parameter for the last year in which a technology can be invested in.
    int mFinalAvailableYear;
    
    //! A list of technology years that were created by interpolations.  This could
    //! be used to avoid writing them back out in toInputXML
    std::vector<int> mInterpolatedTechYears;
    
    bool createAndParseVintage( const xercesc::DOMNode* aNode, const std::string& aTechType );
    
    void interpolateShareWeights( const int aPeriod );
    
    void clearInterpolationRules();
    
    void interpolateVintage( const int aYear, CVintageIterator aPrevTech, CVintageIterator aNextTech );
};

#endif // _TECHNOLOGY_CONTAINER_H_
