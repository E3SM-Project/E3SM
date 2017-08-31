#ifndef _ITECHNOLOGY_CONTAINER_H_
#define _ITECHNOLOGY_CONTAINER_H_
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
 * \file itechnology_container.h
 * \ingroup Objects
 * \brief The ITechnologyContainer interface header file.
 * \author Pralit Patel
 */

#include "util/base/include/iparsable.h"
#include "util/base/include/iround_trippable.h"
#include "util/base/include/time_vector.h"
#include "util/base/include/inamed.h"
#include "util/base/include/ivisitable.h"

// Forward declarations
class ITechnology;
class DependencyFinder;
class IInfo;
class ILandAllocator;
class Demographic;
class InterpolationRule;

/*! 
 * \ingroup Objects
 * \brief This interface allows technologies of the same type (name) to be grouped
 *        and manages parsing and access to the contained technologies.
 * \details The reason to have an iterface is to allow a concrete implmentation as
 *          well as a stub implementation that can seamlessly initialize itself from
 *          a global technology database.
 *          
 * \author Pralit Patel
 */
class ITechnologyContainer : public INamed, public IParsable, public IRoundTrippable, public IVisitable {
    friend class StubTechnologyContainer; // to be able to call clone()
public:
    /*! 
     * \brief Constructor.
     * \details Inlined constructor to avoid compiler problems with abstract
     *          base classes. 
     */
    inline ITechnologyContainer();
    
    /*!
     * \brief Destructor.
     * \details Inlined destructor to avoid compiler problems with abstract base
     *          classes and allow deletion through the base class pointer. 
     */
    inline virtual ~ITechnologyContainer();
    
    /*!
     * \brief Perform any initializations which need to be performed after XMLParse
     *        and before the model begins.
     * \details This will call completeInit on all contained technologies as well
     *          as potenatially performing some initializations for itself such
     *          as interpolation technologies or getting them from the global
     *          technology database.
     * \param aRegionName The region name.
     * \param aSectorName Sector name, also the name of the product.
     * \param aDepDefinder Regional dependency finder.
     * \param aSubsectorInfo Subsector information object.
     * \param aLandAllocator Regional land allocator.
     */
    virtual void completeInit( const std::string& aRegionName, const std::string& aSectorName,
                               const std::string& aSubsectorName, DependencyFinder* aDependencyFinder,
                               const IInfo* aSubsecInfo, ILandAllocator* aLandAllocator ) = 0;
    
    /*!
     * \brief Perform any initializations which need to be performed before the
     *        model begins the given period.
     * \details This will call initCalc on all contained technologies passing
     *          information beteween vintages as necessary.
     * \param aRegionName The region name.
     * \param aSectorName Sector name, also the name of the product.
     * \param aSubsectorInfo Subsector information object.
     * \param aDemographics Regional demographics.
     * \param aPeriod The model period which is about to start.
     */
    virtual void initCalc( const std::string& aRegionName, const std::string& aSectorName,
                           const IInfo* aSubsecInfo, const Demographic* aDemographics, const int aPeriod ) = 0;
    
    /*!
     * \brief Perform any actions after the model has finished running the given
     *        period.
     * \details This will call postCalc on all contained technologies.
     * \param aRegionName The region name.
     * \param aPeriod The model period which just finished.
     */
    virtual void postCalc( const std::string& aRegionName, const int aPeriod ) = 0;
    
    /*!
     * \brief Print any useful debuging information for the given period.
     * \param aPeriod The model period.
     * \param aOut The output stream to print information to.
     * \param aTabs An object to keep track of how many tabs to write.
     */
    virtual void toDebugXML( const int aPeriod, std::ostream& aOut, Tabs* aTabs ) const = 0;
    
    /*!
     * \brief Get what should be the new vintage technology for the given model period.
     * \param aPeriod The model period for which to get the technology.
     * \return The new vintage technology.
     * \warning This method should not be used to iterate over all technologies.
     */
    virtual ITechnology* getNewVintageTechnology( const int aPeriod ) = 0;
    
    /*!
     * \brief Get a const pointer to what should be the new vintage technology for the
     *        given model period.
     * \param aPeriod The model period for which to get the technology.
     * \return The const pointer to the new vintage technology.
     * \warning This method should not be used to iterate over all technologies.
     */
    virtual const ITechnology* getNewVintageTechnology( const int aPeriod ) const = 0;
    
    // Typedef some iterators to abstract away syntax
    typedef std::map<int, ITechnology*>::const_reverse_iterator CTechRangeIterator;
    typedef std::map<int, ITechnology*>::reverse_iterator TechRangeIterator;
    
    /*!
     * \brief Get an iterator which can be used to iterate over all potentially
     *        operating technologies in the given period.
     * \details The begin will be the most current vintage.  This may not be the
     *          the new vintage technology in the case where technologies are no
     *          longer allowed to be invested in at the given period.
     * \param The model period.
     * \return An iterator to the most recent vintage.  This iterator will be
     *         equal to the getVintageEnd() if no technologies are available in
     *         the given period.
     */
    virtual TechRangeIterator getVintageBegin( const int aPeriod ) = 0;
    
    /*!
     * \brief Get an iterator which can be used to iterate over all potentially
     *        operating technologies in the given period.
     * \details The begin will be the most current vintage.  This may not be the
     *          the new vintage technology in the case where technologies are no
     *          longer allowed to be invested in at the given period.
     * \param The model period.
     * \return A const iterator to the most recent vintage.  This iterator will be
     *         equal to the getVintageEnd() if no technologies are available in
     *         the given period.
     */
    virtual CTechRangeIterator getVintageBegin( const int aPeirod ) const = 0;
    
    /*!
     * \brief The end iterator which should be used to compare if the user has
     *        iterated to the end of the technology range.
     * \return The vintage end iterator.
     */
    virtual TechRangeIterator getVintageEnd() = 0;

    /*!
     * \brief The end iterator which should be used to compare if the user has
     *        iterated to the end of the technology range.
     * \return The vintage end const iterator.
     */
    virtual CTechRangeIterator getVintageEnd() const = 0;
    
protected:
    /*!
     * \brief Perform a deep copy of this technology container.
     * \return A pointer to a copy of this technology container.
     */
    virtual ITechnologyContainer* clone() const = 0;
    
    /*!
     * \brief Interpolate a technology to parse the given XML.
     * \details This method will interpolate (if necessary) a technology
     *          in order to parse the given XML.
     * \ param aNode The XML to parse.
     */
    virtual void interpolateAndParse( const xercesc::DOMNode* aNode ) = 0;
};

// Inline function definitions.
ITechnologyContainer::ITechnologyContainer() {
}

ITechnologyContainer::~ITechnologyContainer() {
}

#endif // _ITECHNOLOGY_CONTAINER_H_
