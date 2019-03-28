#ifndef _STUB_TECHNOLOGY_CONTAINER_H_
#define _STUB_TECHNOLOGY_CONTAINER_H_
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
 * \file stub_technology_container.h
 * \ingroup Objects
 * \brief The StubTechnologyContainer class header file.
 * \author Pralit Patel
 */

#include <map>
#include <xercesc/dom/DOMNode.hpp>
#include "technologies/include/itechnology_container.h"

/*! 
 * \ingroup Objects
 * \brief A technology container which indicates that the actual technology conatiner
 *        to use should be pulled from the GlobalTechnologyDatabase.
 * \details Parsing XML into a StubTechnology is allowed and will have the effect
 *          of making adjustments to the local copy of the TechnologyContainer
 *          from the GlobalTechnologyDatabase.  The global technology will get
 *          copied during completeInit and adjustments will be made then by
 *          calling XMLParse on the copied container.  The global technology is
 *          not written back in toInputXML instead only the XML adjustments will
 *          be.  All other ITechnologyContainer calls are merely forwared to the
 *          local copy of the global technology.  The stub technology also allows
 *          interpolating a technology before making adjustments by setting the
 *          allow-interpolate attrubute to true.
 *
 *
 *          <b>XML specification for StubTechnologyContainer</b>
 *          - XML name: -c StubTechnologyContainer::getXMLNameStatic()
 *          - Contained by: Subsector
 *          - Parsing inherited from class: None
 *          - Attributes: name The name of the technologies contained.
 *                        allow-interpolate Whether the XML adjustments with in
 *                        this stub-technology tab may interpolate a technology
 *                        in order to parse then data.
 *          - Elements:
 *              Any XML what so ever.  These elements will be kept around and will
 *              not be parsed until completeInit by the local copy of the global
 *              technology.  This means that parsing error will occur at that time.
 *
 * \author Pralit Patel
 */

class StubTechnologyContainer : public ITechnologyContainer {
public:
    StubTechnologyContainer();
    virtual ~StubTechnologyContainer();
    
    static const std::string& getXMLNameStatic();
    
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
    //! The name of the stub which will be duplicated in each contained
    //! technology however we must keep it around to be able to do the
    //! lookup.
    std::string mName;
    
    //! The cloned global technology to delegate to
    ITechnologyContainer* mTechnology;
    
    //! The vector of XML modifications to make to the global technology
    std::vector<const xercesc::DOMNode*> mXMLAdjustments;
    
    // typdef to help simplify code
    typedef std::vector<const xercesc::DOMNode*>::const_iterator CXMLIterator;
    
    static xercesc::DOMDocument* getDocumentToHoldNodes();
};

#endif // _STUB_TECHNOLOGY_CONTAINER_H_
