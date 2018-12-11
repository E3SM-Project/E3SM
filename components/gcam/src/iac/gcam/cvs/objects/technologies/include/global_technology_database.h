#ifndef _GLOBAL_TECHNOLOGY_DATABASE_H_
#define _GLOBAL_TECHNOLOGY_DATABASE_H_
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
* \file global_technology_database.h
* \ingroup Objects
* \brief The GlobalTechnologyDatabase class header file.
* \author Pralit Patel
*/

#include <xercesc/dom/DOMNode.hpp>
#include "util/base/include/iparsable.h"
#include "util/base/include/iround_trippable.h"

// Forward declarations
class Tabs;
class ITechnologyContainer;

/*! 
 * \ingroup Objects
 * \brief GlobalTechnologyDatabase holds ITechnologyContainers that can be shared
 *       globally.
 * \details This class contains ITechnologyContainers and could but used in conjunction
 *         with SubTechnologyContainers to retrieve and utilize them through the
 *         getTechnology method.  The global technologies are categorized by sector
 *         and subsector names to avoid technology name collisions.
 *
 *          <b>XML specification for GlobalTechnologyDatabase</b>
 *          - XML name: -c GlobalTechnologyDatabase::getXMLNameStatic()
 *          - Contained by: World
 *          - Parsing inherited from class: None
 *          - Attributes: None
 *          - Elements:
 *              - \c location-info
 *                  - Attributes: sector-name, subsector-name
 *                      The sector and subsector names under which this technology
 *                      should be contained.
 *                  - Elements:
 *                      - \c TechnologyContainer::hasTechnologyType() GlobalTechnologyDatabase::mTechnologyList
 *                          The technology container that can be accessed through getTechnology
 *
 * \author Pralit Patel
 */
class GlobalTechnologyDatabase : public IParsable, public IRoundTrippable {
public:
    // Singleton class will only provide a getInstance method and no constructors
    static GlobalTechnologyDatabase* getInstance();
    ~GlobalTechnologyDatabase();
    
    static const std::string& getXMLNameStatic();
    
    const ITechnologyContainer* getTechnology( const std::string& aSectorName,
                                               const std::string& aSubsectorName,
                                               const std::string& aTechnologyName ) const;
    
    // IParsable methods
    virtual bool XMLParse( const xercesc::DOMNode* aNode );
    
    // IRoundTrippable methods
    virtual void toInputXML( std::ostream& aOut, Tabs* aTabs ) const;
    
private:
    //! List of GlobalTechnologies
    std::map<std::pair<std::string, std::string>, std::vector<ITechnologyContainer*> > mTechnologyList;
    
    // some useful iterator typedefs
    typedef std::map<std::pair<std::string, std::string>, std::vector<ITechnologyContainer*> >::const_iterator CTechLocationIterator;
    typedef std::vector<ITechnologyContainer*>::const_iterator CTechListIterator;
    
    // private constructors to enforce a single instance
    GlobalTechnologyDatabase() {}
    
    // intentionally undefined
    GlobalTechnologyDatabase( const GlobalTechnologyDatabase& );
    GlobalTechnologyDatabase& operator=( const GlobalTechnologyDatabase& );
};

#endif // _GLOBAL_TECHNOLOGY_DATABASE_H_

