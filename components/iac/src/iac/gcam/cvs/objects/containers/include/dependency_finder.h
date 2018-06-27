#ifndef _DEPENDENCY_FINDER_H_
#define _DEPENDENCY_FINDER_H_
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
* \file dependency_finder.h
* \ingroup Objects
* \brief The DependencyFinder class header file.
* \author Josh Lurz
*/

#include <vector>
#include <map>
#include <list>
#include <string>

class Marketplace;
class ICycleBreaker;
/*! 
* \ingroup Objects
* \brief This class calculates an ordering of named objects through a scheduling
*        algorithm based on supplied dependencies.
* \details This class can order any type of object that has a name represented by
*          a string. To use this class to create a correct object ordering,
*          the user should first use the addDependency function for each
*          dependency of a object. Once all dependencies have been added,
*          calling createOrdering constructs an internal ordering. This may
*          require a cycle breaker object to resolve dependency cycles.
*          The cycle breaker object must implement the ICycleBreaker interface
*          and should be passed into the constructor.  The cycle breaker should
*          be specific to the type of objects the dependency finder is order.
*          Have a look at SectorCycleBreaker for an example. Additional
*          dependencies can be added once the ordering is created, but
*          createOrdering must be called again. Once createOrdering has
*          completed, the user can then call getOrdering to return an ordered
*          list of objects, starting with the object which has no dependencies.
*          The two key data structures are mDependencyMatrix and mObjectIndices.
*          mDependencyMatrix is a two-dimensional matrix where each cell
*          represents whether there is a dependency from the column object to
*          the row object. mObjectIndices conains a mapping of object name to
*          the index in mDependencyMatrix. Both of these data structures are
*          updated by the addDependency function.
* \author Josh Lurz
*/

class DependencyFinder
{
    friend class SectorCycleBreaker;
public:
    DependencyFinder( ICycleBreaker* aCycleBreaker );
    bool addDependency( const std::string& aObjectName,
                        const std::string& aDependency );
    void createOrdering();
    const std::vector<std::string>& getOrdering() const;
private:
    bool findAndBreakCycle( const std::list<size_t>& aStartPoints );
    bool markCycles( const size_t aCurrVertex, std::list<size_t>& aPath,
                     std::map<std::pair<size_t, size_t>, int>& aEdgesVisited ) const;
    void removeDependency( const size_t aObject, const size_t aDependency );
    const std::string& getNameFromIndex( const size_t aIndex ) const;
    
    // Need the typedef here as addTrackedItem returns a ObjectIndexMap
    // iterator.
    typedef std::map<const std::string, size_t> ObjectIndexMap;   
    ObjectIndexMap::iterator addTrackedItem( const std::string& aItem ); 
    
    //! Mapping of sector name to matrix index to avoid storing all names as
    //! strings within the matrix.
    ObjectIndexMap mObjectIndices;

    //! A square matrix where each cell contains whether the column sector
    //! depends on the row sector.
    std::vector<std::vector<bool> > mDependencyMatrix;
    
    //! The correctly ordered list of sectors.
    std::vector<std::string> mOrdering;

    //! The cycle breaker object
    ICycleBreaker* mCycleBreaker;
};

#endif // _DEPENDENCY_FINDER_H_
