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
* \file dependency_finder.cpp
* \ingroup Objects
* \brief The DependencyFinder class source file.
* \author Josh Lurz
*/

#include "util/base/include/definitions.h"
#include <cassert>
#include <stack>
#include <algorithm>
#include <limits.h>		// INT_MAX and other goodies
#include "containers/include/dependency_finder.h"
#include "util/logger/include/ilogger.h"
#include "containers/include/icycle_breaker.h"

using namespace std;

/*!
 * \brief Constructor.
 * \param aCycleBreaker The cycle breaker object.  This can be null.
 */
DependencyFinder::DependencyFinder( ICycleBreaker* aCycleBreaker ):
mCycleBreaker( aCycleBreaker )
{
}

/*! 
 * \brief Add a dependency for an object.
 * \details This function is used to mark a single dependency, from aObjectName
 *          to aDependency.
 * \param aObjectName Name of the object which has a new dependency.
 * \param aDependency name of the item aObjectName is dependent on.
 * \return Whether the dependency was added to the matrix.
 */
bool DependencyFinder::addDependency( const string& aObjectName,
                                      const string& aDependency ){
    // Check if the object is already in the mapping of object name to matrix index.
    ObjectIndexMap::iterator objectLocation = mObjectIndices.find( aObjectName );
    if( objectLocation == mObjectIndices.end() ){
        // Update the object location iterator after the item is added.
        objectLocation = addTrackedItem( aObjectName );
    }

    // Check if the dependency is already in the mapping.
    ObjectIndexMap::iterator dependencyLocation = mObjectIndices.find( aDependency );
    if( dependencyLocation == mObjectIndices.end() ){
        // Update the dependency location iterator after the item is added.
        dependencyLocation = addTrackedItem( aDependency );
    }

    // The matrix is now setup correctly, add the dependency.
    assert( mDependencyMatrix.size() > objectLocation->second );
    assert( mDependencyMatrix[ objectLocation->second ].size() > dependencyLocation->second );

    // Check if the dependency already exists.
    if( mDependencyMatrix[ objectLocation->second ][ dependencyLocation->second ] ){
        return false;
    }
    // Add the dependency and return that it was a new dependency.
    mDependencyMatrix[ objectLocation->second ][ dependencyLocation->second ] = true;
    return true;
}

/*!
 * \brief Find an ordering of the objects in the dependency finder which orders
 *        each object before each object that depends on it.
 * \details This is referred to as a topological sort. The algorithm is as
 *          follows: Search the adjacency matrix for a vertice, in this
 *          implementation a column in the matrix, with no dependencies. If there
 *          is none, the graph has a cycle and cannot currently be ordered. To
 *          correct this, call findAndBreakCycle to remove a single cycle from
 *          the matrix so the algorithm can continue. Next, add the vertice found
 *          in the previous step to the ordering, and remove it and all
 *          dependencies on it from the adjacency matrix. Start over at the first
 *          step. Repeat this process until there are no nodes left in the graph.
 *          Once this function has been called, the caller may then call
 *          getOrdering to return the ordered vector.
 * \note This function removes nodes by setting the removed flag to true instead
 *       of actually removing the vertice.
 * \pre All relevant objects must have called addDependency() to add their dependencies to the matrix.
 * \todo The control flow of this function could be less confusing.
 * \sa getOrdering
 */
void DependencyFinder::createOrdering() {
    // If there is an existing stored ordering, clear it.
    mOrdering.clear();

    // Create a vector which marks which vertices are removed.
    vector<bool> removed( mDependencyMatrix.size() );

    // Search until the ordering contains all vertices in the matrix.
    while( mOrdering.size() < mDependencyMatrix.size() ){
        unsigned int verticeToRemove = INT_MAX;

        // Search for a vertex with no dependencies.
        for( unsigned int i = 0; i < mDependencyMatrix.size(); ++i ){
            // Only search the vertex if it has not been removed.
            if( removed[ i ] ){
                continue;
            }
            // Check for any dependencies.
            bool depFound = false;
            for( unsigned int j = 0; j < mDependencyMatrix.size(); ++j ){
                // Check if there is a dependency at this location.
                if( !removed[ j ] && mDependencyMatrix[ i ][ j ] ){
                    // Found a dependency so break out of the loop to stop
                    // searching this column and move onto the next.
                    depFound = true;
                    break;
                }
            }
            // If we did not find a dependency, set the index to remove and
            // break the loop so the object can be removed. Otherwise continue
            // searching.
            if( !depFound ){
                verticeToRemove = i;
                break;
            }
        }

        // Check if we found a vertex to remove.
        if( verticeToRemove == INT_MAX ){
            // Since there was no vertex with zero dependencies, this graph has
            // a cycle.
            ILogger& depFinderLog = ILogger::getLogger( "dependency_finder_log" );
            depFinderLog.setLevel( ILogger::DEBUG );
            depFinderLog << "Graph has at least one cycle, attempting to remove." << endl;
            
            // We will start our cycle searches from the vertices that do not
            // have anything dependant on it.  Note that this is pretty much an
            // arbitrary decision.
            list<size_t> startPoints;
            for( size_t col = 0; col < mDependencyMatrix.size(); ++col ) {
                if( !removed[ col ] ) {
                    bool dependsOnCol = false;
                    for( size_t row = 0; row < mDependencyMatrix.size() && !dependsOnCol; ++row ) {
                        if( mDependencyMatrix[ row ][ col ] ) {
                            dependsOnCol = true;
                        }
                    }
                    if( !dependsOnCol ) {
                        startPoints.push_back( col );
                    }
                }
            }
            
            // Find and break cycle does not need to know about the vertices that
            // have already been removed since they are guaranteed to not be part
            // of a cycle.
            if( !findAndBreakCycle( startPoints) ){
                // There was a cycle but it could not be broken.
                ILogger& mainLog = ILogger::getLogger( "main_log" );
                mainLog.setLevel( ILogger::ERROR );
                depFinderLog.setLevel( ILogger::ERROR );
                mainLog << "Failed to remove cycle from dependency matrix. Object ordering failed." << endl;
                depFinderLog << "Failed to remove cycle from dependency matrix. Object ordering failed." << endl;
                mOrdering.clear();
                return;
            }
        }
        else {
            // Add the vertex found to the ordering and remove it from the
            // matrix.
            mOrdering.push_back( getNameFromIndex( verticeToRemove ) );
            removed[ verticeToRemove ] = true;
        }
    }
    
    // For some unknown reason, a blank line can occur in the order.
    // This can cause a failure in the sorting routine. Remove any blank lines to prevent this.
    for ( vector<string>::iterator orderItem = mOrdering.begin(); 
          orderItem != mOrdering.end();  orderItem++ ) {
        if ( *orderItem == "" ) {
            ILogger& depFinderLog = ILogger::getLogger( "dependency_finder_log" );
            depFinderLog.setLevel( ILogger::DEBUG );
            depFinderLog << "Blank item removed from dependency list." << endl;
            mOrdering.erase ( orderItem );
        }
    }

    // Sorting finished, the internal ordering can now be fetched by
    // getOrdering.
}

/*!
 * \brief Get the object ordering.
 * \pre createOrdering has been called.
 * \return The correct object ordering.
 * \sa createOrdering
 */
const vector<string>& DependencyFinder::getOrdering() const {
    // Check if the ordering has been initialized. Print an error if it has
    // not.
    if( mOrdering.empty() ){
        ILogger& depFinderLog = ILogger::getLogger( "dependency_finder_log" );
        depFinderLog.setLevel( ILogger::ERROR );
        depFinderLog << "Returning an empty object ordering." << endl;
    }
    return mOrdering;
}

/*!
 * \brief Remove a dependency from the matrix.
 * \param aObject Object index for which to remove the dependency.
 * \param aDependency Dependency index to remove.
 */
void DependencyFinder::removeDependency( const size_t aObject, const size_t aDependency ){
    // Remove the dependency, or edge.
    assert( mDependencyMatrix.size() > aObject );
    assert( mDependencyMatrix[ aObject ].size() > aDependency );
    mDependencyMatrix[ aObject ][ aDependency ] = false;
}

/*!
 * \brief Add an item which should have its dependencies tracked.
 * \pre The item is not already being tracked.
 * \param aItem Name of the item to track.
 * \return An iterator to the location of the new item.
 */
DependencyFinder::ObjectIndexMap::iterator DependencyFinder::addTrackedItem( const string& aItem ){
    // Add the item to the mapping of item name to index within the matrix.
    const size_t newLocation = mDependencyMatrix.size();

    // Make pair creates a name value pair to insert into the matrix.
    pair<ObjectIndexMap::iterator, bool> newPositionPair = mObjectIndices.insert( make_pair( aItem, newLocation ) );
    
    // Check the precondition that the item does not already exist.
    assert( newPositionPair.second );
    
    // Now add the item to the dependency matrix. Loop through the matrix and
    // add a new position on the end of each row.
    for( unsigned int row = 0; row < mDependencyMatrix.size(); ++row ){
        mDependencyMatrix[ row ].push_back( false );
    }

    // Add a new row to the matrix for the item. Default to not having any
    // dependencies.
    mDependencyMatrix.push_back( vector<bool>( newLocation + 1, false ) );

    // Return an iterator to the position within the index map.
    return newPositionPair.first;
}

/*!
 * \brief Function which locates and breaks a single cycle in the dependency
 *          matrix.
 * \details This function performs a search through the dependency graph starting
 *          at the given start points to find the most commonly traversed edges
 *          involved in a cycle.  See DependencyFinder::markCycles
 *          for more details on the alogorithm.  Once we have finished our searches
 *          we will remove the most traversed dependency.  However since the start
 *          and end points are some what arbitrary we weight the visit counts to
 *          boost out edges from a vertex that combines many highly traversed edges
 *          into fewer edges.  It calls breakCycle to remove the edges between two
 *          nodes in the cycle. Control is then returned to the ordering algorithm,
 *          which may call this function again if multiple cycles exist.
 * \param aStartPoints The vertices to start our cycle searches from.
 * \todo Although this search attempts to remove the optimal edge there still is
 *       no garuntee that it minimizes the number of dependencies broken to remove
 *       all cycles.
 * \return Whether a cycle was found and broken.
 */
bool DependencyFinder::findAndBreakCycle( const list<size_t>& aStartPoints ) {
    if( mCycleBreaker ){
        // Set up a data structure that counts the number of traversals of edges
        // which are part of a cycle.
        map<pair<size_t, size_t>, int> edgeVisits;
        typedef list<size_t>::const_iterator StartIterator;
        for( StartIterator startIter = aStartPoints.begin(); startIter != aStartPoints.end(); ++startIter ) {
            list<size_t> path;
            markCycles( *startIter, path, edgeVisits );
        }
        
        // Create a weighting factor to benefit out-edges from a vertex which
        // combines edges.
        typedef map<pair<size_t, size_t>, int>::const_iterator EdgeVisitIter;
        vector<int> countIn( mDependencyMatrix.size(), 0 );
        vector<int> countOut( mDependencyMatrix.size(), 0 );
        vector<double> countRatio( mDependencyMatrix.size(), 0.0 );
        for( EdgeVisitIter edgeIter = edgeVisits.begin(); edgeIter != edgeVisits.end(); ++edgeIter ) {
            countOut[ (*edgeIter).first.first ] += (*edgeIter).second;
            countIn[ (*edgeIter).first.second ] += (*edgeIter).second;
        }
        for( size_t row = 0; row < mDependencyMatrix.size(); ++row ) {
            if( countIn[ row ] > 0 && countOut[ row ] > 0 ) {
                countRatio[ row ] = static_cast<double>( countIn[ row ] )
                    / static_cast<double>( countOut[ row ] );
            }
        }
        
        // Find the most traversed edge with weighting.
        EdgeVisitIter maxVisits = edgeVisits.end();
        for( EdgeVisitIter edgeIter = edgeVisits.begin(); edgeIter != edgeVisits.end(); ++edgeIter ) {
            if( maxVisits == edgeVisits.end() || ( (*edgeIter).second * countRatio[ (*edgeIter).first.first ] ) 
                >= ( (*maxVisits).second * countRatio[ (*maxVisits).first.first ] ) )
            {
                maxVisits = edgeIter;
            }
        }
        if( maxVisits != edgeVisits.end() ) {
            // break this dependency
            mCycleBreaker->breakCycle( *this, (*maxVisits).first.first, (*maxVisits).first.second );
            return true;
        }
        else {
            // The searches did not find any cycles.
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            ILogger& depFinderLog = ILogger::getLogger( "dependency_finder_log" );
            mainLog.setLevel( ILogger::ERROR );
            depFinderLog.setLevel( ILogger::ERROR );
            mainLog << "Could not find a cycle in the dependency matrix." << endl;
            depFinderLog << "Could not find a cycle in the dependency matrix." << endl;
            return false;
        }
    }

    // A cycle was found, but there is no way to break it so output an error.
    ILogger& depFinderLog = ILogger::getLogger( "dependency_finder_log" );
    depFinderLog.setLevel( ILogger::ERROR );
    depFinderLog << "A cycle was found in the dependency matrix but there is "
                 << "no cycle breaker object." << endl;
    return false;
}

/*
 * \brief Updates the count for the number of times an edge was found to be part
 *        of a cycle.
 * \details Performas a depth first search that stops when it finds an edge which
 *          leads back to a node in the current search path.  It will then increase
 *          a count on that edge to indicate it was part of a cycle.  As it unwinds
 *          the search path if an edge was part of a found cycle it will also increase
 *          that edge's count.
 * \param aCurrVertex The current vertex we are checking.
 * \param aPath The current path through the graph which can be used to check for
 *              if a cycle has been found.
 * \param aEdgesVisited A map to keep track of the number of times an edge has been
 *                      visited while being part of a cycle.
 * \return True if a search from aCurrVertex lead to a cycle.
 */
bool DependencyFinder::markCycles( const size_t aCurrVertex, list<size_t>& aPath,
                                   map<pair<size_t, size_t>, int>& aEdgesVisited ) const
{
    // Check to see if the current vertex is already in the path in which case a
    // cycle is found.
    if( find( aPath.begin(), aPath.end(), aCurrVertex ) != aPath.end() ){
        // Increase the count of this edge and return true that this path was part
        // of a cycle.
        pair<size_t, size_t> currEdge( aPath.back(), aCurrVertex );
        map<pair<size_t, size_t>, int>::iterator edgeIter = aEdgesVisited.find( currEdge );
        if( edgeIter == aEdgesVisited.end() ) {
            aEdgesVisited[ currEdge ] = 1;
        }
        else {
            ++(*edgeIter).second;
        }
        
        // A cycle was found
        return true;
    }
    else {
        // Recursively search all dependencies from this vertex for a cycle.
        aPath.push_back( aCurrVertex );
        bool didFindCycle = false;
        for( size_t potentialDep = 0; potentialDep < mDependencyMatrix.size(); ++potentialDep ) {
            if( mDependencyMatrix[ aCurrVertex ][ potentialDep ] ) {
                if( markCycles( potentialDep, aPath, aEdgesVisited ) ) {
                    didFindCycle = true;
                }
            }
        }
        aPath.pop_back();
        
        // Check if at least one of the out-edges from this vertex was part of a
        // cycle.
        if( didFindCycle && !aPath.empty() ) {
            // This path was part of a cycle so increase the count of the edge
            // which led here just once.
            pair<size_t, size_t> currEdge( aPath.back(), aCurrVertex );
            map<pair<size_t, size_t>, int>::iterator edgeIter = aEdgesVisited.find( currEdge );
            if( edgeIter == aEdgesVisited.end() ) {
                aEdgesVisited[ currEdge ] = 1;
            }
            else {
                ++(*edgeIter).second;
            }
        }
        return didFindCycle;
    }
}

/*!
 * \brief Get the name of an object from the matrix index.
 * \param A matrix index to fetch the name for.
 * \note This function is slow as it performs a reverse map lookup and should not
         be used in sections of code which are called frequently.
 * \return The name of the item associated with the index, NO_NAME if it is not
           found.
 */
const string& DependencyFinder::getNameFromIndex( const size_t aIndex ) const {
    // Search the map linearly as this is the reverse lookup.
    for( ObjectIndexMap::const_iterator item = mObjectIndices.begin(); item != mObjectIndices.end(); ++item ){
        // Check if we found the index we are searching for.
        if( item->second == aIndex ){
            return item->first;
        }
    }
    // The index does not exist.
    const static string& NO_NAME = "NO_NAME";
    return NO_NAME;
}
