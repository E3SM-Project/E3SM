#ifndef _TREE_ITEM_H_
#define _TREE_ITEM_H_
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
* \file tree_item.h
* \ingroup Objects
* \brief Header file for the TreeItem template class.
* \author Josh Lurz
*/

/*!
 * \brief An enum containing the types of searches (BFS/DFS)
 */
enum SearchType {
    /*!
     * \brief Depth-first search.
     */
    eDFS,
    
    /*!
     * \brief Breadth-first search.
     */
    eBFS
};

/*! 
* \ingroup Objects
* \brief This class defines an interface to any leaf or node in a leaf node
*        structure.
* \author Josh Lurz
*/
template<class T> 
class TreeItem {
public:

    typedef typename std::unary_function<const T*, bool> MatchesFunction;

    /*! \brief Get the number of children which the tree item has.
    * \details Returns the number of children which the tree item has. In the
    *          case of a leaf in the tree this will be zero. These children
    *          should be accessible by getChildAt().
    * \return The number of children the tree item has.
    */
    virtual size_t getNumChildren() const = 0;

    /*! \brief Return the child at the index.
    * \details Queries the tree item for the child at a specific index. If this
    *          index is less than zero or greater than or equal to the size,
    *          this will throw an exception.
    * \return The child at the index.
    */
    virtual const T* getChildAt( const size_t aIndex ) const = 0;

    /*! \brief Return the child at the index.
    * \details Queries the tree item for the child at a specific index. If this
    *          index is less than zero or greater than or equal to the size,
    *          this will throw an exception.
    * \return The child at the index.
    */
    virtual T* getChildAt( const size_t aIndex ) = 0;
};

/*!
 * \brief Perform a search of the tree below this item for a tree item with
 *          with the specified Predicate.
 * \details Performs a search of the tree for the item that satisfies the
 *          predicate. This will return null if the tree item is not found.
 * \param aSearchType An enum representing the type of search to be used.
 * \param aNode The node to start search.
 * \param aIsGoal The predicate.
 * \return The first item found below this item that satisfies the predicate.
 */
template<class T, class Predicate>
const T* findItem( SearchType aSearchType, const T* aNode, Predicate aIsGoal ){
    if( aIsGoal( aNode ) ) {
        return aNode;
    }
    
    // breadth first search implies we should check all of our children now before recursing
    if( aSearchType == eBFS ) {
        for( unsigned int child = 0; child < aNode->getNumChildren(); ++child ){
            if( aIsGoal( aNode->getChildAt( child ) ) ) {
                return aNode->getChildAt( child );
            }
        }
    }
    
    // recursively check each subtree
    for( unsigned int child = 0; child < aNode->getNumChildren(); ++child ){
        const T* ret = findItem( aSearchType, aNode->getChildAt( child ), aIsGoal );
        
        // stop searching if something was found in the current subtree
        if( ret ) {
            return ret;
        }
    }
    
    // no matches in this subtree return null
    return 0;
}

/*!
 * \brief Perform a search of the tree below this item for a tree item with
 *          with the specified Predicate.
 * \details Performs a search of the tree for the item that satisfies the
 *          predicate. This will return null if the tree item is not found.
 * \param aSearchType An enum representing the type of search to be used.
 * \param aNode The node to start search.
 * \param aIsGoal The predicate.
 * \return The first item found below this item that satisfies the predicate.
 */
template<class T, class Predicate>
T* findItem( SearchType aSearchType, T* aNode, Predicate aIsGoal ){
    if( aIsGoal( aNode ) ) {
        return aNode;
    }
    
    // breadth first search implies we should check all of our children now before recursing
    if( aSearchType == eBFS ) {
        for( unsigned int child = 0; child < aNode->getNumChildren(); ++child ){
            if( aIsGoal( aNode->getChildAt( child ) ) ) {
                return aNode->getChildAt( child );
            }
        }
    }
    
    // recursively check each subtree
    for( unsigned int child = 0; child < aNode->getNumChildren(); ++child ){
        T* ret = findItem( aSearchType, aNode->getChildAt( child ), aIsGoal );
        
        // stop searching if something was found in the current subtree
        if( ret ) {
            return ret;
        }
    }
    
    // no matches in this subtree return null
    return 0;
}

#endif // _TREE_ITEM_H_
