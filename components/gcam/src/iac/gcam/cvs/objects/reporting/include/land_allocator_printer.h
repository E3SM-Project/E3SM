#ifndef _LAND_ALLOCATOR_PRINTER_H_
#define _LAND_ALLOCATOR_PRINTER_H_
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
* \file land_allocator_printer.h
* \ingroup Objects
* \brief LandAllocatorPrinter class header file.
* \author Jim Naslund
*/

#include <stack>
#include <string>

#include "util/base/include/default_visitor.h"
class ALandAllocatorItem;
/*! 
* \ingroup Objects
* \brief A reporting class which outputs a dot graph of the land allocator for a specified region.
* \details A visitor which can output a dot graph of the land allocator for a region specified by
*          the constructor argument. The graph file must be post-processed by the dot processor to
*          create a viewable graph. The graph printer currently will create nodes for each node in
*          the land allocator, and links between parents and children.  Internal nodes are outputted
*          as circles, leaf nodes are outputted as boxes.
* \author Jim Naslund
*/
class LandAllocatorPrinter : public DefaultVisitor {
public:
    explicit LandAllocatorPrinter( const std::string& aRegionToPrint,
                                   std::ostream& aFile,
                                   const bool aPrintValues,
                                   const bool aPrintSpecificRegion );

    void startVisitRegion( const Region* aRegion, const int aPeriod );
    void endVisitRegion( const Region* aRegion, const int aPeriod );
    void startVisitLandNode( const LandNode * aLandNode, const int aPeriod );
    void endVisitLandNode( const LandNode * aLandNode, const int aPeriod );
    void startVisitLandLeaf( const LandLeaf * aLandLeaf, const int aPeriod );
    void openGraph() const;
    void closeGraph() const;
    static void printToFile( const std::string& aRegion,
                             const std::string& aFileName,
                             const ALandAllocatorItem& aLandAllocatorItem );
private:
    //! The file to which to write.
    std::ostream& mFile;
    
    //! Whether we are printing the current region.
    bool mCorrectRegion;
    
    //! The region for which to print graphs.
    const std::string mRegionToPrint;
    
    //! Stores the parent of the current node
    std::stack<std::string> mParent;

    //! Stores the number of nodes outputted, used ensure node names are unique
    int mNumNodes;

    //! Whether or not to print values on the graph
    bool mPrintValues;

    //! Whether or not to only print a certain region
    bool mPrintSpecificRegion;

    void printNode( const ALandAllocatorItem* aLandItem,
                    const int aPeriod, const bool aIsLeaf ) const;

    void printParentChildRelationship( const ALandAllocatorItem* aLandItem ) const;

    std::string makeNameFromLabel( const std::string& aName ) const;
};

#endif // _LAND_ALLOCATOR_PRINTER_H_
