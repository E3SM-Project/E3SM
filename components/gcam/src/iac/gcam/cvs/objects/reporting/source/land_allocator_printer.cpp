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
 * \file land_allocator_printer.cpp
 * \ingroup Objects
 * \brief The LandAllocatorPrinter class source file.
 *
 * \author Jim Naslund
 */

#include "util/base/include/definitions.h"

#include <iomanip>
#include <boost/lexical_cast.hpp>

#include "reporting/include/land_allocator_printer.h"
#include "util/base/include/util.h"
#include "containers/include/region.h"
#include "land_allocator/include/land_leaf.h"
#include "land_allocator/include/land_node.h"
#include "util/base/include/auto_file.h"

using namespace std;

/*! 
 * \brief Default Constructor
 * \param aRegionToPrint region to print.
 * \param aFile file to print to.
 * \param aPrintValues Where to print land allocation values on the graph.
 * \param aPrintSpecificRegion Whether to print only the region passed in.
 */
LandAllocatorPrinter::LandAllocatorPrinter( const std::string& aRegionToPrint, std::ostream& aFile,
                                            const bool aPrintValues,
                                            const bool aPrintSpecificRegion ):
mFile( aFile ),
mCorrectRegion( false ),
mRegionToPrint( aRegionToPrint ),
mNumNodes( 0 ),
mPrintValues( aPrintValues ),
mPrintSpecificRegion( aPrintSpecificRegion )
{
    // Imbue the output stream with the default locale from the user's machine.
    // This is done so thousands separators will be outputted.
    mFile.imbue( locale( "" ) );
}

/*!
 * \brief Begin visiting a region with the graph printer.
 * \details Opens the graph and prints the header.
 * \param aRegion Region to visit.
 * \param aPeriod Period for which to visit.
 */
void LandAllocatorPrinter::startVisitRegion( const Region* aRegion, const int aPeriod ){
    // Check if this is the region to print.
    if( aRegion->getName() == mRegionToPrint ){
        mCorrectRegion = true;
        // Print the graph header.
        mFile << "digraph " << util::replaceSpaces( aRegion->getName() ) << " {" << endl;
    }
    else {
        // Don't print this region.
        mCorrectRegion = false;
    }
}

/*!
 * \brief End visiting a region with the graph printer.
 * \details Closes the graph.
 * \param aRegion Region to visit.
 * \param aPeriod Period for which to visit.
 */
void LandAllocatorPrinter::endVisitRegion( const Region* aRegion, const int aPeriod ){
    if( mCorrectRegion ){
        // Now close the graph.
        mFile << "}" << endl << endl;
    }
}

/*!
 * \brief Begin visiting a landNode.
 * \details Outputs the node, its parent-child relationship if application and
 * \details pushes its name onto the parent stack.
 * \param aLandNode landNode to visit.
 * \param aPeriod Period for which to visit.
 */
void LandAllocatorPrinter::startVisitLandNode( const LandNode *aLandNode, const int aPeriod ){
    if( mPrintSpecificRegion && !mCorrectRegion ){
        return;
    }

    // Print the node.
    printNode( aLandNode, aPeriod, false );

    if( !mParent.empty() ){
        // Print the parent link.
        printParentChildRelationship( aLandNode );
    }
    mParent.push( util::replaceSpaces( makeNameFromLabel( aLandNode->getName() ) ) );
    mNumNodes++;
}

/*!
 * \brief End visiting a land node.
 * \details Pops its own name off the parent names stack.
 * \param aLandNode node to end visiting.
 * \param aPeriod Period for which to end visiting.
 */
void LandAllocatorPrinter::endVisitLandNode( const LandNode *aLandNode, const int aPeriod ){
    if( mPrintSpecificRegion && !mCorrectRegion ){
        return;
    }
    mParent.pop();
}

/*!
 * \brief Begin visiting a landLeaf.
 * \details Outputs the node, its parent-child relationship (it always has a parent)
 * \param aLandNode landLeaf to visit.
 * \param aPeriod Period for which to visit.
 */
void LandAllocatorPrinter::startVisitLandLeaf( const LandLeaf *aLandLeaf, const int aPeriod ){
    if( mPrintSpecificRegion && !mCorrectRegion ){
        return;
    }
    printNode( aLandLeaf, aPeriod, true );
    printParentChildRelationship( aLandLeaf );
    mNumNodes++;
}

/*!
 * \brief Outputs a node.
 * \details Outputs a node.  Creates a unique name for the node.
 * \param aLandItem the land item to print.
 * \param aPeriod The period.
 * \param aIsLeaf Whether or not the node is a leaf.
 */
void LandAllocatorPrinter::printNode( const ALandAllocatorItem* aLandItem,
                                      const int aPeriod, const bool aIsLeaf ) const{
    string name = aLandItem->getName();
    string nameStripped = util::replaceSpaces( name );
    mFile << "\t" << makeNameFromLabel( name ) << "[label=" << "\"" << nameStripped;
    if( mPrintValues ){
        mFile << "\\n" << setiosflags( ios::fixed ) << setprecision( 0 ) 
              << aLandItem->getLandAllocation( name, aPeriod );
    }
    mFile << "\"";
    if( aIsLeaf ){
        mFile << ", shape=box";
    }
    mFile << "];" << endl;
}

/*!
 * \brief Outputs the parent-child relationship between two nodes.
 * \details Outputs the link from a parent node to its child node.
 * \param aLandItem The child in the relationship.
 */
void LandAllocatorPrinter::printParentChildRelationship( const ALandAllocatorItem* aLandItem ) const{
    mFile << "\t" << mParent.top() << "->" 
          << util::replaceSpaces( makeNameFromLabel( aLandItem->getName() ) ) << ";" << endl;
}

/*!
 * \brief Creates a unique name for the node.
 * \details Prevents duplicate node names from causing problems.
 * \param aName name to make unique.
 * \return unique name for node.
 */
string LandAllocatorPrinter::makeNameFromLabel( const string& aName ) const{
    return aName + boost::lexical_cast<string>( mNumNodes );
}

/*!
 * \brief Opens the graph.
 * \details This method only needs to be called when the graph is being
 *          used for debugging purposes and is being outputted from
 *          inside a land allocator object.
 */
void LandAllocatorPrinter::openGraph() const {
    mFile << "digraph la {" << endl;
}

/*!
 * \brief Closes the graph.
 * \details This method only needs to be called when the graph is being
 *          used for debugging purposes and is being outputted from
 *          inside a land allocator object.
 */
void LandAllocatorPrinter::closeGraph() const {
    mFile << "}" << endl << endl;
}

/*!
 * \brief Static function for printing a land allocator item to a file.
 * \details Static function for printing a land allocator item to a file.  This
 *          function is for debugging purposes and can be used to print a land
 *          allocator to a file at any point during the model's execution.
 * \param aRegion Region to print (not used).
 * \param aFileName Filename to print to.
 * \param aLandAllocatorItem Item to print.
 * \note There is currently no way to control which region this function prints.
 */
void printToFile( const string& aRegion, const string& aFileName,
                  const ALandAllocatorItem& aLandAllocatorItem ){
    const Modeltime* modeltime = scenario->getModeltime();
    for( int period = 0; period < modeltime->getmaxper(); period++ ) {
        AutoOutputFile landAllocatorStream( boost::lexical_cast<string>( period ) + aFileName );
        LandAllocatorPrinter landAllocatorPrinter( aRegion, *landAllocatorStream,
                                                   true, false );
        landAllocatorPrinter.openGraph();
        aLandAllocatorItem.accept( &landAllocatorPrinter, period );
        landAllocatorPrinter.closeGraph();
        landAllocatorPrinter.finish();
    }
}
