#ifndef _GRAPH_PRINTER_H_
#define _GRAPH_PRINTER_H_
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
* \file graph_printer.h
* \ingroup Objects
* \brief GraphPrinter class header file.
* \author Josh Lurz
*/

#include <iosfwd>
#include <sstream>

#include "util/base/include/default_visitor.h"
class Sector;
class Region;
class technology;
/*! 
* \ingroup Objects
* \brief A reporting class which outputs a dot graph for a specified region.
* \details A visitor which can output a dot graph for a region specified by the
*          constructor argument. The graph file must be post-processed by the dot
*          processor to create a viewable graph. The graph printer currently
*          will create nodes for each resource and sector, and links between
*          those based on prices or quantities of consumed inputs.
* \author Josh Lurz
*/
class GraphPrinter : public DefaultVisitor {
public:
    explicit GraphPrinter( const std::string& aRegionToPrint, std::ostream& aFile );
    void startVisitRegion( const Region* aRegion, const int aPeriod );
    void endVisitRegion( const Region* aRegion, const int aPeriod );
    void startVisitResource( const AResource* aResource, const int aPeriod );
    void startVisitSector( const Sector* aSector, const int aPeriod );
	void startVisitFinalDemand( const AFinalDemand* aFinalDemand, const int aPeriod );
	void startVisitTechnology( const Technology* aTechnology, const int aPeriod );
private:
    //! The file to which to write.
    std::ostream& mFile;

    //! Whether we are printing the current region.
    bool mCorrectRegion;

    //! The region for which to print graphs.
    const std::string mRegionToPrint;

    //! The current sector name
    std::string mCurrSectorName;
};

#endif // _GRAPH_PRINTER_H_
