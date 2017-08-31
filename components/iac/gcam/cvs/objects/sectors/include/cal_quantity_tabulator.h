#ifndef _CAL_QUANTITY_TABULATOR_H_
#define _CAL_QUANTITY_TABULATOR_H_
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
* \file cal_quantity_tabulator.h
* \ingroup Objects
* \brief CalQuantityTabulator class header file.
* \author Josh Lurz
*/

#include "util/base/include/default_visitor.h"
#include <map>
#include <string>

/*! 
* \ingroup Objects
* \brief A visitor which tabulates the calibrated and fixed supplies of a
*        region.
* \details This visitor when passed to a Region will tabulate individually the
*          fixed and calibrated supplies of all goods in the region. It will
*          also keep track of whether there are exist variable supplies of each
*          good in the region. This should be used instead of
*          Sector::getCalOutput if secondary outputs need to be included in the
*          sums. The visitor should be used by passing this object as an
*          argument to the accept function, and then access the information
*          through the getSupplyInfo function.
*
* \author Josh Lurz
*/
class CalQuantityTabulator : public DefaultVisitor {
public:
    /*
     * \brief A structure containing calibration information for a single good.
     */
    struct CalInfo {
    public:
        //! Total calibrated supply.
        double mCalQuantity;

        //! Total fixed supply.
        double mFixedQuantity;

        //! Whether all Technologies are calibrated or fixed.
        bool mAllFixed;

        //! Constructor which initializes all variables.
        CalInfo(): mCalQuantity( 0 ), mFixedQuantity( 0 ), mAllFixed( true ){}
    };

    //! Type of the map containing good names and their associated calibration
    //! information.
    typedef std::map<std::string, CalInfo> CalInfoMap;

    CalQuantityTabulator( const std::string& aRegionName );

    // Documentation for visitor methods is inherited.
    virtual void startVisitRegion( const Region* aRegion,
                                   const int aPeriod );
    
    virtual void startVisitResource( const AResource* aResource,
                                     const int aPeriod );

    virtual void startVisitSector( const Sector* aSector,
                                   const int aPeriod );

    virtual void endVisitSector( const Sector* aSector,
                                 const int aPeriod );

    virtual void startVisitSubsector( const Subsector* aSubsector,
                                      const int aPeriod );

    virtual void endVisitSubsector( const Subsector* aSubsector,
                                    const int aPeriod );

    virtual void startVisitTechnology( const Technology* aTechnology,
                                       const int aPeriod );

    virtual void endVisitTechnology( const Technology* aTechnology,
                                     const int aPeriod );
    
    virtual void startVisitOutput( const IOutput* aOutput,
                                   const int aPeriod );

    // Non visitor interface methods.
    void setApplicableSectorType( const std::string& aSectorType );

    const CalInfoMap& getSupplyInfo() const;
private:
    //! Map of output name to a struct containing the amount of calibrated
    //! supply, fixed supply, and whether all technologies were calibrated or
    //! fixed.
    CalInfoMap mCalSupplies;

    //! If set, the type of sector for which to tabulate calibrated and fixed
    //! output. If this is not set(the default), all sector types will be
    //! calibrated.
    std::string mSectorType;

    //! Name of the Region the calibration and fixed outputs are currently being
    //! summed for.
    std::string mCurrentRegionName;

    //! Name of the sector currently being tabulated.
    std::string mCurrentSectorName;

    //! Current primary output calibration or fixed value. This is set to -1 if
    //! the Technology is variable.
    double mCurrentOutput;

    //! State of the current Technology.
    enum CalibrationState {
        //! Variable output.
        eVariable,

        //! Calibrated output.
        eCalibrated,

        //! Fixed output.
        eFixed,

        //! Unknown status.
        eUnknown
    };

    //! Current technology calibration state.
    CalibrationState mTechState;

    //! Current subsector calibration state.
    CalibrationState mSubsectorState;

    //! Whether the visitor is currently within a sector it should tabulate.
    bool mShouldTabulateSector;
};

#endif // _CAL_QUANTITY_TABULATOR_H_
