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

#ifndef _S_CURVE_SHUTDOWN_DECIDER_H_
#define _S_CURVE_SHUTDOWN_DECIDER_H_
#if defined(_MSC_VER)
#pragma once
#endif

/*! 
 * \file s_curve_shutdown_decider.h
 * \ingroup Objects
 * \brief The S_CurveShutdownDecider header file.
 * \author Patrick Luckow
 */
#include "technologies/include/ishutdown_decider.h"

#include <string>
struct ProductionFunctionInfo;
class Tabs;

/*! 
 * \ingroup Objects
 * \brief This object shuts down production for a vintage using a single
 *        s-curve with zero production after the technology lifetime.
 
 // shutdown the production with the function 1/(1+e^(steepness*(years active-halflife))).
 // this prevents a lot of retirement in the early years, and gives the vintage a long tail;
 // All remaining vintage is cut off at the read in lifetime
 
 * \details This object uses a s-curve decay to reduce the production of a
 *          vintage based on a read-in steepness and half life and the number of periods
 *          a vintage has been active. The class calculates a scaling factor which
 *          is used to decrease production as the maximum lifetime of the
 *          technology is approached. The technology will remove the tail of the
 *          production by shutting itself down when it reaches the end of its
 *          lifetime. The coefficient is calculated as 1/(1+e^(steepness*(years active-halflife)))
 *         
 *
 *          <b>XML specification for S_CurveShutdownDecider</b>
 *          - XML name: \c s-curve-shutdown-decider
 *          - Contained by: Technology
 *          - Parsing inherited from class: None
 *          - Attributes: None
 *          - Elements:
 *              - \c shutdown-rate S_CurveShutdownDecider::mShutdownRate
 *     
 * \author Josh Lurz
 */
class S_CurveShutdownDecider: public IShutdownDecider
{
    friend class ShutdownDeciderFactory;
public:
    // IParsedComponent methods.
    virtual S_CurveShutdownDecider* clone() const;

    virtual bool isSameType( const std::string& aType ) const;

    virtual const std::string& getName() const;

    virtual bool XMLParse( const xercesc::DOMNode* aNode );

    virtual void toInputXML( std::ostream& aOut,
                             Tabs* aTabs ) const;

    virtual void toDebugXML( const int aPeriod,
                             std::ostream& aOut,
                             Tabs* aTabs ) const;
    
    // IShutdownDecider methods.
    virtual double calcShutdownCoef( const ProductionFunctionInfo* aFuncInfo,
                                     const double aCalculatedProfits,
                                     const std::string& aRegionName,
                                     const std::string& aSectorName,
                                     const int aInstallationYear,
                                     const int aPeriod ) const;
private:
    //! The steepness of the curve. This rate may be zero
    //! which is the equivalent to not reading in the s-curve shutdown decider.
    double mSteepness;
    
    double mHalfLife;


    S_CurveShutdownDecider();

    static const std::string& getXMLNameStatic();
};

#endif // _S_CURVE_SHUTDOWN_DECIDER_H_
