#ifndef _ACOMPLEX_EMISSIONS_H_
#define _ACOMPLEX_EMISSIONS_H_
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
 * \file acomplex_emissions.h
 * \ingroup Objects
 * \brief AComplexEmissions class header file.
 * \author Jim Naslund
 */

#include "emissions/include/aghg.h"
#include <string>
#include <xercesc/dom/DOMNode.hpp>

class AEmissionsCoef;
class GhgMAC;
class AEmissionsDriver;

/*! 
 * \ingroup Objects
 * \brief An abstract class that represents complex emissions.
 * \author Jim Naslund
 */
class AComplexEmissions: public AGHG {
public:
    virtual ~AComplexEmissions();
    virtual AComplexEmissions* clone() const = 0;

    virtual void copyGHGParameters( const AGHG* prevGHG );

    virtual void initCalc( const std::string& aRegionName,
                           const IInfo* aLocalInfo,
                           const int aPeriod );

    virtual double getGHGValue( const std::string& aRegionName,
                                const std::vector<IInput*>& aInputs,
                                const std::vector<IOutput*>& aOutputs,
                                const ICaptureComponent* aSequestrationDevice,
                                const int aPeriod ) const;

	virtual void calcEmission( const std::string& aRegionName, 
                               const std::vector<IInput*>& aInputs,
                               const std::vector<IOutput*>& aOutputs,
					           const GDP* aGDP,
					           ICaptureComponent* aSequestrationDevice,
                               const int aPeriod );
    
    virtual void doInterpolations( const int aYear, const int aPreviousYear,
                                   const int aNextYear, const AGHG* aPreviousGHG,
                                   const AGHG* aNextGHG );
protected:

    AComplexEmissions();
    AComplexEmissions( const AComplexEmissions& other );
    AComplexEmissions& operator=( const AComplexEmissions& other );
    
    double maxCntrl; //!<  final control fraction for ghg's
    double mGDPMidControl; //!< User inputed variable- represents midpoint of curve for control function
    double mGDPControlRange; //!< User inputed timescale parameter in control function
    double gdpCap; //!< Saved value for GDP per capita. Needed to adjust control.
    double techDiff; //!< technological change parameter- represents percent reduction in gdp0 per year;
    double adjMaxCntrl; //!< multiplier to maxCntrl, keeping current emissions constant
    double multMaxCntrl; //!< multiplier to maxCntrl -- changes current emissions
    double emAdjust; //!< User inputed adjustment to emissions values(0 to 1)
    double finalEmissCoef; //!< user input final emissions factor that is approached asymptotically
    
    //! Global warming potential of the gas.
    double gwp;

    std::auto_ptr<GhgMAC> ghgMac; //!< Marginal Abatement Cost Curve Object
    std::auto_ptr<AEmissionsCoef> mEmissionsCoef; //!< emissions coefficient delegate

    virtual bool XMLDerivedClassParse( const std::string& nodeName, const xercesc::DOMNode* curr );
    virtual void toInputXMLDerived( std::ostream& out, Tabs* tabs ) const = 0;
    virtual void toDebugXMLDerived( const int period, std::ostream& out, Tabs* tabs ) const = 0;
    double adjustControlParameters( const double gdpCap, const double emissDrive, const double macReduction,
                                    const int period );
    virtual void adjustMaxCntrl(const double GDPcap);

    double controlFunction( const double maxCntrlIn, const double aGDPControlRange, const double aGDPMidControl, const double gdpCapIn );
    double calcTechChange( const int period );

private:
    void copy( const AComplexEmissions& other );
    
};

#endif // _ACOMPLEX_EMISSIONS_H_



