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

#ifndef _DEFAULT_TECHNOLOGY_H_
#define _DEFAULT_TECHNOLOGY_H_
#if defined(_MSC_VER)
#pragma once
#endif

/*! 
 * \file default_technology.h
 * \ingroup Objects
 * \brief The DefaultTechnology class header file.
 * \author Josh Lurz
 */

#include <xercesc/dom/DOMNode.hpp>
#include <iosfwd>
#include "technologies/include/technology.h"
class Tabs;

/*! 
 * \ingroup Objects
 * \brief Technology which represents a simple MiniCAM technology.
 * \details The Technology class is where all fuels are either consumed or
 *          transformed. The DefaultTechnology class is based on a MiniCAM-style
 *          logit representation. This class has options for calibration, fixed
 *          output, capture of emissions, and secondary outputs. The class is
 *          used for both supply and demand of goods.
 * \note DefaultTechnology uses implementations of functions shared with other
 *       Technology subclasses. This is done by including implementations of
 *       abstract functions in Technology. DefaultTechnology's implementations
 *       of those abstract functions call the Technology implementations.
 */
class DefaultTechnology: public Technology
{
public:
	DefaultTechnology( const std::string& aName,
                       const int aYear );

	virtual DefaultTechnology* clone() const;

	static const std::string& getXMLNameStatic();
	
    virtual void completeInit( const std::string& aRegionName,
                               const std::string& aSectorName,
                               const std::string& aSubsectorName,
                               DependencyFinder* aDepFinder,
                               const IInfo* aSubsectorIInfo,
                               ILandAllocator* aLandAllocator );

    virtual void initCalc( const std::string& aRegionName,
                           const std::string& aSectorName,
                           const IInfo* aSubsectorInfo,
                           const Demographic* aDemographics,
                           PreviousPeriodInfo& aPrevPeriodInfo,
                           const int aPeriod );

    virtual void postCalc( const std::string& aRegionName,
                           const int aPeriod );	
	
    virtual void production( const std::string& aRegionName,
                             const std::string& aSectorName, 
		                     double aVariableDemand,
                             double aFixedOutputScaleFactor,
                             const GDP* aGDP,
                             const int aPeriod );
    
    virtual void calcCost( const std::string& aRegionName,
                           const std::string& aSectorName,
		                   const int aPeriod );
	
    virtual double calcShare( const std::string& aRegionName,
                              const std::string& aSectorName, 
		                      const GDP* aGDP,
                              const double aLogitExp,
                              const int aPeriod ) const;
protected:

	virtual bool XMLDerivedClassParse( const std::string& aNodeName,
                                       const xercesc::DOMNode* aCurr );

	virtual void toInputXMLDerived( std::ostream& aOut,
                                    Tabs* aTabs ) const;

	virtual void toDebugXMLDerived( const int aPeriod,
                                    std::ostream& aOut,
                                    Tabs* aTabs ) const;

	virtual const std::string& getXMLName() const;
	
    virtual double getTotalInputCost( const std::string& aRegionName,
                                      const std::string& aSectorName,
		                              const int aPeriod ) const;
};

#endif // _DEFAULT_TECHNOLOGY_H_
