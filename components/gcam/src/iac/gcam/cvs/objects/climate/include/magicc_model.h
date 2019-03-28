#ifndef _MAGICC_MODEL_H_
#define _MAGICC_MODEL_H_
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
* \file magicc_model.h
* \ingroup Objects
* \brief The MagiccModel header file.
* \author Josh Lurz
*/

#include <map>
#include <string>
#include <vector>
#include "climate/include/iclimate_model.h"

class Modeltime;
class IVisitor;

/*! 
* \ingroup Objects
* \brief An implementation of the IClimateModel interface using the MAGICC
*        climate module.
* \details The MagiccModel performs climate calculations by passing data to and
*          from the MAGICC Fortran module. No climate calculating code is
*          contained in the C++ MagiccModel code. This wrapper is responsible
*          for reading in a set of default gas emissions for each gas by period,
*          overriding those with values from the model where calculated, and
*          interpolating them into a set of inputs for MAGICC. It then writes
*          those values to a file and calls MAGICC to calculate climate
*          parameters. A subset of those output can then be written by this
*          wrapper to the database and a CSV file.
* \note It is possible to run MAGICC using the Objects framework without running
*       the economic model. This is done by reading in a scenario container with
*       only a modeltime object and an empty world object. It will run off the
*       values in the input_gases.emk file.
* \author Josh Lurz
*/

class MagiccModel: public IClimateModel {
public:
    MagiccModel( const Modeltime* aModeltime );

    virtual void completeInit( const std::string& aScenarioName );
    
    static const std::string& getXMLNameStatic();
    virtual void XMLParse( const xercesc::DOMNode* node );
    virtual void toInputXML( std::ostream& out, Tabs* tabs ) const;
    virtual void toDebugXML( const int period, std::ostream& out, Tabs* tabs ) const;
    
    virtual bool setEmissions( const std::string& aGasName,
                               const int aPeriod,
                               const double aEmission );
	
	virtual bool setLUCEmissions( const std::string& aGasName,
							  const int aYear,
							  const double aEmission );
    
    virtual double getEmissions( const std::string& aGasName,
                                 const int aYear ) const;

    virtual bool runModel();

    virtual double getConcentration( const std::string& aGasName,
                                     const int aYear ) const;

    virtual double getTemperature( const int aYear ) const;
    
    virtual double getForcing( const std::string& aGasName,
                               const int aYear ) const;
    
    virtual double getTotalForcing( const int aYear ) const;

    double getNetTerrestrialUptake( const int aYear ) const;

    double getNetOceanUptake( const int aYear ) const;

    double getNetLandUseChangeEmission( const int aYear ) const;

    virtual int getCarbonModelStartYear() const;

    virtual void printFileOutput() const;
    virtual void printDBOutput() const;
    virtual void accept( IVisitor* aVisitor, const int aPeriod ) const;
    
    static const std::string& getnetDefor80sName();

private:

    bool isValidClimateModelYear( const int aYear ) const;

    int getInputGasIndex( const std::string& aGasName ) const;
    static unsigned int getNumInputGases();
    void readFile();
    void overwriteMAGICCParameters( );
    static int getNumAdditionalGasPoints();

    //! A fixed list of the gases Magicc reads in.
    static const std::string sInputGasNames[];
    
    //! A fixed list of the units for gases Magicc reads in.
    static const std::string sInputGasUnits[];

    //! A map of the gases Magicc can report out.
    std::map<std::string,int> mOutputGasNameMap; 

    //! Return value of getGasIndex if it cannot find the gas.
    static const int INVALID_GAS_NAME = -1;
    
    //! MAGICC critcal start year in gas.emk that must be present
    static const int GAS_EMK_CRIT_YEAR;

    //! Emissions levels by gas and period.
    std::vector<std::vector<double> > mEmissionsByGas;
	
	//! LUC CO2 Emissions by year.
	std::vector<double> mLUCEmissionsByYear;

    //! Name of the scenario.
    std::string mScenarioName;

    //! Name of a GHG input file to use.
    std::string mGHGInputFileName;

    //! A reference to the scenario's modeltime object.
    const Modeltime* mModeltime;

    //! Whether the climate model output is updated.
    bool mIsValid;

    //! Climate Sensitivity.
    double mClimateSensitivity;

    //! Soil Feedback Factor (MAGICC Parameter btSoil)
    double mSoilTempFeedback;

    //! Humus Feedback Factor (MAGICC Parameter btHumus)
    double mHumusTempFeedback;

    //! GPP Feedback Factor (MAGICC Parameter btGPP)
    double mGPPTempFeedback;

    //! 1980s Ocean Uptake (MAGICC Parameter FUSER)
    double mOceanCarbFlux80s;

    //! 1980s net terrestrial Deforestation (MAGICC Parameter DUSER)
    double mNetDeforestCarbFlux80s;

    //! The year the carbon model should start running.
    int mCarbonModelStartYear;
};

#endif // _MAGICC_MODEL_H_
