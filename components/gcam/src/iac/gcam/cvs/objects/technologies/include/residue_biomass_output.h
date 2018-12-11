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

#if !defined( __RESIDUEBIOMASSOUTPUT_H )
#define __RESIDUEBIOMASSOUTPUT_H    // prevent multiple includes

#include "technologies/include/ioutput.h"
#include "util/base/include/value.h"
#include "util/curves/include/cost_curve.h"
#include <vector>

class Curve;

/*!
 * \ingroup objects::biomass
 * \brief A class to output residue biomass supply to the energy
 *        market
 * \details This class contains a set of routines that implement
 *          output of residue biomass supply to an energy market. This object
 *          works for both ag and non-agricultural technologies. If the parameter
 *          mErosCtrl is not specified on input, then this object will work with
 *          non agricultural technologies. NOTE if mErosCtrl is specified for a 
 *          non ag technologies an error will result.
 *
 *
 *   <b>XML specification for ResidueBiomassOutput</b>
 *   - XML name: \c residue-biomass
 *   - Contained by: Technology
 *   - Parsing inherited from class: None.
 *   - Attributes: none
 *   - Elements:
 *   - \c curve-exponent ResidueBiomassOutput::mCostCurve.get/setCurveExponent()
 *   - \c eros-ctrl ResidueBiomassOutput::mErosCtrl
 *   - \c harvest-index ResidueBiomassOutput::mHarvestIndex
 *   - \c mass-conversion ResidueBiomassOutput::mMassConversion
 *   - \c mass-to-energy ResidueBiomassOutput::mMassToEnergy
 *   - \c mid-price ResidueBiomassOutput::mCostCurve.get/setMidprice()
 *
 * \author Kevin Walker
 * \todo Implement XMLDB derived class output so that maximum supply can be written out
 * \todo Figure out way to check for invalid input of mErosCtrl for non ag technology (need to check for existance of land allocator (easy) and existance of specified leaf)
 * \date $ Date $
 * \version $ Revision $
 */
class ResidueBiomassOutput : public IOutput
{
public :
   typedef IOutput parent;

    ResidueBiomassOutput( const std::string& sectorName = std::string() );
    ResidueBiomassOutput( const ResidueBiomassOutput& other );
    virtual ~ResidueBiomassOutput(void);
    ResidueBiomassOutput& operator = ( const ResidueBiomassOutput& other );
    virtual const std::string& getXMLReportingName() const;
    virtual void accept( IVisitor* aVisitor, const int aPeriod ) const;

    virtual IOutput::OutputList calcPhysicalOutput( const double aPrimaryOutput,
                                                    const std::string& aRegionName,
                                                    const ICaptureComponent* aCaptureComponent,
                                                    const int aPeriod ) const;

    virtual ResidueBiomassOutput* clone( void ) const { return new ResidueBiomassOutput( *this ); }

    virtual void completeInit( const std::string& aSectorName, DependencyFinder* aDependencyFinder,
                               const IInfo* aTechInfo, const bool aIsTechOperating );

    virtual double getEmissionsPerOutput( const std::string& aGHGName, const int aPeriod ) const;

    virtual const std::string& getName( ) const { return mName; }

    virtual double getPhysicalOutput( const int aPeriod ) const;

    virtual double getValue( const std::string& aRegionName, const ICaptureComponent* aCaptureComponent,
                             const int aPeriod ) const;

    static const std::string& getXMLNameStatic( );

    virtual void initCalc( const std::string& aRegionName, const std::string& aSectorName, const int aPeriod );

    virtual bool isSameType( const std::string& aType ) const { return getXMLNameStatic().compare( aType ) == 0; }

    virtual void postCalc( const std::string& aRegionName, const int aPeriod );

    virtual void scaleCoefficient( const double aScaler );
    virtual void sendLandAllocator( const ILandAllocator* aLandAllocator, const std::string& aName );
    virtual void setName( const std::string& sectorName ) { mName = sectorName; }
    virtual void setPhysicalOutput( const double aPrimaryOutput, const std::string& aRegionName,
                                    ICaptureComponent* aCaptureComponent, const int aPeriod );

    virtual void setCurrencyOutput( const std::string& aRegionName,  const double aOutput, const int aPeriod ) { }
    virtual double getCurrencyOutput( const int aPeriod ) const { return 0; }
    
    virtual void toDebugXML( const int aPeriod, std::ostream& aOut, Tabs* aTabs ) const;

    virtual void toInputXML( std::ostream& aOut, Tabs* aTabs ) const;
    virtual bool XMLParse( const xercesc::DOMNode* aNode );

    virtual void doInterpolations( const int aYear, const int aPreviousYear,
                                   const int aNextYear, const IOutput* aPreviousInput,
                                   const IOutput* aNextInput );

private :
    const static std::string XML_REPORTING_NAME; //!< tag name for reporting xml db 
    typedef std::vector<Value> value_vector_type;

    //! Physical output by period.
    mutable value_vector_type mPhysicalOutputs;

    /*!
    * Name of the secondary output. Corresponds to a market for this good
    * and a supply sector which supplies this good as its primary output.
    */
    std::string mName;

    //! CO2 emissions coefficient cached from the marketplace.
    Value mCachedCO2Coef;
    
    //! Weak pointer to the land leaf which corresponds to this biomass output
    //! used to save time finding it over and over
    ALandAllocatorItem* mProductLeaf;

    //! The harvest index
    double mHarvestIndex;

    //! The mass per unit area of biomass to be retained to prevent erosion
    double mErosCtrl;

    //! The energy to mass conversion of the crop
    double mMassConversion;
    
    //! Water content of residue
    double mWaterContent;

    //! The mass to energy conversion of the crop residue
    double mMassToEnergy;

    //! Piece-wise linear cost curve 
    std::auto_ptr<Curve> mCostCurve; 

    // These variables are for debugging purposes. Values are written to debug.xml
    // Are made mutable so values can be saved. Are used only for debugging so are not violating const functions.
    //! Mass in crop residue
    mutable double mResMass;

    //! Mass in crop
    mutable double mCropMass;

    //! Residue Available
    mutable double mResAvail;

    //! Mass of biomass to be retained to prevent erosion
    mutable double mMeanErosCtrl;

    //! Max biomass energy supply
    mutable double mMaxBioEnergySupply;

    //! Fraction of max available residue harvested for energy
    mutable double mFractProduced;

};

#endif   // __RESIDUEBIOMASSOUTPUT_H

