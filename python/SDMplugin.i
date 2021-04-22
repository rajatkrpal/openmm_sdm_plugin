%module SDMplugin

%import(module="simtk.openmm") "swig/OpenMMSwigHeaders.i"

%include "std_string.i"

/*
 * The following lines are needed to handle std::vector.
 * Similar lines may be needed for vectors of vectors or
 * for other STL types like maps.
 */

%include "std_vector.i"
namespace std {
  %template(vectord) vector<double>;
  %template(vectori) vector<int>;
};

%{
#include "LangevinIntegratorSDM.h"
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "OpenMMDrude.h"
#include "openmm/RPMDIntegrator.h"
#include "openmm/RPMDMonteCarloBarostat.h"
%}


/*
 * The code below strips all units before the wrapper
 * functions are called. This code also converts numpy
 * arrays to lists.
*/

%pythoncode %{
from simtk import openmm as mm
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import math
%}

%include "SDMUtils.py"

/* strip the units off of all input arguments */
//%pythonprepend %{
//try:
//    args=mm.stripUnits(args)
//except UnboundLocalError:
//    pass
//%}


/*
 * Add units to function inputs.

%pythonappend SDMPlugin::LangevinIntegratorSDM::LangevinIntegratorSDM(double temperature, double frictionCoeff, double stepSize, int ligId, double lamdaId, int atom1, int atom2, double kf, double r0) %{
    val[0] = unit.Quantity(val[0], unit.kelvin)
    val[1] = unit.Quantity(val[1], 1/unit.picosecond)
    val[2] = unit.Quantity(val[2], 1/unit.picosecond)
    val[4] = unit.Quantity(val[4], None)
    val[7] = unit.Quantity(val[7], unit.kilojoule_per_mole/(unit.nanometer*unit.nanometer))
    val[8] = unit.Quantity(val[8], unit.nanometer)
%}
*/



namespace SDMPlugin{

class LangevinIntegratorSDM : public OpenMM::Integrator {
public:
/*
   %apply double INPUT {double temperature};
   %apply double INPUT {double frictionCoeff};
   %apply double INPUT {double stepSize};
   %apply int INPUT {int ligId};
   %apply double INPUT {double lamdaId};
   %apply int INPUT {int atom1};
   %apply int INPUT {int atom2};
   %apply double INPUT {double kf};
   %apply double INPUT {double r0};
*/
   LangevinIntegratorSDM(double temperature, double frictionCoeff, double stepSize,
			   int nParticles);
   void setLambda(double lambdac) ;
   double getLambda() const ;
   double getTemperature() const ;
   void setTemperature(double temp) ;
   double getFriction() const ;
   void setFriction(double coeff) ;
   int getRandomNumberSeed() const ;
   void setRandomNumberSeed(int seed) ;
   double getBindE() const ;
   void setBindE(double be) ;
   double getPotEnergy(void) const ;
   void setPotEnergy(double e) ;
   double getUmax(void) const ;
   void setUmax(double um) ;
   double getAcore(void) const ;
   void setAcore(double a) ;
   double getUbcore(void) const ;
   void setUbcore(double a) ;
   void setBiasMethod(int method) ;
   int getBiasMethod() const ;
   void setSoftCoreMethod(int method) ;
   int getSoftCoreMethod() const ;
   void setGamma(double gammat) ;
   double getGamma() const ;
   void setWBcoeff(double wbcoeff_t) ;
   double getWBcoeff() const ;
   void setW0coeff(double w0coeff_t) ;
   double getW0coeff() const ;
   void setLambda1(double lambda1_t) ;
   double getLambda1() const ;
   void setLambda2(double lambda2_t) ;
   double getLambda2() const ;
   void setAlpha(double alpha_t) ;
   double getAlpha() const ;
   void setU0(double u0_t) ;
   double getU0() const ;
   void setNoneqtmax(double noneq_tmax);
   double getNoneqtmax() const;
   int getNonEquilibrium() const;
   void setNoneqWorkvalue(double noneq_work);
   double getNoneqWorkvalue() const;
   void setlambda1Slope(double ml1);
   double getlambda1Slope() const;
   void setlambda2Slope(double ml2);  
   double getlambda2Slope() const;
   void setu0Slope(double mu0);  
   double getu0Slope() const;
   void setw0Slope(double mw0);  
   double getw0Slope() const;
   void setlambda1intercept(double bl1);  
   double getlambda1intercept() const;
   void setlambda2intercept(double bl2);  
   double getlambda2intercept() const;
   void setu0intercept(double bu0);  
   double getu0intercept() const;
   void setw0intercept(double bw0);
   double getw0intercept() const;
   void setDisplacement(int atom, double dx, double dy, double dz);
   OpenMM::Vec3 getDisplacement(int atom) const;
   virtual void step(int steps) ;
};

}


