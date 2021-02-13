/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2012 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "LangevinIntegratorSDM.h" 
#include "openmm/Context.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/OSRngSeed.h"
#include "openmm/kernels.h"
#include "SDMKernels.h"
#include <string>
#include <vector>
#include <iostream>

using namespace SDMPlugin;
using namespace OpenMM;
using std::string;
using std::vector;

LangevinIntegratorSDM::LangevinIntegratorSDM(double temperature, double frictionCoeff, double stepSize, int nParticles_t) {

      setTemperature(temperature);
      setFriction(frictionCoeff);
      setStepSize(stepSize);
      setConstraintTolerance(1e-5);
      setRandomNumberSeed(osrngseed());

      setUmax(200.0);
      setAcore(1./4.0);
      setUbcore(0.0);


      setSoftCoreMethod(NoSoftCoreMethod);
      setBiasMethod(LinearMethod);
      
      //by default prepare the system in bound state for all methods

      setLambda(1.0);
      
      setGamma(0.0);
      setWBcoeff(1.0);
      setW0coeff(0.0);
      
      setLambda1(1.0);
      setLambda2(1.0);
      setAlpha(1.0);
      setU0(0.0);

      setNonEquilibrium(0);
      work_value = 0.0;

      //set up displacement maps
      nParticles = nParticles_t;
      displ.clear();
      for (int i=0;i<nParticles;i++) displ.push_back(Vec3(0.,0.,0.));
      
}


void LangevinIntegratorSDM::initialize(ContextImpl& contextRef) {
    if (owner != NULL && &contextRef.getOwner() != owner)
        throw OpenMMException("This Integrator is already bound to a context");
    context = &contextRef;

    const System &system = contextRef.getSystem();
    int numforces = system.getNumForces();
    for (int i=0; i<numforces; i++){
      const Force &force = system.getForce(i);
      int group = force.getForceGroup();
      if (!(group == 1 || group == 2)){
	throw OpenMMException("The SDM integrator requires all forces to be either in force group 1 (bonded) or group 2 (nonbonded)");
      }
    }

    owner = &contextRef.getOwner();
    kernel = context->getPlatform().createKernel(IntegrateLangevinStepSDMKernel::Name(), contextRef);
    kernel.getAs<IntegrateLangevinStepSDMKernel>().initialize(contextRef.getSystem(), *this);
}

void LangevinIntegratorSDM::cleanup() {
    kernel = Kernel();
}

vector<string> LangevinIntegratorSDM::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(IntegrateLangevinStepSDMKernel::Name());
    return names;
}

double LangevinIntegratorSDM::computeKineticEnergy() {
    return kernel.getAs<IntegrateLangevinStepSDMKernel>().computeKineticEnergy(*context, *this);
}

/**
 returns the soft-core capped binding energy + derivative
*/
double LangevinIntegratorSDM::SoftCoreF(double u, double umax, double a, double ub, double& fp){
  if(u <= ub){
    fp = 1.;
    return u;
  }
 
  if(softcore_method == NoSoftCoreMethod){
    fp = 1.;
    return u;
  }else if(softcore_method == TanhMethod){
    double x = (u-ub)/umax;
    double t = tanh(x);
    fp = 1. - t*t;
    return umax*t + ub;
  }else if(softcore_method == RationalMethod){
    double gu = (u-ub)/(a*(umax-ub)); //this is y/alpha
    double zeta = 1. + 2.*gu*(gu + 1.) ;
    double zetap = pow( zeta , a );
    double s = 4.*(2.*gu + 1.)/zeta;
    fp = s*zetap/pow(1.+zetap,2);
    return (umax-ub)*(zetap - 1.)/(zetap + 1.) + ub;
  }else{
    throw OpenMMException("Unknown soft core method");
  }
}



void LangevinIntegratorSDM::step(int steps) {
  bool do_energy = true;
  
    for (int i = 0; i < steps; ++i) {

	//evaluate energy and force for state 1, 4 = evaluate force group 2 (non-bonded)
	context->updateContextState();
	double State1Energy = context->calcForcesAndEnergy(true, do_energy, 4);
	//store coordinates and forces of state 1
	kernel.getAs<IntegrateLangevinStepSDMKernel>().SaveState1(*context, *this);

	//make state 2
	kernel.getAs<IntegrateLangevinStepSDMKernel>().MakeState2(*context, *this);
	//evaluate energy and force for state 2, 4 = evaluate force group 2 (non-bonded)
	context->updateContextState();
	double State2Energy = context->calcForcesAndEnergy(true, do_energy, 4);
	//store forces of state 2
	kernel.getAs<IntegrateLangevinStepSDMKernel>().SaveState2(*context, *this);
	
	//restore original state
	kernel.getAs<IntegrateLangevinStepSDMKernel>().RestoreState1(*context, *this);
	//evaluate bonded and restraint forces, 2 = force group 1 (bonded and restraints)
	context->updateContextState();
	double BondAndRestraintEnergy = context->calcForcesAndEnergy(true, do_energy, 2);

	//perform step with hybrid forces
	//std::cout << State1Energy << " " << State2Energy << " " << BondAndRestraintEnergy << std::endl;
        kernel.getAs<IntegrateLangevinStepSDMKernel>().execute(*context, *this,
						State1Energy, State2Energy, BondAndRestraintEnergy);
    }
}
