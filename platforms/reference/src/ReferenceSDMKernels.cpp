/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2013 Stanford University and the Authors.      *
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

#include "ReferenceSDMKernels.h"
#include "openmm/reference/ReferenceConstraints.h"
#include "ReferenceStochasticDynamicsSDM.h"
#include "openmm/Context.h"
#include "openmm/System.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/Integrator.h"
#include "openmm/OpenMMException.h"
#include "openmm/reference/SimTKOpenMMUtilities.h"
#include <cmath>
#include <iostream>
#include <limits>

using namespace SDMPlugin;
using namespace OpenMM;
using namespace std;

static int** allocateIntArray(int length, int width) {
    int** array = new int*[length];
    for (int i = 0; i < length; ++i)
        array[i] = new int[width];
    return array;
}

static RealOpenMM** allocateRealArray(int length, int width) {
    RealOpenMM** array = new RealOpenMM*[length];
    for (int i = 0; i < length; ++i)
        array[i] = new RealOpenMM[width];
    return array;
}

static void disposeIntArray(int** array, int size) {
    if (array) {
        for (int i = 0; i < size; ++i)
            delete[] array[i];
        delete[] array;
    }
}

static void disposeRealArray(RealOpenMM** array, int size) {
    if (array) {
        for (int i = 0; i < size; ++i)
            delete[] array[i];
        delete[] array;
    }
}

static vector<Vec3>& extractPositions(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<Vec3>*) data->positions);
}

static vector<Vec3>& extractVelocities(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<Vec3>*) data->velocities);
}

static vector<Vec3>& extractForces(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<Vec3>*) data->forces);
}

static Vec3& extractBoxSize(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *(Vec3*) data->periodicBoxSize;
}

static ReferenceConstraints& extractConstraints(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *(ReferenceConstraints*) data->constraints;
}

/**
 * Compute the kinetic energy of the system, possibly shifting the velocities in time to account
 * for a leapfrog integrator.
 */
static double computeShiftedKineticEnergy(ContextImpl& context, vector<double>& masses, double timeShift) {
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& velData = extractVelocities(context);
    //vector<Vec3>& forceData = extractForces(context);
    int numParticles = context.getSystem().getNumParticles();
    
    // Compute the shifted velocities.
    //silently ignore time shift (forces are not reliable here because they have not been alchemically hybridized)
    vector<Vec3> shiftedVel(numParticles);
    for (int i = 0; i < numParticles; ++i) {
      shiftedVel[i] = velData[i]; //+forceData[i]*(timeShift/masses[i]);
    }
    
    // Apply constraints to them.
    
    vector<double> inverseMasses(numParticles);
    for (int i = 0; i < numParticles; i++)
        inverseMasses[i] = (masses[i] == 0 ? 0 : 1/masses[i]);
    extractConstraints(context).applyToVelocities(posData, shiftedVel, inverseMasses, 1e-4);
    
    // Compute the kinetic energy.
    
    double energy = 0.0;
    for (int i = 0; i < numParticles; ++i)
        if (masses[i] > 0)
            energy += masses[i]*(shiftedVel[i].dot(shiftedVel[i]));
    return 0.5*energy;
}


ReferenceIntegrateLangevinStepSDMKernel::~ReferenceIntegrateLangevinStepSDMKernel() {
  if (dynamics)
    delete dynamics;
}


void ReferenceIntegrateLangevinStepSDMKernel::initialize(const System& system, const LangevinIntegratorSDM& integrator) {
  int numParticles = system.getNumParticles();
  LigParticle = integrator.getLigParticle();
  int numLigParticles = LigParticle.size();
  masses.resize(numParticles);
  BoundForces.resize(numParticles);
  UnboundForces.resize(numParticles);
  LigandCoordinates.resize(numLigParticles);
 
  for (int i = 0; i < numParticles; ++i)
    masses[i] = static_cast<RealOpenMM>(system.getParticleMass(i));
  SimTKOpenMMUtilities::setRandomNumberSeed((unsigned int) integrator.getRandomNumberSeed()); 
}


void ReferenceIntegrateLangevinStepSDMKernel::SaveBound(ContextImpl& context, const LangevinIntegratorSDM& integrator) {
  int numParticles = context.getSystem().getNumParticles();
  int numLigParticles = LigParticle.size();
  vector<Vec3>& posData = extractPositions(context);
  vector<Vec3>& forceData = extractForces(context);

  for(int p = 0; p<numParticles; p++){
    BoundForces[p][0] = forceData[p][0];
    BoundForces[p][1] = forceData[p][1];
    BoundForces[p][2] = forceData[p][2];
  }
  for(int i = 0; i<numLigParticles; i++){
    int p = LigParticle[i];
    LigandCoordinates[i][0] = posData[p][0];
    LigandCoordinates[i][1] = posData[p][1];
    LigandCoordinates[i][2] = posData[p][2];
  }
}

void ReferenceIntegrateLangevinStepSDMKernel::SaveUnbound(ContextImpl& context, const LangevinIntegratorSDM& integrator) {
  int numParticles = context.getSystem().getNumParticles();
  int numLigParticles = LigParticle.size();
  vector<Vec3>& forceData = extractForces(context);

  for(int p = 0; p<numParticles; p++){
    UnboundForces[p][0] = forceData[p][0];
    UnboundForces[p][1] = forceData[p][1];
    UnboundForces[p][2] = forceData[p][2];
  }
}

void ReferenceIntegrateLangevinStepSDMKernel::RestoreBound(ContextImpl& context, const LangevinIntegratorSDM& integrator) {
  int numLigParticles = LigParticle.size();
  vector<Vec3>& posData = extractPositions(context);

  for(int i = 0; i<numLigParticles; i++){
    int p = LigParticle[i];
    posData[p][0] = LigandCoordinates[i][0];
    posData[p][1] = LigandCoordinates[i][1];
    posData[p][2] = LigandCoordinates[i][2];
  }
}

void ReferenceIntegrateLangevinStepSDMKernel::MakeUnbound(ContextImpl& context, const LangevinIntegratorSDM& integrator) {
  int numLigParticles = LigParticle.size();
  vector<Vec3>& posData = extractPositions(context);
  double displ = 2.38; //DEBUG

  for(int i = 0; i<numLigParticles; i++){
    int p = LigParticle[i];
    posData[p][0] += displ;
    posData[p][1] += displ;
    posData[p][2] += displ;
  }
}


void ReferenceIntegrateLangevinStepSDMKernel::execute(ContextImpl& context, LangevinIntegratorSDM& integrator,
							double BoundEnergy, double UnboundEnergy, double RestraintEnergy){

 double temperature = integrator.getTemperature();
 double friction = integrator.getFriction();
 double stepSize = integrator.getStepSize();
 double lambdac = integrator.getLambda();
 double dlambdac = 0.0;
 double umax = integrator.getUmax();
 double acore = integrator.getAcore();
 double ubcore = integrator.getUbcore();
 double gamma = 0.0;
 double wbcoeff = lambdac;
 double w0coeff = 0.0;
 double lambda1 = lambdac;
 double lambda2 = lambdac;
 double alpha = 1.0;
 double u0 = 0.0;

 //get all slopes(derivatives with respect to lambda for all parameters)
 double m1lambda1 = integrator.getlambda1Slope();
 double m2lambda2 = integrator.getlambda2Slope();
 double mu0 = integrator.getu0Slope();
 double mw0 = integrator.getw0Slope();
 //get all intercepts
 double b1lambda1 = integrator.getlambda1intercept();
 double b2lambda2 = integrator.getlambda2intercept();
 double bu0 = integrator.getu0intercept();
 double bw0 = integrator.getw0intercept();
 if(integrator.getNonEquilibrium() == 1) {
   double Ntmax = integrator.getNoneqtmax();
   lambdac = data.time / Ntmax;
   integrator.setLambda(lambdac);
   //calculate lambda1, lambda2, u0, w0. alpha is constant
   double lambda1 = m1lambda1 * lambdac + b1lambda1;
   double lambda2 = m2lambda2 * lambdac + b2lambda2;
   double u0 = mu0 * lambdac + bu0;
   double w0 = mw0 * lambdac + bw0;
   integrator.setLambda1(lambda1);
   integrator.setLambda2(lambda2);
   integrator.setU0(u0);
   integrator.setW0coeff(w0);
   dlambdac = stepSize / Ntmax;
 }
 
 int method = integrator.getBiasMethod();
 if(method == LangevinIntegratorSDM::QuadraticMethod){
   gamma = integrator.getGamma();
   wbcoeff = integrator.getWBcoeff();
   w0coeff = integrator.getW0coeff();
 } else if( method == LangevinIntegratorSDM::ILogisticMethod){
   lambda1 = integrator.getLambda1();
   lambda2 = integrator.getLambda2();
   alpha = integrator.getAlpha();
   u0 = integrator.getU0();
   w0coeff = integrator.getW0coeff();
 }
 vector<Vec3>& posData = extractPositions(context);
 vector<Vec3>& velData = extractVelocities(context);
 vector<Vec3>& forceData = extractForces(context);
 int numParticles = context.getSystem().getNumParticles();

 //hybrid potential energy
 double fp;
 //double BindE = integrator.SoftCoreF(BoundEnergy - UnboundEnergy, umax, acore, ubcore, fp);
 double BindE = integrator.SoftCoreF(UnboundEnergy - BoundEnergy, umax, acore, ubcore, fp); //DEBUG
 double bfp = 0.0;
 double ebias = 0.0;
 if( method == LangevinIntegratorSDM::QuadraticMethod){
   ebias = 0.5*gamma*BindE*BindE + wbcoeff*BindE + w0coeff;
   bfp = gamma*BindE + wbcoeff;
 }else if( method == LangevinIntegratorSDM::ILogisticMethod){
   double ee = 1.0 + exp(-alpha*(BindE - u0));
   if(alpha > 0){
     ebias = ((lambda2 - lambda1)/alpha) * log(ee);
   }
   ebias += lambda2 * BindE + w0coeff;
   bfp = (lambda2 - lambda1)/ee + lambda1;
 }else{
   ebias = lambdac * BindE;
   bfp = lambdac;
 }
 //double PotEnergy = UnboundEnergy + ebias + RestraintEnergy;
 double PotEnergy = BoundEnergy + ebias + RestraintEnergy;//DEBUG
 integrator.setBindE(BindE);
 integrator.setPotEnergy(PotEnergy);

 if(integrator.getNonEquilibrium() == 1){
   double ee = 1.0 + exp(-alpha*(BindE - u0));
   double dwdl1 = -log(ee)/alpha;
   double dwdl2 = BindE + (log(ee)/alpha);
   double dwdu0 = (lambda2 -lambda1) * exp(-alpha*(BindE - u0))/ee;
   double dwdalpha = 0;
   double dwdw0 = 1;
   
   double dwdlambda = (dwdl1 * m1lambda1) + (dwdl2 * m2lambda2) + (dwdu0 * mu0) + (dwdw0 * mw0);  
   double workValue = integrator.getNoneqWorkvalue();
   //double dworkValue = dlambdac * BindE;
   double dworkValue = dlambdac * dwdlambda;
   integrator.setNoneqWorkvalue(workValue + dworkValue);
 }
 
 //cout << "Bound Energy = " << BoundEnergy << endl;
 //cout << "Unbound Energy = " << UnboundEnergy << endl;
 //cout << "Binding Energy = " << BindE << endl;
 
 // hybrid force, forceData in the r.h.s. at this point holds the restraint forces
 double sp = bfp*fp;
 for (int i = 0 ; i < numParticles; ++i ){
   //forceData[i][0] = sp*BoundForces[i][0]+(1.0-sp)*UnboundForces[i][0]+forceData[i][0] ;
   //forceData[i][1] = sp*BoundForces[i][1]+(1.0-sp)*UnboundForces[i][1]+forceData[i][1] ;
   //forceData[i][2] = sp*BoundForces[i][2]+(1.0-sp)*UnboundForces[i][2]+forceData[i][2] ; 
   //DEBUG
   forceData[i][0] = sp*UnboundForces[i][0]+(1.0-sp)*BoundForces[i][0]+forceData[i][0] ;
   forceData[i][1] = sp*UnboundForces[i][1]+(1.0-sp)*BoundForces[i][1]+forceData[i][1] ;
   forceData[i][2] = sp*UnboundForces[i][2]+(1.0-sp)*BoundForces[i][2]+forceData[i][2] ; 
 }

 if (dynamics == 0 || temperature != prevTemp || friction != prevFriction || stepSize != prevStepSize) {
   // Recreate the computation objects with the new parameters.                                                            

   if (dynamics)
     delete dynamics;
   RealOpenMM tau = static_cast<RealOpenMM>( friction == 0.0 ? 0.0 : 1.0/friction );
   dynamics = new ReferenceStochasticDynamicsSDM(
						   context.getSystem().getNumParticles(),
						   static_cast<RealOpenMM>(stepSize),
						   static_cast<RealOpenMM>(tau),
						   static_cast<RealOpenMM>(temperature) );
   dynamics->setReferenceConstraintAlgorithm(&extractConstraints(context));
   prevTemp = temperature;
   prevFriction = friction;
   prevStepSize = stepSize;
 }

 
 dynamics->update(context.getSystem(), posData, velData, forceData, masses, integrator.getConstraintTolerance());                                                                                                                       

 data.time += stepSize;
 data.stepCount++;
}


double ReferenceIntegrateLangevinStepSDMKernel::computeKineticEnergy(ContextImpl& context, const LangevinIntegratorSDM& integrator) {
  return computeShiftedKineticEnergy(context, masses, 0.0*integrator.getStepSize());//ignore time shift
}

