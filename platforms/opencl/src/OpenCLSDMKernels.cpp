/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2015 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 * -------------------------------------------------------------------------- */

#include "openmm/opencl/OpenCLKernels.h"
#include "openmm/opencl/OpenCLForceInfo.h"
#include "openmm/opencl/OpenCLArray.h"
#include "OpenCLSDMKernels.h"
#include "OpenCLSDMKernelSources.h"
#include "openmm/opencl/OpenCLForceInfo.h"
#include "openmm/LangevinIntegrator.h"
#include "openmm/Context.h"
#include "openmm/internal/OSRngSeed.h"
#include "openmm/opencl/OpenCLBondedUtilities.h"
#include "openmm/opencl/OpenCLExpressionUtilities.h"
#include "openmm/opencl/OpenCLIntegrationUtilities.h"
#include "openmm/opencl/OpenCLNonbondedUtilities.h"
#include <algorithm>
#include <cmath>
#include <set>
#include <iostream>

using namespace SDMPlugin;
using namespace OpenMM;
using namespace std;

static void setPosqCorrectionArg(OpenCLContext& cl, cl::Kernel& kernel, int index) {
  if (cl.getUseMixedPrecision())
    kernel.setArg<cl::Buffer>(index, cl.getPosqCorrection().getDeviceBuffer());
  else
    kernel.setArg<void*>(index, NULL);
}


const float BOLTZMANN = 1.380658e-23f; // (J/K)
const float AVOGADRO = 6.0221367e23f;
const float RGAS = BOLTZMANN*AVOGADRO; // (J/(mol K))
const float BOLTZ = RGAS/1000;         // (kJ/(mol K))






OpenCLIntegrateLangevinStepSDMKernel::~OpenCLIntegrateLangevinStepSDMKernel() {
    if (params != NULL)
        delete params;
    
}

void OpenCLIntegrateLangevinStepSDMKernel::initialize(const System& system, const LangevinIntegratorSDM& integrator) {
    cl.getPlatformData().initializeContexts(system);
    cl.getIntegrationUtilities().initRandomNumberGenerator(integrator.getRandomNumberSeed());
    map<string, string> defines;
    defines["NUM_ATOMS"] = cl.intToString(cl.getNumAtoms());
    defines["PADDED_NUM_ATOMS"] = cl.intToString(cl.getPaddedNumAtoms());
    cl::Program program = cl.createProgram(OpenCLSDMKernelSources::langevin, defines, "");
    kernelSDMForce = cl::Kernel(program, "sdmForce");
    kernel1 = cl::Kernel(program, "integrateLangevinPart1");
    kernel2 = cl::Kernel(program, "integrateLangevinPart2");
    params = new OpenCLArray(cl, 3, cl.getUseDoublePrecision() || cl.getUseMixedPrecision() ? sizeof(cl_double) : sizeof(cl_float), "langevinParams");

    vector<int> lig_particles = integrator.getLigParticle();
    ligsize = lig_particles.size();
    LigParticle = new OpenCLArray(cl, ligsize, sizeof(cl_int), "LigParticle");
    LigParticle->upload(lig_particles);

    //need to worry about mixed/double precisions?
    BoundForces = OpenCLArray::create<mm_float4>(cl, cl.getPaddedNumAtoms(), "BoundForces");
    UnboundForces = OpenCLArray::create<mm_float4>(cl, cl.getPaddedNumAtoms(), "UnboundForces");
    BoundLigandCoordinates = OpenCLArray::create<mm_float4>(cl, ligsize, "BoundLigandCoordinates");

    SaveBoundKernel = cl::Kernel(program, "SaveBound");
    SaveBoundKernel.setArg<cl_int>(0, ligsize);
    SaveBoundKernel.setArg<cl::Buffer>(1, LigParticle->getDeviceBuffer());    
    SaveBoundKernel.setArg<cl::Buffer>(2, cl.getPosq().getDeviceBuffer());
    SaveBoundKernel.setArg<cl::Buffer>(3, cl.getForce().getDeviceBuffer());
    SaveBoundKernel.setArg<cl::Buffer>(4, BoundForces->getDeviceBuffer());
    SaveBoundKernel.setArg<cl::Buffer>(5, BoundLigandCoordinates->getDeviceBuffer());

    SaveUnboundKernel = cl::Kernel(program, "SaveUnbound");
    SaveUnboundKernel.setArg<cl::Buffer>(0, cl.getForce().getDeviceBuffer());
    SaveUnboundKernel.setArg<cl::Buffer>(1, UnboundForces->getDeviceBuffer());

    RestoreBoundKernel = cl::Kernel(program, "RestoreBound");
    RestoreBoundKernel.setArg<cl_int>(0, ligsize);
    RestoreBoundKernel.setArg<cl::Buffer>(1, LigParticle->getDeviceBuffer());    
    RestoreBoundKernel.setArg<cl::Buffer>(2, cl.getPosq().getDeviceBuffer());
    RestoreBoundKernel.setArg<cl::Buffer>(3, BoundLigandCoordinates->getDeviceBuffer());
    
    MakeUnboundKernel = cl::Kernel(program, "MakeUnbound");
    Vec3 displ = integrator.getDisplacement();
    float dx = displ[0]; float dy = displ[1]; float dz = displ[2];
    MakeUnboundKernel.setArg<cl_int>(0, ligsize);
    MakeUnboundKernel.setArg<cl::Buffer>(1, LigParticle->getDeviceBuffer());    
    MakeUnboundKernel.setArg<cl_float>  (2, dx);
    MakeUnboundKernel.setArg<cl_float>  (3, dy);
    MakeUnboundKernel.setArg<cl_float>  (4, dz);
    MakeUnboundKernel.setArg<cl::Buffer>(5, cl.getPosq().getDeviceBuffer());

    //soft core parameters
    umax = integrator.getUmax();
    acore = integrator.getAcore();
    ubcore = integrator.getUbcore();
}


void OpenCLIntegrateLangevinStepSDMKernel::SaveBound(ContextImpl& context, const LangevinIntegratorSDM& integrator) {
  cl.executeKernel(SaveBoundKernel, cl.getNumAtoms());
}

void OpenCLIntegrateLangevinStepSDMKernel::SaveUnbound(ContextImpl& context, const LangevinIntegratorSDM& integrator) {
  cl.executeKernel(SaveUnboundKernel, cl.getNumAtoms());
}



void OpenCLIntegrateLangevinStepSDMKernel::RestoreBound(ContextImpl& context, const LangevinIntegratorSDM& integrator) {
  cl.executeKernel(RestoreBoundKernel, cl.getNumAtoms());
}

void OpenCLIntegrateLangevinStepSDMKernel::MakeUnbound(ContextImpl& context, const LangevinIntegratorSDM& integrator) {
  cl.executeKernel(MakeUnboundKernel, cl.getNumAtoms());
}


void OpenCLIntegrateLangevinStepSDMKernel::execute(ContextImpl& context, LangevinIntegratorSDM& integrator,
						     double BoundEnergy, double UnboundEnergy, double RestraintEnergy) {
    OpenCLIntegrationUtilities& integration = cl.getIntegrationUtilities();
    int numAtoms = cl.getNumAtoms();
    double lambdac = integrator.getLambda();
    double dlambdac = 0.0;
    double gamma = 0.0;
    double wbcoeff = lambdac;
    double w0coeff = 0.0;
    double lambda1 = lambdac;
    double lambda2 = lambdac;
    double u0 = 0.0;
    double alpha = 1.0;
    double stepSize = integrator.getStepSize();
    
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
    if(integrator.getNonEquilibrium() == 1){
      double Ntmax = integrator.getNoneqtmax();
      lambdac = cl.getTime() / Ntmax;
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
    
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;
	kernel1.setArg<cl::Buffer>(0, cl.getVelm().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(1, cl.getForce().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(2, integration.getPosDelta().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(3, params->getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(4, integration.getStepSize().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(5, integration.getRandom().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(0, cl.getPosq().getDeviceBuffer());
        setPosqCorrectionArg(cl, kernel2, 1);
        kernel2.setArg<cl::Buffer>(2, integration.getPosDelta().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(3, cl.getVelm().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(4, integration.getStepSize().getDeviceBuffer());

	kernelSDMForce.setArg<cl::Buffer>(0, cl.getPosq().getDeviceBuffer());
	kernelSDMForce.setArg<cl::Buffer>(1, BoundForces->getDeviceBuffer());
	kernelSDMForce.setArg<cl::Buffer>(2, UnboundForces->getDeviceBuffer());
	kernelSDMForce.setArg<cl::Buffer>(3, cl.getForce().getDeviceBuffer());
    }

    /*DEBUG
    vector<mm_float4> unbound_forces;
    cl.getForce().download(unbound_forces);    
    vector<mm_float4> bound_forces;
    BoundForces->download(bound_forces);
    */
    //hybrid potential energy
    double fp;
    double BindE = integrator.SoftCoreF(BoundEnergy - UnboundEnergy, umax, acore, ubcore, fp);
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
    double PotEnergy = UnboundEnergy + ebias + RestraintEnergy;
    integrator.setPotEnergy(PotEnergy);
    integrator.setBindE(BindE);

    if(integrator.getNonEquilibrium() == 1){
      double ee = 1.0 + exp(-alpha*(BindE - u0));
      double dwdl1 = -log(ee)/alpha;
      double dwdl2 = BindE + (log(ee)/alpha);
      double dwdu0 = (lambda2 -lambda1) * exp(-alpha*(BindE - u0))/ee;
      double dwdalpha = 0;
      double dwdw0 = 1;
   
      double dwdlambda = (dwdl1 * m1lambda1) + (dwdl2 * m2lambda2) + (dwdu0 * mu0) + (dwdw0 * mw0);
      double workValue = integrator.getNoneqWorkvalue();
      double dworkValue = dlambdac * dwdlambda;
      //double dworkValue = dlambdac * BindE;
      integrator.setNoneqWorkvalue(workValue + dworkValue);
    }
    
    //sdm hybrid force
    float sp = bfp*fp;
    kernelSDMForce.setArg<cl_float>(4,sp);
    cl.executeKernel(kernelSDMForce, numAtoms);

    /*DEBUG
    vector<mm_float4> velocities;
    cl.getVelm().download(velocities);
    double max_magnitude = 0;
    int imax = -1;
    for(int i=0; i<velocities.size();i++){
      double magnitude = pow(velocities[i].x,2) + pow(velocities[i].y,2) +  pow(velocities[i].z,2);
      if(magnitude > max_magnitude){
	imax = i;
	max_magnitude = magnitude;
      }
    }
    cout << "VMAX " << imax << " " << max_magnitude << endl;
    */    

    /*DEBUG
    vector<mm_float4> hybrid_forces;
    cl.getForce().download(hybrid_forces);

    vector<mm_float4> velocities;
    cl.getVelm().download(velocities);
    
    double max_magnitude = 0;
    int imax = -1;
    for(int i=0; i<hybrid_forces.size();i++){
      double magnitude = pow(hybrid_forces[i].x,2) + pow(hybrid_forces[i].y,2) +  pow(hybrid_forces[i].z,2);
      if(magnitude > max_magnitude){
	imax = i;
	max_magnitude = magnitude;
      }
      //      cout << "F " << i << " " << bound_forces[i].z << " " << unbound_forces[i].z << " " << hybrid_forces[i].z << " " << velocities[i].z << endl;
      cout << "F " << i << " " << hybrid_forces[i].z << " " << velocities[i].z << endl;
    }
    cout << "FMAX " << imax << " " << hybrid_forces[imax].x << " " << hybrid_forces[imax].y << " " << hybrid_forces[imax].z << endl;
    cout << "lambda= " << fp*lambdac << endl;
    */
    
    double temperature = integrator.getTemperature();
    double friction = integrator.getFriction();
    
    cl.getIntegrationUtilities().setNextStepSize(stepSize);

    /* DEBUG 
    if(cl.getUseDoublePrecision()){
      cout << "Using double precision" << endl;
    }
    if(cl.getUseMixedPrecision()){
      cout << "Using mixed precision" << endl;
    }
    if (!(cl.getUseDoublePrecision() || cl.getUseMixedPrecision())) {
      cout << "Using single precision" << endl;
    }
    */
    
    if (temperature != prevTemp || friction != prevFriction || stepSize != prevStepSize) {
        // Calculate the integration parameters.
        double kT = BOLTZ*temperature;
        double vscale = exp(-stepSize*friction);
        double fscale = (friction == 0 ? stepSize : (1-vscale)/friction);
        double noisescale = sqrt(kT*(1-vscale*vscale));
        if (cl.getUseDoublePrecision() || cl.getUseMixedPrecision()) {
            vector<cl_double> p(params->getSize());
            p[0] = vscale;
            p[1] = fscale;
            p[2] = noisescale;
            params->upload(p);
        }
        else {
            vector<cl_float> p(params->getSize());
            p[0] = (cl_float) vscale;
            p[1] = (cl_float) fscale;
            p[2] = (cl_float) noisescale;
            params->upload(p);
        }
        prevTemp = temperature;
        prevFriction = friction;
        prevStepSize = stepSize;

    }
    
    // Call the first integration kernel.
    kernel1.setArg<cl_uint>(6, integration.prepareRandomNumbers(cl.getPaddedNumAtoms()));
    cl.executeKernel(kernel1, numAtoms);
    
    /*DEBUG
    vector<mm_float4> velocitiesf1;
    cl.getVelm().download(velocitiesf1);
    cout << "VFN=" << velocitiesf1.size() << endl;
    for(int i=0; i<velocitiesf1.size();i++){
      cout << "VF1 " << i << " " << velocitiesf1[i].z << endl;      
    }
    */
    
    // Apply constraints.
    integration.applyConstraints(integrator.getConstraintTolerance());

    // Call the second integration kernel.
    cl.executeKernel(kernel2, numAtoms); 

    //virtual sites
    integration.computeVirtualSites();

    /* DEBUG
    vector<mm_float4> velocitiesf2;
    cl.getVelm().download(velocitiesf2);
    cout << "VFN=" << velocitiesf2.size() << endl;
    for(int i=0; i<velocitiesf2.size();i++){
      cout << "VF2 " << i << " " << velocitiesf2[i].z << endl;      
    }
    */


    
    // Update the time and step count.
    cl.setTime(cl.getTime()+stepSize);
    cl.setStepCount(cl.getStepCount()+1);
    cl.reorderAtoms();
    
    // Reduce UI lag.
#ifdef WIN32
    cl.getQueue().flush();
#endif
}

double OpenCLIntegrateLangevinStepSDMKernel::computeKineticEnergy(ContextImpl& context, const LangevinIntegratorSDM& integrator) {
  return cl.getIntegrationUtilities().computeKineticEnergy(0.0*integrator.getStepSize());//ignore time shift
}
