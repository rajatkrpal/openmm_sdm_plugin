/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2014 Stanford University and the Authors.      *
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

#include "openmm/cuda/CudaKernels.h"
#include "openmm/cuda/CudaArray.h"
#include "CudaSDMKernels.h"
#include "CudaSDMKernelSources.h"
#include "openmm/cuda/CudaForceInfo.h"
#include "openmm/LangevinIntegrator.h"
#include "openmm/internal/OSRngSeed.h"
#include "openmm/Context.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/cuda/CudaBondedUtilities.h"
#include "openmm/cuda/CudaExpressionUtilities.h"
#include "openmm/cuda/CudaIntegrationUtilities.h"
#include "openmm/cuda/CudaNonbondedUtilities.h"
#include <algorithm>
#include <cmath>
#include <set>

using namespace SDMPlugin;
using namespace OpenMM;
using namespace std;


const float BOLTZMANN = 1.380658e-23f; // (J/K)
const float AVOGADRO = 6.0221367e23f;
const float RGAS = BOLTZMANN*AVOGADRO; // (J/(mol K))
const float BOLTZ = RGAS/1000;         // (kJ/(mol K))


CudaIntegrateLangevinStepSDMKernel::~CudaIntegrateLangevinStepSDMKernel() {
    cu.setAsCurrent();
    if (params != NULL)
        delete params;
}

void CudaIntegrateLangevinStepSDMKernel::initialize(const System& system, const LangevinIntegratorSDM& integrator) {
    cu.getPlatformData().initializeContexts(system);
    cu.setAsCurrent();
    cu.getIntegrationUtilities().initRandomNumberGenerator(integrator.getRandomNumberSeed());
    map<string, string> defines;
    CUmodule module = cu.createModule(CudaSDMKernelSources::langevin, defines, "");
    //edit on 7/7/15    
    kernel1 = cu.getKernel(module, "bedamForce");
    kernel2 = cu.getKernel(module, "integrateLangevinPart1");
    kernel3 = cu.getKernel(module, "integrateLangevinPart2");
    kernel4 = cu.getKernel(module, "copyDataToSecondPart");

    params = new CudaArray(cu, 3, cu.getUseDoublePrecision() || cu.getUseMixedPrecision() ? sizeof(double) : sizeof(float), "langevinParams");
    prevStepSize = -1.0;
}

void CudaIntegrateLangevinStepSDMKernel::execute(ContextImpl& context, const LangevinIntegratorSDM& integrator) {
    cu.setAsCurrent();
    CudaIntegrationUtilities& integration = cu.getIntegrationUtilities();
    int numAtoms = cu.getNumAtoms();

    int ligId = integrator.getLigandId();
    double lambdaId = integrator.getLamdaId();
    int atom1 = integrator.getAtom1Number();
    int atom2 = integrator.getAtom2Number()+ligId;
    double kf = integrator.getKf();
    double r0 = integrator.getR0();

    int paddedNumAtoms = cu.getPaddedNumAtoms();
    double temperature = integrator.getTemperature();
    double friction = integrator.getFriction();
    double stepSize = integrator.getStepSize();
    if (temperature != prevTemp || friction != prevFriction || stepSize != prevStepSize) {
        // Calculate the integration parameters.

        double tau = (friction == 0.0 ? 0.0 : 1.0/friction);
        double kT = BOLTZ*temperature;
        double vscale = exp(-stepSize/tau);
        double fscale = (1-vscale)*tau;
        double noisescale = sqrt(2*kT/tau)*sqrt(0.5*(1-vscale*vscale)*tau);
        if (cu.getUseDoublePrecision() || cu.getUseMixedPrecision()) {
            vector<double> p(params->getSize());
            p[0] = vscale;
            p[1] = fscale;
            p[2] = noisescale;
            params->upload(p);
            double2 ss = make_double2(0, stepSize);
            integration.getStepSize().upload(&ss);
        }
        else {
            vector<float> p(params->getSize());
            p[0] = (float) vscale;
            p[1] = (float) fscale;
            p[2] = (float) noisescale;
            params->upload(p);
            float2 ss = make_float2(0, (float) stepSize);
            integration.getStepSize().upload(&ss);
        }
        prevTemp = temperature;
        prevFriction = friction;
        prevStepSize = stepSize;
    }

    //call the two atoms restraint force kernel1
   // void* args1[] = {&numAtoms, &cu.getPosq().getDevicePointer(), &cu.getForce().getDevicePointer(),(void*) &atom1,(void*) &atom2, (void*) &kf, (void*) r0, (void*) &lambdaId)};
    void* args1[] = {&numAtoms, &paddedNumAtoms,&cu.getPosq().getDevicePointer(), &cu.getForce().getDevicePointer(),(void*) &atom1,(void*) &atom2,(void*) &kf, (void*) &r0,(void*) &lambdaId};
    cu.executeKernel(kernel1,args1,numAtoms);


    // Call the first integration kernel.
    
    int randomIndex = integration.prepareRandomNumbers(cu.getPaddedNumAtoms());
    void* args2[] = {&numAtoms, &paddedNumAtoms, &cu.getVelm().getDevicePointer(), &cu.getForce().getDevicePointer(), &integration.getPosDelta().getDevicePointer(),
            &params->getDevicePointer(), &integration.getStepSize().getDevicePointer(), &integration.getRandom().getDevicePointer(), &randomIndex};
    cu.executeKernel(kernel2, args2, numAtoms);

    // Apply constraints.

    integration.applyConstraints(integrator.getConstraintTolerance());

    // Call the second integration kernel.

    CUdeviceptr posCorrection = (cu.getUseMixedPrecision() ? cu.getPosqCorrection().getDevicePointer() : 0);
    void* args3[] = {&numAtoms, &cu.getPosq().getDevicePointer(), &posCorrection, &integration.getPosDelta().getDevicePointer(),
            &cu.getVelm().getDevicePointer(), &integration.getStepSize().getDevicePointer()};
    cu.executeKernel(kernel3, args3, numAtoms);
    integration.computeVirtualSites();

    //call the third coordinate manipulator kernel
    void* args4[] = {&numAtoms, &cu.getPosq().getDevicePointer(),&cu.getVelm().getDevicePointer(),(void*) &ligId};
    cu.executeKernel(kernel4,args4, numAtoms);

    // Update the time and step count.

    cu.setTime(cu.getTime()+stepSize);
    cu.setStepCount(cu.getStepCount()+1);
    cu.reorderAtoms();
}

double CudaIntegrateLangevinStepSDMKernel::computeKineticEnergy(ContextImpl& context, const LangevinIntegratorSDM& integrator) {
  return cu.getIntegrationUtilities().computeKineticEnergy(0.0*integrator.getStepSize());//ignore time shift
}
