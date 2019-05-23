#ifndef OPENCL_SDM_KERNELS_H_
#define OPENCL_SDM_KERNELS_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
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

#include "SDMKernels.h"
#include "openmm/opencl/OpenCLArray.h"
#include "openmm/opencl/OpenCLContext.h"


using namespace OpenMM;
namespace SDMPlugin {

/**
 * This kernel is invoked by LangevinIntegrator to take one time step.
 */
class OpenCLIntegrateLangevinStepSDMKernel : public IntegrateLangevinStepSDMKernel {
public:
 OpenCLIntegrateLangevinStepSDMKernel(std::string name, const OpenMM::Platform& platform, OpenMM::OpenCLContext& cl) : IntegrateLangevinStepSDMKernel(name, platform), cl(cl), hasInitializedKernels(false),
    prevTemp(-1), prevFriction(-1), prevStepSize(-1) {
    }
    ~OpenCLIntegrateLangevinStepSDMKernel();
    /**
     * Initialize the kernel, setting up the particle masses.
     *
     * @param system     the System this kernel will be applied to
     * @param integrator the LangevinIntegrator this kernel will be used for
     */
    void initialize(const OpenMM::System& system, const LangevinIntegratorSDM& integrator);
    /**
     * Execute the kernel.
     *
     * @param context    the context in which to execute this kernel
     * @param integrator the LangevinIntegrator this kernel is being used for
     */
    void execute(OpenMM::ContextImpl& context, LangevinIntegratorSDM& integrator,
		 double BoundEnergy, double UnboundEnergy, double RestraintEnergy);
    /**
     * Compute the kinetic energy.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the LangevinIntegrator this kernel is being used for
     */
    double computeKineticEnergy(OpenMM::ContextImpl& context, const LangevinIntegratorSDM& integrator);

   /* Save ligand coordinates and system forces of bound system */
   void SaveBound(ContextImpl& context, const LangevinIntegratorSDM& integrator);

   /* Save ligand coordinates and system forces of bound system */
   void SaveUnbound(ContextImpl& context, const LangevinIntegratorSDM& integrator);

   /* Restores ligand coordinates at the bound system */
   void RestoreBound(ContextImpl& context, const LangevinIntegratorSDM& integrator);

   /* Displaces ligand so that it is unbound */
   void MakeUnbound(ContextImpl& context, const LangevinIntegratorSDM& integrator);

private:
    OpenCLContext& cl;
    double prevTemp, prevFriction, prevStepSize;
    bool hasInitializedKernels;
    OpenCLArray* params;
    int ligsize;
    OpenCLArray* LigParticle;
    OpenCLArray* BoundForces;
    OpenCLArray* UnboundForces;
    OpenCLArray* BoundLigandCoordinates;
    
    cl::Kernel kernel1, kernel2;
    cl::Kernel kernelSDMForce; //bedam hybrid forces
    cl::Kernel SaveBoundKernel;    //save forces and ligand coordinates of bound state
    cl::Kernel SaveUnboundKernel;    //save forces and ligand coordinates of bound state
    cl::Kernel RestoreBoundKernel; //restores bound coordinates of the ligand
    cl::Kernel MakeUnboundKernel;  //creates unbound state by displacing the ligand

    double umax;
    double acore;
    double ubcore;
};

} // namespace SDMPlugin

#endif /*OPENCL_SDM_KERNELS_H_*/
