#ifndef SDM_KERNELS_H_
#define SDM_KERNELS_H_

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

#include "LangevinIntegratorSDM.h"
#include "openmm/KernelImpl.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include <string>

namespace SDMPlugin {
//namespace OpenMM {
/**
 * This kernel is invoked by LangevinIntegratorSDM to take one time step.
 */

  class IntegrateLangevinStepSDMKernel : public OpenMM::KernelImpl {
public:
    static std::string Name() {
        return "IntegrateLangevinStepSDM";
    }
  IntegrateLangevinStepSDMKernel(std::string name, const OpenMM::Platform& platform) : OpenMM::KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param integrator the LangevinIntegrator this kernel will be used for
     */
    virtual void initialize(const OpenMM::System& system, const LangevinIntegratorSDM& integrator) = 0;
    /**
     * Execute the kernel.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the LangevinIntegrator this kernel is being used for
     */
    virtual void execute(OpenMM::ContextImpl& context, LangevinIntegratorSDM& integrator,
			 double State1Energy, double State2Energy, double RestraintEnergy) = 0;

    
    /**
     * Compute the kinetic energy.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the LangevinIntegrator this kernel is being used for
     */
    virtual double computeKineticEnergy(OpenMM::ContextImpl& context, const LangevinIntegratorSDM& integrator) = 0;
    
    /**
     * Save ligand coordinates of the base state of the system
     */
    //virtual void SaveBaseState(ContextImpl& context, const LangevinIntegratorSDM& integrator) = 0;
   
    /**
     * Save ligand coordinates and system forces of bound system
     */
    virtual void SaveState1(OpenMM::ContextImpl& context, const LangevinIntegratorSDM& integrator) = 0;


    /**
     * Save ligand coordinates and system forces of unbound system
     */
    virtual void SaveState2(OpenMM::ContextImpl& context, const LangevinIntegratorSDM& integrator) = 0;
    
    /**
     * Restores ligand coordinates at the base system
     * 
     */
    virtual void RestoreState1(OpenMM::ContextImpl& context, const LangevinIntegratorSDM& integrator) = 0;

    /**
     * Displaces system so that it is in state 1
     * 
     */
    //virtual void MakeState1(OpenMM::ContextImpl& context, const LangevinIntegratorSDM& integrator) = 0;
    
    /**
     * Displaces system so that it is in state 1
     * 
     */
    virtual void MakeState2(OpenMM::ContextImpl& context, const LangevinIntegratorSDM& integrator) = 0;
     
};


  
  
} // namespace SDMPlugin

#endif /*SDM_KERNELS_H_*/
