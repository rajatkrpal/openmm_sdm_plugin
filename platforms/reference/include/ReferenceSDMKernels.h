#ifndef OPENMM_REFERENCESDMKERNELS_H_
#define OPENMM_REFERENCESDMKERNELS_H_

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

#include "SDMKernels.h"
#include "openmm/reference/ReferencePlatform.h"
#include "openmm/kernels.h"
#include "openmm/reference/SimTKOpenMMRealType.h"
#include "openmm/reference/ReferenceNeighborList.h"
#include "ReferenceStochasticDynamicsSDM.h"

using namespace OpenMM;

//class ReferenceStochasticDynamicsSDM;

namespace SDMPlugin {
//namespace OpenMM {

   /**          
   * This kernel is invoked by LangevinIntegratorSDM to take one time step.
   */
class ReferenceIntegrateLangevinStepSDMKernel : public IntegrateLangevinStepSDMKernel {
public:
  ReferenceIntegrateLangevinStepSDMKernel(std::string name, const Platform& platform, ReferencePlatform::PlatformData& data) : IntegrateLangevinStepSDMKernel(name, platform),
  data(data), dynamics(0) {
  }
  ~ReferenceIntegrateLangevinStepSDMKernel();

   
   /**
   * Initialize the kernel, setting up the particle masses.
   * @param system     the System this kernel will be applied to
   * @param integrator the LangevinIntegrator this kernel will be used for
   */
   void initialize(const System& system, const LangevinIntegratorSDM& integrator);
   /**
    * Execute the kernel.
    * @param context    the context in which to execute this kernel       
    * @param integrator the LangevinIntegrator this kernel is being used for
    */
   void execute(ContextImpl& context, LangevinIntegratorSDM& integrator,
		double BoundEnergy, double UnboundEnergy, double RestraintEnergy);
    /**
     * Compute the kinetic energy.                               
     * @param context    the context in which to execute this kernel
     * @param integrator the LangevinIntegrator this kernel is being used for    
     */
   double computeKineticEnergy(ContextImpl& context, const LangevinIntegratorSDM& integrator);

   /* Save ligand coordinates and system forces of bound system */
   void SaveBound(ContextImpl& context, const LangevinIntegratorSDM& integrator);

   /* Save system forces of unbound system */
   void SaveUnbound(ContextImpl& context, const LangevinIntegratorSDM& integrator);
   
   /* Restores ligand coordinates at the bound system */
   void RestoreBound(ContextImpl& context, const LangevinIntegratorSDM& integrator);

   /* Displaces ligand so that it is unbound */
   void MakeUnbound(ContextImpl& context, const LangevinIntegratorSDM& integrator);


 private:
  ReferencePlatform::PlatformData& data;
  ReferenceStochasticDynamicsSDM* dynamics;
  std::vector<RealOpenMM> masses;
  double prevTemp, prevFriction, prevStepSize;
  std::vector<Vec3> restraintForces;

  std::vector<Vec3> BoundForces;
  std::vector<Vec3> UnboundForces;
  double BoundEnergy, UnboundEnergy;
  std::vector<Vec3> LigandCoordinates;

  vector<int> LigParticle;

  };


} // namespace SDMPlugin

#endif /*OPENMM_REFERENCESDMKERNELS_H_*/
