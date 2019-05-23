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

#include <exception>

#include "CudaSDMKernelFactory.h"
#include "CudaSDMKernels.h"
//#include "openmm/cuda/CudaParallelKernels.h"
//#include "openmm/cuda/CudaPlatform.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/windowsExport.h"

using namespace SDMPlugin;
using namespace OpenMM;

extern "C" OPENMM_EXPORT void registerPlatforms() {
}

extern "C" OPENMM_EXPORT void registerKernelFactories() {
    try {
        Platform& platform = Platform::getPlatformByName("CUDA");
        CudaSDMKernelFactory* factory = new CudaSDMKernelFactory();
        platform.registerKernelFactory(IntegrateLangevinStepSDMKernel::Name(), factory);
    }
    catch (std::exception ex) {
        // Ignore
    }
}

extern "C" OPENMM_EXPORT void registerSDMCudaKernelFactories() {
    try {
        Platform::getPlatformByName("CUDA");
    }
    catch (...) {
        Platform::registerPlatform(new CudaPlatform());
    }
    registerKernelFactories();
}



KernelImpl* CudaSDMKernelFactory::createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const {
    CudaPlatform::PlatformData& data = *static_cast<CudaPlatform::PlatformData*>(context.getPlatformData());
    CudaContext& cu = *static_cast<CudaPlatform::PlatformData*>(context.getPlatformData())->contexts[0];

    //CudaContext& cu = *data.contexts[0];
    
    if (name == IntegrateLangevinStepSDMKernel::Name())
        return new CudaIntegrateLangevinStepSDMKernel(name, platform, cu);
    
    throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '")+name+"'").c_str());
}
