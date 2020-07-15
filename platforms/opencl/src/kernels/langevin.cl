enum {VelScale, ForceScale, NoiseScale, MaxParams};

/**
 * Perform the first step of Langevin integration.
 */

__kernel void integrateLangevinPart1(__global mixed4* restrict velm, __global const real4* restrict force, __global mixed4* restrict posDelta,
        __global const mixed* restrict paramBuffer, __global const mixed2* restrict dt, __global const float4* restrict random, unsigned int randomIndex) {
    mixed vscale = paramBuffer[VelScale];
    mixed fscale = paramBuffer[ForceScale];
    mixed noisescale = paramBuffer[NoiseScale];
    mixed stepSize = dt[0].y;
    int index = get_global_id(0);
    randomIndex += index;
    while (index < NUM_ATOMS) {
        mixed4 velocity = velm[index];
        if (velocity.w != 0.0) {
            mixed sqrtInvMass = sqrt(velocity.w);
            velocity.x = vscale*velocity.x + fscale*velocity.w*force[index].x + noisescale*sqrtInvMass*random[randomIndex].x;
            velocity.y = vscale*velocity.y + fscale*velocity.w*force[index].y + noisescale*sqrtInvMass*random[randomIndex].y;
            velocity.z = vscale*velocity.z + fscale*velocity.w*force[index].z + noisescale*sqrtInvMass*random[randomIndex].z;
            velm[index] = velocity;
            posDelta[index] = stepSize*velocity;
        }
        randomIndex += get_global_size(0);
        index += get_global_size(0);
    }
}

/**
 * Perform the second step of Langevin integration.
 */

__kernel void integrateLangevinPart2(__global real4* restrict posq, __global real4* restrict posqCorrection, __global const mixed4* restrict posDelta, __global mixed4* restrict velm, __global const mixed2* restrict dt) {
#ifdef SUPPORTS_DOUBLE_PRECISION
    double invStepSize = 1.0/dt[0].y;
#else
    float invStepSize = 1.0f/dt[0].y;
    float correction = (1.0f-invStepSize*dt[0].y)/dt[0].y;
#endif
    int index = get_global_id(0);
    while (index < NUM_ATOMS) {
        mixed4 vel = velm[index];
        if (vel.w != 0.0) {
#ifdef USE_MIXED_PRECISION
            real4 pos1 = posq[index];
            real4 pos2 = posqCorrection[index];
            mixed4 pos = (mixed4) (pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z, pos1.w);
#else
            real4 pos = posq[index];
#endif
            mixed4 delta = posDelta[index];
            pos.xyz += delta.xyz;
#ifdef SUPPORTS_DOUBLE_PRECISION
            vel.xyz = convert_mixed4(invStepSize*convert_double4(delta)).xyz;
#else
            vel.xyz = invStepSize*delta.xyz + correction*delta.xyz;
#endif
#ifdef USE_MIXED_PRECISION
            posq[index] = convert_real4(pos);
            posqCorrection[index] = (real4) (pos.x-(real) pos.x, pos.y-(real) pos.y, pos.z-(real) pos.z, 0);
#else
            posq[index] = pos;
#endif
            velm[index] = vel;
        }
        index += get_global_size(0);
    }
}

//lambda-hybridized forces
__kernel void sdmForce(__global const real4* restrict posq,
			 __global real4* restrict force_bound,
			 __global real4* restrict force_unbound,
			 __global real4* restrict force,
			 float lambdac){

    int index = get_global_id(0);
    real lmb  = lambdac;
    real lmb1 = 1.0f - lambdac;
    while (index < NUM_ATOMS) {
      //at this point force[] holds the restraint forces
      //force[index] = lmb*force_bound[index] + lmb1*force_unbound[index] + force[index];
      force[index] = lmb1*force_bound[index] + lmb*force_unbound[index] + force[index];//DEBUG
      index += get_global_size(0);
    }
}


/**
 * restraint force between two groups of atoms (Vsite)
 */
/** disabled
__kernel void CMForceCalcKernel(__global const int*  restrict particle_indexes1,
				__global const int*  restrict particle_indexes2,
				__global const real4* restrict posq, //atomic positions
				__global real4* restrict force,
				int numParticles1, int numParticles2,
				float kf,
				float r0){
				
  uint id = get_global_id(0);
  if(id == 0 && numParticles1 > 0 && numParticles2 > 0){
    real4 rf1 = (real4)(0,0,0,0);
    real4 rf2 = (real4)(0,0,0,0);
    
    real4 cm1 = (real4)(0,0,0,0);
    float vn1 = 1./numParticles1;
    for ( int i = 0; i < numParticles1; i++ ){
      int index = particle_indexes1[i];
      cm1 += vn1*posq[index];
    }
    real4 cm2 = (real4)(0,0,0,0);
    float vn2 = 1./numParticles2;
    for ( int i = 0; i < numParticles2; i++ ){
      int index = particle_indexes2[i];
      cm2 += vn2*posq[index];
    }
    real4 dist = cm2 - cm1;
    float dr = sqrt(dot(dist,dist));
    if(dr>r0) 
    {
      float a = kf*(dr-r0)/dr;
      rf1 = dist *  vn1*a;
      rf2 = dist * -vn2*a;
    }
    
    for ( int i = 0; i < numParticles1; i++ ){
      int index = particle_indexes1[i];
      force[index] += rf1;
    }
    for ( int i = 0; i < numParticles2; i++ ){
      int index = particle_indexes2[i];
      force[index] += rf2;
    }
    
    
  }
}
**/


/**
 * saves bound forces and ligand coordinates
 */
__kernel void SaveBound(int numLigParticles,
			__global   int* restrict LigParticle,
			__global real4* restrict posq,
			__global real4* restrict forces,
			__global real4* restrict BoundForces,
			__global real4* restrict BoundLigandCoordinates) {

  uint index = get_global_id(0);
  while (index < NUM_ATOMS) {
    BoundForces[index] = forces[index];
    index += get_global_size(0);
  }
  
  index = get_global_id(0);
  while(index < numLigParticles){
    int p = LigParticle[index];
    BoundLigandCoordinates[index] = posq[p];
    index += get_global_size(0);
  }

}



/**
 * saves bound forces and ligand coordinates
 */
__kernel void SaveUnbound(__global real4* restrict forces,
			  __global real4* restrict UnboundForces){

  uint index = get_global_id(0);
  while (index < NUM_ATOMS) {
    UnboundForces[index] = forces[index];
    index += get_global_size(0);
  }
}


/**
 * restores bound ligand coordinates
 */
__kernel void RestoreBound(int numLigParticles,
			__global   int* restrict LigParticle,
			__global real4* restrict posq,
			__global real4* restrict BoundLigandCoordinates) {

  uint index = get_global_id(0);
  while (index < numLigParticles) {
    int p = LigParticle[index];
    posq[p] = BoundLigandCoordinates[index];
    index += get_global_size(0);
  }

}

/**
 * displaces ligand to make unbound state
 */
__kernel void MakeUnbound(int numLigParticles,
			  __global   int* restrict LigParticle,
			  float displx, float disply, float displz,
			  __global real4* restrict posq
			  ){
  real4 displacement = (real4)(displx,disply,displz,0);
  uint index = get_global_id(0);
  while (index < numLigParticles) {
    int p = LigParticle[index];
    posq[p] += displacement;
    index += get_global_size(0);
  }

}
