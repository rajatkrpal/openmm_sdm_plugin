#ifndef OPENMM_LANGEVININTEGRATORSDM_H_
#define OPENMM_LANGEVININTEGRATORSDM_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM
 * Langevin Integrator for SDM alchemical calculations                      *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2012 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 * Contributors: Baofeng Zhang, Emilio Gallicchio                             *
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
#include <vector>
#include <iostream>
#include <cmath>
#include "openmm/Integrator.h"
#include "openmm/Kernel.h"
#include "openmm/internal/windowsExport.h"
#include "openmm/Vec3.h"

using namespace OpenMM;
using namespace std;

namespace SDMPlugin {

//namespace OpenMM {
/**
 * This is an Integrator which simulates a System using Langevin dynamics.
 */

class OPENMM_EXPORT LangevinIntegratorSDM : public OpenMM::Integrator { 


public:
    /**
     * Create a LangevinIntegrator.
     * 
     * @param temperature    the temperature of the heat bath (in Kelvin)
     * @param frictionCoeff  the friction coefficient which couples the system to the heat bath (in inverse picoseconds)
     * @param stepSize       the step size with which to integrator the system (in picoseconds)
     *
     * @param LigParticle    particles belonging to the ligand
     *
     * The SDM integrator assumes that the ligand is appropriately restrained by a combination 
     * of restraining forces. In particular, it assumes that a CustomCentroidBondForce
     * restrains the position of the ligand. And that, optionally, a CustomCompoundBondForce 
     * restrains the orientation of the ligand.
     */    
  LangevinIntegratorSDM(double temperature, double frictionCoeff, double stepSize,
			  int nParticles) ;

  /**
   * Set the value the lambda alchemical parameter. 
   * This is used in the linear biasing function (below) and also as a general
   * progress parameter usually ranging from 0 (unbound state) to 1 (bound state) 
   *
   * @param lambdat:     the value of lambda
   */
  void setLambda(double lambdat){
    lambdac = lambdat;
  }
  /**
   * get the value the lambda alchemical parameter
   *
   * @return the value of lambda
   */
  double getLambda() const {
    return lambdac;
  }

  /**
   * The different methods available for the alchemical biasing function w(u), where u(x) is
   * binding energy u(x) = Ubound(x) - Uunbound(x). 
   * The alchemical potential is then
   *    U(x) = Uunbound(x) + w(u(x))
   */
  /* I don't know how to do enums in swig ...
  enum BiasMethod {
    //
    // w(u) = Lambda u + W0
    //
    LinearMethod = 0,

    //
    // w(u) = 1/2 Gamma u^2 + B u + W0
    //
    QuadraticMethod = 1,

    //
    // w(u) = (Lambda2 - Lambda1) ln[1 + exp(-Alpha (u-U0))]/Alpha + Lambda2 u + W0
    //
    // this is the integrated form of the logistic function:
    // wp(u) = (Lambda2 - Lambda1)/(1 + exp(-Alpha (u-U0))) + Lambda1
    //
    ILogisticMethod = 2
  };
  */
  static const int LinearMethod = 0;
  static const int QuadraticMethod = 1;
  static const int ILogisticMethod = 2;
  
  /**
   * Set the biasing method
   */
  void setBiasMethod(int method){
    bias_method = method;
  }
  /**
   * Get the biasing method being used
   */
  int getBiasMethod() const {
    return bias_method;
  }

  /**
   * Set the soft-core method
   * 0 = no soft core: usc = u
   * 1 = tanh soft core: usc = um tanh(u/um)
   * 2 = rational function soft core: usc = um (y^a-1)/(y^a + 1), y = 1 + 2 x + 2 x^2, x = u/(a um) 
   */
  static const int NoSoftCoreMethod = 0;
  static const int TanhMethod = 1;
  static const int RationalMethod = 2;
  void setSoftCoreMethod(int method){
    softcore_method = method;
  }
  /**
   * Get the soft-core method being used
   */
  int getSoftCoreMethod() const {
    return softcore_method;
  }

  /**
   * set the value the gamma parameter of the quadratic bias potential
   *
   * @param gammat:     the value of gamma in 1/kjmol
   */
  void setGamma(double gammat){
    gammac = gammat;
  }
  /**
   * get the value the gamma alchemical parameter
   *
   * @return the value of gamma
   */
  double getGamma() const {
    return gammac;
  }

  /**
   * set the value the linear coefficient of the quadratic bias potential
   *
   * @param wbcoeff_t:   the value of the linear coefficient
   */
  void setWBcoeff(double wbcoeff_t){
    wbcoeff = wbcoeff_t;
  }
  /**
   * get the value the linear coefficient of the quadratic bias potential
   *
   * @return the value of linear coefficient
   */
  double getWBcoeff() const {
    return wbcoeff;
  }


  /**
   * set the value the constant term of the quadratic bias potential
   *
   * @param w0coeff_t:     the value of the constant term in kjmol
   */
  void setW0coeff(double w0coeff_t){
    w0coeff = w0coeff_t;
  }
  /**
   * get the value the constant coefficient of the quadratic bias potential
   *
   * @return the value of the constant coefficient
   */
  double getW0coeff() const {
    return w0coeff;
  }

  /**
   * get/set methods specific for the logistic bias potential. 
   * the other parameters are defined above
   */
  void setLambda1(double lambda1_t) {
    lambda1 = lambda1_t;
  }
  double getLambda1() const {
    return lambda1;
  }
  void setLambda2(double lambda2_t) {
    lambda2 = lambda2_t;
  }
  double getLambda2() const {
    return lambda2;
  }
  void setAlpha(double alpha_t) { //in 1/kjmol
    alpha = alpha_t;
  }
  double getAlpha() const {
    return alpha;
  }
  void setU0(double u0_t) { //in kjmol
    u0 = u0_t;
  }
  double getU0() const {
    return u0;
  }
  
  /**
   * Get the temperature of the heat bath (in Kelvin).
   *
   * @return the temperature of the heat bath, measured in Kelvin
   */
  double getTemperature() const {
    return temperature;
  }
  /**
   * Set the temperature of the heat bath (in Kelvin).
   *
   * @param temp    the temperature of the heat bath, measured in Kelvin
   */
  void setTemperature(double temp) {
    temperature = temp;
  }

  /**
   * Get the friction coefficient which determines how strongly the system is coupled to
   * the heat bath (in inverse ps).
   *
   * @return the friction coefficient, measured in 1/ps
   */
  double getFriction() const {
    return friction;
  }
  /**
   * Set the friction coefficient which determines how strongly the system is coupled to
   * the heat bath (in inverse ps).
   *
   * @param coeff    the friction coefficient, measured in 1/ps
   */
  void setFriction(double coeff) {
    friction = coeff;
  }
  /**
   * Get the random number seed.  See setRandomNumberSeed() for details.
   */
  int getRandomNumberSeed() const {
    return randomNumberSeed;
  }
  /**
   * Set the random number seed.  The precise meaning of this parameter is undefined, and is left up
   * to each Platform to interpret in an appropriate way.  It is guaranteed that if two simulations
   * are run with different random number seeds, the sequence of random forces will be different.  On
   * the other hand, no guarantees are made about the behavior of simulations that use the same seed.
   * In particular, Platforms are permitted to use non-deterministic algorithms which produce different
   * results on successive runs, even if those runs were initialized identically.
   */
  void setRandomNumberSeed(int seed) {
    randomNumberSeed = seed;
  }
    
  /**
   * Advance a simulation through time by taking a series of time steps.
   * 
   * @param steps   the number of time steps to take
   */
  void step(int steps);


  /**
   * Acquire current binding energy
   */
  double getBindE(void) const {
    //cout << "get BindE API " << BindE << endl;//DEBUG
    return BindE;
  }
  /**
   * Set binding energy
   * (this would be called by a computational kernel)
   */
  void setBindE(double be)  {
    //cout << "set BindE API " << be << endl; //DEBUG
    BindE = be; 
  }

  /**
   * Acquire current potential energy
   */
  double getPotEnergy(void) const {
    return PotEnergy;
  }
  /**
   * Set potential energy
   * (this would be called by a computational kernel)
   */
  void setPotEnergy(double e)  {
    PotEnergy = e;
  }

  /**
     The functional form of the binding energy soft core potential is:
     
     up = umax tanh(g(u/umax)) if u  > 0
     up = u                    if u <= 0
     x = u/umax
     g(x) = x/(1 + a (x+b)^(1 - 1/n))
     up never exceeds umax.
     
     defaults:
     umax = 200.0 kJ/mol
     a = 3
     b = 0.1
     n = 12
  */
    
  /**
   * Acquire/set the "umax" soft core parameter
   */
  double getUmax(void) const {
    return Umax;
  }
  void setUmax(double um)  {
    Umax = um;
  }
  /**
   * Acquire/set the "a" soft core parameter
   */
  double getAcore(void) const {
    return Acore;
  }
  void setAcore(double a)  {
    Acore = a;
  }

  /**
   * Acquire/set the "Ub" soft core parameter
   */
  double getUbcore(void) const {
    return Ubcore;
  }
  void setUbcore(double ub)  {
    Ubcore = ub;
  }

  /*
   * soft-core function, u is the argument and fp is the derivative
   */
  double SoftCoreF(double u, double umax, double a, double ub, double& fp);
  
  //set the non-equilibrium max time
  void setNoneqtmax(double noneq_tmax) {
    Ntmax = noneq_tmax;
  }
  double getNoneqtmax() const {
    return Ntmax;
  }

  void setNonEquilibrium( int noneq){
    nonequilibrium = noneq;
  }
  
  int getNonEquilibrium() const {
    return nonequilibrium;
  }

  void setNoneqWorkvalue(double noneq_work){
    work_value = noneq_work;
  }

  double getNoneqWorkvalue() const {
    return work_value;
  }


  void setlambda1Slope(double ml1){
    m1Lambda1 = ml1;
  }
  
  double getlambda1Slope() const {
    return m1Lambda1;
  }

  void setlambda2Slope(double ml2){
    m2Lambda2 = ml2;
  }
  
  double getlambda2Slope() const {
    return m2Lambda2;
  }

  void setu0Slope(double mu0){
    m_u0 = mu0;
  }
  
  double getu0Slope() const {
    return m_u0;
  }

  void setw0Slope(double mw0){
    m_w0 = mw0;
  }
  
  double getw0Slope() const {
    return m_w0;
  }

  void setlambda1intercept(double bl1){
    b1Lambda1 = bl1;
  }
  
  double getlambda1intercept() const {
    return b1Lambda1;
  }

  void setlambda2intercept(double bl2){
    b2Lambda2 = bl2;
  }
  
  double getlambda2intercept() const {
    return b2Lambda2;
  }
  void setu0intercept(double bu0){
    b_u0 = bu0;
  }
  
  double getu0intercept() const {
    return b_u0;
  }

  void setw0intercept(double bw0){
    b_w0 = bw0;
  }
  
  double getw0intercept() const {
    return b_w0;
  }


  void setDisplacement(int atom, double dx, double dy, double dz){
    displ[atom] = Vec3(dx,dy,dz);
  }
  Vec3 getDisplacement(int atom) const {
    return displ[atom]; 
  }

protected:
  /**
   * This will be called by the Context when it is created.  It informs the Integrator
   * of what context it will be integrating, and gives it a chance to do any necessary initialization.
   * It will also get called again if the application calls reinitialize() on the Context.
   */
  void initialize(ContextImpl& context);
  /**
   * This will be called by the Context when it is destroyed to let the Integrator do any necessary
   * cleanup.  It will also get called again if the application calls reinitialize() on the Context.
   */
  void cleanup();
  /**
   * Get the names of all Kernels used by this Integrator.
   */
  std::vector<std::string> getKernelNames();
  
  /**
   * Compute the kinetic energy of the system at the current time.
   */
  double computeKineticEnergy();
    
private:
    double temperature, friction;
    int randomNumberSeed;
    Kernel kernel;

    int bias_method;
    int softcore_method;
    double lambdac, gammac, wbcoeff, w0coeff;
    double lambda1, lambda2, alpha, u0;
    int nParticles; //number of atoms in the system

    //displacement map
    vector<Vec3> displ;

    //soft core parameters
    double Umax;
    double Acore;
    double Ubcore;

    //nonequilibrium max time and work value
    int nonequilibrium;
    double Ntmax, work_value;
    double m1Lambda1, m2Lambda2, m_u0, m_w0, b1Lambda1, b2Lambda2, b_u0, b_w0;

    double PotEnergy;//current alchemical potential energy
    double BindE; // current binding energy
};

} // namespace SDMPlugin

#endif /*OPENMM_LANGEVININTEGRATORSDM_H_*/ 
