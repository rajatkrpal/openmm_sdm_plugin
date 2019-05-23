from __future__ import print_function

from simtk import openmm as mm
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import os, re,time, shutil, math
from simtk.openmm.app.desmonddmsfile import *
from datetime import datetime
from SDMplugin import *

print("Started at: " + str(time.asctime()))

binding_file = 'test.out'
f = open(binding_file, 'w')

lmbd = 0.9
lambda1 = lmbd
lambda2 = lmbd
alpha = 0.0 / kilocalorie_per_mole
u0 = 0.0 * kilocalorie_per_mole
w0coeff = 0.0 * kilocalorie_per_mole

print("lambda = ", lmbd)
print("lambda1 = ", lambda1)
print("lambda2 = ", lambda2)
print("alpha = ", alpha)
print("u0 = ", u0)
print("w0coeff ", w0coeff)

ligfile_input   = 'benzene_lig.dms'
rcptfile_input  = 'bcd_recpt.dms'
ligfile_output  = 'benzene_lig_tmp.dms'
rcptfile_output = 'bcd_recpt_tmp.dms'

#rcptfile_input  = 'mcl1-23_rcpt_2.dms'
#ligfile_input   = 'mcl1-23_lig_2.dms'
#rcptfile_output = 'mcl1-23_rcpt_tmp.dms'
#ligfile_output  = 'mcl1-23_lig_tmp.dms'


shutil.copyfile(rcptfile_input, rcptfile_output)
shutil.copyfile(ligfile_input, ligfile_output)
testDes = DesmondDMSFile([ligfile_output, rcptfile_output])

#system = testDes.createSystem(nonbondedMethod=CutoffNonPeriodic,nonbondedCutoff=1.2*nanometer, OPLS = True, implicitSolvent='AGBNP')
system = testDes.createSystem(nonbondedMethod=CutoffNonPeriodic,nonbondedCutoff=1.2*nanometer, OPLS = True)


#mcl1-23
#natoms_ligand = 43
#lig_atom_restr = 12   #index of ligand atom for Vsite (267CA -- central C of ligand)
#rcpt_atom_restr = 1535   #index of rcpt atom center of Vsite#267CA -- central C of ligand
#kf = 1225.0 #force constant for Vsite in (kj/mol)/nm^2
#r0 = 0.7 #radius of Vsite sphere

#beta-cyclodextrin and benzene
natoms_ligand = 12
lig_atoms = range(natoms_ligand)

rcpt_atom_restr = [6, 72, 61, 50, 39, 28, 17]

lig_atom_restr = lig_atoms
kf = 0.6 * kilocalorie_per_mole/angstrom**2 #force constant for Vsite CM-CM restraint 
r0 = 3.0 * angstrom #radius of Vsite sphere #radius of Vsite sphere

lig_ref_atoms = [0,3,2]
rcpt_ref_atoms = [6,50,28]
for i in range(len(rcpt_atom_restr)):
    rcpt_atom_restr[i] += natoms_ligand
for i in range(len(rcpt_ref_atoms)):
    rcpt_ref_atoms[i] += natoms_ligand


sdm_utils = SDMUtils(system, lig_atoms)
sdm_utils.addRestraintForce(lig_cm_particles = lig_atoms,
                            rcpt_cm_particles = rcpt_atom_restr,
                            kfcm = kf, tolcm = r0,
                            lig_ref_particles = None,
                            rcpt_ref_particles = None,
                            angle_center = 90.0,
                            angletol = 20,
                            kfangle = 0.0,
                            dihedral1center = 90.,
                            kfdihedral1 = 0.0,
                            dihedral2center = 0.,
                            kfdihedral2 = 1.0)

#platform = Platform.getPlatformByName('Reference')
platform = Platform.getPlatformByName('OpenCL')
properties = {}

temperature = 300.0 * kelvin
frictionCoeff = 0.5 / picosecond
MDstepsize = 0.001 * picosecond

integrator = LangevinIntegratorSDM(temperature/kelvin, frictionCoeff/(1/picosecond), MDstepsize/ picosecond, lig_atoms)
integrator.setBiasMethod(sdm_utils.ILogisticMethod)
integrator.setLambda(lmbd)
integrator.setLambda1(lambda1)
integrator.setLambda2(lambda2)
integrator.setAlpha(alpha*kilojoule_per_mole)
integrator.setU0(u0/ kilojoule_per_mole)
integrator.setW0coeff(w0coeff / kilojoule_per_mole) 

simulation = Simulation(testDes.topology, system, integrator,platform, properties)
print("Using platform %s" % simulation.context.getPlatform().getName())
simulation.context.setPositions(testDes.positions)
simulation.context.setVelocities(testDes.velocities)


totalSteps = 3000
nprnt = 100
ntrj = 100
simulation.reporters.append(StateDataReporter(stdout, nprnt, step=True, temperature=True))
simulation.reporters.append(DCDReporter("test.dcd", ntrj))
simulation.reporters.append(PDBReporter("test.pdb", totalSteps))

loops = totalSteps/nprnt
start=datetime.now()
step = 0
for i in range(loops):
    simulation.step(nprnt)
    pot_energy = (integrator.getPotEnergy()*kilojoule_per_mole).value_in_unit(kilocalorie_per_mole)
    bind_energy = (integrator.getBindE()*kilojoule_per_mole).value_in_unit(kilocalorie_per_mole)
    f.write("%f %f %f %f %f %f %f %f\n" % (lmbd, lambda1, lambda2, alpha*kilocalorie_per_mole, u0/kilocalorie_per_mole, w0coeff/kilocalorie_per_mole, pot_energy, bind_energy))
    f.flush()
    step += nprnt
end=datetime.now()

positions = simulation.context.getState(getPositions=True).getPositions()
velocities = simulation.context.getState(getVelocities=True).getVelocities()
testDes.setPositions(positions)
testDes.setVelocities(velocities)
testDes.close()


f.close()
elapsed=end - start
print("MD time="+str(elapsed.seconds+elapsed.microseconds*1e-6)+"s")
