%pythoncode %{

class SDMUtils(object):
    
    def __init__(self, system, ligparticles):
        self.system = system
        self.ligparticles = ligparticles
        self.RestraintControlParameterName = "SDMRestraintControlParameter"
        self.force = None
        self.LinearMethod = 0
        self.QuadraticMethod = 1
        self.ILogisticMethod = 2

        self.NoSoftCoreMethod = 0
        self.TanhSoftCoreMethod = 1
        self.RationalSoftCoreMethod = 2
        
    def _setRcptReferenceParticles(self, rcpt_ref_particles):
        if len(rcpt_ref_particles) != 3:
            raise ValueError("Invalid number of reference particles")
        self.rcpt_ref_particles = rcpt_ref_particles

    def _setLigReferenceParticles(self, lig_ref_particles):
        if len(lig_ref_particles) != 3:
            raise ValueError("Invalid number of reference particles")
        if not all(x in self.ligparticles for x in lig_ref_particles):
            raise ValueError("Invalid ligand atom indexes")
        self.lig_ref_particles = lig_ref_particles

    def getControlParameterName(self):
        return self.RestraintControlParameterName

        
    def addRestraintForce(self,
                          lig_cm_particles = None, rcpt_cm_particles = None,
                          kfcm = 0.0 * kilocalorie_per_mole/angstrom**2,
                          tolcm = 0.0 * angstrom,
                          lig_ref_particles = None, rcpt_ref_particles = None,
                          angle_center = 0.*degrees,
                          kfangle = 0.0 * kilocalorie_per_mole/degrees**2,
                          angletol = 10.0*degrees,
                          dihedral1center = 0.*degrees,
                          kfdihedral1 = 0.0 * kilocalorie_per_mole/degrees**2,
                          dihedral1tol = 10.0*degrees,
                          dihedral2center = 0.*degrees,
                          kfdihedral2 = 0.0 * kilocalorie_per_mole/degrees**2,
                          dihedral2tol = 10.0*degrees):

        if not (lig_cm_particles and rcpt_cm_particles):
            return

        if not (len(lig_cm_particles) > 0 and len(rcpt_cm_particles) > 0):
            return
        
        do_angles = False
        if lig_ref_particles and rcpt_ref_particles:
            if not (len(lig_ref_particles) == 3 and len(rcpt_ref_particles) == 3):
                raise ValueError("Invalid lists of reference atoms")
            do_angles = True
        
        expr = ""
            
        expr += "%s*( (kfcm/2)*step(d12-tolcm)*(d12-tolcm)^2 )                      " % self.RestraintControlParameterName
        
        if do_angles:
            expr += " + %s*( (kfcd0/2)*(step(dm0)*max(0,db0)^2+step(-dm0)*max(0,-da0)^2) )  " % self.RestraintControlParameterName
            expr += " + %s*( (kfcd1/2)*(step(dm1)*max(0,db1)^2+step(-dm1)*max(0,-da1)^2) )  " % self.RestraintControlParameterName
            expr += " + %s*( (kfcd2/2)*(step(dm2)*max(0,db2)^2+step(-dm2)*max(0,-da2)^2) )  " % self.RestraintControlParameterName

        expr += " ; d12 = distance(g1,g2) ; "

        if do_angles:
            expr += "db0 = xb0 - pi*floor(xb0/pi + 0.5)  ; xb0 = theta - b0 ; "
            expr += "da0 = xa0 - pi*floor(xa0/pi + 0.5)  ; xa0 = theta - a0 ; "
            expr += "dm0 = xm0 - pi*floor(xm0/pi + 0.5)  ; xm0 = theta - mid0 ; mid0 = (a0 + b0)/2 ; theta = angle(g3,g6,g7) ; "

            expr += "db1 = xb1 - twopi*floor(xb1/twopi + 0.5)  ; xb1 = phi1 - b1 ; "
            expr += "da1 = xa1 - twopi*floor(xa1/twopi + 0.5)  ; xa1 = phi1 - a1 ; "
            expr += "dm1 = xm1 - twopi*floor(xm1/twopi + 0.5)  ; xm1 = phi1 - mid1 ; mid1 = (a1 + b1)/2 ; phi1 = dihedral(g4,g3,g6,g7) ; "
                
            expr += "db2 = xb2 - twopi*floor(xb2/twopi + 0.5)  ; xb2 = phi2 - b2 ; "
            expr += "da2 = xa2 - twopi*floor(xa2/twopi + 0.5)  ; xa2 = phi2 - a2 ; "
            expr += "dm2 = xm2 - twopi*floor(xm2/twopi + 0.5)  ; xm2 = phi2 - mid2 ; mid2 = (a2 + b2)/2 ; phi2 = dihedral(g3,g6,g7,g8) ; "

        expr += "pi = SDpi ;"
        expr += "twopi = 2*SDpi"
                
        if do_angles:
            self.force =  mm.CustomCentroidBondForce(8,expr)
        else:
            self.force =  mm.CustomCentroidBondForce(2,expr)

        self.system.addForce(self.force)

        self.force.setForceGroup(1) #the restraint force will be evaluated separately
        
        self.force.addGlobalParameter(self.RestraintControlParameterName, 1.0)
        self.force.addGlobalParameter("SDpi", math.pi)

        self.force.addPerBondParameter("kfcm")
        self.force.addPerBondParameter("tolcm")
            
        self.force.addGroup(lig_cm_particles) #g1 CM of lig
        self.force.addGroup(rcpt_cm_particles) #g2 CM of rcpt

        if do_angles:
            
            self.force.addPerBondParameter("kfcd0")
            self.force.addPerBondParameter("a0")
            self.force.addPerBondParameter("b0")

            self.force.addPerBondParameter("kfcd1")
            self.force.addPerBondParameter("a1")
            self.force.addPerBondParameter("b1")

            self.force.addPerBondParameter("kfcd2")
            self.force.addPerBondParameter("a2")
            self.force.addPerBondParameter("b2")
            
            self.force.addGroup([rcpt_ref_particles[0]]) #g3 rcpt ref 0 
            self.force.addGroup([rcpt_ref_particles[1]]) #g4 rcpt ref 1
            self.force.addGroup([rcpt_ref_particles[2]]) #g5 rcpt ref 2
            
            self.force.addGroup([lig_ref_particles[0]]) #g6 lig ref 0 
            self.force.addGroup([lig_ref_particles[1]]) #g7 lig ref 1
            self.force.addGroup([lig_ref_particles[2]]) #g8 lig ref 2

        kfc = kfcm / (kilojoule_per_mole/radians**2)
        tolc = tolcm / nanometer
            
        if do_angles:

            kfcd0 = kfangle / (kilojoule_per_mole/radians**2)
            a0 = (angle_center - angletol)/radians
            b0 = (angle_center + angletol)/radians

            kfcd1 = kfdihedral1 / (kilojoule_per_mole/radians**2)
            a1 = (dihedral1center - dihedral1tol)/radians
            b1 = (dihedral1center + dihedral1tol)/radians
            
            kfcd2 = kfdihedral2 / (kilojoule_per_mole/radians**2)
            a2 = (dihedral2center - dihedral2tol)/radians
            b2 = (dihedral2center + dihedral2tol)/radians

            groups = range(8)
            params = [kfc, tolc, kfcd0, a0, b0, kfcd1, a1, b1, kfcd2, a2, b2]
        else:
            groups = [0,1]
            params = [kfc, tolc]

        self.force.addBond(groups, params)
%}
