%pythoncode %{

class SDMUtils(object):

    def __init__(self, system):
        self.system = system
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
                          dihedral2tol = 10.0*degrees,
                          offset = [0., 0., 0.] * angstrom):

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

        expr += " ; d12 = sqrt((x1 - offx - x2)^2 + (y1 - offy - y2)^2 + (z1 - offz - z2)^2 ) ; "

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
        self.force.addPerBondParameter("offx")
        self.force.addPerBondParameter("offy")
        self.force.addPerBondParameter("offz")

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
        offv = offset / nanometer
        offx = offv[0]
        offy = offv[1]
        offz = offv[2]

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
            params = [kfc, tolc, offx, offy, offz, kfcd0, a0, b0, kfcd1, a1, b1, kfcd2, a2, b2]
        else:
            groups = [0,1]
            params = [kfc, tolc, offx, offy, offz]

        self.force.addBond(groups, params)


    # a force to keep two ligands aligned
    # keeps ref atoms 1 near each other when offset is added
    #    lig2 is assumed to be translated by offset relative to lig1
    #     (lig2pos = lig1pos + offset)
    # keeps theta angle near zero
    # keeps psi angle near zero (assuming theta is near zero)
    def addAlignmentForce(self,
                          liga_ref_particles = None, ligb_ref_particles = None,
                          kfdispl = 0.0 * kilocalorie_per_mole/angstrom**2,
                          ktheta = 0.0 * kilocalorie_per_mole,
                          kpsi = 0.0 * kilocalorie_per_mole,
                          offset = [0., 0., 0.] * angstrom):

        if liga_ref_particles and ligb_ref_particles:
            if not (len(liga_ref_particles) == 3 and len(ligb_ref_particles) == 3):
                raise ValueError("Invalid lists of reference atoms")

        expr  = " (kfdispl/2)*distsq  ; " #displacement potential
        expr += " distsq = (x1 - offx - x2)^2 + "
        expr += "          (y1 - offy - y2)^2 + "
        expr += "          (z1 - offz - z2)^2   " #square distance between b1 and a1 after displacing back b

        displforce =  mm.CustomCompoundBondForce(2, expr);
        self.system.addForce(displforce)
        displforce.addPerBondParameter("kfdispl")
        displforce.addPerBondParameter("offx")
        displforce.addPerBondParameter("offy")
        displforce.addPerBondParameter("offz")
        offv = offset / nanometer
        offx = offv[0]
        offy = offv[1]
        offz = offv[2]
        displforce.addBond([ligb_ref_particles[0],liga_ref_particles[0]] ,
                           [kfdispl/((kilojoule_per_mole/nanometer**2)),
                            offx, offy, offz])
        displforce.setForceGroup(1);


        expr = "(ktheta/2) * (1 - cost) ;" #theta restraint

        expr += "cost = xdn1*xdn2 + ydn1*ydn2 + zdn1*zdn2 ; "
        expr += "xdn1 = xd1/dn1 ; ydn1 = yd1/dn1 ; zdn1 = zd1/dn1 ;"
        expr += "dn1 = sqrt(xd1^2 + yd1^2 + zd1^2 ) ;"
        expr += "xd1 = x2 - x1 ; "
        expr += "yd1 = y2 - y1 ;"
        expr += "zd1 = z2 - z1 ;"
        expr += "xdn2 = xd2/dn2 ; ydn2 = yd2/dn2 ; zdn2 = zd2/dn2 ;"
        expr += "dn2 = sqrt(xd2^2 + yd2^2 + zd2^2 ) ;"
        expr += "xd2 = x4 - x3 ; "
        expr += "yd2 = y4 - y3 ; "
        expr += "zd2 = z4 - z3   "

        thetaforce = mm.CustomCompoundBondForce(4, expr)
        self.system.addForce(thetaforce)
        thetaforce.addPerBondParameter("ktheta")
        thetaforce.addBond([ligb_ref_particles[0],ligb_ref_particles[1],liga_ref_particles[0],liga_ref_particles[1] ] ,
                           [ktheta/kilojoule_per_mole])
        thetaforce.setForceGroup(1);

        expr = "(kpsi/2) * (1 - cosp) ;" #psi restraint

        expr += "cosp = xvn*xwn + yvn*ywn + zvn*zwn ; "
        expr += "xvn = xv/v ; yvn = yv/v; zvn = zv/v ;"
        expr += "v = sqrt(xv^2 + yv^2 + zv^2 ) ;"
        expr += "xv = xd0 - dot01*xdn1 ;"
        expr += "yv = yd0 - dot01*ydn1 ;"
        expr += "zv = zd0 - dot01*zdn1 ;"
        expr += "dot01 =  xd0*xdn1 +  yd0*ydn1 +  zd0*zdn1 ;"
        expr += "xd0 = x3 - x1 ;"
        expr += "yd0 = y3 - y1 ;"
        expr += "zd0 = z3 - z1 ;"
        expr += "xwn = xw/w ; ywn = yw/w; zwn = zw/w ;"
        expr += "w = sqrt(xw^2 + yw^2 + zw^2 ) ;"
        expr += "xw = xd3 - dot31*xdn1 ;"
        expr += "yw = yd3 - dot31*ydn1 ;"
        expr += "zw = zd3 - dot31*zdn1 ;"
        expr += "dot31 =  xd3*xdn1 +  yd3*ydn1 +  zd3*zdn1 ;"
        expr += "xd3 = x5 - x4 ;"
        expr += "yd3 = y5 - y4 ;"
        expr += "zd3 = z5 - z4 ; "
        expr += "xdn1 = xd1/dn1 ; ydn1 = yd1/dn1 ; zdn1 = zd1/dn1 ;"
        expr += "dn1 = sqrt(xd1^2 + yd1^2 + zd1^2 ) ;"
        expr += "xd1 = x2 - x1 ; "
        expr += "yd1 = y2 - y1 ;"
        expr += "zd1 = z2 - z1 "

        psiforce = mm.CustomCompoundBondForce(5, expr)
        self.system.addForce(psiforce)
        psiforce.addPerBondParameter("kpsi")
        psiforce.addBond([ligb_ref_particles[0],ligb_ref_particles[1],ligb_ref_particles[2],
                          liga_ref_particles[0],liga_ref_particles[2] ] ,
                           [0.5*kpsi/kilojoule_per_mole])
        #symmetrize
        psiforce.addBond([liga_ref_particles[0],liga_ref_particles[1],liga_ref_particles[2],
                          ligb_ref_particles[0],ligb_ref_particles[2] ] ,
                           [0.5*kpsi/kilojoule_per_mole])

        psiforce.setForceGroup(1)

%}
