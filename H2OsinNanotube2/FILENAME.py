#################################
########### IMPORTS #############
#################################

# must import argv to pass in parameters such as number of steps
# must import numpy cause.... without it nothing is possible

from sys import argv,exit
import numpy as N
from numpy import zeros, arccos, cos, sin, sqrt, pi, dot, asarray, sign, arctan, exp

#################################
########### SUMMARY #############
#################################

# We a cartesian coordinate space, and euler angle space

# The cartesian coordinate space is used to calculate dipole moment vectors so we can evaluate the potential.

# The potential is a function of the magnitude of the dipole moments, their seperation distance r, and their orientations.

# In the RotOnly case, the orientations are the only variable because the magnitude of the dipole moment doesn't change,
# and the seperation distance doesn't change.

# For our sampling we randomly propose changes in euler angles such that the acceptance ratio is 0.3. This acceptance ratio is a function
# of the average change in angles, which is a function of the Rotational Step (a parameter we pass in using argv).

# Our potential is the dipole dipole interaction potential, and no i'm not going to retype it here cause it will be explictly stated
# when I use it to calculate the potential. You can find it here (except this is technically wrong, its the negative of what is shown here);
# https://chem.libretexts.org/Textbook_Maps/Physical_and_Theoretical_Chemistry_Textbook_Maps/Supplemental_Modules_%28Physical_and_Theoretical_Chemistry%29/Physical_Properties_of_Matter/Atomic_and_Molecular_Properties/Intermolecular_Forces/Specific_Interactions/Dipole-Dipole_Interactions

#################################
######## INITIALIZATION #########
#################################

rhoden = argv[1]
erotden = argv[2]
esqden = argv[3]
rotstep = float(argv[4])
nsteps = int(argv[5])

rhomat = zeros( 23588101, float )
erotmat = zeros( 23588101, float )
esqmat = zeros( 23588101, float )

rhoden = open(rhoden, "r")
erotden = open(erotden, "r")
esqden = open(esqden, "r")

for i in range(23588101):
    rhomat[i] = float(rhoden.readline())
    erotmat[i] = float(erotden.readline())
    esqmat[i] = float(esqden.readline())

k_B = 1.38e-23
temperature = 30.0
P = 4
Nmol = 2
Natom = 3

beta = 1.0/temperature

# The convention of (0,1,2,3) being molecule 1, hydrogen 1, beads 1 through 4, (4,5,6,7) being molecule 1
# hydrogen 2, beads 1 through 4, (8,9,10,11) molecule 1, oxygen, beads 1 through 4, (12,13,14,15) being molecule 2,
# hydrogen 1, beads 1 through 4, (16,17,18,19) being molecule 1 hydrogen 2, beads 1 through 4, (20,21,22,23)
# molecule 2, oxygen, beads 1 through 4 is followed here.

x = zeros( (P*Nmol*Natom,3), float )

xcm = zeros( (P*Nmol , 3), float )
rcm = zeros( (P*Nmol*Natom ,3) , float )
rcm_BFF = zeros( (P*Nmol*Natom ,3) , float )

phi = zeros( (P*Nmol*Natom) , float )
theta = zeros( (P*Nmol*Natom) , float )
chi = zeros( (P*Nmol*Natom) , float )

costheta = zeros( (P*Nmol*Natom) , float )

# x is xcm + rcm
# xcm should always be (0.0, 0.0, 0.0) for rows 0-11, and should always be (0.6, 0.0, 0.0) for rows 12-23, so we set it as such.
# rcm is the vector defined as x - xcm, it descibes the atoms, and the length of this vector is preserved after rotations (i.e det(R) is 1).

for i in range(P*Nmol):
    if i > (P*Nmol/2):
        xcm[i][0] = 0.6

# the water molecule is initialized with its dipole pointing along the z-axis, because of this, the oxygen will have 0.0 for x and y
# coordinates. The water molecule is placed in the xz plane so the hydrogens start with 0.0 for their y coordinates.

for i in range(P*Nmol*Natom):

    if (((i >= 0) and (i < (Natom-2)*P)) or ((i >= P*Natom) and (i < (P*(Natom+1))))):
        rcm_BFF[i][0] = -0.076
        rcm_BFF[i][1] = 0.0
        rcm_BFF[i][2] = -0.0522425412

    if (((i >= P) and (i < (Natom-1)*P)) or ((i >= P*(Natom+1)) and (i < (P*(Natom+2))))):
        rcm_BFF[i][0] = 0.076
        rcm_BFF[i][1] = 0.0
        rcm_BFF[i][2] = -0.0522425412

    if (((i >= 2*P) and (i < Natom*P)) or ((i >= P*(Natom+2)) and (i < (P*(Natom+3))))):
        rcm_BFF[i][0] = 0.0
        rcm_BFF[i][1] = 0.0
        rcm_BFF[i][2] = 0.00658268935

#################################
########## SIMULATION ###########
#################################

# There are two components to the sampling, the ratio of the old and new density matrix values, and exp((-beta/P)*(pot_new - pot_old)).
# Beta is equal to 1/(k_b*T), where k_b is the boltzman constant.

# This value outputs a number from 0 to inf. If this number is greater than a random number between 0 and 1
# then the new position and orientation is accepted.

# We generate new density matrix values using the difference in Euler Angles between corresponding beads.
# Our potential is also a function of the Euler Angles (our only variable) because the orientations are based
# on the Euler Angles, and the potential depends on the orientations.

# We pseodorandomly propose new Euler Angles by multiplying our rotational step (required for acceptance ratio of 0.3) by a
# random number between 0 and 1. This allows for any change in angle (for one simulation step) between 0 and the magnitude of
# our rotational step.

bondlength = 0.05882523055 # Bond Length (nm)
q = (1.8546 / (bondlength*1.0e-9))*3.33564e-30

massO = 15.9993047236
massH = 1.00797597651

massOpercent = massO/(massO + 2*massH)
massHpercent = massH/(massO + 2*massH)

rotmat1 = zeros( (3,3) , float )
rotmat2 = zeros( (3,3) , float )
rotmar = zeros( (3,3) , float )

acceptance_ratio = 0.0
small = 1e-10
average_pot = 0.0
average_erot = 0.0

acceptanceratiofile = open("acceptanceratiofile-"+str(Nmol)+"H20T"+str(temperature)+"P"+str(P)+str(nsteps)+"Steps","w")
potfile = open("potfile-"+str(Nmol)+"H20T"+str(temperature)+"P"+str(P)+str(nsteps)+"Steps","w")
erotfile = open("erotfile-"+str(Nmol)+"H20T"+str(temperature)+"P"+str(P)+str(nsteps)+"Steps","w")

for i in range(nsteps):
    for j in range(P):
        for m in range(Nmol):

            # finding the index

            t0b=j-1
            t1b=j
            t2b=j+1

            if (t0b<0):
                t0b += P

            if (t2b>(P-1)):
                t2b -= P

            t0=m*P+t0b
            t1=m*P+t1b
            t2=m*P+t2b

            rotmat1[0][0] = cos(phi[t1])*cos(theta[t1])*cos(chi[t1]) - sin(phi[t1])*sin(chi[t1])
            rotmat1[0][1] = -cos(phi[t1])*cos(theta[t1])*sin(chi[t1]) - sin(phi[t1])*cos(chi[t1])
            rotmat1[0][2] = cos(phi[t1])*sin(theta[t1])
            rotmat1[1][0] = sin(phi[t1])*cos(theta[t1])*cos(chi[t1]) + cos(phi[t1])*sin(chi[t1])
            rotmat1[1][1] = -sin(phi[t1])*cos(theta[t1])*sin(chi[t1]) + cos(phi[t1])*cos(chi[t1])
            rotmat1[1][2] = sin(phi[t1])*sin(theta[t1])
            rotmat1[2][0] = -sin(theta[t1])*cos(chi[t1])
            rotmat1[2][1] = sin(theta[t1])*sin(chi[t1])
            rotmat1[2][2] = cos(theta[t1])

            rcm = zeros( (P*Nmol*Natom , 3 ), float )

            for z in range(Nmol*Natom):
                for w in range(P):
                    for k in range(3):
                        for l in range(3):
                            rcm[z*P+w][k] += rotmat1[k][l]*rcm_BFF[z*P+w][l]

            for z in range(Nmol*Natom):
                for w in range(P):
                    for k in range(3):
                        if z < Natom:
                            x[z*P+w][k] = rcm[z*P+w][k]
                        else:
                            x[z*P+w][0] = 0.6 + rcm[z*P+w][0]
                            x[z*P+w][1] = rcm[z*P+w][1]
                            x[z*P+w][2] = rcm[z*P+w][2]

            # not really mu, has units [nm] not [C*nm] like an actual dipole momen

            mu_x1 = x[(Natom-1)*P + j][0] - (x[(Natom-2)*P + j][0] + x[(Natom-3)*P + j][0])/2
            mu_y1 = x[(Natom-1)*P + j][1] - (x[(Natom-2)*P + j][1] + x[(Natom-3)*P + j][1])/2
            mu_z1 = x[(Natom-1)*P + j][2] - (x[(Natom-2)*P + j][2] + x[(Natom-3)*P + j][2])/2

            mu_x2 = x[(Natom+2)*P + j][0] - (x[(Natom+1)*P + j][0] + x[(Natom)*P + j][0])/2
            mu_y2 = x[(Natom+2)*P + j][1] - (x[(Natom+1)*P + j][1] + x[(Natom)*P + j][1])/2
            mu_z2 = x[(Natom+2)*P + j][2] - (x[(Natom+1)*P + j][2] + x[(Natom)*P + j][2])/2

            cmvec_x = (massOpercent*(x[(Natom+2)*P + j][0] - x[(Natom-1)*P + j][0])) + massHpercent*((x[(Natom+1)*P + j][0] + x[(Natom)*P + j][0]) - (x[(Natom-2)*P + j][0] + x[(Natom-3)*P + j][0]))
            cmvec_y = (massOpercent*(x[(Natom+2)*P + j][1] - x[(Natom-1)*P + j][1])) + massHpercent*((x[(Natom+1)*P + j][1] + x[(Natom)*P + j][1]) - (x[(Natom-2)*P + j][1] + x[(Natom-3)*P + j][1]))
            cmvec_z = (massOpercent*(x[(Natom+2)*P + j][2] - x[(Natom-1)*P + j][2])) + massHpercent*((x[(Natom+1)*P + j][2] + x[(Natom)*P + j][2]) - (x[(Natom-2)*P + j][2] + x[(Natom-3)*P + j][2]))

            r = sqrt(cmvec_x*cmvec_x + cmvec_y*cmvec_y + cmvec_z*cmvec_z)

            mu1dotmu2 = mu_x1*mu_x2 + mu_y1*mu_y2 + mu_z1*mu_z2
            mu1dotr = mu_x1*cmvec_x + mu_y1*cmvec_y + mu_z1*cmvec_z
            mu2dotr = mu_x2*cmvec_x + mu_y2*cmvec_y + mu_z2*cmvec_z

            pot_old = (q*q/(4.0*pi*k_B*r*r*r*8.854187817e-12))*(mu1dotmu2 - 3.0*mu1dotr*mu2dotr/(r*r))*1e9

            # next we calculate rho_old, which requires first calulcating the difference in Euler Angles between bead t0 and t1,
            # next we use a precalculated density matrix file to determine rho_old1, erot_old1, rho_old2, erot_old2.

            rotmat1[0][0] = cos(phi[t0])*cos(theta[t0])*cos(chi[t0]) - sin(phi[t0])*sin(chi[t0])
            rotmat1[0][1] = -cos(phi[t0])*cos(theta[t0])*sin(chi[t0]) - sin(phi[t0])*cos(chi[t0])
            rotmat1[0][2] = cos(phi[t0])*sin(theta[t0])
            rotmat1[1][0] = sin(phi[t0])*cos(theta[t0])*cos(chi[t0]) + cos(phi[t0])*sin(chi[t0])
            rotmat1[1][1] = -sin(phi[t0])*cos(theta[t0])*sin(chi[t0]) + cos(phi[t0])*cos(chi[t0])
            rotmat1[1][2] = sin(phi[t0])*sin(theta[t0])
            rotmat1[2][0] = -sin(theta[t0])*cos(chi[t0])
            rotmat1[2][1] = sin(theta[t0])*sin(chi[t0])
            rotmat1[2][2] = cos(theta[t0])

            rotmat2[0][0] = cos(phi[t1])*cos(theta[t1])*cos(chi[t1]) - sin(phi[t1])*sin(chi[t1])
            rotmat2[0][1] = -cos(phi[t1])*cos(theta[t1])*sin(chi[t1]) - sin(phi[t1])*cos(chi[t1])
            rotmat2[0][2] = cos(phi[t1])*sin(theta[t1])
            rotmat2[1][0] = sin(phi[t1])*cos(theta[t1])*cos(chi[t1]) + cos(phi[t1])*sin(chi[t1])
            rotmat2[1][1] = -sin(phi[t1])*cos(theta[t1])*sin(chi[t1]) + cos(phi[t1])*cos(chi[t1])
            rotmat2[1][2] = sin(phi[t1])*sin(theta[t1])
            rotmat2[2][0] = -sin(theta[t1])*cos(chi[t1])
            rotmat2[2][1] = sin(theta[t1])*sin(chi[t1])
            rotmat2[2][2] = cos(theta[t1])

            for a in range(3):
                for b in range(3):
                    for c in range(3):
                        rotmar[a][b]=rotmar[a][b]+rotmat1[c][a]*rotmat2[c][b]

            cost = rotmar[2][2]
            if (cost > 1.0):
                cost = 1.0
            elif (cost < -1.0):
                cost = -1.0

            thetar = N.arccos(cost)
            sint = N.sin(thetar)

            if (N.abs(1.0-cost) < small):
                phir = 0.0
                cchi = rotmar[0][0]
                schi = rotmar[1][0]

                if (cchi > 1.0):
                    cchi = 1.0
                elif (cchi < -1.0):
                    cchi = -1.0

                if (schi > 1.0):
                    schi = 1.0
                elif (schi < -1.0):
                    schi = -1.0

                if (schi > 0.0):
                    chir=N.arccos(cchi)
                else:
                    chir=2.0*N.pi-N.arccos(cchi)

            elif (N.abs(1.0+cost) < small):

                phir=0.0
                cchi=rotmar[1][1]
                schi=rotmar[0][1]

                if (cchi > 1.0):
                    cchi = 1.0
                elif (cchi < -1.0):
                    cchi = -1.0

                if (schi > 1.0):
                    schi = 1.0
                elif (schi < -1.0):
                    schi = -1.0

                if (schi > 0.0):
                    chir=N.arccos(cchi)
                else:
                    chir=2.0*N.pi-N.arccos(cchi)
            else:
                cphi = rotmar[0][2]/sint
                sphi = rotmar[1][2]/sint
                cchi = -rotmar[2][0]/sint
                schi = rotmar[2][1]/sint

                if (cphi > 1.0):
                    cphi = 1.0
                elif (cphi < -1.0):
                    cphi = -1.0

                if (sphi > 1.0):
                    sphi = 1.0
                elif (sphi < -1.0):
                    sphi = -1.0

                if (cchi > 1.0):
                    cchi = 1.0
                elif (cchi < -1.0):
                    cchi = -1.0

                if (schi > 1.0):
                    schi = 1.0
                elif (schi < -1.0):
                    schi = -1.0

                if (sphi > 0.0):
                    phir=N.arccos(cphi)
                else:
                    phir=2.0*N.pi-N.arccos(cphi)
                if (schi > 0.0):
                    chir=N.arccos(cchi)
                else:
                    chir=2.0*N.pi-N.arccos(cchi)

            if (phir<0.0):
                phir = 2.0*N.pi+phir
            if (chir<0.0):
                chir = 2.0*N.pi+chir

            phir = phir % (2.0*N.pi)
            chir = chir % (2.0*N.pi)

            if (cost > 1.0):
                cost = 2.0 - cost
            elif (cost < -1.0):
                cost = -2.0 - cost

            thetar = N.arccos(cost)

            phir *= 180.0/N.pi
            thetar *= 180.0/N.pi
            chir *= 180.0/N.pi

            indchi = int(chir)
            indphi = int(phir)
            indthe = int(thetar)

            if indphi == 360:
                indphi == 0

            if indthe == 180:
                indthe == 0

            if indchi == 360:
                indchi == 0

            index_p0=(indthe*361+indphi)*361+indchi # index for point 0
            rho_0 = rhomat[index_p0]
            erot_0 = erotmat[index_p0]
            esq_0 = esqmat[index_p0]

            if indphi == 360:
                indphi == 0

            if indthe == 180:
                indthe == 0

            if indchi == 360:
                indchi == 0

            if (indchi > 360 or indchi < 0):
                print "index of chi is out of range", indchi, chi
                raise()
            elif (indphi > 360 or indphi < 0):
                print "index of phi is out of range", indphi, phi
                raise()
            elif (indthe > 180 or indthe < 0):
                print "index of theta is out of range", indthe, theta
                raise()

            if (indchi != 360):
                index_p1 = (indthe*361+indphi)*361+indchi+1 # next index is just the difference by 1 degree
                # CHI part for each properties
                delchi = rhomat[index_p1]-rho_0
                delch2 = erotmat[index_p1]-erot_0
                delch3 = esqmat[index_p1]-esq_0

            if (indphi != 360):
                index_p1 = (indthe*361+indphi+1)*361+indchi
                # PHI part for each properties
                delphi = rhomat[index_p1]-rho_0
                delph2 = erotmat[index_p1]-erot_0
                delph3 = esqmat[index_p1]-esq_0

            if (indthe != 180):
                index_p1 = ((indthe+1)*361+indphi)*361+indchi
                # THETA part for each properties
                delthe = rhomat[index_p1]-rho_0
                delth2 = erotmat[index_p1]-erot_0
                delth3 = esqmat[index_p1]-esq_0

            # Linear Interpolation for all three angles
            rho_old1 = rho_0+delchi*(chir-N.float(indchi))+delphi*(phir-N.float(indphi))+delthe*(thetar-N.float(indthe))
            erot_old1 = erot_0+delch2*(chir-N.float(indchi))+delph2*(phir-N.float(indphi))+delth2*(thetar-N.float(indthe))
            esq_old1 = esq_0+delch3*(chir-N.float(indchi))+delph3*(phir-N.float(indphi))+delth3*(thetar-N.float(indthe))

            # ok, now we gotta find the rho_old2, erot_old2, esq_old2, its basically the same process

            rotmat1[0][0] = cos(phi[t1])*cos(theta[t1])*cos(chi[t1]) - sin(phi[t1])*sin(chi[t1])
            rotmat1[0][1] = -cos(phi[t1])*cos(theta[t1])*sin(chi[t1]) - sin(phi[t1])*cos(chi[t1])
            rotmat1[0][2] = cos(phi[t1])*sin(theta[t1])
            rotmat1[1][0] = sin(phi[t1])*cos(theta[t1])*cos(chi[t1]) + cos(phi[t1])*sin(chi[t1])
            rotmat1[1][1] = -sin(phi[t1])*cos(theta[t1])*sin(chi[t1]) + cos(phi[t1])*cos(chi[t1])
            rotmat1[1][2] = sin(phi[t1])*sin(theta[t1])
            rotmat1[2][0] = -sin(theta[t1])*cos(chi[t1])
            rotmat1[2][1] = sin(theta[t1])*sin(chi[t1])
            rotmat1[2][2] = cos(theta[t1])

            rotmat2[0][0] = cos(phi[t2])*cos(theta[t2])*cos(chi[t2]) - sin(phi[t2])*sin(chi[t2])
            rotmat2[0][1] = -cos(phi[t2])*cos(theta[t2])*sin(chi[t2]) - sin(phi[t2])*cos(chi[t2])
            rotmat2[0][2] = cos(phi[t2])*sin(theta[t2])
            rotmat2[1][0] = sin(phi[t2])*cos(theta[t2])*cos(chi[t2]) + cos(phi[t2])*sin(chi[t2])
            rotmat2[1][1] = -sin(phi[t2])*cos(theta[t2])*sin(chi[t2]) + cos(phi[t2])*cos(chi[t2])
            rotmat2[1][2] = sin(phi[t2])*sin(theta[t2])
            rotmat2[2][0] = -sin(theta[t2])*cos(chi[t2])
            rotmat2[2][1] = sin(theta[t2])*sin(chi[t2])
            rotmat2[2][2] = cos(theta[t2])

            for a in range(3):
                for b in range(3):
                    for c in range(3):
                        rotmar[a][b]=rotmar[a][b]+rotmat1[c][a]*rotmat2[c][b]

            cost = rotmar[2][2]
            if (cost > 1.0):
                cost = 1.0
            elif (cost < -1.0):
                cost = -1.0

            thetar = N.arccos(cost)
            sint = N.sin(thetar)

            if (N.abs(1.0-cost) < small):
                phir = 0.0
                cchi = rotmar[0][0]
                schi = rotmar[1][0]

                if (cchi > 1.0):
                    cchi = 1.0
                elif (cchi < -1.0):
                    cchi = -1.0

                if (schi > 1.0):
                    schi = 1.0
                elif (schi < -1.0):
                    schi = -1.0

                if (schi > 0.0):
                    chir=N.arccos(cchi)
                else:
                    chir=2.0*N.pi-N.arccos(cchi)

            elif (N.abs(1.0+cost) < small):

                phir=0.0
                cchi=rotmar[1][1]
                schi=rotmar[0][1]

                if (cchi > 1.0):
                    cchi = 1.0
                elif (cchi < -1.0):
                    cchi = -1.0

                if (schi > 1.0):
                    schi = 1.0
                elif (schi < -1.0):
                    schi = -1.0

                if (schi > 0.0):
                    chir=N.arccos(cchi)
                else:
                    chir=2.0*N.pi-N.arccos(cchi)
            else:
                cphi = rotmar[0][2]/sint
                sphi = rotmar[1][2]/sint
                cchi = -rotmar[2][0]/sint
                schi = rotmar[2][1]/sint

                if (cphi > 1.0):
                    cphi = 1.0
                elif (cphi < -1.0):
                    cphi = -1.0

                if (sphi > 1.0):
                    sphi = 1.0
                elif (sphi < -1.0):
                    sphi = -1.0

                if (cchi > 1.0):
                    cchi = 1.0
                elif (cchi < -1.0):
                    cchi = -1.0

                if (schi > 1.0):
                    schi = 1.0
                elif (schi < -1.0):
                    schi = -1.0

                if (sphi > 0.0):
                    phir=N.arccos(cphi)
                else:
                    phir=2.0*N.pi-N.arccos(cphi)
                if (schi > 0.0):
                    chir=N.arccos(cchi)
                else:
                    chir=2.0*N.pi-N.arccos(cchi)

            if (phir<0.0):
                phir = 2.0*N.pi+phir
            if (chir<0.0):
                chir = 2.0*N.pi+chir

            phir = phir % (2.0*N.pi)
            chir = chir % (2.0*N.pi)

            if (cost > 1.0):
                cost = 2.0 - cost
            elif (cost < -1.0):
                cost = -2.0 - cost

            thetar = N.arccos(cost)

            phir *= 180.0/N.pi
            thetar *= 180.0/N.pi
            chir *= 180.0/N.pi

            indchi = int(chir)
            indphi = int(phir)
            indthe = int(thetar)

            if indphi == 360:
                indphi == 0

            if indthe == 180:
                indthe == 0

            if indchi == 360:
                indchi == 0

            index_p0=(indthe*361+indphi)*361+indchi # index for point 0
            rho_0=rhomat[index_p0]
            erot_0=erotmat[index_p0]
            esq_0=esqmat[index_p0]

            if (indchi > 360 or indchi < 0):
                print "index of chi is out of range", indchi, chi
                raise()
            elif (indphi > 360 or indphi < 0):
                print "index of phi is out of range", indphi, phi
                raise()
            elif (indthe > 180 or indthe < 0):
                print "index of theta is out of range", indthe, theta
                raise()

            if (indchi != 360):
                index_p1 = (indthe*361+indphi)*361+indchi+1 # next index is just the difference by 1 degree
                # CHI part for each properties
                delchi = rhomat[index_p1]-rho_0
                delch2 = erotmat[index_p1]-erot_0
                delch3 = esqmat[index_p1]-esq_0

            if (indphi != 360):
                index_p1 = (indthe*361+indphi+1)*361+indchi
                # PHI part for each properties
                delphi = rhomat[index_p1]-rho_0
                delph2 = erotmat[index_p1]-erot_0
                delph3 = esqmat[index_p1]-esq_0

            if (indthe != 180):
                index_p1 = ((indthe+1)*361+indphi)*361+indchi
                # THETA part for each properties
                delthe = rhomat[index_p1]-rho_0
                delth2 = erotmat[index_p1]-erot_0
                delth3 = esqmat[index_p1]-esq_0

            # Linear Interpolation for all three angles
            rho_old2 = rho_0+delchi*(chir-N.float(indchi))+delphi*(phir-N.float(indphi))+delthe*(thetar-N.float(indthe))
            erot_old2 = erot_0+delch2*(chir-N.float(indchi))+delph2*(phir-N.float(indphi))+delth2*(thetar-N.float(indthe))
            esq_old2 = esq_0+delch3*(chir-N.float(indchi))+delph3*(phir-N.float(indphi))+delth3*(thetar-N.float(indthe))


            # concatenate data

            rho_old = rho_old1*rho_old2
            erot_old = (erot_old1 + erot_old2)/2.0


            # now we save our cartesian and euler angle coordinates before we propose new euler angles

            rcm_old = zeros((Nmol*Natom*P,3),float)

            for z in range(Nmol*Natom):

                rcm_old[z*P + j][0] = rcm[z*P + j][0]
                rcm_old[z*P + j][1] = rcm[z*P + j][1]
                rcm_old[z*P + j][2] = rcm[z*P + j][2]

            for z in range(Nmol*Natom):
                for w in range(P):
                    costheta[z*P + w] = cos(theta[z*P + w])

            phi_old = phi[t1]
            theta_old = theta[t1]
            costheta_old = costheta[t1]
            chi_old = chi[t1]


            # now that we have saved our values we can go ahead and propose new ones

            phi[t1] += 2.0*N.pi*rotstep*(N.random.random()-0.5)
            costheta[t1] += rotstep*(N.random.random()-0.5)
            chi[t1] += 2.0*N.pi*rotstep*(N.random.random()-0.5)

            if (phi[t1] < 0.0):
                phi[t1] += 2.0*N.pi

            if (chi[t1] < 0.0):
                chi[t1] += 2.0*N.pi

            phi[t1] = phi[t1] % (2.0*N.pi)
            chi[t1] = chi[t1] % (2.0*N.pi)

            if (costheta[t1] > 1.0):
                costheta[t1] = 2.0 - costheta[t1]

            elif (costheta[t1] < -1.0):
                costheta[t1] = -2.0 - costheta[t1]

            theta[t1] = N.arccos(costheta[t1])

            # Now we use these proposed values to calculate basically all of that stuff again but we'll label things new instead of old

            rotmat1[0][0] = cos(phi[t1])*cos(theta[t1])*cos(chi[t1]) - sin(phi[t1])*sin(chi[t1])
            rotmat1[0][1] = -cos(phi[t1])*cos(theta[t1])*sin(chi[t1]) - sin(phi[t1])*cos(chi[t1])
            rotmat1[0][2] = cos(phi[t1])*sin(theta[t1])
            rotmat1[1][0] = sin(phi[t1])*cos(theta[t1])*cos(chi[t1]) + cos(phi[t1])*sin(chi[t1])
            rotmat1[1][1] = -sin(phi[t1])*cos(theta[t1])*sin(chi[t1]) + cos(phi[t1])*cos(chi[t1])
            rotmat1[1][2] = sin(phi[t1])*sin(theta[t1])
            rotmat1[2][0] = -sin(theta[t1])*cos(chi[t1])
            rotmat1[2][1] = sin(theta[t1])*sin(chi[t1])
            rotmat1[2][2] = cos(theta[t1])

            rcm = zeros( (P*Nmol*Natom , 3 ), float )

            for z in range(Nmol*Natom):
                for w in range(P):
                    for k in range(3):
                        for l in range(3):
                            rcm[z*P+w][k] += rotmat1[k][l]*rcm_BFF[z*P+w][l]

            for z in range(Nmol*Natom):
                for w in range(P):
                    for k in range(3):
                        if z < Natom:
                            x[z*P+w][k] = rcm[z*P+w][k]
                        else:
                            x[z*P+w][0] = 0.6 + rcm[z*P+w][0]
                            x[z*P+w][1] = rcm[z*P+w][1]
                            x[z*P+w][2] = rcm[z*P+w][2]

            mu_x1 = x[(Natom-1)*P + j][0] - (x[(Natom-2)*P + j][0] + x[(Natom-3)*P + j][0])/2
            mu_y1 = x[(Natom-1)*P + j][1] - (x[(Natom-2)*P + j][1] + x[(Natom-3)*P + j][1])/2
            mu_z1 = x[(Natom-1)*P + j][2] - (x[(Natom-2)*P + j][2] + x[(Natom-3)*P + j][2])/2

            mu_x2 = x[(Natom+2)*P + j][0] - (x[(Natom+1)*P + j][0] + x[(Natom)*P + j][0])/2
            mu_y2 = x[(Natom+2)*P + j][1] - (x[(Natom+1)*P + j][1] + x[(Natom)*P + j][1])/2
            mu_z2 = x[(Natom+2)*P + j][2] - (x[(Natom+1)*P + j][2] + x[(Natom)*P + j][2])/2

            cmvec_x = (massOpercent*(x[(Natom+2)*P + j][0] - x[(Natom-1)*P + j][0])) + massHpercent*((x[(Natom+1)*P + j][0] + x[(Natom)*P + j][0]) - (x[(Natom-2)*P + j][0] + x[(Natom-3)*P + j][0]))
            cmvec_y = (massOpercent*(x[(Natom+2)*P + j][1] - x[(Natom-1)*P + j][1])) + massHpercent*((x[(Natom+1)*P + j][1] + x[(Natom)*P + j][1]) - (x[(Natom-2)*P + j][1] + x[(Natom-3)*P + j][1]))
            cmvec_z = (massOpercent*(x[(Natom+2)*P + j][2] - x[(Natom-1)*P + j][2])) + massHpercent*((x[(Natom+1)*P + j][2] + x[(Natom)*P + j][2]) - (x[(Natom-2)*P + j][2] + x[(Natom-3)*P + j][2]))

            r = sqrt(cmvec_x*cmvec_x + cmvec_y*cmvec_y + cmvec_z*cmvec_z)

            mu1dotmu2 = mu_x1*mu_x2 + mu_y1*mu_y2 + mu_z1*mu_z2
            mu1dotr = mu_x1*cmvec_x + mu_y1*cmvec_y + mu_z1*cmvec_z
            mu2dotr = mu_x2*cmvec_x + mu_y2*cmvec_y + mu_z2*cmvec_z

            pot_new = (q*q/(4.0*pi*k_B*r*r*r*8.854187817e-12))*(mu1dotmu2 - 3.0*mu1dotr*mu2dotr/(r*r))*1e9

            # next we calculate rho_new, which requires first calulcating the difference in Euler Angles between bead t0 and t1,
            # next we use a precalculated density matrix file to determine rho_new1, erot_new1, rho_new2, erot_new2.

            rotmat1[0][0] = cos(phi[t0])*cos(theta[t0])*cos(chi[t0]) - sin(phi[t0])*sin(chi[t0])
            rotmat1[0][1] = -cos(phi[t0])*cos(theta[t0])*sin(chi[t0]) - sin(phi[t0])*cos(chi[t0])
            rotmat1[0][2] = cos(phi[t0])*sin(theta[t0])
            rotmat1[1][0] = sin(phi[t0])*cos(theta[t0])*cos(chi[t0]) + cos(phi[t0])*sin(chi[t0])
            rotmat1[1][1] = -sin(phi[t0])*cos(theta[t0])*sin(chi[t0]) + cos(phi[t0])*cos(chi[t0])
            rotmat1[1][2] = sin(phi[t0])*sin(theta[t0])
            rotmat1[2][0] = -sin(theta[t0])*cos(chi[t0])
            rotmat1[2][1] = sin(theta[t0])*sin(chi[t0])
            rotmat1[2][2] = cos(theta[t0])

            rotmat2[0][0] = cos(phi[t1])*cos(theta[t1])*cos(chi[t1]) - sin(phi[t1])*sin(chi[t1])
            rotmat2[0][1] = -cos(phi[t1])*cos(theta[t1])*sin(chi[t1]) - sin(phi[t1])*cos(chi[t1])
            rotmat2[0][2] = cos(phi[t1])*sin(theta[t1])
            rotmat2[1][0] = sin(phi[t1])*cos(theta[t1])*cos(chi[t1]) + cos(phi[t1])*sin(chi[t1])
            rotmat2[1][1] = -sin(phi[t1])*cos(theta[t1])*sin(chi[t1]) + cos(phi[t1])*cos(chi[t1])
            rotmat2[1][2] = sin(phi[t1])*sin(theta[t1])
            rotmat2[2][0] = -sin(theta[t1])*cos(chi[t1])
            rotmat2[2][1] = sin(theta[t1])*sin(chi[t1])
            rotmat2[2][2] = cos(theta[t1])

            for a in range(3):
                for b in range(3):
                    for c in range(3):
                        rotmar[a][b]=rotmar[a][b]+rotmat1[c][a]*rotmat2[c][b]

            cost = rotmar[2][2]
            if (cost > 1.0):
                cost = 1.0
            elif (cost < -1.0):
                cost = -1.0

            thetar = N.arccos(cost)
            sint = N.sin(thetar)

            if (N.abs(1.0-cost) < small):
                phir = 0.0
                cchi = rotmar[0][0]
                schi = rotmar[1][0]

                if (cchi > 1.0):
                    cchi = 1.0
                elif (cchi < -1.0):
                    cchi = -1.0

                if (schi > 1.0):
                    schi = 1.0
                elif (schi < -1.0):
                    schi = -1.0

                if (schi > 0.0):
                    chir=N.arccos(cchi)
                else:
                    chir=2.0*N.pi-N.arccos(cchi)

            elif (N.abs(1.0+cost) < small):

                phir=0.0
                cchi=rotmar[1][1]
                schi=rotmar[0][1]

                if (cchi > 1.0):
                    cchi = 1.0
                elif (cchi < -1.0):
                    cchi = -1.0

                if (schi > 1.0):
                    schi = 1.0
                elif (schi < -1.0):
                    schi = -1.0

                if (schi > 0.0):
                    chir=N.arccos(cchi)
                else:
                    chir=2.0*N.pi-N.arccos(cchi)
            else:
                cphi = rotmar[0][2]/sint
                sphi = rotmar[1][2]/sint
                cchi = -rotmar[2][0]/sint
                schi = rotmar[2][1]/sint

                if (cphi > 1.0):
                    cphi = 1.0
                elif (cphi < -1.0):
                    cphi = -1.0

                if (sphi > 1.0):
                    sphi = 1.0
                elif (sphi < -1.0):
                    sphi = -1.0

                if (cchi > 1.0):
                    cchi = 1.0
                elif (cchi < -1.0):
                    cchi = -1.0

                if (schi > 1.0):
                    schi = 1.0
                elif (schi < -1.0):
                    schi = -1.0

                if (sphi > 0.0):
                    phir=N.arccos(cphi)
                else:
                    phir=2.0*N.pi-N.arccos(cphi)
                if (schi > 0.0):
                    chir=N.arccos(cchi)
                else:
                    chir=2.0*N.pi-N.arccos(cchi)

            if (phir<0.0):
                phir = 2.0*N.pi+phir
            if (chir<0.0):
                chir = 2.0*N.pi+chir

            phir = phir % (2.0*N.pi)
            chir = chir % (2.0*N.pi)

            if (cost > 1.0):
                cost = 2.0 - cost
            elif (cost < -1.0):
                cost = -2.0 - cost

            thetar = N.arccos(cost)

            phir *= 180.0/N.pi
            thetar *= 180.0/N.pi
            chir *= 180.0/N.pi

            indchi = int(chir)
            indphi = int(phir)
            indthe = int(thetar)

            index_p0=(indthe*361+indphi)*361+indchi # index for point 0
            rho_0 = rhomat[index_p0]
            erot_0 = erotmat[index_p0]
            esq_0 = esqmat[index_p0]

            if indphi == 360:
                indphi == 0

            if indthe == 180:
                indthe == 0

            if indchi == 360:
                indchi == 0

            if (indchi > 360 or indchi < 0):
                print "index of chi is out of range", indchi, chi
                raise()
            elif (indphi > 360 or indphi < 0):
                print "index of phi is out of range", indphi, phi
                raise()
            elif (indthe > 180 or indthe < 0):
                print "index of theta is out of range", indthe, theta
                raise()

            if (indchi != 360):
                index_p1 = (indthe*361+indphi)*361+indchi+1 # next index is just the difference by 1 degree
                # CHI part for each properties
                delchi = rhomat[index_p1]-rho_0
                delch2 = erotmat[index_p1]-erot_0
                delch3 = esqmat[index_p1]-esq_0

            if (indphi != 360):
                index_p1 = (indthe*361+indphi+1)*361+indchi
                # PHI part for each properties
                delphi = rhomat[index_p1]-rho_0
                delph2 = erotmat[index_p1]-erot_0
                delph3 = esqmat[index_p1]-esq_0

            if (indthe != 180):
                index_p1 = ((indthe+1)*361+indphi)*361+indchi
                # THETA part for each properties
                delthe = rhomat[index_p1]-rho_0
                delth2 = erotmat[index_p1]-erot_0
                delth3 = esqmat[index_p1]-esq_0

            # Linear Interpolation for all three angles
            rho_new1 = rho_0+delchi*(chir-N.float(indchi))+delphi*(phir-N.float(indphi))+delthe*(thetar-N.float(indthe))
            erot_new1 = erot_0+delch2*(chir-N.float(indchi))+delph2*(phir-N.float(indphi))+delth2*(thetar-N.float(indthe))
            esq_new1 = esq_0+delch3*(chir-N.float(indchi))+delph3*(phir-N.float(indphi))+delth3*(thetar-N.float(indthe))

            # ok, now we gotta find the rho_old2, erot_old2, esq_old2, its basically the same process

            rotmat1[0][0] = cos(phi[t1])*cos(theta[t1])*cos(chi[t1]) - sin(phi[t1])*sin(chi[t1])
            rotmat1[0][1] = -cos(phi[t1])*cos(theta[t1])*sin(chi[t1]) - sin(phi[t1])*cos(chi[t1])
            rotmat1[0][2] = cos(phi[t1])*sin(theta[t1])
            rotmat1[1][0] = sin(phi[t1])*cos(theta[t1])*cos(chi[t1]) + cos(phi[t1])*sin(chi[t1])
            rotmat1[1][1] = -sin(phi[t1])*cos(theta[t1])*sin(chi[t1]) + cos(phi[t1])*cos(chi[t1])
            rotmat1[1][2] = sin(phi[t1])*sin(theta[t1])
            rotmat1[2][0] = -sin(theta[t1])*cos(chi[t1])
            rotmat1[2][1] = sin(theta[t1])*sin(chi[t1])
            rotmat1[2][2] = cos(theta[t1])

            rotmat2[0][0] = cos(phi[t2])*cos(theta[t2])*cos(chi[t2]) - sin(phi[t2])*sin(chi[t2])
            rotmat2[0][1] = -cos(phi[t2])*cos(theta[t2])*sin(chi[t2]) - sin(phi[t2])*cos(chi[t2])
            rotmat2[0][2] = cos(phi[t2])*sin(theta[t2])
            rotmat2[1][0] = sin(phi[t2])*cos(theta[t2])*cos(chi[t2]) + cos(phi[t2])*sin(chi[t2])
            rotmat2[1][1] = -sin(phi[t2])*cos(theta[t2])*sin(chi[t2]) + cos(phi[t2])*cos(chi[t2])
            rotmat2[1][2] = sin(phi[t2])*sin(theta[t2])
            rotmat2[2][0] = -sin(theta[t2])*cos(chi[t2])
            rotmat2[2][1] = sin(theta[t2])*sin(chi[t2])
            rotmat2[2][2] = cos(theta[t2])

            for a in range(3):
                for b in range(3):
                    for c in range(3):
                        rotmar[a][b]=rotmar[a][b]+rotmat1[c][a]*rotmat2[c][b]

            cost = rotmar[2][2]
            if (cost > 1.0):
                cost = 1.0
            elif (cost < -1.0):
                cost = -1.0

            thetar = N.arccos(cost)
            sint = N.sin(thetar)

            if (N.abs(1.0-cost) < small):
                phir = 0.0
                cchi = rotmar[0][0]
                schi = rotmar[1][0]

                if (cchi > 1.0):
                    cchi = 1.0
                elif (cchi < -1.0):
                    cchi = -1.0

                if (schi > 1.0):
                    schi = 1.0
                elif (schi < -1.0):
                    schi = -1.0

                if (schi > 0.0):
                    chir=N.arccos(cchi)
                else:
                    chir=2.0*N.pi-N.arccos(cchi)

            elif (N.abs(1.0+cost) < small):

                phir=0.0
                cchi=rotmar[1][1]
                schi=rotmar[0][1]

                if (cchi > 1.0):
                    cchi = 1.0
                elif (cchi < -1.0):
                    cchi = -1.0

                if (schi > 1.0):
                    schi = 1.0
                elif (schi < -1.0):
                    schi = -1.0

                if (schi > 0.0):
                    chir=N.arccos(cchi)
                else:
                    chir=2.0*N.pi-N.arccos(cchi)
            else:
                cphi = rotmar[0][2]/sint
                sphi = rotmar[1][2]/sint
                cchi = -rotmar[2][0]/sint
                schi = rotmar[2][1]/sint

                if (cphi > 1.0):
                    cphi = 1.0
                elif (cphi < -1.0):
                    cphi = -1.0

                if (sphi > 1.0):
                    sphi = 1.0
                elif (sphi < -1.0):
                    sphi = -1.0

                if (cchi > 1.0):
                    cchi = 1.0
                elif (cchi < -1.0):
                    cchi = -1.0

                if (schi > 1.0):
                    schi = 1.0
                elif (schi < -1.0):
                    schi = -1.0

                if (sphi > 0.0):
                    phir=N.arccos(cphi)
                else:
                    phir=2.0*N.pi-N.arccos(cphi)
                if (schi > 0.0):
                    chir=N.arccos(cchi)
                else:
                    chir=2.0*N.pi-N.arccos(cchi)

            if (phir<0.0):
                phir = 2.0*N.pi+phir
            if (chir<0.0):
                chir = 2.0*N.pi+chir

            phir = phir % (2.0*N.pi)
            chir = chir % (2.0*N.pi)

            if (cost > 1.0):
                cost = 2.0 - cost
            elif (cost < -1.0):
                cost = -2.0 - cost

            thetar = N.arccos(cost)

            phir *= 180.0/N.pi
            thetar *= 180.0/N.pi
            chir *= 180.0/N.pi

            indchi = int(chir)
            indphi = int(phir)
            indthe = int(thetar)

            index_p0=(indthe*361+indphi)*361+indchi # index for point 0
            rho_0=rhomat[index_p0]
            erot_0=erotmat[index_p0]
            esq_0=esqmat[index_p0]

            if indphi == 360:
                indphi == 0

            if indthe == 180:
                indthe == 0

            if indchi == 360:
                indchi == 0

            if (indchi > 360 or indchi < 0):
                print "index of chi is out of range", indchi, chi
                raise()
            elif (indphi > 360 or indphi < 0):
                print "index of phi is out of range", indphi, phi
                raise()
            elif (indthe > 180 or indthe < 0):
                print "index of theta is out of range", indthe, theta
                raise()

            if (indchi != 360):
                index_p1 = (indthe*361+indphi)*361+indchi+1 # next index is just the difference by 1 degree
                # CHI part for each properties
                delchi = rhomat[index_p1]-rho_0
                delch2 = erotmat[index_p1]-erot_0
                delch3 = esqmat[index_p1]-esq_0

            if (indphi != 360):
                index_p1 = (indthe*361+indphi+1)*361+indchi
                # PHI part for each properties
                delphi = rhomat[index_p1]-rho_0
                delph2 = erotmat[index_p1]-erot_0
                delph3 = esqmat[index_p1]-esq_0

            if (indthe != 180):
                index_p1 = ((indthe+1)*361+indphi)*361+indchi
                # THETA part for each properties
                delthe = rhomat[index_p1]-rho_0
                delth2 = erotmat[index_p1]-erot_0
                delth3 = esqmat[index_p1]-esq_0

            # Linear Interpolation for all three angles
            rho_new2 = rho_0+delchi*(chir-N.float(indchi))+delphi*(phir-N.float(indphi))+delthe*(thetar-N.float(indthe))
            erot_new2 = erot_0+delch2*(chir-N.float(indchi))+delph2*(phir-N.float(indphi))+delth2*(thetar-N.float(indthe))
            esq_new2 = esq_0+delch3*(chir-N.float(indchi))+delph3*(phir-N.float(indphi))+delth3*(thetar-N.float(indthe))


            # concatenate data

            rho_new = rho_new1*rho_new2
            erot_new = (erot_new1 + erot_new2)/2.0


            # sampling time!

            rd = (rho_new/rho_old)*exp((-(beta/P)*(pot_new - pot_old)))

            k = (i*P + j)*Nmol + m
            start_at = 10*Nmol*P

            if (rd > N.random.random()):

                acceptance_ratio += 1.0/(nsteps*Nmol*P)

                if k == start_at: # start averaging after some time steps (concept known as equilibration)
                    temp_val = pot_new
                    temp_val2 = erot_new
                    average_pot = temp_val
                    average_erot = temp_val2

                elif k > start_at:
                    new_val = pot_new
                    new_val2 = erot_new
                    weight1 = float(k-start_at)/float(k-start_at+1)
                    weight2 = float(1)/float(k-start_at+1)
                    temp_val = average_pot*weight1 + new_val*weight2
                    temp_val2 = average_erot*weight1 + new_val2*weight2
                    average_pot = temp_val
                    average_erot = temp_val2

            else:  # if we don't accept we must keep our original values

                if k == start_at:
                    temp_val = pot_old
                    temp_val2 = erot_old
                    average_pot = temp_val
                    average_erot = temp_val2

                elif k > start_at:
                    new_val = pot_old
                    new_val2 = erot_old
                    weight1 = float(k-start_at)/float(k-start_at+1)
                    weight2 = float(1)/float(k-start_at+1)
                    temp_val = average_pot*weight1 + new_val*weight2
                    temp_val2 = average_erot*weight1 + new_val2*weight2
                    average_pot = temp_val
                    average_erot = temp_val2

                for z in range(Nmol*Natom):

                    rcm[z*P + j][0] = rcm_old[z*P + j][0]
                    rcm[z*P + j][1] = rcm_old[z*P + j][1]
                    rcm[z*P + j][2] = rcm_old[z*P + j][2]

                phi[t1] = phi_old
                costheta[t1] = costheta_old
                theta[t1] = theta_old
                chi[t1] = chi_old


print acceptance_ratio, "acceptance_ratio"
print average_pot, "average_pot"
print average_erot, "average_erot"

potfile.write(str(average_pot)+"\n")
erotfile.write(str(average_erot)+"\n")
acceptanceratiofile.write(str(acceptance_ratio)+"\n")
