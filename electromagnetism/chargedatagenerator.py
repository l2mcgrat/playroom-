# Charges making a ring

from sys import argv, exit
import numpy as np
from numpy import pi, sin, cos

numberofcharges = int(argv[1])     # many typically (100,1000)

# ringfieldgeneratior

ringradius = float(argv[2])        # in nm

ring = open("RING"+str(ringradius)+"nmfield"+str(numberofcharges)+"charges","w")

for i in range(numberofcharges):
    for j in range(2):

        theta = i*2*np.pi/float(numberofcharges)

        if j == 0:
            ringchargeposition = ringradius*cos(theta)

        if j == 1:
            ringchargeposition = ringradius*sin(theta)

        ring.write(str(ringchargeposition)+"\n")

ring.close()

length = float(argv[3]) # in nm

line = open("LINE"+str(length)+"nmfield"+str(numberofcharges)+"charges","w")

for i in range(numberofcharges):
    for j in range(2):

        if (j == 0):
            linechargeposition = 0.0

        if (j == 1):
            linechargeposition = length*i/(numberofcharges - 1)

        line.write(str(linechargeposition)+"\n")

line.close()

chargevalue = 1.6e-19

charges = open("CHARGES"+str(numberofcharges)+"file","w")

for i in range(numberofcharges):

    charges.write(str(chargevalue)+"\n")

charges.close()
