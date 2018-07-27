# Charges making a ring

from sys import argv, exit
import numpy as np
from numpy import pi, sin, cos

numberofcharges = int(argv[1])     # many typically (100,1000)

locations = open("locations"+"of"+str(numberofcharges)+"pointcharges","w")

for i in range(numberofcharges):
    for j in range(3):

        if (j == 0):
            pointchargeposition = (np.random.random() - 0.5)*3.6e-9

        if (j == 1):
            pointchargeposition = (np.random.random() - 0.5)*3.6e-9

        if (j == 2):
            pointchargeposition = (np.random.random() - 0.5)*3.6e-9

        locations.write(str(pointchargeposition)+"\n")

locations.close()

chargevalue = 1.6e-19

charges = open("CHARGES"+str(numberofcharges)+"file","w")

for i in range(numberofcharges):

    charges.write(str(chargevalue)+"\n")

charges.close()
