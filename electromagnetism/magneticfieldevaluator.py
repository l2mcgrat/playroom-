import matplotlib.pyplot as plt
import numpy as np
from numpy import zeros, exp, pi, sqrt
from sys import argv

positions = argv[1]

numberofcharges = int(positions[positions.find("d")+1:positions.find("c")])
length = float(positions[positions.find("E")+1:positions.find("n")])

xcharges = np.zeros( numberofcharges , float )
ycharges = np.zeros( numberofcharges , float )

# opening file (read only)

positions = open(positions, "r")

for i in range(numberofcharges):

    position = float(positions.readline())
    xcharges[i] = position*1.0e-9

    position = float(positions.readline())
    ycharges[i] = position*1.0e-9

positions.close()

stationairycharges = argv[2]

# store values of charges

charge = np.zeros( numberofcharges , float )

stationairycharges = open(stationairycharges, "r")

for i in range(numberofcharges):

    stationairycharge = float(stationairycharges.readline())
    charge[i] = stationairycharge

stationairycharges.close()

# initialize position of test charge

xposition = 1.0e-9  # nm
yposition = 5.0e-5  # nm

E_field = np.zeros( 2 , float  )

for i in range(numberofcharges):
    for j in range(2):

        if (j == 0):
            E_field[j] += charge[i]*(xposition - xcharges[i])/(4*pi*(8.85e-12)*(( (xposition - xcharges[i])**2.0 + (yposition - ycharges[i])**2.0)**(3.0/2.0)))

        if (j == 1):
            E_field[j] += charge[i]*(yposition - ycharges[i])/(4*pi*(8.85e-12)*(( (xposition - xcharges[i])**2.0 + (yposition - ycharges[i])**2.0)**(3.0/2.0)))

print "Ex = ", E_field[0], " N/C .", " Charge density = ", (1.6e-19)*(numberofcharges/length)*1.0e9, " C/m"
print "Ey = ", E_field[1], " N/C .", " Charge density = ", (1.6e-19)*(numberofcharges/length)*1.0e9, " C/m"

