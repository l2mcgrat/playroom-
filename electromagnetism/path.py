import matplotlib.pyplot as plt
import numpy as np
from numpy import zeros, exp, pi, sqrt
from sys import argv

# argv[1] is the file containing the positions of all the charges

# its called GEOMETRYfieldfile

# the x values are all stored in the odd rows of the file

# the y values are all stored in the even rows of the file

positions = argv[1]

# importing number of charges for looping purposes

numberofcharges = int(positions[positions.find("d")+1:positions.find("c")])

# initializing arrays that store positions

xcharges = np.zeros( numberofcharges , float )
ycharges = np.zeros( numberofcharges , float )

# opening file (read only)

positions = open(positions, "r")

for i in range(numberofcharges):

    position = float(positions.readline())
    xcharges[i] = position

    position = float(positions.readline())
    ycharges[i] = position

positions.close()

# obtain values of charges

stationairycharges = argv[2]
movingcharge = 1.6e-19

# store values of charges

charge = np.zeros( numberofcharges , float )

stationairycharges = open(stationairycharges, "r")

for i in range(numberofcharges):

    stationairycharge = float(stationairycharges.readline())
    charge[i] = stationairycharge

stationairycharges.close()

# import number of simulation steps

numberofsteps = 10

# import time step

dt = 1e-12

# import mass of charges

m = 9.11e-31

# initialize positions of charge(s) array

xposition = np.zeros( numberofsteps , float )
yposition = np.zeros( numberofsteps , float )

# obtain initial position of particles

position[0] = 9.0
yposition[0] = 7.0

# initialize velocities of charge(s) array

xvelocity = np.zeros( numberofsteps , float )
yvelocity = np.zeros( numberofsteps , float )

# obtain initial velocities of particles

xvelocity[0] = 0.0
yvelocity[0] = 0.0

############## SIMULATION DETAILS ###############

# now we are ready to simulate the motion of the charges

# we have the initial conditions in the form of the position(s) and velocitie(s) of the particle

# we also have a set of charges that determine our potential energy U(x,y)

# we will first update velocites using knowledge of the gradient of the potential

# we will then update the positions using knowledge of the velocity

################## SIMULATION ###################

for i in range((numberofsteps-1)):
    for j in range(numberofcharges):

        xvelocity[i+1] += (-1.0)*xposition[i]*( dt / m )*movingcharge*charge[j]/((4.0*pi*8.854187e-12)*( ( sqrt( ( xposition[i] - xcharges[j] )**2 + ( yposition[i] - ycharges[j] )**2 ) )**3 ))*1e9
        yvelocity[i+1] += (-1.0)*yposition[i]*( dt / m )*movingcharge*charge[j]/((4.0*pi*8.854187e-12)*( ( sqrt( ( xposition[i] - xcharges[j] )**2 + ( yposition[i] - ycharges[j] )**2 ) )**3 ))*1e9

#    print xvelocity[i+1], yvelocity[i+1]

    xposition[i+1] = xposition[i] + xvelocity[i+1]*dt
    yposition[i+1] = yposition[i] + yvelocity[i+1]*dt

    print xposition[i+1], yposition[i+1]

# lol, yup thats actually it

############### END SIMULATION #################

# plotting positions of charges

plt.plot(xcharges, ycharges, 'r-')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Positions of Charges')
plt.show()

# plotting position over time

plt.plot(xposition, yposition, 'r-')
plt.xlabel('x (nm)')
plt.ylabel('y (nm)')
plt.title('Charges Path')
plt.show()
