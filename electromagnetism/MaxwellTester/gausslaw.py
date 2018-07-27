import matplotlib.pyplot as plt
import numpy as np
from numpy import zeros, exp, pi, sqrt
from sys import argv

positions = argv[1] # obtains positions of charges

numberofcharges = int(positions[positions.find("f")+1:positions.find("p")])

xcharges = np.zeros( numberofcharges , float )
ycharges = np.zeros( numberofcharges , float )
zcharges = np.zeros( numberofcharges , float )

# opening position file (read only)

positions = open(positions, "r")

for i in range(numberofcharges):

    position = float(positions.readline())
    xcharges[i] = position*1.0e-9

    position = float(positions.readline())
    ycharges[i] = position*1.0e-9

    position = float(positions.readline())
    zcharges[i] = position*1.0e-9

positions.close()

stationairycharges = argv[2]

# store values of charges

charge = np.zeros( numberofcharges , float )

stationairycharges = open(stationairycharges, "r")

for i in range(numberofcharges):

    stationairycharge = float(stationairycharges.readline())
    charge[i] = stationairycharge

stationairycharges.close()

surface = argv[3] # obtain file which contains the surface information

numberofsurfacepoints = int(surface[surface.find("m")+1:surface.find("S")]) # extracts number of points from the filename

radius = (float(surface[surface.find("R")+1:surface.find("n")]))*1.0e-9 # extracts the radius from the filename

surfaceelement = 4*pi*(radius**2)/numberofsurfacepoints # finds the value of dS

surfacepoint = np.zeros( (numberofsurfacepoints, 3) , float )

surface = open(surface, "r")

for i in range(numberofsurfacepoints):
    for j in range(3):

        if (j == 0):
            surfacelocation = float(surface.readline())
            surfacepoint[i][j] = surfacelocation*1.0e-9
        if (j == 1):
            surfacelocation = float(surface.readline())
            surfacepoint[i][j] = surfacelocation*1.0e-9
        if (j == 2):
            surfacelocation = float(surface.readline())
            surfacepoint[i][j] = surfacelocation*1.0e-9

surface.close()

# must find the normal at all points in order to evaluate the dot product

normal = np.zeros( (numberofsurfacepoints, 3) , float )

for i in range(numberofsurfacepoints):
    for j in range(3):
        normal[i][j] = surfacepoint[i][j]/( sqrt( ( surfacepoint[i][0] )**2.0 + ( surfacepoint[i][1] )**2.0 + ( surfacepoint[i][2] )**2.0 ) )

 # integration of the electricfield around a surface S.

integral = 0.0

electricfield = np.zeros( 3, float )

for i in range(numberofsurfacepoints):

    electricfield[0] = 0.0
    electricfield[1] = 0.0
    electricfield[2] = 0.0

    for k in range(numberofcharges):
        for j in range(3):

            electricfield[j] += charge[k]*(surfacepoint[i][j] - xcharges[k])/(4*pi*(8.85e-12)*( ( (surfacepoint[i][0] - xcharges[k])**2.0 + (surfacepoint[i][1] - ycharges[k])**2.0 + (surfacepoint[i][2] - zcharges[k])**2.0)**(3.0/2.0)))

    integral += ( electricfield[0]*normal[i][0] + electricfield[1]*normal[i][1] + electricfield[2]*normal[i][2]  )*surfaceelement

print integral, numberofcharges*(1.6e-19)/(8.85e-12)

