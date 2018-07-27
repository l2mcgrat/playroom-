from sys import argv, exit
import numpy as np
from numpy import pi, sin, cos, sqrt

numberofsurfacepoints = int(argv[1])  # must have an integer sqrt

# surfacefieldgeneratior

surfaceradius = float(argv[2])        # in nm

R = open("R"+str(surfaceradius)+"nm"+str(numberofsurfacepoints+2)+"SURFACEPOINTS","w")

for i in range(int(sqrt(numberofsurfacepoints))):
    for k in range(int(sqrt(numberofsurfacepoints))):
        for j in range(3):
            phi = i*2*pi/float(sqrt(numberofsurfacepoints))
            theta = (k + 1)*pi/float(((sqrt(numberofsurfacepoints)) + 1))

            if j == 0:
                surfaceposition = surfaceradius*sin(theta)*cos(phi)

            if j == 1:
                surfaceposition = surfaceradius*sin(theta)*sin(phi)

            if j == 2:
                surfaceposition = surfaceradius*cos(theta)

            R.write(str(surfaceposition)+"\n")

# we included all but the north and south polls with the surface file thus far

for i in range(2):
    for j in range(3):
        if j == 0:
            surfaceposition = 0.0

        if j == 1:
            surfaceposition = 0.0

        if j == 2:
            surfaceposition = ((-1.0)**(i+1))*surfaceradius

        R.write(str(surfaceposition)+"\n")

R.close()

