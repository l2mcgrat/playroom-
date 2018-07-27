import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from sys import argv
import numpy as np
from numpy import zeros,sqrt,mean

positions = argv[1]
numberofcharges = int(argv[2])
numberofdimensions = 2

x = zeros( numberofcharges , float )
y = zeros( numberofcharges , float )

positions = open(positions, "r")

for i in range(numberofcharges):
    for j in range(numberofdimensions):
        if j == 0:
            x[i] = float(positions.readline())
        if j == 1:
            y[i] = float(positions.readline())

positions.close()

plt.plot(x, y, 'ro')
plt.xlabel('x(nm)')
plt.ylabel('y(nm)')
plt.title('field of charges')
plt.show()
