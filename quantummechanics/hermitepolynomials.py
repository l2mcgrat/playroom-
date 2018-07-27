import numpy as np
from sys import argv

def hermitepolynomial( quantum_number , x ):

    hermite = np.zeros( 11 , float )

    hermite[0] = 1.0
    hermite[1] = 2.0*x
    hermite[2] = 4.0*(x**2) - 2.0
    hermite[3] = 8.0*(x**3) - 12.0*x
    hermite[4] = 16.0*(x**4) - 48.0*(x**2) + 12.0
    hermite[5] = 32.0*(x**5) - 160.0*(x**3) + 120.0*x
    hermite[6] = 64.0*(x**6) - 480.0*(x**4) + 720.0*(x**2) - 120.0
    hermite[7] = 128.0*(x**7) - 1344.0*(x**5) + 3360.0*(x**3) - 1680.0
    hermite[8] = 256.0*(x**8) - 3584.0*(x**6) + 13440.0*(x**4) - 13440.0*(x**2) + 1680.0
    hermite[9] = 512.0*(x**9) - 9216.0*(x**7) + 48384.0*(x**5) - 80640.0*(x**3) + 30240.0*x
    hermite[10] = 1024.0*(x**10) - 23040.0*(x**8) + 161280.0*(x**6) - 403200.0*(x**4) + 302400.0*(x**2)

    return hermite[quantum_number]


def innerproduct( quantum_number1, quantum_number2, x):

    product = hermitepolynomial( quantum_number1, x )*hermitepolynomial( quantum_number2, x )

    return product 


def approximateintegral( quantum_number1, quantum_number2 ):

    integral = 0.0
    howslowdoyouwantit = int(argv[1])
    bound = int(argv[2]) 
    xstuff = np.linspace( (-1.0*bound), bound, howslowdoyouwantit )
    xvalues = np.zeros( (howslowdoyouwantit - 1) , float )
    yvalues = np.zeros( (howslowdoyouwantit - 1) , float )
 
    for i in range(howslowdoyouwantit - 1):
        xvalues[i] = (xstuff[i] + xstuff[i+1])/2

    for i in range(howslowdoyouwantit - 1):
        yvalues[i] = innerproduct( quantum_number1, quantum_number2, xvalues[i])
        integral += yvalues[i]*xvalues[i]

    return integral


def lolsomethingswrong():

    matrix = np.zeros( (11,11) , float )

    for i in range(11):
        for j in range(11):
            matrix[i][j] = approximateintegral(i, j)
            if matrix[i][j] < 1e-5:
                matrix[i][j] = 0.0
            else:
                matrix[i][j] = 1.0

    print matrix 

    return

lolsomethingswrong()
