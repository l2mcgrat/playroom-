import numpy as N
from numpy import zeros, exp                                                                        
from sys import argv

DataFile = N.loadtxt(argv[1])

DataSize = int(argv[2])

DataValue = zeros( (DataSize, 9) , float )
Symbol = zeros( (DataSize, 9) , str )

for i in range(DataSize):
    for j in range(9):
        DataValue[i][j] = DataFile[i,j]
        if DataValue[i][j] == 1.0:
            Symbol[i][j] = 'X'
        if DataValue[i][j] == 0.0:
            Symbol[i][j] = 'O'
        if DataValue[i][j] == -1.0:
            Symbol[i][j] = ' '

for i in range(DataSize):

    print ' GAME',(i+1)
    print ''
    print '###################'
    print '##',Symbol[i][0],'###',Symbol[i][1],'###',Symbol[i][2],'##'
    print '###################'
    print '##',Symbol[i][3],'###',Symbol[i][4],'###',Symbol[i][5],'##'
    print '###################'
    print '##',Symbol[i][6],'###',Symbol[i][7],'###',Symbol[i][8],'##'
    print '###################'
    print ''
