import numpy as N
from numpy import zeros, exp                                                                        
from sys import argv

####################################################
######### IMPORTING GAME DATA AND LABELS ###########
####################################################

DataSize = int(argv[2])                                      #  Number of Examples 

DataSets = zeros( (DataSize, 9) , float )
Labels = zeros( (DataSize, 3) , float )
DataFile  = N.loadtxt(argv[1])                     #  Importing Data

for i in range(DataSize):                          #  Loading in Training Data
    for j in range(9):
        DataSets[i][j] = DataFile[i][j] 

for i in range(DataSize):                          #  Loading in Labels
    for j in range(3):
        Labels[i][j] = DataFile[i+DataSize][j]     

####################################################
############### TRAINING ALGORITHM #################
####################################################

NeuralNet = zeros( (3,9) , float)                  #  Initializing our small as fuck Neural Net
reps = int(argv[3])
DecentRate = 0.001  # delta wij                    #  Magnitude of updates to weights 

for i in range(3):                                 #  Initialization of NeuralNet using random number between -0.5 and 0.5
    for j in range(9):
        NeuralNet[i][j] = (N.random.random() - 0.5) 

for l in range(reps):                              #  Repeats training examples
    for n in range(DataSize):                      #  Loops through Training Data

        Output = zeros( 3 , float )

        for i in range(3):                         #  Feeding data through the Neural Network
            for j in range(9):
                Output[i] += NeuralNet[i][j]*DataSets[n][j]

        for i in range(3):                         #  Gradient Decent
            for j in range(9):
                NeuralNet[i][j] -= 2*(Output[i] - Labels[n][i])*DataSets[n][j]*DecentRate

####################################################
############### TESTING NEURAL NET #################
####################################################

Assessment_File = open("AssessmentFile_N"+str(DataSize),"w")   

Result =  'idk', 0 
Answer =  'idk', 0 
Score = 0.0

for d in range(DataSize):                          #  Tests all examples 
    Output = zeros( 3, float )
    for i in range(3):                             #  Feeding data through the Neural Network
        for j in range(9):
            Output[i] += NeuralNet[i][j]*DataSets[d][j]

# print NeuralNet

    ghjk = -1000.0

    if Output[0] > ghjk:                          #  Reads the output from NeuralNet and decides on its answer
        Answer = 'Win', 1
        ghjk = Output[0]

    if Output[1] > ghjk:
        Answer = 'Tie', 2
        ghjk = Output[1]

    if Output[2] > ghjk:
        Answer = 'Loss', 3

    asdf = 0.0
    
    if Labels[d][0] > asdf:                       #  Reads each label to see what the result of the game was
        Result = 'Win', 1
        asdf = Labels[d][0]

    if Labels[d][1] > asdf:
        Result = 'Tie', 2
        asdf = Labels[d][1]

    if Labels[d][2] > asdf:
        Result = 'Loss', 3

    Assessment_File.write("Neural Net thinks game "+str(d+1)+" ended in a "+str(Answer[0])+". "+"In reality, game "+str(d+1)+" ended in a "+str(Result[0])+"\n")

    if Answer[1] == Result[1]:                    #  Tallies up how many times it got it right
        Score += 1.0/float(DataSize)

Score *= 100.0
reaction = ''

if Score <= 50.0:                                 #  Adds a reaction comment based on the score :P 
    reaction = ' ! :o'

if (Score >= 50.0) and (Score <= 80.0):
    reaction = ' . :|'

if (Score >= 80.0) and (Score <= 94.0):
    reaction = ' ! :)'

if (Score >= 94.0):
    reaction = ' !!! :D'

Assessment_File.write("NeuralNet scored "+str(Score)+"%"+str(reaction))
Assessment_File.close()
