from numpy import zeros, exp 
from sys import argv

#  NEURAL NET PLAYROOM 1

#  Here I will attempt to create a Neural Net for Interpreting Numbers using inputs as pixels on a screen (6 x 3). 

#  It will take in the pixels as Neurons and use one 10x18 matrix to create a set of 10 output Neurons that each correspond 
#  to the likelyhood that its one of the digits. 

#  The pixels will be either 1 or 0 (blank) and will look like as follow: 

#   1 1 1      1        1 1    1 1 1    1   1    1 1 1    1 1 1    1 1 1    1 1 1    1 1 1
#   1   1    1 1      1   1        1    1   1    1        1        1   1    1   1    1   1
#   1   1      1          1    1 1 1    1 1 1    1 1 1    1 1 1    1   1    1 1 1    1 1 1
#   1   1      1        1          1        1        1    1   1        1    1   1        1
#   1   1      1      1            1        1        1    1   1        1    1   1        1
#   1 1 1    1 1 1    1 1 1    1 1 1        1    1 1 1    1 1 1        1    1 1 1        1 

#  Each of the "pixels" will be stored in 6x3 matricies which I will now initialize

Zero = zeros( (6,3) , float )
One = zeros( (6,3) , float )
Two = zeros( (6,3) , float )
Three = zeros( (6,3) , float )
Four = zeros( (6,3) , float )
Five = zeros( (6,3) , float )
Six = zeros( (6,3) , float )
Seven = zeros( (6,3) , float )
Eight = zeros( (6,3) , float )
Nine = zeros( (6,3) , float )

#  Filling Zero

for i in range(6):
    Zero[i][0] = 1.0
    Zero[i][2] = 1.0

Zero[0][1] = 1.0
Zero[5][1] = 1.0

#  Filling One

for i in range(6):
    One[i][1] = 1.0

One[1][0] = 1.0
One[5][0] = 1.0
One[5][2] = 1.0

#  Filling Two

for i in range(3):
    Two[i][2] = 1.0
    Two[5][i] = 1.0

Two[1][0] = 1.0
Two[4][0] = 1.0
Two[0][1] = 1.0
Two[3][1] = 1.0

#  Filling Three

for i in range(6):
    Three[i][2] = 1.0

for i in range(3):
    Three[0][i] = 1.0
    Three[2][i] = 1.0
    Three[5][i] = 1.0

#  Filling Four

for i in range(6):
    Four[i][2] = 1.0

for i in range(3):
    Four[i][0] = 1.0
    Four[2][i] = 1.0

#  Filling Five

for i in range(3):
    Five[0][i] = 1.0
    Five[2][i] = 1.0
    Five[5][i] = 1.0
    Five[i][0] = 1.0

Five[3][2] = 1.0
Five[4][2] = 1.0

#  Filling Six

for i in range(6):
    Six[i][0] = 1.0

for i in range(3):
    Six[0][i] = 1.0
    Six[2][i] = 1.0
    Six[5][i] = 1.0

Six[3][2] = 1.0
Six[4][2] = 1.0

#  Filling Seven

for i in range(6):
    Seven[i][2] = 1.0

for i in range(3):
    Seven[0][i] = 1.0
    Seven[i][0] = 1.0

#  Filling Eight

for i in range(6):
    Eight[i][0] = 1.0
    Eight[i][2] = 1.0

for i in range(3):
    Eight[0][i] = 1.0
    Eight[2][i] = 1.0
    Eight[5][i] = 1.0

#  Filling Nine

for i in range(6):
    Nine[i][2] = 1.0

for i in range(3):
    Nine[i][0] = 1.0
    Nine[0][i] = 1.0
    Nine[2][i] = 1.0

# Now I will set up print statements to check if the values were inputted as desired 

#print ' '
#print ' '

# CHECK Zero

#print  Zero[0][0], Zero[0][1], Zero[0][2]
#print  Zero[1][0], Zero[1][1], Zero[1][2]
#print  Zero[2][0], Zero[2][1], Zero[2][2]
#print  Zero[3][0], Zero[3][1], Zero[3][2]
#print  Zero[4][0], Zero[4][1], Zero[4][2]
#print  Zero[5][0], Zero[5][1], Zero[5][2]

#print ' '
#print ' '

# CHECK One

#print  One[0][0], One[0][1], One[0][2]
#print  One[1][0], One[1][1], One[1][2]
#print  One[2][0], One[2][1], One[2][2]
#print  One[3][0], One[3][1], One[3][2]
#print  One[4][0], One[4][1], One[4][2]
#print  One[5][0], One[5][1], One[5][2]

#print ' '
#print ' '

# CHECK Two

#print  Two[0][0], Two[0][1], Two[0][2]
#print  Two[1][0], Two[1][1], Two[1][2]
#print  Two[2][0], Two[2][1], Two[2][2]
#print  Two[3][0], Two[3][1], Two[3][2]
#print  Two[4][0], Two[4][1], Two[4][2]
#print  Two[5][0], Two[5][1], Two[5][2]

#print ' '
#print ' '

# CHECK Three

#print  Three[0][0], Three[0][1], Three[0][2]
#print  Three[1][0], Three[1][1], Three[1][2]
#print  Three[2][0], Three[2][1], Three[2][2]
#print  Three[3][0], Three[3][1], Three[3][2]
#print  Three[4][0], Three[4][1], Three[4][2]
#print  Three[5][0], Three[5][1], Three[5][2]

#print ' '
#print ' '

# CHECK Four

#print  Four[0][0], Four[0][1], Four[0][2]
#print  Four[1][0], Four[1][1], Four[1][2]
#print  Four[2][0], Four[2][1], Four[2][2]
#print  Four[3][0], Four[3][1], Four[3][2]
#print  Four[4][0], Four[4][1], Four[4][2]
#print  Four[5][0], Four[5][1], Four[5][2]

#print ' '
#print ' '

# CHECK Five

#print  Five[0][0], Five[0][1], Five[0][2]  
#print  Five[1][0], Five[1][1], Five[1][2]
#print  Five[2][0], Five[2][1], Five[2][2]
#print  Five[3][0], Five[3][1], Five[3][2]
#print  Five[4][0], Five[4][1], Five[4][2]
#print  Five[5][0], Five[5][1], Five[5][2]

#print ' '
#print ' '

# CHECK Six

#print  Six[0][0], Six[0][1], Six[0][2]
#print  Six[1][0], Six[1][1], Six[1][2]
#print  Six[2][0], Six[2][1], Six[2][2]
#print  Six[3][0], Six[3][1], Six[3][2]
#print  Six[4][0], Six[4][1], Six[4][2]
#print  Six[5][0], Six[5][1], Six[5][2]

#print ' '
#print ' '

# CHECK Seven

#print  Seven[0][0], Seven[0][1], Seven[0][2]  
#print  Seven[1][0], Seven[1][1], Seven[1][2]
#print  Seven[2][0], Seven[2][1], Seven[2][2]
#print  Seven[3][0], Seven[3][1], Seven[3][2]
#print  Seven[4][0], Seven[4][1], Seven[4][2]
#print  Seven[5][0], Seven[5][1], Seven[5][2]

#print ' '
#print ' '

# CHECK Eight

#print  Eight[0][0], Eight[0][1], Eight[0][2]  
#print  Eight[1][0], Eight[1][1], Eight[1][2]
#print  Eight[2][0], Eight[2][1], Eight[2][2]
#print  Eight[3][0], Eight[3][1], Eight[3][2]
#print  Eight[4][0], Eight[4][1], Eight[4][2]
#print  Eight[5][0], Eight[5][1], Eight[5][2]

#print ' '
#print ' '

# CHECK Nine 

#print  Nine[0][0], Nine[0][1], Nine[0][2]  
#print  Nine[1][0], Nine[1][1], Nine[1][2]
#print  Nine[2][0], Nine[2][1], Nine[2][2]
#print  Nine[3][0], Nine[3][1], Nine[3][2]
#print  Nine[4][0], Nine[4][1], Nine[4][2]
#print  Nine[5][0], Nine[5][1], Nine[5][2]

#print ' '
#print ' '

##############################################

# Set up numebrs for Primative Neural Network

Number = zeros( (10,18), float )
    
for i in range(3):
    Number[0][i] = Zero[0][i]
    Number[0][i+3] = Zero[1][i]
    Number[0][i+6] = Zero[2][i]
    Number[0][i+9] = Zero[3][i]
    Number[0][i+12] = Zero[4][i]
    Number[0][i+15] = Zero[5][i]

#print Number0

for i in range(3):
    Number[1][i] = One[0][i]
    Number[1][i+3] = One[1][i]
    Number[1][i+6] = One[2][i]
    Number[1][i+9] = One[3][i]
    Number[1][i+12] = One[4][i]
    Number[1][i+15] = One[5][i]

#print Number1

for i in range(3):
    Number[2][i] = Two[0][i]
    Number[2][i+3] = Two[1][i]
    Number[2][i+6] = Two[2][i]
    Number[2][i+9] = Two[3][i]
    Number[2][i+12] = Two[4][i]
    Number[2][i+15] = Two[5][i]

#print Number2

for i in range(3):
    Number[3][i] = Three[0][i]
    Number[3][i+3] = Three[1][i]
    Number[3][i+6] = Three[2][i]
    Number[3][i+9] = Three[3][i]
    Number[3][i+12] = Three[4][i]
    Number[3][i+15] = Three[5][i]

#print Number3

for i in range(3):
    Number[4][i] = Four[0][i]
    Number[4][i+3] = Four[1][i]
    Number[4][i+6] = Four[2][i]
    Number[4][i+9] = Four[3][i]
    Number[4][i+12] = Four[4][i]
    Number[4][i+15] = Four[5][i]

#print Number4

for i in range(3):
    Number[5][i] = Five[0][i]
    Number[5][i+3] = Five[1][i]
    Number[5][i+6] = Five[2][i]
    Number[5][i+9] = Five[3][i]
    Number[5][i+12] = Five[4][i]
    Number[5][i+15] = Five[5][i]

#print Number5

for i in range(3):
    Number[6][i] = Six[0][i]
    Number[6][i+3] = Six[1][i]
    Number[6][i+6] = Six[2][i]
    Number[6][i+9] = Six[3][i]
    Number[6][i+12] = Six[4][i]
    Number[6][i+15] = Six[5][i]

#print Number6

for i in range(3):
    Number[7][i] = Seven[0][i]
    Number[7][i+3] = Seven[1][i]
    Number[7][i+6] = Seven[2][i]
    Number[7][i+9] = Seven[3][i]
    Number[7][i+12] = Seven[4][i]
    Number[7][i+15] = Seven[5][i]

#print Number7

for i in range(3):
    Number[8][i] = Eight[0][i]
    Number[8][i+3] = Eight[1][i]
    Number[8][i+6] = Eight[2][i]
    Number[8][i+9] = Eight[3][i]
    Number[8][i+12] = Eight[4][i]
    Number[8][i+15] = Eight[5][i]

#print Number8

for i in range(3):
    Number[9][i] = Nine[0][i]
    Number[9][i+3] = Nine[1][i]
    Number[9][i+6] = Nine[2][i]
    Number[9][i+9] = Nine[3][i]
    Number[9][i+12] = Nine[4][i]
    Number[9][i+15] = Nine[5][i]

#print Number9


############################################################
############################################################

True_Value = zeros( (10,10) , float )

for i in range(10):
    for j in range(10):
        if (i == j):
            True_Value[i][j] = 1.0

############################################################
############################################################

Brain = zeros( (10,18), float )
Output = zeros( 10, float )
Input = zeros( 18, float )
Sigma_Output = zeros( 10, float )
Answer = 0
Decent_Rate = 0.001

#######################################
######  Iteration 2: Primitive  #######
#######################################

for i in range(10):
    for j in range(18):
        Brain[i][j] = 1.0 / 18.0

#####################
#### Train Brain ####
#####################

# Number of repeated training sessions == reps 

reps = 1000

for r in range(reps):
    for k in range(10):
        for i in range(10):
            for j in range(18):
                Output[i] += Brain[i][j]*Number[k][j] 

        for i in range(10):
            Sigma_Output[i] = 1 / ( 1 + exp(-Output[i]))

        for i in range(10):
            for j in range(18):
                Brain[i][j] -= 2*(Output[i] - True_Value[k][i])*(exp(-Output[i])/((1 + exp(-Output[i]))**2))*Number[k][j]*Decent_Rate

        offset = zeros( 10 , float ) 

        for i in range(10):
            for j in range(18):
                if offset[i] > Brain[i][j]:
                    offset[i] = Brain[i][j] 
           
        for i in range(10):   
            offset[i] = (-1)*offset[i] 
            for j in range(18):
                Brain[i][j] += offset[i]

        Z = zeros( 10 , float )

        for i in range(10):
            for j in range(18):
                Z[i] += Brain[i][j] 

        for i in range(10):
            for j in range(18):
                Brain[i][j] /= Z[i]

        for i in range(10):
            Output[i] = 0.0

# print Brain

#####################
## Verifying Brain ##
#####################

Best_Value = -0.01
Correct_Answer = int(argv[1])

for i in range(10):
    for j in range(18):
        Output[i] += Brain[i][j]*Number[Correct_Answer][j]

print Output

for i in range(10):
    if (Output[i] > Best_Value):
        Best_Value = Output[i]
        Answer = i
    else: 
        None

print 'The Brain Guesses the number is: ', Answer 
print 'The Real Number is:', Correct_Answer 
