from numpy import zeros 
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

# Set up Most Primative Neural Network

Number0 = zeros( 18, float) 
Number1 = zeros( 18, float)
Number2 = zeros( 18, float)
Number3 = zeros( 18, float)
Number4 = zeros( 18, float)
Number5 = zeros( 18, float)
Number6 = zeros( 18, float)
Number7 = zeros( 18, float)
Number8 = zeros( 18, float)
Number9 = zeros( 18, float)

for i in range(3):
    Number0[i] = Zero[0][i]
    Number0[i+3] = Zero[1][i]
    Number0[i+6] = Zero[2][i]
    Number0[i+9] = Zero[3][i]
    Number0[i+12] = Zero[4][i]
    Number0[i+15] = Zero[5][i]

#print Number0

for i in range(3):
    Number1[i] = One[0][i]
    Number1[i+3] = One[1][i]
    Number1[i+6] = One[2][i]
    Number1[i+9] = One[3][i]
    Number1[i+12] = One[4][i]
    Number1[i+15] = One[5][i]

#print Number1

for i in range(3):
    Number2[i] = Two[0][i]
    Number2[i+3] = Two[1][i]
    Number2[i+6] = Two[2][i]
    Number2[i+9] = Two[3][i]
    Number2[i+12] = Two[4][i]
    Number2[i+15] = Two[5][i]

#print Number2

for i in range(3):
    Number3[i] = Three[0][i]
    Number3[i+3] = Three[1][i]
    Number3[i+6] = Three[2][i]
    Number3[i+9] = Three[3][i]
    Number3[i+12] = Three[4][i]
    Number3[i+15] = Three[5][i]

#print Number3

for i in range(3):
    Number4[i] = Four[0][i]
    Number4[i+3] = Four[1][i]
    Number4[i+6] = Four[2][i]
    Number4[i+9] = Four[3][i]
    Number4[i+12] = Four[4][i]
    Number4[i+15] = Four[5][i]

#print Number4

for i in range(3):
    Number5[i] = Five[0][i]
    Number5[i+3] = Five[1][i]
    Number5[i+6] = Five[2][i]
    Number5[i+9] = Five[3][i]
    Number5[i+12] = Five[4][i]
    Number5[i+15] = Five[5][i]

#print Number5

for i in range(3):
    Number6[i] = Six[0][i]
    Number6[i+3] = Six[1][i]
    Number6[i+6] = Six[2][i]
    Number6[i+9] = Six[3][i]
    Number6[i+12] = Six[4][i]
    Number6[i+15] = Six[5][i]

#print Number6

for i in range(3):
    Number7[i] = Seven[0][i]
    Number7[i+3] = Seven[1][i]
    Number7[i+6] = Seven[2][i]
    Number7[i+9] = Seven[3][i]
    Number7[i+12] = Seven[4][i]
    Number7[i+15] = Seven[5][i]

#print Number7

for i in range(3):
    Number8[i] = Eight[0][i]
    Number8[i+3] = Eight[1][i]
    Number8[i+6] = Eight[2][i]
    Number8[i+9] = Eight[3][i]
    Number8[i+12] = Eight[4][i]
    Number8[i+15] = Eight[5][i]

#print Number8

for i in range(3):
    Number9[i] = Nine[0][i]
    Number9[i+3] = Nine[1][i]
    Number9[i+6] = Nine[2][i]
    Number9[i+9] = Nine[3][i]
    Number9[i+12] = Nine[4][i]
    Number9[i+15] = Nine[5][i]

#print Number9


############################################################
############################################################

Output0 = zeros( 10, float ) 
Output1 = zeros( 10, float )
Output2 = zeros( 10, float )
Output3 = zeros( 10, float )
Output4 = zeros( 10, float )
Output5 = zeros( 10, float )
Output6 = zeros( 10, float )
Output7 = zeros( 10, float )
Output8 = zeros( 10, float )
Output9 = zeros( 10, float )

Output0[0] = 1.0
Output1[1] = 1.0
Output2[2] = 1.0
Output3[3] = 1.0
Output4[4] = 1.0
Output5[5] = 1.0
Output6[6] = 1.0
Output7[7] = 1.0
Output8[8] = 1.0
Output9[9] = 1.0

############################################################
############################################################

Brain = zeros( (10,18), float )
Output = zeros( 10, float )
Input = zeros( 18, float )
Answer = 0
Best_Value = -0.01
Decent_Rate = 0.01

#######################################
####  Iteration 1: Most Primitive  ####
#######################################

for i in range(10):
    for j in range(18):
        Brain[i][j] = 1.0 / 18.0

#####################
#### Train Brain ####
#####################

# Number of repeated training sessions == reps 

reps = 1

for n in range(reps):

    # Train for 0

    for i in range(10):
        for j in range(18):
            Output[i] += Brain[i][j]*Number0[j] 

    for i in range(10):
        for j in range(18):
            Brain[i][j] -= 2*(Output[i] - Output0[i])*Number0[j]*Decent_Rate

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

    print Brain

    # Train for 1

    for i in range(10):
        for j in range(18):
            Output[i] += Brain[i][j]*Number1[j] 

    for i in range(10):
        for j in range(18):
            Brain[i][j] -= 2*(Output[i] - Output1[i])*Number1[j]*Decent_Rate

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

    # Train for 2

    for i in range(10):
        for j in range(18):
            Output[i] += Brain[i][j]*Number2[j] 

    for i in range(10):
        for j in range(18):
            Brain[i][j] -= 2*(Output[i] - Output2[i])*Number2[j]*Decent_Rate

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

    # Train for 3

    for i in range(10):
        for j in range(18):
            Output[i] += Brain[i][j]*Number3[j] 

    for i in range(10):
        for j in range(18):
            Brain[i][j] -= 2*(Output[i] - Output3[i])*Number3[j]*Decent_Rate

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

    # Train for for

    for i in range(10):
        for j in range(18):
            Output[i] += Brain[i][j]*Number4[j] 

    for i in range(10):
        for j in range(18):
            Brain[i][j] -= 2*(Output[i] - Output4[i])*Number4[j]*Decent_Rate

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

    # Train for 5

    for i in range(10):
        for j in range(18):
            Output[i] += Brain[i][j]*Number5[j] 

    for i in range(10):
        for j in range(18):
            Brain[i][j] -= 2*(Output[i] - Output5[i])*Number5[j]*Decent_Rate

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

    # Train for 6

    for i in range(10):
        for j in range(18):
            Output[i] += Brain[i][j]*Number6[j] 

    for i in range(10):
        for j in range(18):
            Brain[i][j] -= 2*(Output[i] - Output6[i])*Number6[j]*Decent_Rate

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

    # Train for 7

    for i in range(10):
        for j in range(18):
            Output[i] += Brain[i][j]*Number7[j] 

    for i in range(10):
        for j in range(18):
            Brain[i][j] -= 2*(Output[i] - Output7[i])*Number7[j]*Decent_Rate

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

    # Train for 8

    for i in range(10):
        for j in range(18):
            Output[i] += Brain[i][j]*Number8[j] 

    for i in range(10):
        for j in range(18):
            Brain[i][j] -= 2*(Output[i] - Output8[i])*Number8[j]*Decent_Rate

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

    # Train for 9

    for i in range(10):
        for j in range(18):
            Output[i] += Brain[i][j]*Number9[j] 

    for i in range(10):
        for j in range(18):
            Brain[i][j] -= 2*(Output[i] - Output9[i])*Number9[j]*Decent_Rate

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

# print Brain

#####################
## Verifying Brain ##
#####################

Correct_Answer = argv[1]

for i in range(10):
    for j in range(18):
        Output[i] += Brain[i][j]*Number2[j]

print Output

for i in range(10):
    if (Output[i] > Best_Value):
        Best_Value = Output[i]
        Answer = i
    else: 
        None

print 'The Brain Guesses the number is: ', Answer 
print 'The Real Number is:', Correct_Answer 
