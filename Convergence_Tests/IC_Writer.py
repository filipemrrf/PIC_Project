"""
 " @file IC_Writer.py
 " @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 " @brief Writes the initial conditions for the main executable to run
 " @version 1.0
 " @date 2023-05-24
 " 
 " @copyright Copyright (c) 2023
 " 
"""

from numpy import cos, exp, sin, pi

def Constant(x):
    return 0

def Cos(x):
    # Sets the amplitude of the sin
    A = 1

    return A*cos(2*pi*x)

def Gauss(x):
    # Sets the standard deviation and the amplitude of the gaussian
    A = 0.002
    std_dev = 0.5
    c = 1/(std_dev*std_dev)

    return A*exp(-0.5*c*x*x)+1

def Log_deriv_Gauss(x):
    # Sets the standard deviation and the amplitude of the gaussian
    A = 0.002
    std_dev = 0.5
    c = 1/(std_dev*std_dev)

    return -(A*c*x*exp(-0.5*c*x*x))/(1+exp(-0.5*c*x*x))

def Sin(x):
    # Sets the amplitude of the sin
    A = 1

    return A*sin(2*pi*x)

# Sets the required parameters for the script to run
filename = "0"
function = Constant
size = 10
point_density = 50

# Initializes the number of points
NPoints = int(point_density*size)

# Loops through all the files that will be written
for i in range(0, 5):
    # Opens the file that will be written
    FILE = open(("Initial_Conditions/IC-" + filename + "-" + str(NPoints) + "p.dat"), "w")

    # Writtes the initial conditions for the equation
    for j in range(NPoints):
        # Devines an auxiliary variable x
        x = j*size/NPoints

        # Writes the IC to the file
        FILE.write(str(x) + " " + str(function(x)) + "\n")

    # Closes the file
    FILE.close()

    # Updates how many points the next initial conditions will have
    NPoints *= 2