"""
 " @file IC_Writer.py
 " @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 " @brief Writes the initial conditions for the main executable to run
 " @version 2.0
 " @date 2024-11-14
 " 
 " @copyright Copyright (c) 2024
 " 
"""

from numpy import cos, exp, sin, pi, sqrt

def Constant(x):
    return 0

def Cos(x):
    # Sets the amplitude of the sin
    A = 1

    return A*cos(2*pi*x)

def Gauss(x):
    # Sets the standard deviation and the amplitude of the gaussian
    A = 1
    std_dev = 0.1
    c = 1/(std_dev*std_dev)

    return A*exp(-0.5*c*x*x)

def Gauss_Deriv(x):
    # Sets the standard deviation and the amplitude of the gaussian
    A = 1
    std_dev = 0.1
    c = 1/(std_dev*std_dev)

    return -A*c*x*exp(-0.5*c*x*x)*pow(Omega(x), 2)/L(x)

def Log_deriv_Gauss(x):
    # Sets the standard deviation and the amplitude of the gaussian
    A = 0.002
    std_dev = 0.5
    c = 1/(std_dev*std_dev)

    return -(A*c*x*exp(-0.5*c*x*x))/(1+A*exp(-0.5*c*x*x))

def Sin(x):
    # Sets the amplitude of the sin
    A = 1

    return A*sin(2*pi*x)

def Omega(x):
    S = 1

    return (1 - x*x/(S*S))/2

def L(x):
    S = 1

    return (1 + x*x/(S*S))/2

def H(x):
    S = 1

    return (2*S*x)/(S*S + x*x)

def A(x):
    S = 1

    return -(S*S + x*x)/(2*S*S)

def B(x):
    S = 1

    return -8*(S**6)*(x**2+S**2)/pow(x**4 + S**4 + (x**2)*(4*(S**4) - 2*(S**2)), 2)

# Sets the required parameters for the script to run
filename = "Gauss"
function = Gauss
size = 1
point_density = 100

# Initializes the number of points
NPoints = int(point_density*size)

# Loops through all the files that will be written
for i in range(0, 5):
    # Opens the file that will be written
    FILE = open(("Initial_Conditions/IC-" + filename + "-" + str(NPoints+1) + "p.dat"), "w")

    # Writtes the initial conditions for the equation
    for j in range(NPoints+1):
        # Devines an auxiliary variable x
        x = j*size/NPoints

        # Writes the IC to the file
        FILE.write(str(x) + " " + str(function(x)) + "\n")

    # Closes the file
    FILE.close()

    # Updates how many points the next initial conditions will have
    NPoints *= 2