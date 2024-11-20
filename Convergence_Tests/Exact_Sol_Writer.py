"""
 " @file Exact_Sol_Writer.py
 " @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 " @brief Writes the exact solution for the exact convergence tests
 " @version 1.0
 " @date 2024-11-14
 " 
 " @copyright Copyright (c) 2024
 " 
"""

from numpy import exp, sin, pi, sqrt

# Exact solution to simple wave equation
def sin(t,x):
    return sin(2*pi*(x-t))

def compact_wave_Psi(t,x):
    A = 1
    std = 0.1
    c = 1/(std*std)

    S = 1

    return 0.5*(A*exp(-0.5*c*(t + (2*S*S*x)/(-S*S + x*x) + sqrt((S*S*(S*S+x*x)*(S*S+x*x))/((S*S-x*x)*(S*S-x*x))))**2) + \
            A*exp(-0.5*c*(t + 2*x/(1-x*x/(S*S)) + sqrt(S*S + 4*x*x/((-1+x*x/(S*S))*(1-x*x/(S*S)))))**2))


# Sets the required parameters for the script to run
filename = "compact_wave_Psi"
function = compact_wave_Psi

size = 2
x0 = -1
point_density = 50
dx = size/point_density

T = 5
cfl = 0.25
dt = cfl*dx

# Initializes the number of points
NPoints = int(point_density*size)

# Loops through all the files that will be written
for i in range(0, 5):
    # Opens the file that will be written
    FILE = open(("Exact_Solutions/Exact-" + filename + "-" + str(NPoints+1) + "p.dat"), "w")

    # Writtes the initial conditions for the equation
    for time in range(T):
        # Calculates the time variable
        t = time*dt
        FILE.write("\"Time = " + str(t) + "\n")

        for space in range(NPoints+1):
            # Calculates the space variable
            x = x0 + space*dx

            # Writes the IC to the file
            FILE.write(str(x) + " " + str(function(t,x)) + "\n")

    # Closes the file
    FILE.close()

    # Updates how many points the next initial conditions will have
    NPoints *= 2