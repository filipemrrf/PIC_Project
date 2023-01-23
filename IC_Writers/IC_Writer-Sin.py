"""
 " @file IC_Writer-Sin.py
 " @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 " @brief Writes the initial conditions for the main executable to run
 " @version 2.1
 " @date 2023-01-13
 " 
 " @copyright Copyright (c) 2023
 " 
"""

import numpy

# Initializes the number of points
NPoints = 50

for i in range(0, 5):
    # Opens the file that will be written
    FILE = open(("Data/IC-Sin-" + str(NPoints) + "p.dat"), "w")

    # Writtes the number of equations in the ODE systemspace step, the number of points to the file and the parameters for the equation
    FILE.write("#NEq: 2\n")
    FILE.write("#step_x: " + str(1/NPoints) + "\n")
    FILE.write("#NPoints: " + str(NPoints+1) + "\n")
    FILE.write("#pars: 3 1.0 2.0 0.01\n")

    # Writtes the initial conditions for the equation
    for j in range(NPoints):
        FILE.write(str(numpy.sin(2*numpy.pi*j/NPoints)))
        FILE.write(" ")
        FILE.write(str(2*numpy.pi*numpy.cos(2*numpy.pi*j/NPoints)))
        FILE.write("\n")

    FILE.write("0.0 " + str(2*numpy.pi))

    # Closes the file
    FILE.close()

    # Updates how many points the next initial conditions will have
    NPoints *= 2