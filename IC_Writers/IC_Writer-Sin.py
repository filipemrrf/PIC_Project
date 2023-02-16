"""
 " @file IC_Writer-Sin.py
 " @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 " @brief Writes the initial conditions for the main executable to run
 " @version 2.2
 " @date 2023-02-16
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