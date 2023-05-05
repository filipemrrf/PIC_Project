"""
 " @file IC_Writer-Sin.py
 " @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 " @brief Writes the initial conditions for the main executable to run
 " @version 3.0
 " @date 2023-03-28
 " 
 " @copyright Copyright (c) 2023
 " 
"""

import numpy

# Initializes the number of points
NPoints = 50

for i in range(0, 5):
    # Opens the file that will be written
    FILE_PHI = open(("Initial_Conditions/IC-Sin(2pi.x)-" + str(NPoints) + "p.dat"), "w")
    FILE_PI = open(("Initial_Conditions/IC-2pi.Cos(2pi.x)-" + str(NPoints) + "p.dat"), "w")

    # Writtes the initial conditions for the equation
    for j in range(NPoints):
        FILE_PHI.write(str(j/NPoints) + " " + str(numpy.sin(2*numpy.pi*(j/NPoints))) + "\n")
        FILE_PI.write(str(j/NPoints) + " " + str(2*numpy.pi*numpy.cos(2*numpy.pi*(j/NPoints))) + "\n")

    # Closes the file
    FILE_PHI.close()
    FILE_PI.close()

    # Updates how many points the next initial conditions will have
    NPoints *= 2