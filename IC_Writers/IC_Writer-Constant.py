"""
 " @file IC_Writer-Constant.py
 " @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 " @brief Writes the initial conditions for the main executable to run
 " @version 1.0
 " @date 2023-04-24
 " 
 " @copyright Copyright (c) 2023
 " 
"""

from numpy import exp

# Sets the standard deviation of the gaussian
std_dev = 0.1
c = 1/(std_dev*std_dev)

# Sets the full size of the space and the space of interest
size = 1

# Initializes the number of points
NPoints = int(50*size)

# Loops through all the files that will be written
for i in range(0, 5):
    # Opens the file that will be written
    FILE = open(("Initial_Conditions/IC-1-" + str(NPoints) + "p.dat"), "w")

    # Writtes the initial conditions for the equation
    for j in range(NPoints):
        # Devines an auxiliary variable x
        x = j*size/NPoints

        # Writes the IC to the file
        FILE.write(str(x) + " " + str(1) + "\n")

    # Closes the file
    FILE.close()

    # Updates how many points the next initial conditions will have
    NPoints *= 2