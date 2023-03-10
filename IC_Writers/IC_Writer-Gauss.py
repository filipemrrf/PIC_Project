"""
 " @file IC_Writer-Gauss.py
 " @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 " @brief Writes the initial conditions for the main executable to run
 " @version 2.0
 " @date 2023-02-27
 " 
 " @copyright Copyright (c) 2023
 " 
"""

from numpy import exp

# Sets the standard deviation of the gaussian
std_dev = 0.1
c = 1/(std_dev*std_dev)

# Sets the full size of the space and the space of interest
full_size = 30
relevant_size = 5
size_ratio = full_size/relevant_size

# Initializes the number of points
NPoints = int(50*size_ratio)

# Loops through all the files that will be written
for i in range(0, 5):
    # Opens the file that will be written
    FILE = open(("Data/IC-Gauss-" + str(NPoints) + "p.dat"), "w")

    # Writtes the initial conditions for the equation
    for j in range(NPoints-1):
        # Devines an auxiliary variable x
        x = j*relevant_size/NPoints

        # Writes the IC to the file
        FILE.write(str(exp(-0.5*c*x*x)))
        FILE.write(" ")
        FILE.write(str(-x*c*exp(-0.5*c*x*x)))
        FILE.write("\n")

    FILE.write("0.0 0.0")

    # Closes the file
    FILE.close()

    # Updates how many points the next initial conditions will have
    NPoints *= 2