"""
 " @file IC_Writer-Gauss.py
 " @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 " @brief Writes the initial conditions for the main executable to run
 " @version 3.0
 " @date 2023-04-18
 " 
 " @copyright Copyright (c) 2023
 " 
"""

from numpy import exp

# Sets the standard deviation of the gaussian
std_dev = 0.1
c = 1/(std_dev*std_dev)

# Sets the full size of the space and the space of interest
size = 5

# Initializes the number of points
NPoints = int(50*size)

# Loops through all the files that will be written
for i in range(0, 5):
    # Opens the file that will be written
    FILE_Phi = open(("Initial_Conditions/IC-exp(-50x^2)-" + str(NPoints) + "p.dat"), "w")
    FILE_Pi = open(("Initial_Conditions/IC-(-100).exp(-50x^2)-" + str(NPoints) + "p.dat"), "w")

    # Writtes the initial conditions for the equation
    for j in range(NPoints):
        # Devines an auxiliary variable x
        x = j*size/NPoints

        # Writes the IC to the file
        FILE_Phi.write(str(x) + " " + str(exp(-0.5*c*x*x)) + "\n")
        FILE_Pi.write(str(x) + " " + str(-c*exp(-0.5*c*x*x)) + "\n")

    # Closes the file
    FILE_Phi.close()
    FILE_Pi.close()

    # Updates how many points the next initial conditions will have
    NPoints *= 2