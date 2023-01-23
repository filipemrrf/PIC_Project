"""
 " @file Get_Data.py
 " @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 " @brief Uses the main executable to get data
 " @version 1.0
 " @date 2023-01-09
 " 
 " @copyright Copyright (c) 2023
 " 
"""

import os

# Declaration of the variables to control the data aquisition process
IC = ""
EQ = ""
T = 0
W = 1


# Displays the information the aquisition process has started
os.system("echo Starting data aquisition process")

# Initializes the number of points
NPoints = 50

# Loops through all the number of points that the equation will be solved for
for i in range(0, 5):
    # Displays the information the equation is being solved for a specific number of points
    os.system("echo Solving the equation for initial conditions with " + str(NPoints) + " points")

    # Solves the equation with the parameters given
    os.system("./test -IC Data/IC-" + IC + "-" + str(NPoints) + "p.dat -EQ " + EQ + " -FN Data/" + EQ + "-" + str(NPoints) + "p.dat -T " + str(T) + " -W " + str(W))

    os.system("echo \n")

    # Updates how many points the next initial conditions will have and how many timesteps will be written to disk
    NPoints *= 2
    W *= 2

# Displays the information the aquisition process has ended
os.system("echo Data aquisition process finished")