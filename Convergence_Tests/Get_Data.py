"""
 " @file Get_Data.py
 " @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 " @brief Uses the main executable to get data
 " @version 2.0
 " @date 2023-01-16
 " 
 " @copyright Copyright (c) 2023
 " 
"""

import os

# Declaration of the variables to control the data aquisition process
EQ = ""
params = "0"
IC = ""
T = 0
W = 0

spherical = False

# Displays the information the aquisition process has started
os.system("echo Starting data aquisition process")

# Initializes the number of points
NPoints = 50

# Loops through all the number of points that the equation will be solved for
for i in range(0, 5):
    # Displays the information the equation is being solved for a specific number of points
    os.system("echo Solving the equation for initial conditions with " + str(NPoints) + " points")


    # Writes the parameters file for the run
    # Opens the file
    FILE = open("Data/Parameters-" + EQ + "-" + str(NPoints) + "p.txt", "w")

    # Writes the data to the file
    FILE.write("#Eq: " + EQ + "\n")
    FILE.write("#params: " + params + "\n")
    FILE.write("\n")
    FILE.write("#IC: Data/IC-" + IC + "-" + str(NPoints) + "p.dat\n")
    FILE.write("#step_x: " + str(1/NPoints) + "\n")
    FILE.write("#NPoints: " + str(NPoints+1) + "\n")
    FILE.write("\n")
    FILE.write("#FN: Data/" + EQ + "-" + str(NPoints) +"p.dat\n")
    FILE.write("#T: " + str(T) + "\n")
    FILE.write("#W: " + str(W) + "\n")

    # Closes the file
    FILE.close()

    # Solves the equation with the parameters given
    os.system("./main Data/Parameters-" + EQ + "-" + str(NPoints) + "p.txt")

    os.system("echo \n")

    # Updates how many points the next initial conditions will have and how many timesteps will be written to disk
    NPoints *= 2
    W *= 2

# Displays the information the aquisition process has ended
os.system("echo Data aquisition process finished")