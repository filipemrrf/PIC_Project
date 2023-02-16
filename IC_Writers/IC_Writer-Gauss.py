"""
 " @file IC_Writer-Gauss.py
 " @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 " @brief Writes the initial conditions for the main executable to run
 " @version 1.0
 " @date 2023-01-16
 " 
 " @copyright Copyright (c) 2023
 " 
"""

from numpy import exp

# Initializes the number of points
NPoints = 50

for i in range(0, 5):
    # Opens the file that will be written
    FILE = open(("Data/IC-Gauss-" + str(NPoints) + "p.dat"), "w")

    # Writtes the number of equations in the ODE systemspace step, the number of points to the file and the parameters for the equation
    FILE.write("#NEq: 2\n")
    FILE.write("#step_x: " + str(1/NPoints) + "\n")
    FILE.write("#NPoints: " + str(NPoints) + "\n")
    FILE.write("#pars: 1 1.0\n")

    # Writtes the initial conditions for the equation
    for j in range(NPoints):
        FILE.write(str(exp(-(j/NPoints)*(j/NPoints)/0.03)))
        FILE.write(" ")
        FILE.write(str(-0.06*(j/NPoints)*exp(-(j/NPoints)*(j/NPoints)*0.03)))
        FILE.write("\n")

    FILE.write("0.0 0.0")

    # Closes the file
    FILE.close()

    # Updates how many points the next initial conditions will have
    NPoints *= 2