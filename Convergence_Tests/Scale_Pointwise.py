"""
 " @file Scale_Pointwise.py
 " @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 " @brief Scales pointwise convergence of the data
 " @version 2.0
 " @date 2023-11-14
 " 
 " @copyright Copyright (c) 2024
 " 
"""

import sys

# Reads the command line arguments
for i in range(1, len(sys.argv)):
    # Chooses the directory where the results of the comparison will be saved
    if sys.argv[i] == "-DIR":
        directory = sys.argv[i+1]

    # Chooses the number of points for the files
    if sys.argv[i] == "-NP":
        NPoints = int(sys.argv[i+1])

    # Chooses the number of runs
    if sys.argv[i] == "-Runs":
        Runs = int(sys.argv[i+1])

    # Defines the scale factor
    if sys.argv[i] == "-S":
        Scale = int(sys.argv[i+1])


# Loops through the several files
for i in range(1, Runs-1):
    # Opens the file to be read and to be written
    IN = open(directory + str(2*NPoints+1) + "p," + str(4*NPoints+1) + "p-Point_Comparison.dat", "r")
    OUT = open(directory + str(2*NPoints+1) + "p," + str(4*NPoints+1) + "p-Point_Comparison-Scaled.dat", "w")

    # Reads the original file and writes the rescaled one
    for l in IN:
        if (l[0] == "\"") or (l == "\n"):
            OUT.write(l)
        else:
            aux = l.split()
            OUT.write(aux[0] + " " + str(pow(Scale, i)*float(aux[1])) + "\n")

    # Closes the files
    IN.close()
    OUT.close()

    # Updates the number of points
    NPoints *= 2