"""
 " @file Convergence_Test.py
 " @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 " @brief Compares the solutions given and checks their conversion
 " @version 3.0
 " @date 2023-05-17
 " 
 " @copyright Copyright (c) 2022
 " 
"""

import sys
import matplotlib.pyplot as plt
import numpy as np

# Reads the file and stores its information in Data
def Read_File(filename, time, space, intensity):
    IN = open(filename, "r")

    auxv = []

    for l in IN:
        if l[0] == "\"":
            aux = l.split()
            time.append(float(aux[2]))
        elif l != "\n":
            aux = l.split()
            space.append(float(aux[0]))
            auxv.append(float(aux[1]))
        else:
            intensity.append(auxv.copy())
            auxv.clear()

    IN.close()

# Declares the variables needed for the plot
time = []
space = []
intensity = []

# Reads the input file
Read_File(sys.argv[1], time, space, intensity)

# Reads the command line arguments to define variables needed for the comparison
directory = sys.argv[2]

# Create the intensity plot
plt.imshow(intensity, cmap='RdYlBu', interpolation='nearest', aspect='auto', extent=[min(space), max(space), max(time), min(time)])

# Add colorbar for reference
plt.colorbar()

# Add labels
plt.ylabel("Time")
plt.xlabel("Space")

# Saves the figure
plt.savefig(directory + "Intensity.png")