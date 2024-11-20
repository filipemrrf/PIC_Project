"""
 " @file Intensity_Plot.py
 " @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 " @brief Makes an intensity plot of the intensity of the solutions
 " @version 2.0
 " @date 2024-11-19
 " 
 " @copyright Copyright (c) 2023
 " 
"""

import argparse
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# Reads the file and stores its information in Data
def Read_File(filename, time, space, intensity, Tmax):
    IN = open(filename, "r")

    auxv = []

    for l in IN:
        if l[0] == "\"":
            aux = l.split()
            if Tmax != None:
                if float(aux[2]) > Tmax:
                    break
            time.append(float(aux[2]))
        elif l != "\n":
            aux = l.split()
            space.append(float(aux[0]))
            auxv.append(float(aux[1]))
        else:
            intensity.append(auxv.copy())
            auxv.clear()

    IN.close()

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description='Process the arguments to get the data')

# Define the arguments
parser.add_argument('solution', type=str, help='The solution to plot')
parser.add_argument('--Dir', type=str, default='', help='Directory to save the plots (default: \'.\')')
parser.add_argument('--Tmax', type=float, help='The max time to plot')

# Parse the arguments
args = parser.parse_args()

filename = args.solution # The name of the file to read
directory = args.Dir # The directory to save the plots
Tmax = args.Tmax # The max time to plot

# Declares the variables needed for the plot
time = []
space = []
intensity = []

# Reads the input file
Read_File(filename, time, space, intensity, Tmax)

# Inverts the order of the intensity list
intensity = intensity[::-1]

# Defines the colormap
colors = [(0, 'red'), (0.5, 'white'), (1, 'blue')]  # Blue to white to red
custom_cmap = LinearSegmentedColormap.from_list('custom_cmap', colors)

# Create the intensity plot
plt.imshow(intensity, cmap=custom_cmap, interpolation='nearest', aspect='auto', extent=[min(space), max(space), min(time), max(time)])

# Add colorbar for reference
plt.colorbar()

# Add labels
plt.ylabel("Time")
plt.xlabel("Space")

# Saves the figure
plt.savefig(directory + "Intensity.png")