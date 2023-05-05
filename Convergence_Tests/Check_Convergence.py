"""
 " @file Check_Convergence.py
 " @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 " @brief Checks the convergence of the data
 " @version 5.0
 " @date 2023-04-01
 " 
 " @copyright Copyright (c) 2023
 " 
"""

import os
import matplotlib.pyplot as plt

# Declaration of the variables to control the convergence analysis
EQ = "spherical_wave"
field = "Pi"
folder = "Results/spherical_wave-4th_order/"
Space_Size = 5
S = 16

Exact_Sol = [False, "Sin"]
Spherical = True


# Displays the information the convergence analysis has started
os.system("echo Beginning convergence analysis\n")

# Initializes the number of points
NPoints = 500

# Creates a folder for all of the comparison data to be saved
os.system("mkdir " + folder + "Convergence-" + field)

# Decides which of the convergence scripts will be run
if Spherical:
    script = "Convergence_Tests/Convergence_Test-Spherical.py"
else:
    script = "Convergence_Tests/Convergence_Test.py"

# Loops through all the comparisons that need to be computed
for i in range(3):
    # Displays the information the convergence check is running for a specific set of points
    os.system("echo Checking convergence between the solutions with " + str(int(NPoints/2)) + ", " + str(NPoints) + " and " + str(NPoints*2))

    # Does the convergence check with all the given parameters
    os.system("python3 " + script + " " + \
            folder + str(int(NPoints/2)) + "p/" + EQ + "_" + field + "-" + str(int(NPoints/2)) + "p.dat " + \
                folder + str(NPoints) + "p/" + EQ + "_" + field + "-" + str(NPoints) + "p.dat " + \
                    folder + str(2*NPoints) + "p/" + EQ + "_" + field + "-" + str(2*NPoints) + "p.dat " + \
                        "-SX " + str(2*Space_Size/NPoints) + " " + str(Space_Size/NPoints) + " " + str(Space_Size/(2*NPoints)) + " "\
                            "-NP " + str(int(NPoints/2)) + " " + str(NPoints) + " " + str(2*NPoints) + " " + \
                                "-DIR " + (folder + "Convergence-" + field + "/ ") + "-S " + str(S))

    os.system("echo \n")

    # Checks if the solutions are to be compared with an exact solution
    if Exact_Sol[0]:
        # If it's the first iteration of the loop, there is one extra comparison to be made
        if i == 0:
            # Creates a folder for all of the comparison data to be saved
            os.system("mkdir " + folder + "Convergence-" + field + "-Exact_Sol/ " + EQ)

            # Displays the information the convergence check with the exact solution is running for a specific set of points
            os.system("echo Checking convergence between the solutions with " + str(int(NPoints/2)) + " and " + str(NPoints) + " and the exact solution")

            # Does the convergence check with all the given parameters
            os.system("python3 Convergence_Tests/Sol_Comparer-" + Exact_Sol[1] + ".py " + \
                folder + str(int(NPoints/2)) + "p/" + EQ + "_" + field + "-" + str(int(NPoints/2)) + "p.dat " + \
                    folder + str(NPoints) + "p/" + EQ + "_" + field + "-" + str(NPoints) + "p.dat " + \
                        "-SX "  + str(2*Space_Size/NPoints) + " " + str(Space_Size/NPoints) + " " + " "\
                            "-NP " + str(int(NPoints/2)) + " " + str(NPoints) + " " + \
                                "-DIR " + (folder + "Convergence-" + field + "-Exact_Sol/ ") + "-S " + str(S))

            os.system("echo \n")
        

        # Displays the information the convergence check with the exact solution is running for a specific set of points
        os.system("echo Checking convergence between the solutions with " + str(NPoints) + " and " + str(2*NPoints) + " and the exact solution")

        # Does the convergence check with all the given parameters
        os.system("python3 Convergence_Tests/Sol_Comparer-" + Exact_Sol[1] + ".py " + \
            folder + str(NPoints) + "p/" + EQ + "_" + field + "-" + str(NPoints) + "p.dat " + \
                folder + str(2*NPoints) + "p/" + EQ + "_" + field + "-" + str(2*NPoints) + "p.dat " + \
                    "-SX "  + str(Space_Size/NPoints) + " " + str(Space_Size/(2*NPoints)) + " " + \
                        "-NP " + str(NPoints) + " " + str(2*NPoints) + " " + \
                            "-DIR " + (folder + "Convergence-" + field + "-Exact_Sol/ ") + "-S " + str(S))

        os.system("echo \n")
        
    # Updates how many points the next initial conditions will have
    NPoints *= 2


# Displays the information that the norm comparison between all files is being plotted
os.system("echo Plotting the norm comparison with all the data")

# Re-initializes the number of points
NPoints = 250

# Declares a list for the plot's legend
legend = []

# Loops through all the comparisons between numerical solutions
for i in range(3):
    # Opens the comparison file
    IN = open(folder + "Convergence-" + field + "/" + str(NPoints) + "p," + str(4*NPoints) + "p-Norm_Comparison.dat", "r")

    # Declares lists to hold the plot's data
    time = []
    norm = []

    # Loops through the file
    for l in IN:
        # Fills the lists 
        if l[0] != "#":
            aux = l.split(" ")

            time.append(float(aux[0]))
            norm.append(float(aux[1]))

    # Draws the graph
    plt.plot(time.copy(), norm.copy())
    legend.append(str(NPoints) + "p," + str(4*NPoints) + "p Norm Comparison")
    time.clear()
    norm.clear()

    # Closes the file
    IN.close()

    # Updates the number of points
    NPoints *= 2


# Saves the graph
plt.title(EQ + " Norm Convergence Comparison")
plt.ylim([0,5])
plt.xlabel("Time")
plt.ylabel("Log2(Norm Convergence Comparison)")
plt.legend(legend)
plt.savefig(folder + "Convergence-" + field + "/Norm_Convergence_Full_Comparison.png")


if Exact_Sol[0]:
    # Clears the graph so a new one can be drawn and resets the number of points
    plt.clf()
    legend.clear()
    NPoints = 50

    # Loops through all the comparisons with the exact solutions
    for i in range(4):
        # Opens the comparison file
        IN = open(folder + "Convergence-" + field + "-Exact_Sol/" + str(NPoints) + "p," + str(2*NPoints) + "p-Exact_Sol_Norm_Comparison.dat", "r")

        # Declares lists to hold the plot's data
        time = []
        norm = []

        # Loops through the file
        for l in IN:
            # Fills the lists 
            if l[0] != "#":
                aux = l.split(" ")

                time.append(float(aux[0]))
                norm.append(float(aux[1]))

        # Draws the graph
        plt.plot(time.copy(), norm.copy())
        legend.append(str(NPoints) + "p,Exact_Sol," + str(2*NPoints) + "p Norm Comparison")
        time.clear()
        norm.clear()

        # Closes the file
        IN.close()

        # Updates the number of points
        NPoints *= 2

    # Saves the graph
    plt.title(EQ + " Norm Convergence Comparison")
    plt.ylim([0,5])
    plt.xlabel("Time")
    plt.ylabel("Log2(Norm Convergence Comparison)")
    plt.legend(legend)
    plt.savefig(folder + "Convergence-" + field + "-Exact_Sol/Norm_Convergence_Full_Comparison.png")

os.system("echo \n")

# Displays the information the convergence analysis has finished
os.system("echo Convergence analysis finished")