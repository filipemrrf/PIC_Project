"""
 " @file Check_Convergence.py
 " @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 " @brief Checks the convergence of the data
 " @version 4.0
 " @date 2023-02-16
 " 
 " @copyright Copyright (c) 2023
 " 
"""

import os
import matplotlib.pyplot as plt

# Declaration of the variables to control the convergence analysis
EQ = ""
Space_Size = 0
S = 0
Exact_Sol = [False, ""]
Spherical = False


# Displays the information the convergence analysis has started
os.system("echo Beginning convergence analysis\n")

# Initializes the number of points
NPoints = 100

# Creates a folder for all of the comparison data to be saved
os.system("mkdir Results/Convergence-" + EQ)

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
            "Data/" + EQ + "-" + str(int(NPoints/2)) + "p.dat Data/" + EQ + "-" + str(NPoints) + "p.dat Data/" + EQ + "-" + str(NPoints*2) + "p.dat " + \
                "-SX " + str(2*Space_Size/NPoints) + " " + str(Space_Size/NPoints) + " " + str(Space_Size/(2*NPoints)) + " "\
                    "-NP " + str(int(NPoints/2)) + " " + str(NPoints) + " " + str(2*NPoints) + " " + \
                        "-DIR " + ("Results/Convergence-" + EQ + "/ ") + "-S " + str(S))

    os.system("echo \n")

    # Checks if the solutions are to be compared with an exact solution
    if Exact_Sol[0]:
        # If it's the first iteration of the loop, there is one extra comparison to be made
        if i == 0:
            # Creates a folder for all of the comparison data to be saved
            os.system("mkdir Results/Convergence_Exact_Sol-" + EQ)

            # Displays the information the convergence check with the exact solution is running for a specific set of points
            os.system("echo Checking convergence between the solutions with " + str(int(NPoints/2)) + " and " + str(NPoints) + " and the exact solution")

            # Does the convergence check with all the given parameters
            os.system("python3 Convergence_Tests/Sol_Comparer-" + Exact_Sol[1] + ".py " + \
                "Data/" + EQ + "-" + str(int(NPoints/2)) + "p.dat Data/" + EQ + "-" + str(NPoints) + "p.dat " + \
                    "-SX "  + str(2*Space_Size/NPoints) + " " + str(Space_Size/NPoints) + " " + " "\
                        "-NP " + str(int(NPoints/2)) + " " + str(NPoints) + " " + \
                            "-DIR " + ("Results/Convergence_Exact_Sol-" + EQ + "/ ") + "-S " + str(S))

            os.system("echo \n")
        

        # Displays the information the convergence check with the exact solution is running for a specific set of points
        os.system("echo Checking convergence between the solutions with " + str(NPoints) + " and " + str(2*NPoints) + " and the exact solution")

        # Does the convergence check with all the given parameters
        os.system("python3 Convergence_Tests/Sol_Comparer-" + Exact_Sol[1] + ".py " + \
            "Data/" + EQ + "-" + str(NPoints) + "p.dat Data/" + EQ + "-" + str(2*NPoints) + "p.dat " + \
                "-SX "  + str(Space_Size/NPoints) + " " + str(Space_Size/(2*NPoints)) + " " + \
                    "-NP " + str(NPoints) + " " + str(2*NPoints) + " " + \
                        "-DIR " + ("Results/Convergence_Exact_Sol-" + EQ + "/ ") + "-S " + str(S))

        os.system("echo \n")
        
    # Updates how many points the next initial conditions will have
    NPoints *= 2


# Displays the information that the norm comparison between all files is being plotted
os.system("echo Plotting the norm comparison with all the data")

# Re-initializes the number of points
NPoints = 50

# Declares a list for the plot's legend
legend = []

# Loops through all the comparisons between numerical solutions
for i in range(3):
    # Opens the comparison file
    IN = open("Results/Convergence-" + EQ + "/" + str(NPoints) + "p," + str(4*NPoints) + "p-Norm_Comparison.dat", "r")

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
plt.savefig("Results/Convergence-" + EQ + "/Norm_Convergence_Full_Comparison.png")


if Exact_Sol[0]:
    # Clears the graph so a new one can be drawn and resets the number of points
    plt.clf()
    legend.clear()
    NPoints = 50

    # Loops through all the comparisons with the exact solutions
    for i in range(4):
        # Opens the comparison file
        IN = open("Results/Convergence_Exact_Sol-" + EQ + "/" + str(NPoints) + "p," + str(2*NPoints) + "p-Exact_Sol_Norm_Comparison.dat", "r")

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
    plt.savefig("Results/Convergence_Exact_Sol-" + EQ + "/Norm_Convergence_Full_Comparison.png")

os.system("echo \n")

# Displays the information the convergence analysis has finished
os.system("echo Convergence analysis finished")