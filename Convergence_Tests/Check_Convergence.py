"""
 " @file Check_Convergence.py
 " @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 " @brief Checks the convergence of the data
 " @version 5.1
 " @date 2023-04-01
 " 
 " @copyright Copyright (c) 2023
 " 
"""

import os
import matplotlib.pyplot as plt

# Does the convergence check with all the given parameters
def Convergence_Test(EQ, folder, field, NPoints, Spherical):
    os.system("python3 Convergence_Tests/Convergence_Test.py" + " " + \
        folder + str(NPoints) + "p/" + EQ + "_" + field + "-" + str(NPoints) + "p.dat " + \
            folder + str(2*NPoints) + "p/" + EQ + "_" + field + "-" + str(2*NPoints) + "p.dat " + \
                folder + str(4*NPoints) + "p/" + EQ + "_" + field + "-" + str(4*NPoints) + "p.dat " + \
                    folder + str(8*NPoints) + "p/" + EQ + "_" + field + "-" + str(8*NPoints) + "p.dat " + \
                        folder + str(16*NPoints) + "p/" + EQ + "_" + field + "-" + str(16*NPoints) + "p.dat " + \
                            "-SX " + str(Space_Size/NPoints) + " " + \
                                "-NP " + str(NPoints) + " " + \
                                    "-DIR " + (folder + "Convergence-" + field + "/ ") + "-S " + str(Spherical) + " "\
                                        "-CFL " + str(0.25))
    
# Does the convergence test with the exact solution
def Exact_Convergence(EQ, folder, field, NPoints, Space_Size, Spherical):
    os.system("python3 Convergence_Tests/Convergence_Test-Exact_Sol.py" + " " + \
        folder + str(NPoints) + "p/" + EQ + "_" + field + "-" + str(NPoints) + "p.dat " + \
            folder + str(2*NPoints) + "p/" + EQ + "_" + field + "-" + str(2*NPoints) + "p.dat " + \
                folder + str(4*NPoints) + "p/" + EQ + "_" + field + "-" + str(4*NPoints) + "p.dat " + \
                    folder + str(8*NPoints) + "p/" + EQ + "_" + field + "-" + str(8*NPoints) + "p.dat " + \
                        folder + str(16*NPoints) + "p/" + EQ + "_" + field + "-" + str(16*NPoints) + "p.dat " + \
                            "-SX " + str(Space_Size/NPoints) + " " + \
                                "-NP " + str(NPoints) + " " + \
                                    "-L " + str(Space_Size) + " "\
                                        "-DIR " + (folder + "Convergence-" + field + "/ ") + "-S " + str(Spherical) + " "\
                                            "-CFL " + str(0.25))

# Scales the pointwise conversion to the base resolution
def Scale_Pointwise(folder, field, NPoints, Scale):
    os.system("python3 " + "Convergence_Tests/Scale_Pointwise.py" + " " + \
        "-NP " + str(NPoints) + " " + "-DIR " + (folder + "Convergence-" + field + "/ ") + "-S " + str(Scale))

# Chooses the equation to analyse
EQ = "adm_evolution"
Acc = 2

# Declaration of the variables to control the cmunionvergence analysis
if(EQ == "simple_wave"):
    if(Acc == 2):
        folder = "Results/simple_wave-2nd_order/"
        S = 4

    if(Acc == 4):
        folder = "Results/simple_wave-4th_order/"
        S = 16

    Fields = ["Phi", "Pi"]
    Space_Size = 1
    
    Exact_Sol = [True, ""]
    Spherical = "0"

if(EQ == "non_linear_simple_wave"):
    if(Acc == 2):
        folder = "Results/non_linear_simple_wave-2nd_order/"
        S = 4

    if(Acc == 4):
        folder = "Results/non_linear_simple_wave-4th_order/"
        S = 16

    Fields = ["Phi", "Pi"]
    Space_Size = 1

    Exact_Sol = [False, ""]
    Spherical = "0"

if(EQ == "spherical_wave"):
    if(Acc == 2):
        folder = "Results/spherical_wave-2nd_order/"
        S = 4

    if(Acc == 4):
        folder = "Results/spherical_wave-4th_order/"
        S = 16

    Fields = ["Phi", "Pi"]
    Space_Size = 5

    Exact_Sol = [False, ""]
    Spherical = "1"

if(EQ == "non_linear_spherical_wave"):
    if(Acc == 2):
        folder = "Results/non_linear_spherical_wave-2nd_order/"
        S = 4

    if(Acc == 4):
        folder = "Results/non_linear_spherical_wave-4th_order/"
        S = 16

    Fields = ["Phi", "Pi"]
    Space_Size = 5

    Exact_Sol = [False, ""]
    Spherical = "1"

if(EQ == "adm_evolution"):
    folder = "Results/adm_evolution-2nd_order/"
    S = 4

    Fields = ["A", "KA", "B", "KB", "lambda"]
    Space_Size = 10

    Exact_Sol = [False, ""]
    Spherical = "1"

# Displays the information the convergence analysis has started
os.system("echo Beginning convergence analysis\n\n")

# Loops through the various fields
for field in Fields:
    # Initializes the number of points
    Point_Density = 50
    NPoints = Point_Density*Space_Size

    # Creates a folder for all of the comparison data to be saved
    os.system("mkdir " + folder + "Convergence-" + field)

    # Displays the information the convergence check is running for a specific set of points
    os.system("echo Checking convergence between the solutions for " + field + "\n")
    Convergence_Test(EQ, folder, field, NPoints, Spherical)

    # Scales the Pointwise convergence tests
    Scale_Pointwise(folder, field, NPoints, Acc*Acc)

    # Makes legend for the plot
    legend = ["low resolution", str(Acc) + " x medium resolution", str(pow(Acc,2)) + " x high resolution", str(pow(Acc,3)) + " x higher resolution"]

    # Loops through all the comparisons between numerical solutions
    for i in range(4):
        filename = folder + "Convergence-" + field + "/" + str(NPoints) + "p," + str(2*NPoints) + "p-Point_Comparison"
        if i!=0:
            filename += "-Scaled"
        filename += ".dat"
        # Opens the comparison file
        IN = open(filename, "r")

        # Declares lists to store the data
        time = []
        pointwise = []

        # Loops through the file
        for l in IN:
            # Saves the time to the respective list
            if l[0] == "\"":
                aux = l.split()
                time.append(float(aux[2]))
            # Saves the pointwise error to the respective list
            else:
                if l[0] != "\n":
                    aux = l.split()
                    if float(aux[0]) == 0.0:
                        pointwise.append(float(aux[1]))

        # Draws the graph
        plt.plot(time.copy(), pointwise.copy())
        time.clear()
        pointwise.clear()

        # Closes the file
        IN.close()

        # Updates the number of points
        NPoints *= 2

    # Saves the graph
    plt.xlabel("Time")
    plt.ylabel("Pointwise Convergence Comparison")
    plt.legend(legend)
    plt.savefig(folder + "Convergence-" + field + "/Pointwise_Convergence_Full_Comparison.png")

    # Clears the graph
    plt.figure().clear()
    plt.close()
    plt.cla()
    plt.clf()
    legend.clear()


# Checks if the solutions are to be compared with an exact solution
if Exact_Sol[0]:
    # Loops through the various fields
    for field in Fields:
        # Initializes the number of points
        Point_Density = 50
        NPoints = Point_Density*Space_Size

        # Creates a folder for all of the comparison data to be saved
        os.system("mkdir " + folder + "Convergence-" + field)

        # Displays the information the convergence check is running for a specific set of points
        os.system("echo Checking convergence with the exact solution for " + field + "\n")
        Exact_Convergence(EQ, folder, field, NPoints, Space_Size, Spherical)

        # Scales the Pointwise convergence tests
        Scale_Pointwise(folder, field, NPoints, Acc*Acc)

        # Makes legend for the plot
        legend = ["low resolution", str(Acc) + " x medium resolution", str(pow(Acc,2)) + " x high resolution", str(pow(Acc,3)) + " x higher resolution", str(pow(Acc,4)) + " x highest resolution"]

        # Loops through all the comparisons between numerical solutions
        for i in range(4):
            filename = folder + "Convergence-" + field + "/" + str(NPoints) + ",Exact_Solution-Point_Comparison"
            if i!=0:
                filename += "-Scaled"
            filename += ".dat"
            # Opens the comparison file
            IN = open(filename, "r")

            # Declares lists to store the data
            time = []
            pointwise = []

            # Loops through the file
            for l in IN:
                # Saves the time to the respective list
                if l[0] == "\"":
                    aux = l.split()
                    time.append(float(aux[2]))
                # Saves the pointwise error to the respective list
                else:
                    if l[0] != "\n":
                        aux = l.split()
                        if float(aux[0]) == 0.0:
                            pointwise.append(float(aux[1]))

            # Draws the graph
            plt.plot(time.copy(), pointwise.copy())
            time.clear()
            pointwise.clear()

            # Closes the file
            IN.close()

            # Updates the number of points
            NPoints *= 2

        # Saves the graph
        plt.xlabel("Time")
        plt.ylabel("Pointwise Convergence Comparison")
        plt.legend(legend)
        plt.savefig(folder + "Convergence-" + field + "/Pointwise_Convergence_Full_Comparison.png")

        # Clears the graph
        plt.figure().clear()
        plt.close()
        plt.cla()
        plt.clf()
        legend.clear()

# Displays the information the convergence analysis has finished
os.system("echo Convergence analysis finished")