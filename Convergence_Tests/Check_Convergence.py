"""
 " @file Check_Convergence.py
 " @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 " @brief Checks the convergence of the data
 " @version 6.0
 " @date 2024-11-19
 " 
 " @copyright Copyright (c) 2024
 " 
"""

import os
import argparse
import matplotlib.pyplot as plt

# Does the convergence check with all the given parameters
def Self_Convergence_Test(folder, field, NPoints, Space_Size, Runs, Spherical):
    files = ""

    N = NPoints

    for i in range(Runs):
        files += folder + str(N+1) + "p/" + field + ".dat "
        N *= 2

    os.system("python3 Convergence_Tests/Self_Convergence_Test.py" + " " + files + "--N " + str(NPoints) + " " + "--SX " + str(Space_Size/NPoints) +\
                " --Runs " + str(Runs) + " --Dir " + (folder + "Self_Convergence-" + field + "/ ") + "--S " + str(Spherical) + " --CFL " + str(0.25))

# Does the convergence test with the exact solution
def Exact_Convergence(folder, field, NPoints, Space_Size, Runs, Spherical):
    os.system("python3 Convergence_Tests/Convergence_Test-Exact_Sol.py" + " " + folder + str(NPoints+1) + "p/" + field + ".dat " + \
        folder + str(2*NPoints+1) + "p/" + field + ".dat " + folder + str(4*NPoints+1) + "p/" + field + ".dat " + \
            folder + str(8*NPoints+1) + "p/" + field + ".dat " + folder + str(16*NPoints+1) + "p/" + field + ".dat " +  \
                "--N " + str(NPoints) + " " + "--SX " + str(Space_Size/NPoints) + " " + "--Runs " + str(Runs) + " " \
                    "--Dir " + (folder + "Exact_Sol_Convergence-" + field + "/ ") + "--S " + str(Spherical) + " --CFL " + str(0.25)) + \
                        " --Exact " + Exact_Sol[1]

# Scales the pointwise conversion to the base resolution
def Scale_Pointwise(folder, field, NPoints, Runs, Scale):
    os.system("python3 " + "Convergence_Tests/Scale_Pointwise.py" + " " + \
        "-NP " + str(NPoints) + " -Runs " + str(Runs) + " -DIR " + (folder + "Self_Convergence-" + field + "/ ") + "-S " + str(Scale))

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description='Process the arguments to get the data')

# Define the arguments
parser.add_argument('equation', type=str, help='The equation to analyze')
parser.add_argument('--Acc', type=int, default=2, help='Accuracy order of the solution (default: 2)')
parser.add_argument('--N', type=int, default=100, help='Number of points (default: 100)')
parser.add_argument('--Size', type=float, default=1, help='Dimensions of the space (default: 1)')
parser.add_argument('--Runs', type=int, default=5, help='Number of different resolutions to check convergence for (default: 5)')

# Parse the arguments
args = parser.parse_args()

EQ = args.equation # Chooses the equation to analyse
Acc = args.Acc # Sets the accuracy order of the solution
N = args.N # Sets the initial number of points
Space_Size = args.Size # Sets the size of the space
Runs = args.Runs # Sets the number of different resolutions to check convergence for

# Declaration of the variables to control the cmunionvergence analysis
if(EQ == "simple_wave"):
    if(Acc == 2):
        folder = "Results/simple_wave-2nd_order/"
        S = 4

    if(Acc == 4):
        folder = "Results/simple_wave-4th_order/"
        S = 16

    Fields = ["Phi", "Pi"]
    
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

    Exact_Sol = [False, ""]
    Spherical = "1"

if(EQ == "spherical_reduced_wave"):
    if(Acc == 2):
        folder = "Results/spherical_reduced_wave-2nd_order/"
        S = 4

    Fields = ["Psi", "Phi", "Pi"]

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

    Exact_Sol = [False, ""]
    Spherical = "1"

if(EQ == "adm_evolution"):
    folder = "Results/adm_evolution-2nd_order-long/"
    S = 4

    Fields = ["A", "DA", "KA", "B", "DB", "KB", "lambda"]

    Exact_Sol = [False, ""]
    Spherical = "1"

if(EQ == "compact_wave_equation"):
    if(Acc == 2):
        folder = "Results/compact_wave_equation-2nd_order/"
        S = 4

    Fields = ["Psi", "Phi", "Pi", "Constraint"]
    
    Exact_Sol = [False, ""]
    Spherical = "0"

if(EQ == "spherical_compact_wave_equation"):
    if(Acc == 2):
        folder = "Results/spherical_compact_wave_equation-2nd_order/"
        S = 4

    Fields = ["Psi", "Phi", "Pi", "Constraint"]
    
    Exact_Sol = [False, ""]
    Spherical = "1"

if(EQ == "power_non_linear_spherical_compact_wave_equation"):
    if(Acc == 2):
        folder = "Results/power_non_linear_spherical_compact_wave_equation-2nd_order/"
        S = 4

    Fields = ["Psi", "Phi", "Pi", "Constraint"]
    
    Exact_Sol = [False, ""]
    Spherical = "1"

# Displays the information the convergence analysis has started
os.system("echo Beginning convergence analysis\n\n")

# Loops through the various fields
for field in Fields:
    # Creates a folder for all of the comparison data to be saved
    os.system("mkdir " + folder + "Self_Convergence-" + field)

    # Resets the number of points
    NPoints = N

    # Displays the information the convergence check is running for a specific set of points
    os.system("echo Checking convergence between the solutions for " + field + "\n")
    Self_Convergence_Test(folder, field, NPoints, Space_Size, Runs, Spherical)

    # Scales the Pointwise convergence tests
    Scale_Pointwise(folder, field, NPoints, Runs, Acc*Acc)

    # Makes legend for the plot
    legend = []

    # Loops through all the comparisons between numerical solutions
    for i in range(Runs-1):
        filename = folder + "Self_Convergence-" + field + "/" + str(NPoints+1) + "p," + str(2*NPoints+1) + "p-Point_Comparison"
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
                    if float(aux[0]) == 1.0:
                        pointwise.append(float(aux[1]))

        # Draws the graph
        plt.plot(time.copy(), pointwise.copy())
        legend.append("|" + str(NPoints+1) + "p - " + str(2*NPoints+1) + "p|")
        time.clear()
        pointwise.clear()

        # Closes the file
        IN.close()

        # Updates the number of points
        NPoints *= 2

    # Saves the graph
    plt.xlabel("Time")
    plt.ylabel("Pointwise Convergence Comparison")
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.legend(legend)
    plt.savefig(folder + "Self_Convergence-" + field + "/Pointwise_Convergence_Full_Comparison.png")

    # Clears the graph
    plt.figure().clear()
    plt.close()
    plt.cla()
    plt.clf()
    legend.clear()


# Checks if the solutions are to be compared with an exact solution <---Check this is right
if Exact_Sol[0]:
    # Loops through the various fields
    for field in Fields:
        # Creates a folder for all of the comparison data to be saved
        os.system("mkdir " + folder + "Exact_Sol_Convergence-" + field)

        # Resets the number of points
        NPoints = N

        # Displays the information the convergence check is running for a specific set of points
        os.system("echo Checking convergence with the exact solution for " + field + "\n")
        Exact_Convergence(folder, field, NPoints, Space_Size, Runs, Spherical)

        # Scales the Pointwise convergence tests
        Scale_Pointwise(folder, field, NPoints, Acc*Acc)

        # Makes legend for the plot
        legend = ["low resolution", str(Acc) + " x medium resolution", str(pow(Acc,2)) + " x high resolution", str(pow(Acc,3)) + " x higher resolution", str(pow(Acc,4)) + " x highest resolution"]

        # Loops through all the comparisons between numerical solutions
        for i in range(Runs):
            filename = folder + "Exact_Sol_Convergence-" + field + "/" + str(NPoints) + ",Exact_Solution-Point_Comparison"
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
        plt.savefig(folder + "Exact_Sol_Convergence-" + field + "/Pointwise_Convergence_Full_Comparison.png")

        # Clears the graph
        plt.figure().clear()
        plt.close()
        plt.cla()
        plt.clf()
        legend.clear()

# Creates the intensity plot of the solution
os.system("echo Creating the intensity plot of the solution\n")
os.system("python3 Convergence_Tests/Intensity_Plot.py " + folder + str(NPoints+1) + "p/" + Fields[0] + ".dat " + "--Dir " + folder)

# Displays the information the convergence analysis has finished
os.system("echo Convergence analysis finished\n")