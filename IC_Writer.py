import numpy

# Opens the file that will be written
FILE = open("IC.dat", "w")

# Writtes the number of equations in the ODE systemspace step and the number of points to the file
FILE.write("#NEq: 2\n")
FILE.write("#step_x: 0.02\n")
FILE.write("#NPoints: 51\n")

# Writtes the initial conditions for the equation
for i in range(50):
    FILE.write(str(numpy.sin(2*numpy.pi*i/50)))
    FILE.write(" ")
    FILE.write(str(2*numpy.pi*numpy.cos(2*numpy.pi*i/50)))
    FILE.write("\n")

# Closes the file
FILE.close()