import numpy as np
import matplotlib.pyplot as plt

#data_files = ["Results/Incoming_Waves/spherical_compact_wave_equation-2nd_order-1/IC-IncomingPsi-1601p.dat", "Results/Incoming_Waves/spherical_compact_wave_equation-2nd_order-2/IC-IncomingPsi-1601p.dat", "Results/Incoming_Waves/spherical_compact_wave_equation-2nd_order-3/IC-IncomingPsi-1601p.dat", "Results/Incoming_Waves/spherical_compact_wave_equation-2nd_order-4/IC-IncomingPsi-1601p.dat", "Results/Incoming_Waves/spherical_compact_wave_equation-2nd_order-5/IC-IncomingPsi-1601p.dat", "Results/Incoming_Waves/spherical_compact_wave_equation-2nd_order-6/IC-IncomingPsi-1601p.dat", "Results/Incoming_Waves/spherical_compact_wave_equation-2nd_order-7/IC-IncomingPsi-1601p.dat", "Results/Incoming_Waves/spherical_compact_wave_equation-2nd_order-8/IC-IncomingPsi-1601p.dat", "Results/Incoming_Waves/spherical_compact_wave_equation-2nd_order-9/IC-IncomingPsi-1601p.dat", "Results/Incoming_Waves/spherical_compact_wave_equation-2nd_order-10/IC-IncomingPsi-1601p.dat"]

data_files = ["Results/Incoming_Waves/error1.dat", "Results/Incoming_Waves/error2.dat", "Results/Incoming_Waves/error3.dat", "Results/Incoming_Waves/error4.dat", "Results/Incoming_Waves/error5.dat", "Results/Incoming_Waves/error6.dat", "Results/Incoming_Waves/error7.dat", "Results/Incoming_Waves/error8.dat", "Results/Incoming_Waves/error9.dat", "Results/Incoming_Waves/error10.dat"]

legend = ["$x_0 = 1$", "$x_0 = 2$", "$x_0 = 3$", "$x_0 = 4$", "$x_0 = 5$", "$x_0 = 6$", "$x_0 = 7$", "$x_0 = 8$", "$x_0 = 9$", "$x_0 = 10$"]


# Load data
for i in range(len(data_files)):
    FILE = open(data_files[i], 'r')

    # Declares lists to store the data
    x = []
    y = []

    # Loops through the file
    for l in FILE:
        if l[0] != "\n":
            aux = l.split()
            x.append(float(aux[0]))
            y.append(float(aux[1]))

    # Draws the graph
    plt.plot(x.copy(), y.copy())

    x.clear()
    y.clear()

    # Closes the file
    FILE.close()

plt.xlabel("Time")
plt.ylabel("Normalized Norm of the Error")
plt.xlim([0,10])
plt.ylim([0,6])
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.legend(legend)
plt.savefig("Report-PIC2/Images/Incoming_Error.png")
