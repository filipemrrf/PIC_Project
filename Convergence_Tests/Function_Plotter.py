import matplotlib.pyplot as plt
import numpy as np

# Define the function you want to plot
def f1(r):
    return -2*r**2/((1+r**2)*(-1 + 1/np.sqrt(1 + 4*r**2/((-1+r**2)**2))))

def f2(r):
    return (-1 + r**2)**2/((1+r**2)*(-2-(((-1 + r**2)**2) * (-1 + 1/np.sqrt(1 + 4*r**2/((-1+r**2)**2))))/2*r**2))

# Generate x values (domain of the function)
x = np.linspace(0, 1, 500)  # Generates 500 points between -10 and 10

# Compute y values (range of the function)
y1 = f1(x)
y2 = f2(x)

# Plot the function
plt.plot(x, y1, label='Outgoing')
plt.plot(x, y2, label='Incoming')

# Add labels and title
plt.xlabel('Space')
plt.ylabel('Characteristic Speed')

plt.xlim([0,1])


# Show the grid
plt.grid(True)

# Add a legend
plt.legend()

# Display the plot
plt.savefig("Report-PIC2/Images/Good_Speeds.png")