import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

"""
Notes:

The growth factor can likely just be derived from the doubling rate, a and b (which are alpha and beta) should 
both be dependent on the concentrations of antibiotic and substrate, choose how to define this, must also be a 
maximum value for these, choose where to define this value - could change between states?
"""


def grow_phase(N, t, CA, CS0, a_max, b_max, K1, K2, G, kmaxs, kmaxp, K_A):
    """
    Function of ODEs to model the growth phase of the biofilm, measuring the number of susceptible
    and percepter cells
    :param N: Number of persister and susceptible cells
    :param CA: Concentration of antibiotic (always 0 when growing, gL^-1)
    :param CS0: Concentration of substrate at the start (gL^-1)
    :param a_max: Maximum susceptible to persister switching rate (h^-1)
    :param b_max: Maximum persister to susceptible switching rate (h^-1)
    :param G: Growth rate of the bacteria (h^-1)
    :param kmaxs: Maximum killing rate of susceptible cells (h^-1)
    :param kmaxp: Maximum killing rate of persister cells (h^-1)
    :param K_A: Half saturation constant so when CA = MIC, killing rate = mu_max = 1.25 h^-1
    :return: y: vector of the two differential equations
    """
    Ns, Np = N
    eat_rate = 0.0001  # define the eating rate of the bacteria (how Ns affects CS)

    # Conc. of substrate depends on constant rate and the number of susceptible cells
    CS = CS0 - (Ns * eat_rate)

    # Calculate a and b depending on the concentration of antibiotic and substrate
    a = (a_max * (1 - (CS/(CS + K1)))) + (a_max * (CA/(CA + K2)))
    b = 0.5 * (0.1 * b_max * (CS/(CS + K1))) + (b_max * (1 - (CA/(CA + K2))))  # 0.1 to match the default NetLogo model
    dNsdt = G*Ns + b*Np - a*Ns - kmaxs * CA/(CA + K_A)  # ODE for the rate of change of susceptible cells
    dNpdt = a*Ns - b*Np - kmaxp * CA/(CA + K_A)  # ODE for the rate of change of persister cells

    y = [dNsdt, dNpdt]

    return y


def main():
    """
    Main function where parameter values are defined and the process is run
    """
    Ns0, Np0 = 10, 0  # Initial number of each type of cell
    N = [Ns0, Np0]
    CS0 = 0.4  # Bulk concentration of substrate
    G = 1.25  # Growing rate of the bacteria (Need to equate to mu_max)
    K1, K2 = 0.0035, 0.000005  # Half saturation rate for the substrate and antibiotic respectively
    t = np.linspace(0, 8, 36000)  # Define the time range to solve over
    a_max, b_max = 1, 1
    kmaxs, kmaxp = 10, 0.1  # Define the killing rates of both types of bacteria
    K_A = 6.4  # Define the half-saturation constant for the killing rates of susceptible and persister cells

    sol = odeint(grow_phase, N, t, args=(0, CS0, a_max, b_max, K1, K2, G, kmaxs, kmaxp, K_A))

    plt.plot(t, sol[:, 0], 'b', label='Ns')
    plt.plot(t, sol[:, 1], 'g', label='Np')
    plt.legend(loc='best')
    plt.xlabel('t')
    plt.show()


if __name__ == '__main__':
    main()


