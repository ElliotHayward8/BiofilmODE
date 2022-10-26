import matplotlib as plt
import numpy as np
from scipy.integrate import odeint

"""
Notes:

The growth factor can likely just be derived from the doubling rate, a and b (which are alpha and beta) should 
both be dependent on the concentrations of antibiotic and substrate, choose how to define this, must also be a 
maximum value for these, choose where to define this value - could change between states?
"""

def grow_phase(Ns, Np, CA, CS0, a_max, b_max, K1, K2):
    """
    Function of ODEs to model the growth phase of the biofilm, measuring the number of susceptible
    and percepter cells
    :param Ns: Number of susceptible cells
    :param Np: Number of percepter cells
    :param CA: Concentration of antibiotic (at the moment always 0)
    :param CS0: Concentration of substrate at the start
    :param a_max: Maximum susceptible to persister switching rate
    :param b_max: Maximum persister to susceptible switching rate
    :return: y: vector of the two differential equations
    """

    eat_rate = 1 # define the eating rate of the bacteria (find value)

    # Conc. of substrate depends on constant rate and the number of susceptible cells
    CS = CS0 - (Ns * eat_rate)

    # Calculate a and b depending on the concentration of antibiotic and substrate
    a = (a_max * (1 - (CS/(CS + K1)))) + (a_max * (CA/(CA + K2)))
    b = (b_max * (CS/(CS + K1))) + (b_max * (1 - (CA/(CA + K2))))

    dNsdt = G*Ns + b*Np # ODE for the rate of change of susceptible cells
    dNpdt = a*Ns # ODE for the rate of change of persister cells

    y = [dNsdt, dNpdt]
    return y



def main():
    Ns0 = 10  # Initial number of susceptible cells
    Np0 = 0  # Initial number of percepter cells
    CS0 = 0.4 # Bulk concentration of substrate
    t = np.linspace(0, 36000, 36000) # Define the time range to solve over


if __name__ == '__main__':
    main()


