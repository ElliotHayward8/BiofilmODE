import matplotlib as plt
import numpy as np
from scipy.integrate import odeint

def grow_phase(Ns, Np, CA, CS):
    """
    Function of ODEs to model the growth phase of the biofilm, measuring the number of susceptible
    and percepter cells
    :param Ns: Number of susceptible cells
    :param Np: Number of percepter cells
    :param CA: Concentration of antibiotic
    :param CS: Concentration of substrate
    :return: y: vector of the two differential equations
    """

    # Calculate a and b depending on the concentration of antibiotic and substrate
    a =
    b =

    dNsdt = Growth + b*Np
    dNpdt = a*Ns

    y = [dNsdt, dNpdt]
    return y



def main():
    Ns0 = 10  # Initial number of susceptible cells
    Np0 = 0  # Initial number of percepter cells


if __name__ == '__main__':
    main()


