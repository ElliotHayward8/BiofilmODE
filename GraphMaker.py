#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.integrate import odeint
import csv




def main():
    """
        Main function where parameter values are defined and the process is run
        """
    MIC = 0.000005  # Minimal Inhibitory Concentration of the antibiotic
    CA = 1000 * MIC  # Define the concentration of antibiotic present during the treatment phase
    Ns0, Np0 = 10, 0  # Initial number of each type of cell
    N = [Ns0, Np0]
    CS0 = 0.4  # Bulk concentration of substrate
    """
    Change G to reflect slower growth dependent on NetLogo model compared to this model
    """
    G = 1.278  # Growing rate of the bacteria (Need to equate to mu_max), include the substrate concentration
    K1, K2 = 0.0035, 0.000005  # Half saturation constant for the substrate and antibiotic respectively
    a_max, b_max = 1, 1  # Typical values - (1, 1)
    kmaxs, kmaxp = 10, 0.1  # Define the killing rates of both types of bacteria (10, 0.1)
    K_A = 6.4 * MIC  # Define the half-saturation constant for the killing rates of susceptible and persister cells
    K_S = 0.0035  # Define the half-saturation constant for the substrate (0.0035)






if __name__ == '__main__':
    main()