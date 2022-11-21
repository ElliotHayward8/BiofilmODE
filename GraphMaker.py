#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.integrate import odeint
import csv
from BiofilmODE import constant_switch, constant_treatment


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
    a_max, b_max, a_min, b_min = 1, 1, 0.1, 0.1  # Typical values - (1, 1)
    kmaxs, kmaxp = 10, 0.1  # Define the killing rates of both types of bacteria (10, 0.1)
    K_A = 6.4 * MIC  # Define the half-saturation constant for the killing rates of susceptible and persister cells
    K_S = 0.0035  # Define the half-saturation constant for the substrate (0.0035)

    full_sol1, t_total1, end_time1 = constant_switch(N, 5, 0.5, 0.5, CA, CS0, a_min, b_min, K1, K2, G, kmaxs,
                                                     kmaxp, K_A, K_S)

    con_sol1, con_t1, con_end1 = constant_treatment(N, 5, 2, CA, CS0, a_min, b_min, K1, K2, G, kmaxs, kmaxp, K_A, K_S)

    full_sol2, t_total2, end_time2 = constant_switch(N, 5, 0.5, 0.5, CA, CS0, a_min, b_max, K1, K2, G, kmaxs,
                                                     kmaxp, K_A, K_S)

    con_sol2, con_t2, con_end2 = constant_treatment(N, 5, 2, CA, CS0, a_min, b_max, K1, K2, G, kmaxs, kmaxp, K_A, K_S)

    full_sol3, t_total3, end_time3 = constant_switch(N, 5, 0.5, 0.5, CA, CS0, a_min, b_min, K1, K2, G, kmaxs,
                                                     kmaxp, K_A, K_S)

    con_sol3, con_t3, con_end3 = constant_treatment(N, 5, 2, CA, CS0, a_max, b_max, K1, K2, G, kmaxs, kmaxp, K_A, K_S)

    full_sol4, t_total4, end_time4 = constant_switch(N, 5, 0.5, 0.5, CA, CS0, a_min, b_min, K1, K2, G, kmaxs,
                                                     kmaxp, K_A, K_S)

    con_sol4, con_t4, con_end4 = constant_treatment(N, 5, 2, CA, CS0, a_max, b_max, K1, K2, G, kmaxs, kmaxp, K_A, K_S)

    fig, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8)) = plt.subplots(2, 4)

    ax1.set_xlabel('Time (mins)'), ax2.set_xlabel('Time (mins)'), ax3.set_xlabel('Time (mins)')
    ax4.set_xlabel('Time (mins)'), ax5.set_xlabel('Time (mins)'), ax6.set_xlabel('Time (mins)')
    ax7.set_xlabel('Time (mins)'), ax8.set_xlabel('Time (mins)'), ax1.set_ylabel('Number of cells')
    ax2.set_ylabel('Number of cells'), ax3.set_ylabel('Number of cells'), ax4.set_ylabel('Number of cells')
    ax5.set_ylabel('Number of cells'), ax6.set_ylabel('Number of cells'), ax6.set_ylabel('Number of cells')
    ax8.set_ylabel('Number of cells')
    ax1.legend(loc=1), ax2.legend(loc=1), ax3.legend(loc=1), ax4.legend(loc=1), ax5.legend(loc=1), ax6.legend(loc=1)
    ax7.legend(loc=1), ax8.legend(loc=1)
    fig.set_size_inches(18.5, 10.5)
    fig.tight_layout(pad=5)
    plt.savefig("AlphaBetaComp.pdf", bbox_inches='tight', pad_inches=0)
    plt.show()


if __name__ == '__main__':
    main()
