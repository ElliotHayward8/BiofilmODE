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

    full_sol1, t_total1, end_time1, T_time1 = constant_switch(N, 5, 0.5, 0.5, CA, CS0, a_min, b_min, K1, K2, G, kmaxs,
                                                              kmaxp, K_A, K_S)

    con_sol1, con_t1, con_end1, con_T_time1 = constant_treatment(N, 5, 2, CA, CS0, a_min, b_min, K1, K2, G, kmaxs,
                                                                 kmaxp, K_A, K_S)

    full_sol2, t_total2, end_time2, T_time2 = constant_switch(N, 5, 0.5, 0.5, CA, CS0, a_min, b_max, K1, K2, G, kmaxs,
                                                              kmaxp, K_A, K_S)

    con_sol2, con_t2, con_end2, con_T_time2 = constant_treatment(N, 5, 2, CA, CS0, a_min, b_max, K1, K2, G, kmaxs,
                                                                 kmaxp, K_A, K_S)

    full_sol3, t_total3, end_time3, T_time3 = constant_switch(N, 5, 0.5, 0.5, CA, CS0, a_max, b_min, K1, K2, G, kmaxs,
                                                              kmaxp, K_A, K_S)

    con_sol3, con_t3, con_end3, con_T_time3 = constant_treatment(N, 5, 2, CA, CS0, a_max, b_min, K1, K2, G, kmaxs,
                                                                 kmaxp, K_A, K_S)

    full_sol4, t_total4, end_time4, T_time4 = constant_switch(N, 5, 0.5, 0.5, CA, CS0, a_max, b_max, K1, K2, G, kmaxs,
                                                              kmaxp, K_A, K_S)

    con_sol4, con_t4, con_end4, con_T_time4 = constant_treatment(N, 5, 2, CA, CS0, a_max, b_max, K1, K2, G, kmaxs,
                                                                 kmaxp, K_A, K_S)

    # Convert the time values from hours into minutes
    con_t4, con_t3, con_t2, con_t1 = con_t4*60, con_t3*60, con_t2*60, con_t1*60
    t_total1, t_total2, t_total3, t_total4 = t_total1*60, t_total2*60, t_total3*60, t_total4*60
    T_time1, T_time2, T_time3, T_time4 = T_time1*60, T_time2*60, T_time3*60, T_time4*60
    con_T_time1, con_T_time2, con_T_time3, con_T_time4 = con_T_time1*60, con_T_time2*60, con_T_time3*60, con_T_time4*60

    text_properties = {'fontsize': 16}
    leg_properties = {'size': 13}
    top_y, bottom_y = 6000, -15

    fig, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8)) = plt.subplots(2, 4)
    ax1.set_ylim(bottom=bottom_y, top=top_y), ax2.set_ylim(bottom=bottom_y, top=top_y)
    ax3.set_ylim(bottom=bottom_y, top=top_y), ax4.set_ylim(bottom=bottom_y, top=top_y)
    ax5.set_ylim(bottom=bottom_y, top=top_y), ax6.set_ylim(bottom=bottom_y, top=top_y)
    ax7.set_ylim(bottom=bottom_y, top=top_y), ax8.set_ylim(bottom=bottom_y, top=top_y)
    ax1.text(1200, 3500, 'T-time = ' + str(int(T_time1)), fontdict=text_properties)
    ax2.text(350, 3500, 'T-time = ' + str(int(T_time2)), fontdict=text_properties)
    ax3.text(1500, 3500, 'T-time = ' + str(int(T_time3)), fontdict=text_properties)
    ax4.text(400, 3500, 'T-time = ' + str(int(T_time4)), fontdict=text_properties)
    ax5.text(1175, 3500, 'T-time = ' + str(int(con_T_time1)), fontdict=text_properties)
    ax6.text(500, 3500, 'T-time = ' + str(int(con_T_time2)), fontdict=text_properties)
    ax7.text(1750, 3500, 'T-time = ' + str(int(con_T_time3)), fontdict=text_properties)
    ax8.text(550, 3500, 'T-time = ' + str(int(con_T_time4)), fontdict=text_properties)
    ax1.plot(t_total1, full_sol1[:, 0], 'b', label='$N_s$'), ax1.plot(t_total1, full_sol1[:, 1], 'r', label='$N_p$')
    ax5.plot(con_t1, con_sol1[:, 0], 'b', label='$N_s$'), ax5.plot(con_t1, con_sol1[:, 1], 'r', label='$N_p$')
    ax2.plot(t_total2, full_sol2[:, 0], 'b', label='$N_s$'), ax2.plot(t_total2, full_sol2[:, 1], 'r', label='$N_p$')
    ax6.plot(con_t2, con_sol2[:, 0], 'b', label='$N_s$'), ax6.plot(con_t2, con_sol2[:, 1], 'r', label='$N_p$')
    ax3.plot(t_total3, full_sol3[:, 0], 'b', label='$N_s$'), ax3.plot(t_total3, full_sol3[:, 1], 'r', label='$N_p$')
    ax7.plot(con_t3, con_sol3[:, 0], 'b', label='$N_s$'), ax7.plot(con_t3, con_sol3[:, 1], 'r', label='$N_p$')
    ax4.plot(t_total4, full_sol4[:, 0], 'b', label='$N_s$'), ax4.plot(t_total4, full_sol4[:, 1], 'r', label='$N_p$')
    ax8.plot(con_t4, con_sol4[:, 0], 'b', label='$N_s$'), ax8.plot(con_t4, con_sol4[:, 1], 'r', label='$N_p$')

    ax1.set_xlabel('Time (mins)', fontsize=13), ax2.set_xlabel('Time (mins)', fontsize=13)
    ax3.set_xlabel('Time (mins)', fontsize=13), ax4.set_xlabel('Time (mins)', fontsize=13)
    ax5.set_xlabel('Time (mins)', fontsize=13), ax6.set_xlabel('Time (mins)', fontsize=13)
    ax7.set_xlabel('Time (mins)', fontsize=13), ax8.set_xlabel('Time (mins)', fontsize=13)
    ax1.set_ylabel('Number of cells', fontsize=13), ax2.set_ylabel('Number of cells', fontsize=13)
    ax3.set_ylabel('Number of cells', fontsize=13), ax4.set_ylabel('Number of cells', fontsize=13)
    ax5.set_ylabel('Number of cells', fontsize=13), ax6.set_ylabel('Number of cells', fontsize=13)
    ax7.set_ylabel('Number of cells', fontsize=13), ax8.set_ylabel('Number of cells', fontsize=13)

    ax1.legend(loc=1, prop=leg_properties), ax2.legend(loc=1, prop=leg_properties)
    ax3.legend(loc=1, prop=leg_properties), ax4.legend(loc=1, prop=leg_properties)
    ax5.legend(loc=1, prop=leg_properties), ax6.legend(loc=1, prop=leg_properties)
    ax7.legend(loc=1, prop=leg_properties), ax8.legend(loc=1, prop=leg_properties)
    fig.set_size_inches(14, 9)
    fig.tight_layout(pad=2)
    plt.savefig("AlphaBetaComp.pdf", bbox_inches='tight', pad_inches=0)
    plt.show()

    fig, (ax1, ax2) = plt.subplots(1, 2)

    ax1.plot(t_total4, full_sol4[:, 0], 'b', label='No. of susceptibles')
    ax1.plot(t_total4, full_sol4[:, 1], 'r', label='No. of persisters')
    ax2.plot(con_t4, con_sol4[:, 0], 'b', label='No. of susceptibles')
    ax2.plot(con_t4, con_sol4[:, 1], 'r', label='No. of persisters')
    # print(end_time4*60, con_end4*60)
    ax2.text(800, 3500, '$t_T$ = ' + str(int(con_T_time4)), fontdict=text_properties)
    ax1.text(600, 3500, '$t_T$ = ' + str(int(T_time4)), fontdict=text_properties)
    ax1.tick_params(axis='both', labelsize=16), ax2.tick_params(axis='both', labelsize=16)
    ax1.set_title('Periodic Treatment', size=16)
    ax2.set_title('Continuous Treatment', size=16)
    ax1.legend(loc=1, prop=leg_properties), ax2.legend(loc=1, prop=leg_properties)
    ax1.set_xlabel('Time (mins)', fontsize=16), ax2.set_xlabel('Time (mins)', fontsize=16)
    ax1.set_ylabel('Number of cells', fontsize=16), ax2.set_ylabel('Number of cells', fontsize=16)
    fig.set_size_inches(12, 6)
    fig.tight_layout(pad=2)
    plt.savefig("2figcomp.pdf")
    plt.show()


if __name__ == '__main__':
    main()
