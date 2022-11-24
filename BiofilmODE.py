#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.integrate import odeint
import csv

"""
Notes:

The growth factor can likely just be derived from the doubling rate, a and b (which are alpha and beta) should 
both be dependent on the concentrations of antibiotic and substrate, choose how to define this, must also be a 
maximum value for these, choose where to define this value - could change between states?
"""


def grow_phase(N, t, CA, CS0, a_max, b_max, K1, K2, G, kmaxs, kmaxp, K_A, K_S):
    """
    Function of ODEs to model the growth phase of the biofilm, measuring the number of susceptible
    and percepter cells
    :param N: Number of persister and susceptible cells
    :param CA: Concentration of antibiotic (always 0 when growing, gL^-1)
    :param CS0: Concentration of substrate at the start (gL^-1)
    :param a_max: Maximum susceptible to persister switching rate (h^-1)
    :param b_max: Maximum persister to susceptible switching rate (h^-1)
    :param G: Growth rate of the bacteria (h^-1)
    :param K1: Half-saturation constant for substrate-dependent switching
    :param K2: Half-saturation constant for substrate-dependent switching
    :param kmaxs: Maximum killing rate of susceptible cells (h^-1)
    :param kmaxp: Maximum killing rate of persister cells (h^-1)
    :param K_A: Half-saturation constant so when CA = MIC, killing rate = mu_max = 1.25 h^-1
    :param K_S: Half-saturation constant for the substrate S
    :return: y: vector of the two differential equations
    """
    Ns, Np = N
    eat_rate = 0.00001  # define the eating rate of the bacteria (how Ns affects CS)

    # Concentration of substrate depends on constant rate and the number of susceptible cells
    CS = CS0 - (Ns * eat_rate)

    # Growth rate depends on the concentration of substrate
    G = G * (CS / (CS + K_S))

    # Calculate a and b depending on the concentration of antibiotic and substrate

    a = (a_max * (1 - (CS / (CS + K1)))) + (a_max * (CA / (CA + K2)))
    b = 0.5 * (b_max * (CS / (CS + K1))) + (b_max * (1 - (CA / (CA + K2))))
    """
    elif switch_type == 2:  # Substrate dependent switching
        a = a_max * (1 - (CS / (CS + K1)))
        b = 0.1 * b_max * (CS / (CS + K1))
    elif switch_type == 3:  # Antibiotic dependent switching
        a = a_max * (CA / (CA + K2))
        b = b_max * (1 - (CA / (CA + K2)))
    else:
        raise ValueError('switch_type must be 1, 2 or 3')
    """
    dNsdt = G * Ns + b * Np - a * Ns - Ns * (kmaxs * (CA / (CA + K_A)))  # Rate of change of susceptible cells
    dNpdt = a * Ns - b * Np - Np * (kmaxp * (CA / (CA + K_A)))  # Rate of change of persister cells

    if CA != 0:
        if dNsdt > 0:
            print(G * Ns, b * Np, a * Ns, kmaxs * (CA / (CA + K_A)))
    return [dNsdt, dNpdt]


def run_biofilm(N, T, CA, CS0, a_max, b_max, K1, K2, G, kmaxs, kmaxp, K_A, K_S):
    """
    Function to run the ODEs for time T with the inputted parameters
    :param N: Number of persister and susceptible cells
    :param CA: Concentration of antibiotic (always 0 when growing, gL^-1)
    :param CS0: Concentration of substrate at the start (gL^-1)
    :param a_max: Maximum susceptible to persister switching rate (h^-1)
    :param b_max: Maximum persister to susceptible switching rate (h^-1)
    :param G: Growth rate of the bacteria (h^-1)
    :param K1: Half-saturation constant for substrate-dependent switching
    :param K2: Half-saturation constant for substrate-dependent switching
    :param kmaxs: Maximum killing rate of susceptible cells (h^-1)
    :param kmaxp: Maximum killing rate of persister cells (h^-1)
    :param K_A: Half-saturation constant so when CA = MIC, killing rate = mu_max = 1.25 h^-1
    :param K_S: Half-saturation constant for the substrate S
    :return: y: vector of the two differential equations
    """
    t_run = np.linspace(0, T, T * 3600)

    sol = odeint(grow_phase, N, t_run, args=(CA, CS0, a_max, b_max, K1, K2, G, kmaxs, kmaxp, K_A, K_S))

    return sol


def constant_treatment(N, T_grow, T1, CA, CS0, a_max, b_max, K1, K2, G, kmaxs, kmaxp, K_A, K_S):
    """
    A biofilm treatment strategy where the antibiotic is applied constantly
    :param T_grow: The time the biofilm grows before it is first treated
    :param T1: Length of each iteration
    """
    end_val = 0.5
    t_grow = np.linspace(0, T_grow, int(3600 * T_grow))
    t_treat = np.linspace(0, T1, int(3600 * T1))

    sol = odeint(grow_phase, N, t_grow, args=(0, CS0, a_max, b_max, K1, K2, G, kmaxs, kmaxp, K_A, K_S))
    full_sol, t_total, treat_time = sol, t_grow, 0
    final_N = [full_sol[-1, 0], full_sol[-1, 1]]

    total_time = T_grow
    check = 0
    while sum(final_N) > end_val and total_time < 250:
        treat_sol = odeint(grow_phase, final_N, t_treat, args=(CA, CS0, a_max, b_max, K1, K2, G, kmaxs, kmaxp,
                                                               K_A, K_S))
        final_N = [treat_sol[-1, 0], treat_sol[-1, 1]]
        for i in range(len(treat_sol)):
            if check == 0:
                if (treat_sol[i, 0] + treat_sol[i, 1]) < end_val:
                    check, treat_sol, t_treat, end_time = 1, treat_sol[:i+1], t_treat[:i+1], t_grow[i] + total_time
                    treat_time += t_grow[i]

        full_sol, t_total = np.concatenate((full_sol, treat_sol)), np.append(t_total, t_treat + total_time)
        total_time = round(total_time + T1, 1)
        if check == 0:
            treat_time += T1

    return full_sol, t_total, end_time, treat_time


def constant_switch(N, T_grow, T1, T2, CA, CS0, a_max, b_max, K1, K2, G, kmaxs, kmaxp, K_A, K_S):
    """
    A biofilm treatment strategy where the antibiotic is applied and removed for constant time periods (keep T">0.1)
    :param T_grow: The time the biofilm grows before it is first treated
    :param T1: Length of each treatment
    :param T2: Length between treatments
    """
    end_val, t_grow, t_treat = 0.5, np.linspace(0, T_grow, int(3600 * T_grow)), np.linspace(0, T1, int(3600 * T1))
    t_regrow = np.linspace(0, T2, int(3600 * T2))

    sol = odeint(grow_phase, N, t_grow, args=(0, CS0, a_max, b_max, K1, K2, G, kmaxs, kmaxp, K_A, K_S))
    full_sol, t_total, treat_time = sol, t_grow, 0
    final_N = [full_sol[-1, 0], full_sol[-1, 1]]

    total_time, check = T_grow, 0
    while sum(final_N) > end_val and total_time < 150:
        treat_sol = odeint(grow_phase, final_N, t_treat, args=(CA, CS0, a_max, b_max, K1, K2, G, kmaxs, kmaxp,
                                                               K_A, K_S))
        final_N = [treat_sol[-1, 0], treat_sol[-1, 1]]
        for i in range(len(treat_sol)):
            if check == 0:
                if (treat_sol[i, 0] + treat_sol[i, 1]) < end_val:
                    check = 1
                    treat_sol, t_treat, end_time = treat_sol[:i+1], t_treat[:i+1], t_grow[i] + total_time
                    treat_time += t_grow[i]

        full_sol, t_total = np.concatenate((full_sol, treat_sol)), np.append(t_total, t_treat + total_time)

        if check == 0:
            treat_time += T1
            total_time += T1
            regrow_sol = odeint(grow_phase, final_N, t_regrow, args=(0, CS0, a_max, b_max, K1, K2, G, kmaxs, kmaxp,
                                                                     K_A, K_S))
            final_N = [regrow_sol[-1, 0], regrow_sol[-1, 1]]
            full_sol, t_total = np.concatenate((full_sol, regrow_sol)), np.append(t_total, t_regrow + total_time)
            total_time += T2

    # Ensures end_time has a value
    if total_time >= 150:
        end_time = total_time

    return full_sol, t_total, end_time, treat_time


def simple_constant_switch(N, T_grow, T1, T2, CA, CS0, a_max, b_max, K1, K2, G, kmaxs, kmaxp, K_A, K_S):
    """
    A biofilm treatment strategy where the antibiotic is applied and removed for constant time periods (keep T">0.1)
    :param T_grow: The time the biofilm grows before it is first treated
    :param T1: Length of each treatment
    :param T2: Length between treatments
    """
    end_val, t_grow, t_treat = 0.5, np.linspace(0, T_grow, int(3600 * T_grow)), np.linspace(0, T1, int(3600 * T1))
    t_regrow = np.linspace(0, T2, int(3600 * T2))

    sol = odeint(grow_phase, N, t_grow, args=(0, CS0, a_max, b_max, K1, K2, G, kmaxs, kmaxp, K_A, K_S))
    treat_total = 0
    final_N = [sol[-1, 0], sol[-1, 1]]

    total_time, check = T_grow, 0
    while sum(final_N) > end_val and treat_total < 25:
        treat_sol = odeint(grow_phase, final_N, t_treat, args=(CA, CS0, a_max, b_max, K1, K2, G, kmaxs, kmaxp,
                                                               K_A, K_S))
        treat_total += T1
        final_N = [treat_sol[-1, 0], treat_sol[-1, 1]]
        for i in range(len(treat_sol)):
            if check == 0:
                if (treat_sol[i, 0] + treat_sol[i, 1]) < end_val:
                    check = 1
                    end_time = t_grow[i] + total_time
                    treat_total += t_treat[i] - T1

        if check == 0:
            total_time += T1
            total_time += T2
            regrow_sol = odeint(grow_phase, final_N, t_regrow, args=(0, CS0, a_max, b_max, K1, K2, G, kmaxs, kmaxp,
                                                                     K_A, K_S))
            final_N = [regrow_sol[-1, 0], regrow_sol[-1, 1]]

    # Ensures end_time has a value
    if treat_total >= 25:
        end_time = total_time

    return end_time, treat_total


def param_scan(N, T_grow, T1s, T2s, CA, CS0, a_maxs, b_maxs, K1, K2, G, kmaxs, kmaxp, K_A, K_S):
    """
    A function to scan through potential values to see which treatment times are best for differing a and b values
    :param T1s: List of treatment times used
    :param T2s: List of recovery times used
    :param a_maxs: List of a_max values used
    :param b_maxs: List of b_max values used
    :return: For each combination of a_max and b_max it provides all the treatment times for each T1-T2 combination
    """
    a_max_results, b_max_results = [], []
    T1_results, T2_results, end_time_results = [], [], []

    for a_max in a_maxs:
        for b_max in b_maxs:
            print('a = ' + str(a_max) + '        b = ' + str(b_max))  # For tracking the progress of the scan
            for T1 in T1s:
                for T2 in T2s:
                    full_sol, t_total, end_time = constant_switch(N, T_grow, T1, T2, CA, CS0, a_max, b_max, K1, K2, G,
                                                                  kmaxs, kmaxp, K_A, K_S)

                    if T1_results:
                        a_max_results.append(a_max), b_max_results.append(b_max), T1_results.append(T1)
                        T2_results.append(T2), end_time_results.append(end_time)
                    else:
                        a_max_results, b_max_results, T1_results, T2_results = [a_max], [b_max], [T1], [T2]
                        end_time_results = [end_time]

    return a_max_results, b_max_results, T1_results, T2_results, end_time_results


def best_param_scan(N, T_grow, T1s, T2s, CA, CS0, a_maxs, b_maxs, K1, K2, G, kmaxs, kmaxp, K_A, K_S):
    """
    A function to scan through potential values to see which treatment times are best for differing a and b values
    :param T1s: List of treatment times used
    :param T2s: List of recovery times used
    :param a_maxs: List of a_max values used
    :param b_maxs: List of b_max values used
    :return: For each combination of a_max and b_max it the optimal T1-T2 combination
    """
    a_max_vals, b_max_vals = [], []
    best_T1, best_T2, best_end_time, best_treat_time = [], [], [], []
    for a_max in a_maxs:
        for b_max in b_maxs:
            print('a = ' + str(a_max) + '        b = ' + str(b_max))  # For tracking the progress of the scan
            a_max_vals.append(a_max), b_max_vals.append(b_max)
            best_T1.append(0), best_T2.append(0), best_end_time.append(500), best_treat_time.append(500)
            for T1 in T1s:
                for T2 in T2s:
                    end_time, treat_total = simple_constant_switch(N, T_grow, T1, T2, CA, CS0, a_max, b_max, K1, K2, G,
                                                                   kmaxs, kmaxp, K_A, K_S)
                    # Should I be tracking treat_total or end_time????
                    if treat_total < best_treat_time[-1]:
                        best_treat_time[-1], best_end_time[-1], best_T1[-1], best_T2[-1] = treat_total, end_time, T1, T2

    fieldnames = ['a_max', 'b_max', 'T1', 'T2', 't_total', 'treat_total']
    with open('output.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(i for i in fieldnames)
        for i in range(len(a_max_vals)):
            writer.writerow([a_max_vals[i], b_max_vals[i], best_T1[i], best_T2[i], best_end_time[i],
                             best_treat_time[i]])

    return a_max_vals, b_max_vals, best_T1, best_T2, best_end_time


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

    """
    grow_hours, treat_hours, recov_hours = 8, 2, 5  # Define the number of hours for each stage of simulation
    t_grow, t_treat = np.linspace(0, grow_hours, grow_hours * 3600), np.linspace(0, treat_hours, treat_hours * 3600)
    t_recov = np.linspace(0, recov_hours, recov_hours * 3600)
    
    sol = odeint(grow_phase, N, t_grow, args=(0, CS0, a_max, b_max, K1, K2, G, kmaxs, kmaxp, K_A, K_S))
    final_N = [sol[-1, 0], sol[-1, 1]]

    sol2 = odeint(grow_phase, final_N, t_treat, args=(CA, CS0, a_max, b_max, K1, K2, G, kmaxs, kmaxp, K_A, K_S))
    final_N = [sol2[-1, 0], sol2[-1, 1]]

    sol3 = odeint(grow_phase, final_N, t_recov, args=(0, CS0, a_max, b_max, K1, K2, G, kmaxs, kmaxp, K_A, K_S))

    # Combine the solutions to produce one graph
    full_sol = np.concatenate((sol, sol2, sol3))
    t_total = np.concatenate((t_grow, t_treat + grow_hours, t_recov + grow_hours + treat_hours))

    plt.plot(t_total, full_sol[:, 0], 'b', label='Ns'), plt.plot(t_total, full_sol[:, 1], 'g', label='Np')
    plt.legend(loc='best'), plt.xlabel('t'), plt.ylabel('Number of cells')
    plt.grid()
    plt.show()
    """

    # Define parameters for the parameter scan
    # a_maxs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    # b_maxs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    b_s = np.linspace(1.0, 0.1, 40)
    a_s = np.linspace(1.0, 0.1, 40)  # np.linspace(0.1, 1.0, 2)
    # T1s = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    # T2s = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    T1s, T2s = np.linspace(1/12, 3, 36), np.linspace(1/12, 3, 36)

    """
    a_max_results, b_max_results, T1_results, T2_results, end_time_results = param_scan(N, 5, T1s, T2s, CA, CS0, a_s,
                                                                                        b_s, K1, K2, G, kmaxs, kmaxp,
                                                                                        K_A, K_S)
    """

    a_max_results, b_max_results, T1_results, T2_results, end_time_results = best_param_scan(N, 5, T1s, T2s, CA, CS0,
                                                                                             a_s, b_s, K1, K2, G,
                                                                                             kmaxs, kmaxp, K_A, K_S)

    # T1_results, T2_results = np.asarray(T1_results), np.asarray(T2_results)
    # T1_results, T2_results = np.reshape(T1_results, (len(b_s), len(a_s))).T, np.reshape(T2_results,
    #                                                                                     (len(b_s), len(a_s))).T
    #
    # fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(8, 4))
    # tkx, x0 = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], -0.5
    # tky = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
    # x_len, y_len = len(a_s) + x0, len(b_s) + x0
    # im1, im2 = ax1.imshow(T1_results*60, cmap='viridis'), ax2.imshow(T2_results*60, cmap='viridis')
    #
    # ax1.set_xticks(np.linspace(x0, x_len, 10), labels=tkx), ax2.set_xticks(np.linspace(x0, x_len, 10), labels=tkx)
    # ax1.set_yticks(np.linspace(x0, y_len, 10), labels=tky), ax2.set_yticks(np.linspace(x0, y_len, 10), labels=tky)
    # ax1.set_xticks(np.arange(len(a_s))-0.5, minor=True), ax2.set_xticks(np.arange(len(a_s))-0.5, minor=True)
    # ax1.set_yticks(np.arange(len(b_s))-0.5, minor=True), ax2.set_yticks(np.arange(len(b_s))-0.5, minor=True)
    # ax1.grid(which='minor', color='k', linestyle='-', linewidth=0.3), ax2.grid(which='minor', color='k', linestyle='-',
    #                                                                            linewidth=0.3)
    # ax1.tick_params(axis='both', which='both', length=0), ax2.tick_params(axis='both', which='both', length=0)
    #
    # divider1, divider2 = make_axes_locatable(ax1), make_axes_locatable(ax2)
    # cax1, cax2 = divider1.append_axes("right", size="5%", pad=0.2), divider2.append_axes("right", size="5%", pad=0.5)
    # cbar1, cbar2 = plt.colorbar(im1, cax=cax1), plt.colorbar(im2, cax=cax2)
    # cbar1.set_label(label='Optimised treatment duration (min)', size=10)
    # cbar2.set_label(label='Optimised treatment duration (min)', size=10)
    # ax1.set_xlabel('$a_{max}$'), ax2.set_xlabel('$a_{max}$'), ax1.set_ylabel('$b_{max}$'), ax2.set_ylabel('$b_{max}$')
    # fig.tight_layout(pad=5)
    # # plt.savefig("ParameterScan.pdf", bbox_inches='tight', pad_inches=0)
    # plt.show()
    """
    full_sol, t_total, end_time1 = constant_switch(N, 5, 0.1, 2, CA, CS0, a_max, b_max, K1, K2, G, kmaxs,
                                                   kmaxp, K_A, K_S)

    end_time2, treat_total = simple_constant_switch(N, 5, 0.1, 2, CA, CS0, a_max, b_max, K1, K2, G, kmaxs, kmaxp, K_A,
                                                    K_S)
    """
    #
    # con_sol, con_t, con_end = constant_treatment(N, 5, 2, CA, CS0, a_max, b_max, K1, K2, G, kmaxs, kmaxp, K_A, K_S)
    #
    # plt.plot(con_t, con_sol[:, 0], 'r', label='Constant Ns'), plt.plot(con_t, con_sol[:, 1], 'm', label='Constant Np')
    # plt.plot(t_total, full_sol[:, 0], 'b', label='Ns'), plt.plot(t_total, full_sol[:, 1], 'g', label='Np')
    # plt.legend(loc='best'), plt.xlabel('Time'), plt.ylabel('Number of cells')
    # plt.axvline(x=end_time, color='c')
    # plt.axvline(x=con_end, color='y')
    # plt.show()


if __name__ == '__main__':
    main()
