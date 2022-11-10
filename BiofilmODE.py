import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

"""
Notes:

The growth factor can likely just be derived from the doubling rate, a and b (which are alpha and beta) should 
both be dependent on the concentrations of antibiotic and substrate, choose how to define this, must also be a 
maximum value for these, choose where to define this value - could change between states?
"""


def grow_phase(N, t, CA, CS0, a_max, b_max, K1, K2, G, kmaxs, kmaxp, K_A, K_S, switch_type):
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
    eat_rate = 0.00003  # define the eating rate of the bacteria (how Ns affects CS)

    # Concentration of substrate depends on constant rate and the number of susceptible cells
    CS = CS0 - (Ns * eat_rate)

    # Growth rate depends on the concentration of substrate
    G = G * (CS / (CS + K_S))

    # Calculate a and b depending on the concentration of antibiotic and substrate
    if switch_type == 1:  # Combination dependent switching
        a = (a_max * (1 - (CS / (CS + K1)))) + (a_max * (CA / (CA + K2)))
        b = 0.5 * (0.1 * b_max * (CS / (CS + K1))) + (
                    b_max * (1 - (CA / (CA + K2))))  # 0.1 is to match the default model
    elif switch_type == 2:  # Substrate dependent switching
        a = a_max * (1 - (CS / (CS + K1)))
        b = 0.1 * b_max * (CS / (CS + K1))
    elif switch_type == 3:  # Antibiotic dependent switching
        a = a_max * (CA / (CA + K2))
        b = b_max * (1 - (CA / (CA + K2)))
    else:
        raise ValueError('switch_type must be 1, 2 or 3')

    dNsdt = G * Ns + b * Np - a * Ns - Ns * (
                kmaxs * (CA / (CA + K_A)))  # ODE for the rate of change of susceptible cells
    dNpdt = a * Ns - b * Np - Np * (kmaxp * (CA / (CA + K_A)))  # ODE for the rate of change of persister cells

    y = [dNsdt, dNpdt]
    if CA != 0:
        if dNsdt > 0:
            print(G * Ns, b * Np, a * Ns, kmaxs * (CA / (CA + K_A)))
    return y


def run_biofilm(N, T, CA, CS0, a_max, b_max, K1, K2, G, kmaxs, kmaxp, K_A, K_S, switch_type):
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

    sol = odeint(grow_phase, N, t_run, args=(CA, CS0, a_max, b_max, K1, K2, G, kmaxs, kmaxp, K_A, K_S, switch_type))

    return sol


def consant_treatment(N, T_grow, T1, CA, CS0, a_max, b_max, K1, K2, G, kmaxs, kmaxp, K_A, K_S, switch_type):
    """
    A biofilm treatment strategy where the antibiotic is applied constantly
    :param T_grow: The time the biofilm grows before it is first treated
    :param T1: Length of each iteration
    """
    end_val = 0.5
    t_grow = np.linspace(0, T_grow, int(3600 * T_grow))
    t_treat = np.linspace(0, T1, int(3600 * T1))

    sol = odeint(grow_phase, N, t_grow, args=(0, CS0, a_max, b_max, K1, K2, G, kmaxs, kmaxp, K_A, K_S, switch_type))
    full_sol, t_total = sol, t_grow
    final_N = [full_sol[-1, 0], full_sol[-1, 1]]

    total_time = T_grow
    check = 0
    while sum(final_N) > end_val:
        treat_sol = odeint(grow_phase, final_N, t_treat, args=(CA, CS0, a_max, b_max, K1, K2, G, kmaxs, kmaxp,
                                                               K_A, K_S, switch_type))

        final_N = [treat_sol[-1, 0], treat_sol[-1, 1]]
        for i in range(len(treat_sol)):
            if check == 0:
                if (treat_sol[i, 0] + treat_sol[i, 1]) < end_val:
                    check = 1
                    print(i, t_total[i])
                    treat_sol, t_treat = treat_sol[:i+1], t_treat[:i+1]
                    end_time = t_total[i] + total_time

        full_sol, t_total = np.concatenate((full_sol, treat_sol)), np.append(t_total, t_treat + total_time)
        total_time = round(total_time + T1, 1)

    return full_sol, t_total, end_time


def constant_switch(N, T_grow, T1, T2, CA, CS0, a_max, b_max, K1, K2, G, kmaxs, kmaxp, K_A, K_S, switch_type):
    """
    A biofilm treatment strategy where the antibiotic is applied and removed for constant time periods (keep T">0.1)
    :param T_grow: The time the biofilm grows before it is first treated
    :param T1: Length of each treatment
    :param T2: Length between treatments
    """
    end_val = 0.5
    t_grow = np.linspace(0, T_grow, int(3600 * T_grow))
    t_treat = np.linspace(0, T1, int(3600 * T1))
    t_regrow = np.linspace(0, T2, int(3600 * T2))

    sol = odeint(grow_phase, N, t_grow, args=(0, CS0, a_max, b_max, K1, K2, G, kmaxs, kmaxp, K_A, K_S, switch_type))
    full_sol, t_total = sol, t_grow
    final_N = [full_sol[-1, 0], full_sol[-1, 1]]

    total_time = T_grow
    check = 0
    while sum(final_N) > end_val:
        treat_sol = odeint(grow_phase, final_N, t_treat, args=(CA, CS0, a_max, b_max, K1, K2, G, kmaxs, kmaxp,
                                                               K_A, K_S, switch_type))
        final_N = [treat_sol[-1, 0], treat_sol[-1, 1]]
        for i in range(len(treat_sol)):
            if check == 0:
                if (treat_sol[i, 0] + treat_sol[i, 1]) < end_val:
                    check = 1
                    print(i, t_total[i])
                    treat_sol, t_treat = treat_sol[:i+1], t_treat[:i+1]
                    end_time = t_total[i] + total_time

        full_sol, t_total = np.concatenate((full_sol, treat_sol)), np.append(t_total, t_treat + total_time)
        total_time = round(total_time + T1, 1)
        regrow_sol = odeint(grow_phase, final_N, t_regrow, args=(0, CS0, a_max, b_max, K1, K2, G, kmaxs, kmaxp,
                                                                 K_A, K_S, switch_type))

        final_N = [regrow_sol[-1, 0], regrow_sol[-1, 1]]
        if check == 0:
            full_sol, t_total = np.concatenate((full_sol, regrow_sol)), np.append(t_total, t_regrow + total_time)
            total_time = round(total_time + T2, 1)

    return full_sol, t_total, end_time


def main():
    """
    Main function where parameter values are defined and the process is run
    """
    MIC = 0.000005  # Minimal Inhibitory Concentration of the antibiotic
    CA = 1000 * MIC  # Define the concentration of antibiotic present during the treatment phase
    Ns0, Np0 = 10, 0  # Initial number of each type of cell
    N = [Ns0, Np0]
    CS0 = 0.4  # Bulk concentration of substrate
    G = 1.278  # Growing rate of the bacteria (Need to equate to mu_max), include the substrate concentration
    K1, K2 = 0.0035, 0.000005  # Half saturation constant for the substrate and antibiotic respectively
    switch = 1  # Define the switching type (1 = combo, 2 = substrate, 3 = antibiotic)
    a_max, b_max = 1, 1  # (1, 1)
    kmaxs, kmaxp = 10, 0.1  # Define the killing rates of both types of bacteria (10, 0.1)
    K_A = 6.4 * MIC  # Define the half-saturation constant for the killing rates of susceptible and persister cells
    K_S = 0.0035  # Define the half-saturation constant for the substrate (0.0035)

    """
    grow_hours, treat_hours, recov_hours = 8, 2, 5  # Define the number of hours for each stage of simulation
    t_grow, t_treat = np.linspace(0, grow_hours, grow_hours * 3600), np.linspace(0, treat_hours, treat_hours * 3600)
    t_recov = np.linspace(0, recov_hours, recov_hours * 3600)
    
    sol = odeint(grow_phase, N, t_grow, args=(0, CS0, a_max, b_max, K1, K2, G, kmaxs, kmaxp, K_A, K_S, switch))
    final_N = [sol[-1, 0], sol[-1, 1]]

    sol2 = odeint(grow_phase, final_N, t_treat, args=(CA, CS0, a_max, b_max, K1, K2, G, kmaxs, kmaxp, K_A, K_S, switch))
    final_N = [sol2[-1, 0], sol2[-1, 1]]

    sol3 = odeint(grow_phase, final_N, t_recov, args=(0, CS0, a_max, b_max, K1, K2, G, kmaxs, kmaxp, K_A, K_S, switch))

    # Combine the solutions to produce one graph
    full_sol = np.concatenate((sol, sol2, sol3))
    t_total = np.concatenate((t_grow, t_treat + grow_hours, t_recov + grow_hours + treat_hours))

    plt.plot(t_total, full_sol[:, 0], 'b', label='Ns'), plt.plot(t_total, full_sol[:, 1], 'g', label='Np')
    plt.legend(loc='best'), plt.xlabel('t'), plt.ylabel('Number of cells')
    # plt.grid()
    plt.show()
    """

    full_sol, t_total, end_time = constant_switch(N, 8, 2, 0.4, CA, CS0, a_max, b_max, K1, K2, G, kmaxs, kmaxp, K_A, K_S, switch)

    plt.plot(t_total, full_sol[:, 0], 'b', label='Ns'), plt.plot(t_total, full_sol[:, 1], 'g', label='Np')
    plt.legend(loc='best'), plt.xlabel('Time'), plt.ylabel('Number of cells')
    plt.axvline(x=end_time, color='r')
    plt.show()
    print(end_time)


if __name__ == '__main__':
    main()
