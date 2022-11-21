import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.integrate import odeint
import csv


def main():
    """
    Main function where parameter values are defined and the process is run
    """

    a_max_results, b_max_results, T1_results, T2_results, end_time_results, treat_times = [], [], [], [], [], []

    file = open('output.csv')
    csvreader = csv.reader(file)
    header = []
    header = next(csvreader)

    rows = []
    for row in csvreader:
        rows.append(row)

    for row in rows:
        a_max_results.append(float(row[0])), b_max_results.append(float(row[1])), T1_results.append(float(row[2]))
        T2_results.append(float(row[3])), end_time_results.append(float(row[4])), treat_times.append(float(row[5]))

    a_s, b_s = set(a_max_results), set(b_max_results)

    T1_results, T2_results = np.asarray(T1_results), np.asarray(T2_results)
    T1_results, T2_results = np.reshape(T1_results, (len(b_s), len(a_s))).T, np.reshape(T2_results,
                                                                                        (len(b_s), len(a_s))).T

    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(8, 4))
    tkx, x0 = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], -0.5
    tky = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
    x_len, y_len = len(a_s) + x0, len(b_s) + x0

    im1, im2 = ax1.imshow(T1_results*60, cmap='viridis'), ax2.imshow(T2_results*60, cmap='viridis')

    ax1.set_xticks(np.linspace(x0, x_len, 10), labels=tkx), ax2.set_xticks(np.linspace(x0, x_len, 10), labels=tkx)
    ax1.set_yticks(np.linspace(x0, y_len, 10), labels=tky), ax2.set_yticks(np.linspace(x0, y_len, 10), labels=tky)
    # Set ticks to surround each plotted box
    ax1.set_xticks(np.arange(len(a_s))-0.5, minor=True), ax2.set_xticks(np.arange(len(a_s))-0.5, minor=True)
    ax1.set_yticks(np.arange(len(b_s))-0.5, minor=True), ax2.set_yticks(np.arange(len(b_s))-0.5, minor=True)
    # Plot a box around each value
    ax1.grid(which='minor', color='k', linestyle='-', linewidth=0.3), ax2.grid(which='minor', color='k', linestyle='-',
                                                                               linewidth=0.3)
    ax1.tick_params(axis='both', which='both', length=0), ax2.tick_params(axis='both', which='both', length=0)

    divider1, divider2 = make_axes_locatable(ax1), make_axes_locatable(ax2)
    cax1, cax2 = divider1.append_axes("right", size="5%", pad=0.2), divider2.append_axes("right", size="5%", pad=0.5)
    cbar1, cbar2 = plt.colorbar(im1, cax=cax1), plt.colorbar(im2, cax=cax2)
    cbar1.set_label(label='Optimised treatment duration (min)', size=10)
    cbar2.set_label(label='Optimised treatment duration (min)', size=10)
    ax1.set_xlabel('$a_{max}$'), ax2.set_xlabel('$a_{max}$'), ax1.set_ylabel('$b_{max}$'), ax2.set_ylabel('$b_{max}$')
    fig.tight_layout(pad=5)
    plt.savefig("ParameterScanTrial.pdf", bbox_inches='tight', pad_inches=0)
    plt.show()


if __name__ == '__main__':
    main()
