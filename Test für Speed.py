import random
from alive_progress import alive_bar
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import norm
from astropy import modeling
from numba import jit

@jit
def RandoMando(steps, start=(0, 0)):
    position = start
    position_set = {position}
    position_list = [position]
    msd_cumulative_list = [0]
    cumulative_sum = 0

    for step in range(1, steps + 1):
        dir_check = True
        tries = 0
        while dir_check:
            new_position = position
            direction = random.randint(1, 4)
            match direction:
                case 1:
                    new_position = (position[0], position[1] + 1)
                case 2:
                    new_position = (position[0], position[1] - 1)
                case 3:
                    new_position = (position[0] + 1, position[1])
                case 4:
                    new_position = (position[0] - 1, position[1])

            if selfavoiding_random_walk:
                if new_position not in position_set:
                    position = new_position
                    position_set.add(position)
                    position_list.append(position)
                    dir_check = False
            else:
                position = new_position
                position_list.append(position)
                dir_check = False

            tries += 1
            if tries > 4:
                return None, None

        squared_distance = position[0] ** 2 + position[1] ** 2
        cumulative_sum += squared_distance
        msd_cumulative = cumulative_sum / step
        msd_cumulative_list.append(msd_cumulative)

    return position_list, msd_cumulative_list

def RandoMandoDangoLaengo(steps):
    while True:
        positions, msd_cumulative_values = RandoMando(steps)
        if positions and msd_cumulative_values:
            return msd_cumulative_values

def nen_fit(tau, a, alpha):
    return a * tau ** alpha

def multiple_iterations(steps, iterations):
    all_msd_values = []

    with alive_bar(iterations, force_tty=True) as bar:
        for i in range(iterations):
            msd_values = RandoMandoDangoLaengo(steps)
            all_msd_values.append(msd_values)
            bar()

    avg_msd = [sum(step_values) / len(step_values) for step_values in zip(*all_msd_values)]

    tau = np.arange(1, len(avg_msd) + 1)
    msd = np.array(avg_msd)

    # noinspection PyTupleAssignmentBalance
    popt, pcov = curve_fit(nen_fit, tau, msd, p0=(1, 0.5))

    msd_fit = nen_fit(tau, *popt)

    plt.figure(figsize=(12, 6))

    for index, msd_values in enumerate(all_msd_values):
        if index < 1000:
            plt.plot(range(1, len(msd_values) + 1), msd_values, alpha=0.3, linewidth=1, color='gray')

    plt.plot(tau, msd_fit, color='blue', linestyle='--', linewidth=2,
             label=f'Fit: a={popt[0]:.4f}, α={popt[1]:.4f}')

    plt.xlabel("Schritte")
    plt.ylabel("Mittlere quadratische Entfernung")
    plt.title(f"Mittlere quadratische Entfernung für {iterations} Iterationen und {steps} Schritte")
    plt.legend()
    plt.grid()
    plt.show()

    return avg_msd, all_msd_values, popt

selfavoiding_random_walk = False
steps = 100
iterations = 1000000
average_msd, all_msd, fit_params = multiple_iterations(steps, iterations)