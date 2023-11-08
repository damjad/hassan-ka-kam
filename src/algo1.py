import csv

import numpy as np
import matplotlib.pyplot as plt

data_points_c = {}
data_points_a_c = {}


def eq_34__E_c(__f_co):
    return 4700 * np.sqrt(__f_co)


def eq_04__sig_cc(__sig_cu, __sig_l):
    return __sig_cu * (1 + 3.5 * (__sig_l / __sig_cu))


def eq_05__eq_07__eps_cc(__sig_l, __sig_cu, __eps_c0):
    ratio = __sig_l / __sig_cu

    if ratio <= 0.1:
        power = 1.2
    else:
        power = 1

    return __eps_c0 * (1 + 17.5 * np.float_power(ratio, power))


def eq_06__eps_c0(__sig_cu):
    return 0.00076 + np.sqrt((0.626 * __sig_cu - 4.33) * np.float_power(10, -7))


def eq_08__sig_c(__sig_cc, __eps, __eps_c, __eps_cc, __r):
    return __sig_cc * ((__eps / __eps_cc) * __r) / (__r - 1 + np.float_power(__eps_c / __eps_cc, __r))


def eq_09__r(__E_c, __sig, __eps_cc):
    return __E_c / (__E_c - (__sig / __eps_cc))


def eq_11__eps_a_c(__eps_u_c, __eps_cc, __eps_cu):
    return __eps_u_c + (__eps_cc - __eps_cu) * __eps_u_c / __eps_cu


def eq_12__sig_a_c(__sig_c, __sig_cc, __sig_cu):
    return __sig_c - (__sig_cc - __sig_cu)


def eq_15__d(__sig_a_c, __sig_cu):
    return 1 - (__sig_a_c / __sig_cu)


def eq_24__eps_p_c(__eps_c, __sig_c, __E_c, __sig_l, __v_c):
    return __sig_c / __E_c - 2 * sig_l / __E_c * __v_c


def delta_sig_l(__sig_l, __sig_cu):
    ratio = __sig_l / __sig_cu
    if ratio <= 0.1:
        return 2.5 / 100.0 * __sig_cu
    else:
        return 10.0 / 100.0 * __sig_cu


def delta_eps_c(__eps_c, __eps_c0):
    # TODO: fix me, assuming it is 1% of eps_c0
    if not __eps_c:
        return 0.001 * __eps_c0
    else:
        return 0.01 * __eps_c


def run(eps_c0, sig_cu, sig_l, sig_lmax, E_c, eps_cu, eps_cmax, v_c):
    if sig_l > sig_lmax:
        return  # TODO: return something

    sig_cc = eq_04__sig_cc(sig_cu, sig_l)
    eps_cc = eq_05__eq_07__eps_cc(sig_l, sig_cu, eps_c0)

    # TODO: assumption sig=sig_cc
    r = eq_09__r(E_c, sig_cc, eps_cc)

    prev_sig_l = sig_l
    sig_l = sig_l + delta_sig_l(sig_l, sig_cu)
    eps_c = 0
    sig_a_c = None
    eps_a_c = None

    while True:
        if eps_c > eps_cmax:
            return run(eps_c0, sig_cu, sig_l, sig_lmax, E_c, eps_cu, eps_cmax, v_c)

        sig_c = eq_08__sig_c(sig_cc, eps_c, eps_c, eps_cc, r)
        sig_u_c = sig_c  # TODO: assumption. fix me
        eps_u_c = eps_c  # TODO: assumption. fix me. maybe it's correct. because in flowchart and equation 10. they are used interchangebly

        if sig_c > sig_cu:
            # TODO: hack to save the curve from going above sig_cu
            sig_c = sig_cu

        if eps_c <= eps_cu:
            sig_a_c = sig_u_c
            eps_a_c = eq_11__eps_a_c(eps_u_c, eps_cc, eps_cu)
            eps_p_c = eq_24__eps_p_c(eps_c, sig_c, E_c, sig_l, v_c)
        else:
            if eps_c > eps_cc:
                sig_a_c = eq_12__sig_a_c(sig_c, sig_cc, sig_cu)
                eps_a_c = eps_c
                eps_p_c = eq_24__eps_p_c(eps_c, sig_c, E_c, sig_l, v_c)
                d = eq_15__d(sig_a_c, sig_cu)

        eps_c = eps_c + delta_eps_c(eps_c, eps_c0)

        # capture data points
        if prev_sig_l not in data_points_c:
            data_points_c[prev_sig_l] = []
        data_points_c[prev_sig_l].append((eps_c, sig_c))

        if prev_sig_l not in data_points_a_c:
            data_points_a_c[prev_sig_l] = []
        data_points_a_c[prev_sig_l].append((eps_a_c, sig_a_c))


def plot_crvs(data_points, title):
    for sig_l, points in data_points.items():
        x_vals, y_vals = zip(*points)
        plt.plot(x_vals, y_vals, label=sig_l)

    plt.legend()
    plt.xlabel('Eps')
    plt.ylabel('Sig')
    plt.title(title)
    plt.show()


def take_inputs(sig_cu, sig_lmax, eps_cmax):
    # inputs

    if not sig_cu:
        inp = input('Please enter sig_cu: ')
        if not inp:
            raise ValueError('Please provide sig_cu')

        sig_cu = np.float64(inp)
    inp = None
    # inp = input('Please enter E_c: ')
    if inp:
        E_c = np.float64(inp)
    else:
        # TODO: assumption
        f_co = sig_cu
        E_c = eq_34__E_c(f_co)

    # inp = input('Please enter eps_c0: ')
    if inp:
        eps_c0 = np.float64(inp)
    else:
        eps_c0 = eq_06__eps_c0(sig_cu)

    eps_cu = eps_c0

    if not sig_lmax:
        inp = input('Please enter sig_lmax: ')
        if not inp:
            raise ValueError('Please provide sig_lmax')
        sig_lmax = np.float64(inp)

    if not eps_cmax:
        inp = input('Please enter eps_cmax: ')
        if not inp:
            raise ValueError('Please provide eps_cmax')
        eps_cmax = np.float64(inp)

    return eps_c0, sig_cu, sig_lmax, E_c, eps_cu, eps_cmax


def plot_curves():
    plot_crvs(data_points_c, 'Confined uniaxial stress-strain curves')
    plot_crvs(data_points_a_c, 'Adjusted uniaxial stress-strain curves.')


def write_file(data_points, name_prefix):
    # Write data points to CSV files
    for sig_l, points in data_points.items():
        filename = f'{name_prefix}_{sig_l}.csv'
        with open(filename, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['x', 'y'])  # Write header
            for x, y in points:
                writer.writerow([x, y])  # Write data points


def write_to_files():
    write_file(data_points_c, 'curves_')
    write_file(data_points_a_c, 'adjusted_curves_')


if __name__ == '__main__':
    # initialize
    sig_l = 0

    # assumptions
    v_c = 0.2
    sig_cu = 30.0
    sig_lmax = 2
    eps_cmax = 0.006

    eps_c0, sig_cu, sig_lmax, E_c, eps_cu, eps_cmax = take_inputs(sig_cu, sig_lmax, eps_cmax)

    run(eps_c0, sig_cu, sig_l, sig_lmax, E_c, eps_cu, eps_cmax, v_c)
    plot_curves()
    write_to_files()
