import numpy as np
import matplotlib.pyplot as plt

import matplotlib.pylab as pylab
params = {'legend.fontsize': 16,
          'legend.title_fontsize': 16,
          'axes.labelsize': 18,
          'axes.titlesize': 18,
          'xtick.labelsize': 18,
          'ytick.labelsize': 18}
pylab.rcParams.update(params)

plt.rcParams.update({
    "text.usetex": True,
    "text.latex.preamble": [r'\usepackage{amsmath}']
})


a = 6371.22e3

def computeViscosityParameter(order, n, t):
    return 1. / t  * np.power(n * (n + 1.) / (a * a), -order / 2)

def computeDamping(order, n, dt, nu, discrete = True):
    if discrete:
        return 1. / (1. + dt * nu * np.power(n * (n + 1) / (a * a), order / 2))
    else:
        return np.exp(-dt * nu * np.power(n * (n + 1) / (a * a), order / 2))


markers = ['o', 'x', 's', '^']
linestyles = ['-', '--', ':']
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

orders = [2,4,6]
modes = np.arange(1, 8192, 1)
input_times = np.array([.25, .5, 1, 2, 4, 8, 16, 32, 64, 128]) * 60. * 60.

###################################################
## plot coefficients in function of damping time ##
###################################################
fig, ax = plt.subplots()

labels = []
styles = []

labels2 = []
styles2 = []

imode = 0
for mode in [32, 128]:
    coefs_time = {}
    facs_time = {}
    iorder = 0
    ## in function of damping time
    for order in orders:
        coefs_time[order] = []
        facs_time[order] = []
        for time in input_times:

            fac = np.power( mode * (mode + 1) / (a * a), order)
            coef = 1. / fac / time

            coef = computeViscosityParameter(order, mode, time)

            coefs_time[order].append(abs(coef))
            ##facs_time[order].append(abs(fac))

        linestyle = linestyles[imode]
        marker = markers[iorder]
        color = colors[iorder]

        plt.plot(input_times / 3600., coefs_time[order], marker = marker, linestyle = linestyle, color = color)

        iorder += 1

        if imode == 0:
            styles2.append(plt.Line2D((0, 1), (0, 0), color=color, marker = marker, linestyle='-'));
            labels2.append(r'$q = {}$'.format(order))


    styles.append(plt.Line2D((0, 1), (0, 0), color='black', linestyle=linestyle));
    labels.append(r'$M = {}$'.format(mode))
    imode += 1


legend = ax.legend(styles, labels, loc = 4);
legend2 = plt.legend(styles2, labels2, loc = 1);
ax.add_artist(legend)
ax.add_artist(legend2);

plt.xscale("log")
plt.yscale("log")
plt.xlabel(r'$\tau \mathrm{(h)}$')
plt.ylabel(r'$\nu (\mathrm{m}^q\mathrm{s}^{-1})$')
###plt.legend(loc = 1)
plt.ylim(bottom = 1e2)
plt.tight_layout()
plt.savefig('coefs_order_time.pdf')
plt.close()


########################################################
## plot damping factors in function of the wavenumber ##
########################################################
fig, ax = plt.subplots()

labels = {}
styles = {}

nb_dt = 1
dt = 120.
discrete = True
iorder = 0
for order in orders:
    damping_factor = {}

    if order == 2:
        coefs = [1e5, 1e6, 1e7]
    elif order == 4:
        coefs = [1e15, 1e16, 1e17]
    elif order == 6:
        coefs = [1e25, 1e26, 1e27]
    inu = 0
    for coef in coefs:
        damping_factor[order] = []
        for mode in modes:

            d = computeDamping(order, mode, dt, coef, discrete = discrete)

            d = np.power(d, nb_dt)
            damping_factor[order].append(d)

        linestyle = linestyles[inu]
        marker = markers[iorder]
        color = colors[iorder]

        plt.plot(modes, damping_factor[order], marker = marker, linestyle = linestyle, color = color, markevery = .1)

        if inu == 0:
            styles[iorder] = []
            labels[iorder] = []

        styles[iorder].append(plt.Line2D((0, 1), (0, 0), color=color, marker = marker, linestyle=linestyle));
        labels[iorder].append(r'$q = {}$'.format(order) + "; " + r'$\nu = 10^{{{}}}$'.format(int(np.log10(coef))))

        inu += 1

    iorder += 1

for iorder in range(len(orders)):
    if iorder == 0:
        loc = 2
    elif iorder == 1:
        loc = 1
    elif iorder == 2:
        loc = 4
    legend = ax.legend(styles[iorder], labels[iorder], loc = loc)
    ax.add_artist(legend)

plt.xscale("log")
plt.ylim(0.995, 1.0025)
plt.xlabel(r'$n$')
if discrete:
    if nb_dt == 1:
        plt.ylabel(r'$\hat{b}_{n,\Delta t}$')
    else:
        plt.ylabel(r'$\hat{{{0}}}_{{{1}}}^{{{2}}}$'.format("b", "n,\Delta t", nb_dt))
plt.tight_layout()
plt.savefig('damping_factor2.pdf')
plt.close()

plt.show()

