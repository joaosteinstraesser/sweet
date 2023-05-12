import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import sys
import math

matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
matplotlib.rc('legend',fontsize=12)
matplotlib.rc('axes',labelsize=14)

small_stability_region = 1e-12

def plotContourSeveralSchemes(y, Ak, methods, case, parareal = False, figname_parareal = "", levels_dict = 1, invert_scale = False, multiple_pi = True, linewidth = 1.5, main_style = 'color', several_y = False):

    dirname = "combined_contours";

    if not os.path.isdir(dirname):
        os.makedirs(dirname);

    extent = [range_min_x, range_max_x, range_min_y, range_max_y];

    if case == "":
        label_x = r'$\Re(y)$';
        label_y = r'$\Im(y)$';
    elif case == "RR":
        label_x = r'$\Re(x)$';
        label_y = r'$\Re(y)$';
    elif case == "RI":
        label_x = r'$\Re(x)$';
        label_y = r'$\Im(y)$';
    elif case == "IR":
        label_x = r'$\Im(x)$';
        label_y = r'$\Re(y)$';
    elif case == "II":
        label_x = r'$\Im(x)$';
        label_y = r'$\Im(y)$';
    elif case == "R0":
        label_x = r'$Re(\xi_{\mathbf{N}})$';
        label_y = r'$Im(\xi_{\mathbf{N}})$';
        if invert_scale:
            label_x = r'$Re\left(\frac{\xi_{\mathbf{N}}}{m_c}\right)$';
            label_y = r'$Im\left(\frac{\xi_{\mathbf{N}}}{m_c}\right)$';

    y_value = "";
    y_title = "";
    if case == "R0":
        if multiple_pi:
            y_value = "_y" + str(y / np.pi) + "pi";
            y_title = " for " + r'$\xi_{\mathbf{L}}$' + r'$ = {}\pi i $'.format(y / np.pi)
        else:
            y_value = "_y" + str(round(y / 2.5e-4)) + "xi";
            y_title = " for " + r'$\xi_{\mathbf{L}}$' + r'$ = {} $'.format(round(y / 2.5e-4)) + r'$\tilde{\xi}_{\mathbf{L}}$'
        print(y, y.imag, y / np.pi)

    ## plot for all kdxs only contour
    fig, ax = plt.subplots(1, 1);
    linestyles_all = ["-", "--", "-.", ":"]
    colors_all = plt.rcParams['axes.prop_cycle'].by_key()['color'];
    icolor = 0;
    styles = [];
    labels = [];


    ax.axhline(y=0, linewidth = .5, linestyle = "--",  color='k')
    ax.axvline(x=0, linewidth = .5, linestyle = "--",  color='k')

    if levels_dict == 1:
        for method in methods:
            linestyles = [linestyles_all[icolor % len(linestyles_all)]];
            colors = [colors_all[icolor % len(colors_all)]];
            #####im = ax.imshow(Ak[method], cmap = plt.cm.jet, extent = extent, origin = "lower");
            #####fig.colorbar(im, ax = ax)
            CS = ax.contour(Ak[method], [1.0], linestyles = linestyles, colors = colors, extent = extent, origin = "lower", linewidths = linewidth);
            label = method;
            ##CS.collections[0].set_label(label);


            styles.append(plt.Line2D((0, 1), (0, 0), color=colors_all[icolor], linestyle=linestyles_all[icolor%4]));
            labels.append(label);

            icolor += 1;
    elif levels_dict == 2:
        styles2 = [];
        labels2 = [];
        for method in Ak.keys():
            if main_style == 'color':
                colors = [colors_all[icolor % len(colors_all)]];
            elif main_style == 'linestyle':
                linestyles = [linestyles_all[icolor % len(linestyles_all)]];
            ilinestyle = 0;
            for method2 in Ak[method].keys():
                if main_style == 'color':
                    linestyles = [linestyles_all[ilinestyle % len(linestyles_all)]];
                elif main_style == 'linestyle':
                    colors = [colors_all[ilinestyle % len(colors_all)]];
                CS = ax.contour(Ak[method][method2], [1.0], linestyles = linestyles, colors = colors, extent = extent, origin = "lower", linewidths = linewidth);
                label = method + "; " + method2;

                if icolor == 0:
                    if main_style == 'color':
                        styles2.append(plt.Line2D((0, 1), (0, 0), color='black', linestyle=linestyles_all[ilinestyle]));
                    elif main_style == 'linestyle':
                        styles2.append(plt.Line2D((0, 1), (0, 0), color=colors_all[ilinestyle], linestyle="-"));
                    labels2.append(method2);

                ilinestyle += 1;

            if main_style == 'color':
                styles.append(plt.Line2D((0, 1), (0, 0), color=colors_all[icolor], linestyle="-"));
            elif main_style == 'linestyle':
                styles.append(plt.Line2D((0, 1), (0, 0), color='black', linestyle=linestyles_all[icolor]));
            labels.append(method);

            icolor += 1;


    legend = ax.legend(styles, labels, loc = 1);

    if levels_dict == 2:
        legend2 = plt.legend(styles2, labels2, loc = 4);
        ax.add_artist(legend)
        ax.add_artist(legend2);

    ##ax.legend();
    if not several_y:
        ax.set_title("Stability region for " + y_title);

    ax.set_xlabel(label_x);
    ax.set_ylabel(label_y);

    figname = dirname + "/stability_contour";
    if not parareal:
        for method in methods:
            figname += "_" + method;
        if several_y:
            figname += "_several_y"
    else:
        figname += figname_parareal;
    if not several_y:
        figname += "_case" + y_value;
    fig.savefig(figname + ".pdf");

    plt.close()


def ups(n, z):

    small = 1e-10;

    idx = np.argwhere(np.abs(z) < small);

    if type(z) == int or np.linalg.norm(z) < small:
        return 1. / 6.;
    else:
        if n == 1:
            out = (-4. - z + np.exp(z) * (4 - 3. * z + z * z) ) / (z * z * z);
        elif n == 2:
            out = (2. + z + np.exp(z) * (-2. + z) ) / (z * z * z);
        elif n == 3:
            out = (-4. - 3. * z - z * z + np.exp(z) * (4 - z) ) / (z * z * z);

    if idx.size > 0:
        out[idx] = 1. / 6.;

    return out;

def phin(n, z):

    small = 1e-10;

    if n == 0:
        return np.exp(z);

    idx = np.argwhere(np.abs(z) < small);
    idx2 = np.argwhere(np.abs(z) >= small);

    np.seterr(divide='ignore');

    if type(z) == int or np.linalg.norm(z) < small:
        if np.abs(z) < small:
            if n == 1:
                return 1.;
            elif n == 2:
                return .5;

    out = 1. / z * (phin(n - 1, z) - phin( n - 1, 0));

    if idx.size > 0:
        if n == 1:
            out[idx] = 1.;
        elif n == 2:
            out[idx] = .5;


    return out

def phin_modified(n, z, alpha):

    small = 1e-10;


    A_cn = (2 + z) / (2 - z);
    theta_cn = np.angle(A_cn);
    zm = alpha * z + (1. - alpha) * 1j * theta_cn;

    if n == 0:
        return np.exp(zm);

    idx = np.argwhere(np.abs(z) < small);
    idx2 = np.argwhere(np.abs(z) >= small);

    np.seterr(divide='ignore');

    if type(z) == int or np.linalg.norm(z) < small:
        if np.abs(z) < small:
            if n == 1:
                return 1.;
            elif n == 2:
                return .5;

    out = 1. / z * (phin_modified(n - 1, z, alpha) - phin_modified( n - 1, 0, alpha));

    if idx.size > 0:
        if n == 1:
            out[idx] = 1.;
        elif n == 2:
            out[idx] = .5;


    return out


def psin(n, z):

    s = 0;
    for i in range(1, n):
        s += phin(i, -z);
    return np.power(-1, n + 1) * phin(n, -z) + s;


def combination(m, n):

    return math.comb(m, n);


def computeStabilityRegion(x, y, kdx, alpha, method):

    if method == "IMEX":

        #####return (2. + y) / (2. - y) *  ( (1. + x / 2.) ** 2 );
        return ( (4. + y ) / ( 4. - y )  ) ** 2 * .5 * x * (2. + x);

    elif method == "SL-SI-SETTLS":

        eikdx = np.exp(-1j * kdx);

        a = 1. - .5 * y;
        b = - (  .5 * y * eikdx + x * eikdx + .5 * x +  eikdx  );
        c = .5 * x * eikdx;

        delta = b * b - 4. * a * c;
        sqrt = np.sqrt(delta.astype(complex));
        Ak_plus =  (- b + sqrt ) / (2. * a);
        Ak_minus = ( -b - sqrt ) / (2. * a);

        return [Ak_plus, Ak_minus];

    elif method == "SL-ETD1RK":

        ey = np.exp(y);
        emy = np.exp(-y);
        eikdx = np.exp(-1j * kdx)

        if type(y) is float:
            if not y == 0:
                computed = ey * eikdx * ( 1. - x / y * (emy - 1) )
            else:
                computed = eikdx * (1. + x);

        else:
            computed = ey * eikdx * ( 1. - x / y * (emy - 1)  )
        return computed;


    elif method == "SL-ETD2RK":

        eikdx = np.exp(-1j * kdx)
        r1 = phin(0, y) * eikdx * ( 1 + x * psin(1, y));
        computed = r1 + x * phin(0, y) * psin(2, y) * (r1 - eikdx);

        return computed;

    elif method == "SL-ETD2RK-bis":

        eikdx = np.exp(-1j * kdx)
        r1 = phin(0, y) * eikdx * ( 1 + x * psin(1, y));
        computed = phin(0, y) * ( eikdx * ( 1 + .5 * x * psin(1, y)) + .5 * x * psin(1, y) * r1 );

        return computed;

    elif method == "CN":

        return (2. * x + y + 2.) / (2. - y);

    elif method == "ETD1RK":

        computed = phin(0, y) + x * phin(1, y);

        return computed;

    elif method == "ETD2RK":

        r1 = phin(0, y) + x * phin(1, y);
        computed = r1 + x * phin(2, y) * (r1 - 1);

        return computed;

    elif method == "ETD3RK":

        An = phin(0, .5 * y)      + .5 * x * phin(1, .5 * y);
        Bn = phin(0, .5 * y)      + .5 * x * phin(1, .5 * y) * (2. * An - 1);

        R0 = 1.;
        R1 = x;
        R2 = x * An;
        R3 = x * Bn;

        computed = phin(0, y) * R0 + ups(1, y) * R1 + 4. * ups(2, y) * R2 + ups(3, y) * R3;

        return computed;

    elif method == "ETD4RK":

        An = phin(0, .5 * y)      + .5 * x * phin(1, .5 * y);
        Bn = phin(0, .5 * y)      + .5 * x * phin(1, .5 * y) * An;
        Cn = phin(0, .5 * y) * An + .5 * x * phin(1, .5 * y) * (2.* Bn - 1.);

        R0 = 1.;
        R1 = x;
        R2 = x * (An + Bn);
        R3 = x * Cn;

        computed = phin(0, y) * R0 + ups(1, y) * R1 + 2. * ups(2, y) * R2 + ups(3, y) * R3;

        return computed;

    elif method == "ETD2RK-bis":

        r1 = phin(0, y) + x * phin(1, y);
        computed = phin(0, y) + x * phin(1, y) * .5 * (r1 + 1.);

        return computed;



def computeStabilityRegionLinearProblem(y_real, y_imag, kdxs, alphas, method):

    imag_range_line = imag_range[:, 0];

    Ak = {};
    Ak_all = np.ones_like(y_real);
    A_relative = {};
    theta = {};

    theta_relative_x = {};
    theta_relative_y = {};

    y = y_real + 1j * y_imag;

    u_exact = np.exp(y);
    A_exact = np.abs(u_exact);
    theta_exact = np.angle(u_exact);

    if method == "1997":

        for kdx in kdxs:

            Ak_plus, Ak_minus = computeStabilityRegion(x, y, kdx, 0, method);

            A_relative[kdx] = Ak_plus / A_exact;
            theta[kdx] = np.angle(Ak_plus);
            theta_relative[kdx] = theta[kdx] / theta_exact;

            Ak_plus = np.abs(Ak_plus) <= 1. ;
            Ak_minus = np.abs(Ak_minus) <= 1. ;

            Ak[kdx] = Ak_plus * Ak_minus;
            Ak_all = Ak_all * Ak[kdx];


    elif method == "SL-SI-SETTLS":

        for kdx in kdxs:

            Ak_plus, Ak_minus = computeStabilityRegion(x, y, kdx, 0, method);

            A_relative[kdx] = Ak_plus / A_exact;
            theta[kdx] = np.angle(Ak_plus);
            theta_relative[kdx] = theta[kdx] / theta_exact;

            Ak_plus = np.abs(Ak_plus) <= 1. ;
            Ak_minus = np.abs(Ak_minus) <= 1. ;

            Ak[kdx] = Ak_plus * Ak_minus;
            Ak_all = Ak_all * Ak[kdx]


    plotStabilityTheta(y_real, y_imag, kdxs, alphas, Ak, Ak_all, theta, theta_exact, method, "")


def computeAndPlotStabilityRegionsParareal(x, y, kdxs, kdxs_extended, alphas, tsm_fine, tsm_coarse, case, fixed_ys, niter, time_slice, coarsening_factor, plot = False):
    ## x: real range
    ## y: imag range
    ## kdx: values of k*s (SL schemes)
    ## kdxs_extended: more values of k*s in order to compute intersection of stability regions
    ## alphas: weights in combined methods (e.g. CN-ETD2RK)
    ## tsm_fine: fine timestepping_method
    ## tsm_coarse: coarse_timestepping method
    ## case: RR, RI, IR, II or R0
        ## RR: real part of lambda_L and real part of lambda_N
        ## R0: x = x + i * y = lambda_N; fixed y = lambda_L
    ## fixed_ys (float): fixed value for lambda_L (should contain a single value)
    ## niter: parareal iteration
    ## time_slice: coarse time slice
    ## coarsening factor: Dt / dt


    ## The stability functions of tsm_fine and tsm_coarse are computed using computeStabilityRegionNonlinearProblem
    ## This function receives the same arguments, expect those specific for Parareal
    ## In the case "R0", x is set to x + i * y
    ## and y is set to the single value contained in fixed_ys
    ## The original value of y is not further used


    pass;


def computeStabilityRegionNonlinearProblem(x, y, kdxs, kdxs_extended, alphas, method, case, fixed_ys, plot = True, get_contour = True, parareal = False, niter = None, time_slice = None, coarsening_factor = None, tsm_fine = None, tsm_coarse = None, nrelax = 0, nlevels = 2, invert_scale = False):

    Ak = {};
    Ak_all = np.ones_like(x);
    A_relative = {};
    theta = {};
    theta_relative = {};

    x_orig = x;
    y_orig = y;

    if case == "RR":
        x = x_orig;
        y = y_orig;
    elif case == "RI":
        x = x_orig;
        y = 1j * y_orig;
    elif case == "IR":
        x = 1j * x_orig;
        y = y_orig;
    elif case == "II":
        x = 1j * x_orig;
        y = 1j * y_orig;
    elif case == "R0":
        x = x_orig + 1j * y_orig;
        ##y = 0. * y_orig;
        ##y = 1j * 0. *  np.pi;

    u_exact = np.exp(y) + np.exp(x);
    A_exact = np.abs(u_exact);

    u_exact_x = np.exp(1j * x_line)
    theta_exact_x = np.angle(u_exact_x);

    methods = ["CN", "AB2AM2", "ETD1RK", "ETD2RK", "ETD3RK", "ETD4RK", "ETD2RK-bis", "ETD2RK-backward", "ETD2RK-alternative", "ETD2RK-Leapfrog", "ETD2RK-B", "ETD1RK-imp", "ETD1RK-ter", "ETD1RK-quater", "ETD2RK-bis", "Genilson-old", "Genilson-new", "IMEX", "SI-SETTLS"];
    methods_kdx = ["SL-EXP-SETTLS", "SL-ETD1RK", "SL-ETD2RK", "SL-ETD2RK-bis", "SL-SI-SETTLS", "SL-alternative", "SL-ETD1RK-SETTLS", "SL-ETD2RK-SETTLS"];
    methods_alpha = ["ETD2RK-CN", "ETD2RK-CN-argument-L","ETD2RK-CN-argument-L-NL","ETD2RK-CN-argument-NL", "ETD2RK-CN-modified-L", "ETD2RK-CN-modified-L-NL", "ETD2RK-CN-modified-NL"];
    methods_kdx_alpha = ["SL-ETD2RK-SL-SI-SETTLS"];
    methods_quadratic = ["SL-EXP-SETTLS", "AB2AM2", "ETD2RK-alternativa", "ETD2RK-Leapfrog", "SL-ETD2RK-SL-SI-SETTLS", "SL-SI-SETTLS", "SL-ETD1RK-SETTLS", "SL-ETD2RK-SETTLS", "SI-SETTLS"];

    if case == "R0":
        ys = np.array(fixed_ys) * 1j * np.pi;
    else:
        ys = [y];

    for y in ys:
        if case == "R0":
            y_orig = y;

        if parareal:

            method = "parareal_" + tsm_fine + "_" + tsm_coarse;

            if tsm_fine in methods_kdx:
                kdxs_fine = kdxs;
            else:
                kdxs_fine = [kdxs[0]];

            if tsm_coarse in methods_kdx:
                kdxs_coarse = kdxs;
            else:
                kdxs_coarse = [kdxs[0]];

            Ak = 1;

            for kdx_fine in kdxs_fine:
                for kdx_coarse in kdxs_coarse:

                    Ak_fine = np.empty(2, dtype = np.ndarray);
                    Ak_fine_power = np.empty(2, dtype = np.ndarray);
                    Ak_coarse_power = {};

                    if not invert_scale:
                        if nlevels == 2:

                            if tsm_fine in methods_quadratic:
                                Ak_fine[0], Ak_fine[1] = computeStabilityRegion(x / coarsening_factor, y / coarsening_factor, kdx_fine, 0, tsm_fine);
                            else:
                                Ak_fine[0] = computeStabilityRegion(x / coarsening_factor, y / coarsening_factor, kdx_fine, 0, tsm_fine);
                                Ak_fine[1] = Ak_fine[0];

                            for i in range(2):
                                Ak_fine[i] = Ak_fine[i] ** coarsening_factor;
                                Ak_fine_power[i] = Ak_fine[i] ** nrelax;

                            Ak_coarse = np.empty(2, dtype = np.ndarray);
                            if tsm_coarse in methods_quadratic:
                                Ak_coarse[0], Ak_coarse[1] = computeStabilityRegion(x, y, kdx_coarse, 0, tsm_coarse);
                            else:
                                Ak_coarse[0] = computeStabilityRegion(x, y, kdx_coarse, 0, tsm_coarse);
                                Ak_coarse[1] = Ak_coarse[0];

                        else:

                            cfactor_fine = coarsening_factor ** (nlevels - 1);
                            if tsm_fine in methods_quadratic:
                                Ak_fine[0], Ak_fine[1] = computeStabilityRegion(x / cfactor_fine, y / cfactor_fine, kdx_fine, 0, tsm_fine);
                            else:
                                Ak_fine[0] = computeStabilityRegion(x / cfactor_fine, y / cfactor_fine, kdx_fine, 0, tsm_fine);
                                Ak_fine[1] = Ak_fine[0];

                            for i in range(2):
                                Ak_fine[i] = Ak_fine[i] ** cfactor_fine;
                                Ak_fine_power[i] = Ak_fine[i] ** nrelax;

                            ####Ak_coarse = np.empty(2, dtype = np.ndarray);
                            ####if tsm_coarse in methods_quadratic:
                            ####    Ak_coarse[0], Ak_coarse[1] = computeStabilityRegion(x, y, kdx_coarse, 0, tsm_coarse);
                            ####else:
                            ####    Ak_coarse[0] = computeStabilityRegion(x, y, kdx_coarse, 0, tsm_coarse);
                            ####    Ak_coarse[1] = Ak_coarse[0];



                            for l in range(1, nlevels):
                                Ak_coarse_power[l] = np.empty(2, dtype = np.ndarray);
                                cfactor_coarse = coarsening_factor ** (nlevels - 1 - l);
                                if tsm_coarse in methods_quadratic:
                                    Ak_coarse_power[l][0], Ak_coarse_power[l][1] = computeStabilityRegion(x / cfactor_coarse, y / cfactor_coarse, kdx_coarse, 0, tsm_coarse);
                                else:
                                    Ak_coarse_power[l][0] = computeStabilityRegion(x / cfactor_coarse, y / cfactor_coarse, kdx_coarse, 0, tsm_coarse);
                                    Ak_coarse_power[l][1] = Ak_coarse_power[l][0];

                                for i in range(2):
                                    Ak_coarse_power[l][i] = np.power(Ak_coarse_power[l][i], cfactor_coarse );
                    #######else:
                    #######    if tsm_fine in methods_quadratic:
                    #######        Ak_fine[0], Ak_fine[1] = computeStabilityRegion(x, y, kdx_fine, 0, tsm_fine);
                    #######    else:
                    #######        Ak_fine[0] = computeStabilityRegion(x, y, kdx_fine, 0, tsm_fine);
                    #######        Ak_fine[1] = Ak_fine[0];

                    #######    for i in range(2):
                    #######        Ak_fine[i] = Ak_fine[i]### ** coarsening_factor;
                    #######        Ak_fine_power[i] = Ak_fine[i] ** nrelax;

                    #######    Ak_coarse = np.empty(2, dtype = np.ndarray);
                    #######    if tsm_coarse in methods_quadratic:
                    #######        Ak_coarse[0], Ak_coarse[1]  = computeStabilityRegion(x * coarsening_factor, y * coarsening_factor, kdx_coarse, 0, tsm_coarse);
                    #######    else:
                    #######        Ak_coarse[0] = computeStabilityRegion(x * coarsening_factor, y * coarsening_factor, kdx_coarse, 0, tsm_coarse);
                    #######        Ak_coarse[1] = Ak_coarse[0];

                    #######    for i in range(2):
                    #######        Ak_coarse[i] = Ak_coarse[i] ** ( 1. / coarsening_factor);

                    #######    if nlevels > 2:
                    #######        for l in range(1, nlevels):
                    #######            Ak_coarse_power[l] = np.empty(2, dtype = np.ndarray);
                    #######            for i in range(2):
                    #######                Ak_coarse_power[l][i] = Ak_coarse[i] ** ( (nlevels - 1 - l) * coarsening_factor );

                    Aks = np.empty(4, dtype = np.ndarray);
                    for iF in range(2):
                        for iC in range(2):
                            idx = 2 * iF + iC;
                            Aks[idx] = 0;
                            #######if nrelax == 0:
                            #######    for i in range(niter + 1):
                            #######        Aks[idx] += combination(time_slice, i) * np.power(Ak_fine[iF] - Ak_coarse[iC], i) * np.power(Ak_coarse[iC], time_slice - i);
                            #######elif nrelax == 1:
                            #######    for i in range((niter // 2) + 1):
                            #######        Aks[idx] += combination(time_slice - i, i) * np.power(Ak_fine[iF] * (Ak_fine[iF] - Ak_coarse[iC]), i) * np.power(Ak_coarse[iC], time_slice - 2 * i);
                            if nlevels == 2:
                                for i in range((niter // (nrelax + 1)) + 1):
                                    Aks[idx] += combination(time_slice - nrelax * i, i) * np.power(Ak_fine_power[iF] * (Ak_fine[iF] - Ak_coarse[iC]), i) * np.power(Ak_coarse[iC], time_slice - (nrelax + 1) * i);
                            else:
                                ###S = ( Ak_fine_power[iF] **  (nlevels - 1) ) * (Ak_fine[iF] ** (nlevels - 1) - Ak_coarse_power[1][iC]);
                                S = Ak_fine_power[iF]  * (Ak_fine[iF]  - Ak_coarse_power[1][iC]);
                                for l in range(1, nlevels - 1):
                                    S += (Ak_coarse_power[l][iC] ** nrelax) * ( Ak_coarse_power[l][iC] - Ak_coarse_power[l + 1][iC] );
                                for i in range((niter // (nrelax + 1)) + 1):
                                    Aks[idx] += combination(time_slice - nrelax * i, i) * np.power(S, i) * np.power(Ak_coarse_power[nlevels - 1][iC], time_slice - (nrelax + 1) * i);
                                ##print(np.linalg.norm(Ak_coarse_power[nlevels - 1][iC] - Ak_coarse[iC]))
                            ###A_relative = np.abs(Ak) / np.abs(A_exact);
                            A_relative = np.abs(Aks[idx]) / np.abs(A_exact);
                            Aks[idx] = np.abs(Aks[idx]) <= 1 + small_stability_region;

                    for i in range(4):
                        Ak *= Aks[i];

            theta_relative_x = None
            theta_relative_y = None

        elif method in methods:

            if not method in methods_quadratic:
                Ak = computeStabilityRegion(x, y, 0, 0, method);
                A_relative = np.abs(Ak) / np.abs(A_exact);
                if get_contour:
                    Ak = np.abs(Ak) <= 1 + small_stability_region;

                if plot_dispersion:
                    Ak_x = computeStabilityRegion(1j * x_line, 0., 0, 0, method);
                    Ak_y = computeStabilityRegion(0., 1j * x_line, 0, 0, method);
            else:
                Ak_plus, Ak_minus = computeStabilityRegion(x, y, 0, 0, method);
                A_relative = np.abs(Ak_plus) / np.abs(A_exact);
                Ak_plus = np.abs(Ak_plus) <= 1. + small_stability_region;
                Ak_minus = np.abs(Ak_minus) <= 1. + small_stability_region;
                Ak = Ak_plus * Ak_minus;

                if plot_dispersion:
                    Ak_x = computeStabilityRegion(1j * x_line, 0., 0, 0, method)[0];
                    Ak_y = computeStabilityRegion(0., 1j * x_line, 0, 0, method)[0];

            if plot_dispersion:
                theta_taylor_x = computeTaylorDispersion(x_line, 0., 0, 0, method);
                theta_x = np.angle(Ak_x)
                theta_relative_x  = {};
                theta_relative_x["computed"] = theta_x / theta_exact_x;
                theta_relative_x["taylor"] = theta_taylor_x / theta_exact_x;

                theta_taylor_y = computeTaylorDispersion(0., x_line, 0, 0, method);
                theta_y = np.angle(Ak_y);
                theta_relative_y = {};
                theta_relative_y["computed"] = theta_y / theta_exact_x;
                theta_relative_y["taylor"] = theta_taylor_y / theta_exact_x;

        elif method in methods_kdx:
            theta_relative_x  = {};
            theta_relative_y  = {};

            for kdx in kdxs_extended:

                if not method in methods_quadratic:

                    Ak[kdx] = computeStabilityRegion(x, y, kdx, 0, method);
                    A_relative[kdx] = Ak[kdx] / A_exact;
                    Ak[kdx] = np.abs(Ak[kdx]) <= 1 + small_stability_region
                    Ak_all = Ak_all * Ak[kdx];

                    if plot_dispersion:
                        Ak_x = computeStabilityRegion(1j * x_line, 0., kdx, 0, method);
                        Ak_y = computeStabilityRegion(0., 1j * x_line, kdx, 0, method);

                else:
                    Ak_plus, Ak_minus = computeStabilityRegion(x, y, kdx, 0, method);
                    A_relative[kdx] = Ak_plus / A_exact;
                    Ak_plus = np.abs(Ak_plus) <= 1. + small_stability_region;
                    Ak_minus = np.abs(Ak_minus) <= 1. + small_stability_region;
                    Ak[kdx] = Ak_plus * Ak_minus;
                    Ak_all = Ak_all * Ak[kdx]

                    if plot_dispersion:
                        Ak_x = computeStabilityRegion(1j * x_line, 0., kdx, 0, method)[0];
                        Ak_y = computeStabilityRegion(0., 1j * x_line, kdx, 0, method)[0];


                if plot_dispersion:
                    u_exact_x = np.exp(1j * x_line + 1j * kdx);
                    theta_exact_x = np.angle(u_exact_x);

                    theta_taylor_x = computeTaylorDispersion(x_line, 0., kdx, 0, method);
                    theta_x = np.angle(Ak_x)
                    theta_relative_x[kdx] = {};
                    theta_taylor_x = theta_x;
                    theta_relative_x[kdx]["computed"] = theta_x / theta_exact_x;
                    theta_relative_x[kdx]["taylor"] = theta_taylor_x / theta_exact_x;

                    theta_taylor_y = computeTaylorDispersion(0., x_line, kdx, 0, method);
                    theta_y = np.angle(Ak_y);
                    theta_taylor_y = theta_y;
                    theta_relative_y[kdx] = {};
                    theta_relative_y[kdx]["computed"] = theta_y / theta_exact_x;
                    theta_relative_y[kdx]["taylor"] = theta_taylor_y / theta_exact_x;

            if not plot:
                Ak = Ak_all;

        elif method in methods_alpha:

            theta_relative_x  = {};
            theta_relative_y  = {};

            for alpha in alphas:

                if not method in methods_quadratic:

                    Ak[alpha] = computeStabilityRegion(x, y, 0, alpha, method);
                    A_relative[alpha] = np.abs(Ak[alpha]) / np.abs(A_exact);
                    Ak[alpha] = np.abs(Ak[alpha]) <= 1 + small_stability_region;
                    Ak_all = Ak_all * Ak[alpha];

                    if plot_dispersion:
                        Ak_x = computeStabilityRegion(1j * x_line, 0., 0, alpha, method);
                        Ak_y = computeStabilityRegion(0., 1j * x_line, 0, alpha, method);

                else:

                    Ak_plus, Ak_minus = computeStabilityRegion(x, y, 0, alpha, method);
                    A_relative[alpha] = Ak_plus / A_exact;
                    Ak_plus = np.abs(Ak_plus) <= 1. + small_stability_region;
                    Ak_minus = np.abs(Ak_minus) <= 1. + small_stability_region;
                    Ak[alpha] = Ak_plus * Ak_minus;
                    Ak_all = Ak_all * Ak[alpha]

                    if plot_dispersion:
                        Ak_x = computeStabilityRegion(1j * x_line, 0., 0, alpha, method)[0];
                        Ak_y = computeStabilityRegion(0., 1j * x_line, 0, alpha, method)[0];

                if plot_dispersion:
                    theta_taylor_x = computeTaylorDispersion(x_line, 0., 0, alpha, method);
                    theta_x = np.angle(Ak_x)
                    theta_relative_x[alpha] = {};
                    theta_relative_x[alpha]["computed"] = theta_x / theta_exact_x;
                    theta_relative_x[alpha]["taylor"] = theta_taylor_x / theta_exact_x;

                    theta_taylor_y = computeTaylorDispersion(0., x_line, 0, alpha, method);
                    theta_y = np.angle(Ak_y);
                    theta_relative_y[alpha] = {};
                    theta_relative_y[alpha]["computed"] = theta_y / theta_exact_x;
                    theta_relative_y[alpha]["taylor"] = theta_taylor_y / theta_exact_x;



        elif method in methods_kdx_alpha:

            theta_relative_x  = {};
            theta_relative_y  = {};
            Ak_all = {};

            for kdx in kdxs_extended:

                theta_relative_x[kdx] = {};
                theta_relative_y[kdx] = {};
                Ak[kdx] = {};
                Ak_all[kdx] = np.ones_like(x, dtype = float);

                for alpha in alphas:
                    print (case, kdx, alpha);

                    Ak_plus, Ak_minus = computeStabilityRegion(x, y, kdx, alpha, method);
                    A_relative[kdx] = Ak_plus / A_exact;
                    Ak_plus = np.abs(Ak_plus) <= 1. + small_stability_region;
                    Ak_minus = np.abs(Ak_minus) <= 1. + small_stability_region;
                    Ak[kdx][alpha] = Ak_plus * Ak_minus;
                    Ak_all[kdx] = Ak_all[kdx] * Ak[kdx][alpha]

                    if plot_dispersion:

                        Ak_x = computeStabilityRegion(1j * x_line, 0., kdx, alpha, method)[0];
                        theta_taylor_x = computeTaylorDispersion(x_line, 0., kdx, alpha, method);
                        theta_x = np.angle(Ak_x)
                        theta_relative_x[kdx][alpha] = {};
                        theta_relative_x[kdx][alpha]["computed"] = theta_x / theta_exact_x;
                        theta_relative_x[kdx][alpha]["taylor"] = theta_taylor_x / theta_exact_x;

                        Ak_y = computeStabilityRegion(0., 1j * x_line, kdx, alpha, method)[0];
                        theta_taylor_y = computeTaylorDispersion(0., x_line, kdx, alpha, method);
                        theta_y = np.angle(Ak_y);
                        theta_relative_y[kdx][alpha] = {};
                        theta_relative_y[kdx][alpha]["computed"] = theta_y / theta_exact_x;
                        theta_relative_y[kdx][alpha]["taylor"] = theta_taylor_y / theta_exact_x;

            else:
                sys.exit("Unknown method " + method);

        if not plot_dispersion:

            theta_relative_x = None;
            theta_relative_y = None;

        if plot:
            if not parareal:
                plotStabilityTheta(x_orig, y_orig, kdxs, alphas, Ak, Ak_all, A_relative, theta_relative_x, theta_relative_y, method, case)
            else:
                plotStabilityTheta(x_orig, y_orig, empty_array, alphas, Ak, Ak_all, A_relative, theta_relative_x, theta_relative_y, method, case)


    return Ak;


## discretization size in each plot direction
N = 1000;

## axes_limits
range_min_x = -3;
range_max_x = 3;
range_min_y = -6.3;
range_max_y = 1.5;

###lim_x_line = 1e-3;
###x_line = np.linspace(-lim_x_line, lim_x_line, N);

orig_range_x = np.linspace(range_min_x, range_max_x, N);
orig_range_y = np.linspace(range_min_y, range_max_y, N);
real_range = orig_range_x;
imag_range = orig_range_y;

real_range, imag_range = np.meshgrid(real_range, imag_range);

cases = ["R0"];

## range kappa * dx (for computing stability regions of SL schemes)
kdxs = np.array([0., .5, 1., 1.5, 2.]) * np.pi;
kdxs_extended = np.append(kdxs, np.linspace(0, 2. * np.pi, 20));

empty_array = np.array([]);

## fixed values of xi_L
fixed_ys = np.array([0, 5e3, 1e4, 2.5e4]) * 2.5e-4 / np.pi;
multiple_pi = False

plot_dispersion = False;

linewidth = 1.

fine_tsms = ["IMEX"];
schemes = ["IMEX", "SL-SI-SETTLS"]

###############################################################
## plot stability regions of serial schemes for all y_values ##
###############################################################
for case in cases:
    Ak = {};
    for scheme in schemes:
        Ak[scheme] = {}
        for y_value in fixed_ys:
            if y_value == 0:
                key_str = r'$0$'
            else:
                key_str = r'$\xi_{\mathbf{L}} = $' +  r'${}$'.format(round(y_value * np.pi / 2.5e-4)) + r'$\tilde{\xi}_{\mathbf{L}}$'
            Ak[scheme][key_str] = computeStabilityRegionNonlinearProblem(real_range, imag_range, kdxs, kdxs_extended, empty_array, scheme, case, [y_value], plot = False);
    plotContourSeveralSchemes(y_value * np.pi, Ak, list(Ak.keys()), case, parareal = False, levels_dict = 2, multiple_pi = multiple_pi, linewidth = linewidth, main_style = 'linestyle', several_y = True);


###########################################################
## Plot parareal stability (nrelax = 0) along iterations ##
###########################################################
for nrelax in [0]:
    for fine_tsm in fine_tsms:
        for y_value in fixed_ys:
            for case in cases:
                Ak = {};
                for coarse_tsm in schemes:
                    Ak[coarse_tsm] = {};
                    for niter in [0, 1, 5, 10]:
                        for coarsening_factor in [2]:
                            for time_slice in [100]:
                                    filename = "_parareal_F" + fine_tsm + "_nrelax" + str(nrelax) + "__iterations";
                                    if ischeme == 1:
                                        filename += "__" + ischeme;
                                    Ak[coarse_tsm][r'$k = {}$'.format(niter)] = computeStabilityRegionNonlinearProblem(real_range, imag_range, kdxs_extended, empty_array, empty_array, filename, case, [y_value], plot = False, parareal = True, niter = niter, time_slice = time_slice, coarsening_factor = coarsening_factor, tsm_fine = fine_tsm, tsm_coarse = coarse_tsm, nrelax = nrelax);
                plotContourSeveralSchemes(y_value * np.pi, Ak, list(Ak.keys()), case, parareal = True, figname_parareal = filename, levels_dict = 2, multiple_pi = multiple_pi, linewidth = linewidth, main_style = 'linestyle');
################################
################################
################################    ## ín function of spatial coarsening
################################    for nrelax in [0]:
################################        for fine_tsm in fine_tsms:
################################            for y_value in fixed_ys:
################################                for case in cases:
################################                    Ak = {};
################################                    for coarse_tsm in schemes:
################################                        Ak[coarse_tsm] = {};
################################                        for niter in [5]:
################################                            for coarsening_factor in [2, 4, 8, 16]:
################################                                for time_slice in [100]:
################################                                        filename = "_parareal_F" + fine_tsm + "_nrelax" + str(nrelax) + "__coarsening";
################################                                        if ischeme == 1:
################################                                            filename += "__" + ischeme;
################################                                        Ak[coarse_tsm][r'$m_c = {}$'.format(coarsening_factor)] = computeStabilityRegionNonlinearProblem(real_range, imag_range, kdxs_extended, empty_array, empty_array, filename, case, [y_value], plot = False, parareal = True, niter = niter, time_slice = time_slice, coarsening_factor = coarsening_factor, tsm_fine = fine_tsm, tsm_coarse = coarse_tsm, nrelax = nrelax);
########################                    plotContourSeveralSchemes(y_value * np.pi, Ak, list(Ak.keys()), case, parareal = True, figname_parareal = filename, levels_dict = 2);
########################
########################
########################
########################    ## ín function of spatial coarsening, plot in function of dt instead of Dt
########################    for nrelax in [0]:
########################        for fine_tsm in fine_tsms:
########################            for y_value in fixed_ys:
########################                for case in cases:
########################                    Ak = {};
########################                    for coarse_tsm in schemes:
########################                        Ak[coarse_tsm] = {};
########################                        for niter in [5]:
########################                            ##for coarsening_factor in [1, 2, 4, 8]:
########################                            for coarsening_factor in [1/2, 1/4, 1/8, 1/16]:
########################                                for time_slice in [100]:
########################                                        filename = "_parareal_F" + fine_tsm + "_nrelax" + str(nrelax) + "__coarsening_cfactor";
########################                                        if ischeme == 1:
########################                                            filename += "__" + ischeme;
########################
########################                                        N_aux = N;
########################                                        orig_range_aux = np.linspace(range_min, range_max, N_aux);
########################                                        real_range_aux = orig_range_aux;
########################                                        imag_range_aux = orig_range_aux;
########################                                        real_range_aux, imag_range_aux = np.meshgrid(real_range_aux, imag_range_aux);
########################
########################                                        #########N2 = 1000;
########################                                        #########N_aux = N2 * coarsening_factor;
########################                                        ############range_min = -8.;
########################                                        ############range_max = 8.;
########################                                        #########range_min_aux = -6.3 * coarsening_factor;
########################                                        #########range_max_aux = 6.3 * coarsening_factor;
########################
########################                                        #########orig_range_aux = np.linspace(range_min_aux, range_max_aux, N_aux + 1);
########################                                        #########real_range_aux = orig_range_aux;
########################                                        #########imag_range_aux = orig_range_aux;
########################
########################                                        #########print(coarsening_factor, real_range_aux[0], real_range_aux[-1], len(real_range_aux))
########################
########################                                        #########real_range_aux = real_range_aux[N_aux//2 - N2//2: N_aux//2 + N2//2]
########################                                        #########imag_range_aux = imag_range_aux[N_aux//2 - N2//2: N_aux//2 + N2//2]
########################
########################                                        #########print(N2, N_aux)
########################                                        #########print(coarsening_factor, real_range_aux[0], real_range_aux[-1])
########################
########################                                        #########real_range_aux, imag_range_aux = np.meshgrid(real_range_aux, imag_range_aux);
########################
########################
########################                                        ###Ak[coarse_tsm][r'$m_c = {}$'.format(coarsening_factor)] = computeStabilityRegionNonlinearProblem(real_range_aux, imag_range_aux, kdxs, empty_array, empty_array, filename, case, [y_value], plot = False, parareal = True, niter = niter, time_slice = time_slice, coarsening_factor = coarsening_factor, tsm_fine = fine_tsm, tsm_coarse = coarse_tsm, nrelax = nrelax, invert_scale = True);
########################                                        Ak[coarse_tsm][r'$m_c = {}$'.format(coarsening_factor)] = computeStabilityRegionNonlinearProblem(real_range_aux, imag_range_aux, kdxs, empty_array, empty_array, filename, case, [y_value], plot = False, parareal = True, niter = niter, time_slice = time_slice, coarsening_factor = coarsening_factor, tsm_fine = fine_tsm, tsm_coarse = coarse_tsm, nrelax = nrelax, invert_scale = False);
########################
########################                    orig_range_min = range_min;
########################                    orig_range_max = range_max;
########################                    range_min = range_min;
########################                    range_max = range_max;
########################                    plotContourSeveralSchemes(y_value * np.pi, Ak, list(Ak.keys()), case, parareal = True, figname_parareal = filename, levels_dict = 2, invert_scale = True);
########################                    range_min = orig_range_min;
########################                    range_max = orig_range_max;
########################
########################
########################
########################
########################
########################
########################    ## in function of timeslice
########################    for nrelax in [0]:
########################        for fine_tsm in fine_tsms:
########################            for y_value in fixed_ys:
########################                for case in cases:
########################                    Ak = {};
########################                    for coarse_tsm in schemes:
########################                        Ak[coarse_tsm] = {};
########################                        for niter in [5]:
########################                            for coarsening_factor in [2]:
########################                                for time_slice in [10, 20, 50, 100]:
########################                                        filename = "_parareal_F" + fine_tsm +  "_nrelax" + str(nrelax) + "__timeslice";
########################                                        if ischeme == 1:
########################                                            filename += "__" + ischeme;
########################                                        Ak[coarse_tsm][r'$n = {}$'.format(time_slice)] = computeStabilityRegionNonlinearProblem(real_range, imag_range, kdxs_extended, empty_array, empty_array, filename, case, [y_value], plot = False, parareal = True, niter = niter, time_slice = time_slice, coarsening_factor = coarsening_factor, tsm_fine = fine_tsm, tsm_coarse = coarse_tsm, nrelax = nrelax);
########################                    plotContourSeveralSchemes(y_value * np.pi, Ak, list(Ak.keys()), case, parareal = True, figname_parareal = filename, levels_dict = 2);

## in function of nrelax
for fine_tsm in fine_tsms:
    for y_value in fixed_ys:
        for case in cases:
            Ak = {};
            for coarse_tsm in ["IMEX", "ETD2RK"]:
                Ak[coarse_tsm] = {}
                for nrelax in [0, 1, 2, 3]:
                    for niter in [5]:
                        for coarsening_factor in [2]:
                            for time_slice in [100]:
                                    filename = "_parareal_F" + fine_tsm + "__nrelax";
                                    Ak[coarse_tsm][r'$N_{{relax}} = {}$'.format(nrelax)] = computeStabilityRegionNonlinearProblem(real_range, imag_range, kdxs_extended, empty_array, empty_array, filename, case, [y_value], plot = False, parareal = True, niter = niter, time_slice = time_slice, coarsening_factor = coarsening_factor, tsm_fine = fine_tsm, tsm_coarse = coarse_tsm, nrelax = nrelax);
            plotContourSeveralSchemes(y_value * np.pi, Ak, list(Ak.keys()), case, parareal = True, figname_parareal = filename, levels_dict = 2, multiple_pi = multiple_pi, linewidth = linewidth, main_style = 'linestyle');


############################ in function of iteration for n coarse schemes
##################################coarse_tsms = ["IMEX", "SL-SI-SETTLS", "ETD2RK", "SL-ETD2RK"]
##########################coarse_tsms = ["IMEX", "SL-SI-SETTLS"]
##########################for fine_tsm in fine_tsms:
##########################    for y_value in fixed_ys:
##########################        for case in cases:
##########################            Ak = {};
##########################            for i in range(1,len(coarse_tsms) + 1):
##########################                print("CCCCCCCC", coarse_tsms[:i])
##########################                for coarse_tsm in coarse_tsms[:i]:
##########################                    Ak[coarse_tsm] = {}
##########################                    for nrelax in [0]:
##########################                        filename = "_parareal_F" + fine_tsm + "_nrelax" + str(nrelax) + "__iterations__C" + str(i);
##########################                        for niter in [0, 1, 5,10]:
##########################                            for coarsening_factor in [2]:
##########################                                for time_slice in [100]:
##########################                                        Ak[coarse_tsm][r'$k = {}$'.format(niter)] = computeStabilityRegionNonlinearProblem(real_range, imag_range, kdxs_extended, empty_array, empty_array, filename, case, [y_value], plot = False, parareal = True, niter = niter, time_slice = time_slice, coarsening_factor = coarsening_factor, tsm_fine = fine_tsm, tsm_coarse = coarse_tsm, nrelax = nrelax);
##########################                plotContourSeveralSchemes(y_value * np.pi, Ak, list(Ak.keys()), case, parareal = True, figname_parareal = filename, levels_dict = 2, multiple_pi = multiple_pi);
##########################
##########################
############################ in function of nrelax for n coarse schemes
##############################coarse_tsms = ["IMEX", "SL-SI-SETTLS", "ETD2RK", "SL-ETD2RK"]
##########################coarse_tsms = ["IMEX", "SL-SI-SETTLS"]
##########################for fine_tsm in fine_tsms:
##########################    for y_value in fixed_ys:
##########################        for case in cases:
##########################            Ak = {};
##########################            for i in range(1,len(coarse_tsms) + 1):
##########################                print("CCCCCCCC", coarse_tsms[:i])
##########################                filename = "_parareal_F" + fine_tsm + "__nrelax__C" + str(i);
##########################                for coarse_tsm in coarse_tsms[:i]:
##########################                    Ak[coarse_tsm] = {}
##########################                    for nrelax in [0, 1, 2, 3]:
##########################                        for niter in [2]:
##########################                            for coarsening_factor in [2]:
##########################                                for time_slice in [100]:
##########################                                        Ak[coarse_tsm][r'$N_{{relax}} = {}$'.format(nrelax)] = computeStabilityRegionNonlinearProblem(real_range, imag_range, kdxs_extended, empty_array, empty_array, filename, case, [y_value], plot = False, parareal = True, niter = niter, time_slice = time_slice, coarsening_factor = coarsening_factor, tsm_fine = fine_tsm, tsm_coarse = coarse_tsm, nrelax = nrelax);
##########################                plotContourSeveralSchemes(y_value * np.pi, Ak, list(Ak.keys()), case, parareal = True, figname_parareal = filename, levels_dict = 2, multiple_pi = multiple_pi);


## in function of nrelax only for ETD2RK
for fine_tsm in fine_tsms:
    for y_value in fixed_ys:
        for case in cases:
            Ak = {};
            for coarse_tsm in ["ETD2RK"]:
                for nrelax in [-1, 0, 1, 3, 4, 5, 10]:
                    for niter in [5]:
                        for coarsening_factor in [2]:
                            for time_slice in [100]:
                                    filename = "_parareal_F" + fine_tsm + "__" + coarse_tsm + "_only" + "__nrelax";
                                    if nrelax < 0:
                                        Ak[coarse_tsm] = computeStabilityRegionNonlinearProblem(real_range, imag_range, kdxs, kdxs_extended, empty_array, coarse_tsm, case, [y_value], plot = False);
                                    else:
                                        Ak["MGRIT; " + r'$N_{{relax}} = {}$'.format(nrelax)] = computeStabilityRegionNonlinearProblem(real_range, imag_range, kdxs_extended, empty_array, empty_array, filename, case, [y_value], plot = True, parareal = True, niter = niter, time_slice = time_slice, coarsening_factor = coarsening_factor, tsm_fine = fine_tsm, tsm_coarse = coarse_tsm, nrelax = nrelax);
                ###Ak[coarse_tsm] = {}
            plotContourSeveralSchemes(y_value * np.pi, Ak, list(Ak.keys()), case, parareal = True, figname_parareal = filename, levels_dict = 1, multiple_pi = multiple_pi, linewidth = linewidth, main_style = 'linestyle');
###################
##################### in function of nrelax and nlevels
###################for fine_tsm in fine_tsms:
###################    for y_value in fixed_ys:
###################        for case in cases:
###################            for nrelax in [0, 1, 2, 3, 4, 5]:
###################                Ak = {};
###################                for coarse_tsm in ["IMEX", "SL-SI-SETTLS", "ETD2RK", "SL-ETD2RK"]:
###################                    Ak[coarse_tsm] = {}
###################                    for nlevels in [2, 3, 4, 5]:
###################                        for niter in [5]:
###################                            for coarsening_factor in [2]:
###################                                for time_slice in [100]:
###################                                        filename = "_parareal_F" + fine_tsm + "_nrelax" + str(nrelax) + "__nlevels";
###################                                        Ak[coarse_tsm][r'$N_{{levels}} = {}$'.format(nlevels)] = computeStabilityRegionNonlinearProblem(real_range, imag_range, kdxs, empty_array, empty_array, filename, case, [y_value], plot = False, parareal = True, niter = niter, time_slice = time_slice, coarsening_factor = coarsening_factor, tsm_fine = fine_tsm, tsm_coarse = coarse_tsm, nrelax = nrelax, nlevels = nlevels);
###################                plotContourSeveralSchemes(y_value * np.pi, Ak, list(Ak.keys()), case, parareal = True, figname_parareal = filename, levels_dict = 2);


##########for nrelax in [0, 1]:
##########    for y_value in fixed_ys:
##########        for case in cases:
##########            Ak = {};
##########            for fine_tsm in ["IMEX"]:
##########                for coarse_tsm in ["IMEX", "ETD2RK", "SL-ETD2RK", "SL-SI-SETTLS"]:
##########                    for niter in [5]:
##########                        for coarsening_factor in [2]:
##########                            for time_slice in [10, 20, 50, 100]:
##########                                    filename = "_parareal_" + "_nrelax" + str(nrelax) + "__tsm_coarse";
##########                                    Ak[coarse_tsm] = computeStabilityRegionNonlinearProblem(real_range, imag_range, kdxs, empty_array, empty_array, filename, case, [y_value], plot = False, parareal = True, niter = niter, time_slice = time_slice, coarsening_factor = coarsening_factor, tsm_fine = fine_tsm, tsm_coarse = coarse_tsm, nrelax = nrelax);
##########            plotContourSeveralSchemes(y_value * np.pi, Ak, list(Ak.keys()), case, parareal = True, figname_parareal = filename);


###for scheme in schemes_alpha:
###    for case in cases:
###        computeStabilityRegionNonlinearProblem(real_range, imag_range, empty_array, empty_array, alphas, scheme, case, fixed_ys = fixed_ys);


###for scheme in schemes_kdx:
###    for case in cases:
###        Ak = computeStabilityRegionNonlinearProblem(real_range, imag_range, kdxs, kdxs_extended, empty_array, scheme, case, fixed_ys);
###        print(Ak)

###for y_value in fixed_ys:
###    for case in cases:
###        Ak = {};
###        for scheme in schemes_kdx:
###            Ak[scheme] = computeStabilityRegionNonlinearProblem(real_range, imag_range, kdxs, kdxs, empty_array, scheme, case, [y_value], plot = False)[kdxs[1]];
###            print(case, scheme, y_value, Ak[scheme])
###        plotContourSeveralSchemes(y_value * np.pi, Ak, schemes_kdx, case);


##for scheme in schemes_kdx_alpha:
##    for case in cases:
##        computeStabilityRegionNonlinearProblem(real_range, imag_range, kdxs, kdxs_extended, alphas, scheme, case);

##computeStabilityRegionLinearProblem(real_range, imag_range, kdxs, "1997");
##computeStabilityRegionLinearProblem(real_range, imag_range, kdxs, empty_array, "SETTLS");
##computeStabilityRegionNonlinearProblem(real_range, imag_range, kdxs, "SL-EXP-SETTLS", "RR");
##computeStabilityRegionNonlinearProblem(real_range, imag_range, kdxs, "SL-ETD2RK", "RR");
##computeStabilityRegionNonlinearProblem(real_range, imag_range, np.array([]), "CN", "RR");
##computeStabilityRegionNonlinearProblem(real_range, imag_range, np.array([]), "AB2AM2", "RR");
##computeStabilityRegionNonlinearProblem(real_range, imag_range, np.array([]), "CN", "RI");
##computeStabilityRegionNonlinearProblem(real_range, imag_range, np.array([]), "AB2AM2", "RI");
##computeStabilityRegionNonlinearProblem(real_range, imag_range, np.array([]), "CN", "IR");
##computeStabilityRegionNonlinearProblem(real_range, imag_range, np.array([]), "AB2AM2", "IR");
##computeStabilityRegionNonlinearProblem(real_range, imag_range, np.array([]), "CN", "II");
##computeStabilityRegionNonlinearProblem(real_range, imag_range, np.array([]), "AB2AM2", "II");

##computeStabilityRegionNonlinearProblem(real_range, imag_range, kdxs, "SL-EXP-SETTLS", "RI");
##computeStabilityRegionNonlinearProblem(real_range, imag_range, kdxs, "SL-ETD2RK", "RI");
##computeStabilityRegionNonlinearProblem(real_range, imag_range, kdxs, "SL-EXP-SETTLS", "IR");
##computeStabilityRegionNonlinearProblem(real_range, imag_range, kdxs, "SL-ETD2RK", "IR");
##computeStabilityRegionNonlinearProblem(real_range, imag_range, kdxs, "SL-EXP-SETTLS", "II");
##computeStabilityRegionNonlinearProblem(real_range, imag_range, kdxs, "SL-ETD2RK", "II");

## plot 1D, few kdxs
##computeStabilityRegionLinearProblem(0. * orig_range, orig_range, kdxs, "1997");
##computeStabilityRegionLinearProblem(0. * orig_range, orig_range, kdxs, "SETTLS");
##computeStabilityRegionNonlinearProblem(orig_range, orig_range, kdxs, "SL-EXP-SETTLS", "II");
##computeStabilityRegionNonlinearProblem(orig_range, orig_range, kdxs, "SL-ETD2RK", "II");
##computeStabilityRegionLinearProblem(0. * orig_range, orig_range, np.array([]), "CN");
##computeStabilityRegionLinearProblem(0. * orig_range, orig_range, np.array([]), "AB2AM2");


## plot 2D, several kdxs
kdxs = np.linspace(0., 4., N) * np.pi;
##computeStabilityRegionLinearProblem(real_range, imag_range, kdxs, "1997");
##computeStabilityRegionLinearProblem(real_range, imag_range, kdxs, "SETTLS");
##computeStabilityRegionNonlinearProblem(real_range, imag_range, kdxs, "SL-EXP-SETTLS", "RR");
##computeStabilityRegionNonlinearProblem(real_range, imag_range, kdxs, "SL-ETD2RK", "RR");
##computeStabilityRegionNonlinearProblem(real_range, imag_range, np.array([]), "CN", "RR");
##computeStabilityRegionNonlinearProblem(real_range, imag_range, np.array([]), "AB2AM2", "RR");

##computeStabilityRegionNonlinearProblem(real_range, imag_range, kdxs, "SL-EXP-SETTLS", "RI");
##computeStabilityRegionNonlinearProblem(real_range, imag_range, kdxs, "SL-ETD2RK", "RI");
##computeStabilityRegionNonlinearProblem(real_range, imag_range, kdxs, "SL-EXP-SETTLS", "IR");
##computeStabilityRegionNonlinearProblem(real_range, imag_range, kdxs, "SL-ETD2RK", "IR");
##computeStabilityRegionNonlinearProblem(real_range, imag_range, kdxs, "SL-EXP-SETTLS", "II");
##computeStabilityRegionNonlinearProblem(real_range, imag_range, kdxs, "SL-ETD2RK", "II");

