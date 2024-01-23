import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as op
from astropy.table import QTable
from astropy.modeling import models, fitting
from astropy.stats import sigma_clip

from utils.fit_full2dcor import lnlike_correlated
from utils.plot_cov_ellipses import draw_ellipses

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--datasets",
        help="give the datasets to plot",
        nargs="+",
        default=["val04", "fit07", "gor09"],
        choices=["val04", "gor03_smc", "gor03_lmc", "fit07", "gor09", "gor24_smc",
                 "gor24_smc_nobump", "gor24_smc_bump", "gor24_smc_flat", "gor24_smc_lowebv"],
    )
    parser.add_argument("--sprops", help="sample properties", action="store_true")
    parser.add_argument("--spropsebv", help="sample properties versus ebv", action="store_true")
    parser.add_argument("--spropsav", help="sample properties versus av", action="store_true")
    parser.add_argument("--gdprops", help="N(HI)/E(B-V) properties", action="store_true")
    parser.add_argument("--fm90main", help="only plot the main FM90 parameters", action="store_true")
    parser.add_argument("--ebv", help="plot FM90 versus E(B-V)", action="store_true")
    parser.add_argument("--av", help="plot FM90 versus A(V)", action="store_true")
    parser.add_argument("--nhi", help="plot FM90 versus N(HI)", action="store_true")
    parser.add_argument("--rv", help="plot FM90 versus R(V)", action="store_true")
    parser.add_argument("--irv", help="plot FM90 versus 1/R(V)", action="store_true")
    parser.add_argument("--nouncs", help="do not plot uncs", action="store_true")
    parser.add_argument("--ebvcut", help="only plot data equal or above E(B-V) value",
                        type=float, default=0.0)
    parser.add_argument("--showstats", help="print summary stats for each dataset", action="store_true")    
    parser.add_argument("--fit", help="Fit lines for some plots", action="store_true")
    parser.add_argument("--paper", help="portrait format", action="store_true")
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    # plot types, colors and alphas
    ptypes = {
        "val04": ("k.", 0.1, "MW: VCG04"),
        "gor03_smc": ("m<", 0.25, "SMC: G03"),
        "gor03_lmc": (("tab:orange", ">"), 0.25, "LMC: G03"),
        "fit07": ("k+", 0.1, "MW: FM07"),
        "gor09": ("kD", 0.25, "MW: GCC09"),
        "gor24_smc": ("b>", 0.5, "SMC: G24"),
        "gor24_smc_nobump": ("bo", 0.5, "SMC: Weak/absent 2175 A bump"),
        "gor24_smc_bump": ("rP", 0.5, "SMC: Significant 2175 A bump"),
        "gor24_smc_flat": ("cs", 0.5, "SMC: Flat"),
        "gor24_smc_lowebv": (("tab:brown", "v"), 0.5, r"SMC: $E(B-V)_\mathrm{SMC} < 0.1$"),
    }

    # get the data to plot
    allnames = []
    alldata = []
    summarystats = {}
    for cset in args.datasets:
        fname = f"data/{cset}_ensemble_params.dat"
        allnames.append(cset)
        tdata = QTable.read(fname, format="ascii.ipac")

        # now add data if missing and derivable from expected columns
        if "B3" not in tdata.colnames:
            tdata["B3"] = tdata["C3"] / (tdata["gamma"] ** 2)
            if "C3_unc" in tdata.colnames:
                tdata["B3_unc"] = np.absolute(tdata["B3"]) * np.sqrt(tdata["C3_unc"] ** 2 +  2.0 * (tdata["gamma_unc"].value ** 2))
            # temp fix until LMC and MW GCC09 can be refit with B3
            # C3 and gamma strongly correlated
            if args.fit:
                tdata["B3_unc"] *= 0.2

        if "IRV" not in tdata.colnames:
            tdata["IRV"] = 1. / tdata["RV"]
            tdata["IRV_unc"] = tdata["IRV"] * tdata["RV_unc"] / tdata["RV"]

        # divide by 10^21 to make easier to undertand and fit numbers
        if "NHI" in tdata.colnames:
            tdata["NHI"] /= 1e21
            tdata["NHI_unc"] /= 1e21

        if ("NHI" in tdata.colnames) & ("NHI_EBV" not in tdata.colnames):
            tdata["NHI_EBV"] = tdata["NHI"] / tdata["EBV"]
            tdata["NHI_EBV_unc"] = (tdata["NHI_unc"] / tdata["NHI"]) **2 + (tdata["EBV_unc"] / tdata["EBV"]) **2
            tdata["NHI_EBV_unc"] = tdata["NHI_EBV"] * np.sqrt(tdata["NHI_EBV_unc"])

        if ("NHI" in tdata.colnames) & ("NHI_AV" not in tdata.colnames):
            tdata["NHI_AV"] = tdata["NHI"] / tdata["AV"]
            tdata["NHI_AV_unc"] = (tdata["NHI_unc"] / tdata["NHI"]) **2 + (tdata["AV_unc"] / tdata["AV"]) **2
            tdata["NHI_AV_unc"] = tdata["NHI_AV"] * np.sqrt(tdata["NHI_AV_unc"])

        if args.ebvcut > 0.0:
            tdata = tdata[tdata["EBV"].value >= args.ebvcut]

        alldata.append(tdata)

        summarystats[cset] = {}
        if args.showstats:
            print(f"Summary statistics for {cset}")
        for ccol in tdata.colnames:
            if ("name" not in ccol.lower()) & ("unc" not in ccol.lower()):
                ave = np.average(tdata[ccol].data)
                std = np.std(tdata[ccol].data)
                summarystats[cset][ccol] = (ave, std)
                if args.showstats:
                    print(f"{ccol}: ave = {ave} +/- {std}")

    # make the plots
    fontsize = 14
    font = {"size": fontsize}
    plt.rc("font", **font)
    plt.rc("lines", linewidth=2)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("ytick.major", width=2)

    if args.paper:
        fsize = (12, 14)
        nrows = 3
        ncols = 2
        pi = [0, 1, 2, 3, 4, 5]
    else:
        fsize = (18, 10)
        nrows = 2
        ncols = 3
        pi = [0, 1, 3, 4, 2, 5]

    # default values
    yplabels = ["$C_1$ = UV intercept", "$C_2$ = UV slope", "$B_3$ = bump amplitude",
                "$C_4$ = FUV rise amplitude", "$x_o$ = bump center", r"$\gamma$ = bump width"]
    yptags = ["C1", "C2", "B3", "C4", "x0", "gamma"]
    fitlines = [False] * (nrows * ncols)
    show_gd = None
    if args.sprops:
        ostr = "sprops"
        fsize = (12, 10)
        nrows = 2
        ncols = 2
        pi = [0, 1, 2, 3]
        xplabels = ["$E(B-V)$", "$A(V)$", "$E(B-V)$", "$A(V)$"]
        xptags = ["EBV", "AV", "EBV", "AV"]
        yplabels = ["$R(V)$", "$R(V)$", "$N(HI)/E(B-V)$ [$10^{21}$]", "$N(HI)/A(V)$ [$10^{21}$]"]
        yptags = ["RV", "RV", "NHI_EBV", "NHI_AV"]
    elif args.spropsebv:
        ostr = "sprops_ebv"
        fsize = (14, 6)
        nrows = 1
        ncols = 2
        pi = [0, 1]
        xplabels = ["$E(B-V)$", "$E(B-V)$"]
        xptags = ["EBV", "EBV"]
        yplabels = ["$R(V)$", "$N(HI)$ [$10^{21}$]"]
        yptags = ["RV", "NHI"]
        show_gd = [False, True]
    elif args.spropsav:
        ostr = "sprops_av"
        fsize = (14, 6)
        nrows = 1
        ncols = 2
        pi = [0, 1]
        xplabels = ["$A(V)$", "$A(V)$"]
        xptags = ["AV", "AV"]
        yplabels = ["$R(V)$", "$N(HI)$ [$10^{21}$]"]
        yptags = ["RV", "NHI"]
        show_gd = [False, True]
    elif args.gdprops:
        ostr = "gdprops"
        fsize = (12, 10)
        nrows = 2
        ncols = 2
        pi = [0, 1, 2, 3]
        xplabels = ["$A(V)$", "$C_2$ = UV slope", "$B_3$ = bump amplitude", "$C_4$ = FUV rise amplitude"]
        xptags = ["AV", "C2", "B3", "C4"]
        yplabels = ["$N(HI)/A(V)$ [$10^{21}$]", "$N(HI)/A(V)$ [$10^{21}$]",
                    "$N(HI)/A(V)$ [$10^{21}$]", "$N(HI)/A(V)$ [$10^{21}$]"]
        yptags = ["NHI_AV", "NHI_AV", "NHI_AV", "NHI_AV"]
        fitlines = [False, True, True, True]
    elif args.ebv:
        ostr = "ebv"
        npts = len(yplabels)
        xplabels = ["$E(B-V)$"] * npts
        xptags = ["EBV"] * npts
    elif args.av:
        ostr = "av"
        npts = len(yplabels)
        xplabels = ["$A(V)$"] * npts
        xptags = ["AV"] * npts
    elif args.nhi:
        ostr = "nhi"
        npts = len(yplabels)
        xplabels = ["$N(HI)$"] * npts
        xptags = ["NHI"] * npts
    elif args.rv:
        ostr = "rv"
        npts = len(yplabels)
        xplabels = ["$R(V)$"] * npts
        xptags = ["RV"] * npts
    elif args.irv:
        ostr = "irv"
        npts = len(yplabels)
        xplabels = ["1/$R(V)$"] * npts
        xptags = ["IRV"] * npts
    elif args.fm90main:
        ostr = "fm90main"
        fsize = (12, 10)
        nrows = 2
        ncols = 2
        pi = [0, 1, 2, 3]
        fitlines = [True] * nrows * ncols
        xplabels = ["$C_2$ = UV slope", "$C_2$ = UV slope", "$C_2$ = UV slope", "$C_4$ = FUV rise amplitude"]
        xptags = ["C2", "C2", "C2", "C4"]
        yplabels = ["$C_1$ = UV intercept", "$B_3$ = bump amplitude", "$C_4$ = FUV rise amplitude",
                    "$B_3$ = bump amplitude"]
        yptags = ["C1", "B3", "C4", "B3"]
    else:  # plot fm90 vs fm90
        ostr = "fm90"
        xplabels = ["$C_2$ = UV slope", "$C_2$ = UV slope", "$C_2$ = UV slope", "$C_2$ = UV slope",
                    "$B_3$ = bump amplitude", "$x_o$ = bump center [$\mu$m$^{-1}$]"]
        xptags = ["C2", "C2", "C2", "C2", "B3", "x0"]
        yplabels = ["$C_1$ = UV intercept", "$B_3$ = bump amplitude",
                    "$C_4$ = FUV rise amplitude", "$x_0$ = bump center [$\mu$m$^{-1}$]",
                    r"$\gamma$ = bump width [$\mu$m$^{-1}$]",
                    "$\gamma$ = bump width [$\mu$m$^{-1}$]"]
        yptags = ["C1", "B3", "C4", "x0", "gamma", "gamma"]

    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=fsize)

    for i in range(nrows * ncols):
        xvals = []
        xvals_unc = []
        yvals = []
        yvals_unc = []
        ebv = []
        ebv_unc = []
        for cname, cdata in zip(allnames, alldata):
            ptype, palpha, clabel = ptypes[cname]
            if i > 0:
                plabel = None
            else:
                plabel = clabel

            xdata = cdata[xptags[i]].data
            ydata = cdata[yptags[i]].data

            # for fitting
            xvals = np.concatenate((xvals, xdata))
            yvals = np.concatenate((yvals, ydata))

            # check if uncertainties are included
            if f"{xptags[i]}_unc" in cdata.colnames:
                xdata_unc = cdata[f"{xptags[i]}_unc"].data
                xvals_unc = np.concatenate((xvals_unc, xdata_unc))
            else:
                xdata_unc = None
            if f"{yptags[i]}_unc" in cdata.colnames:
                ydata_unc = cdata[f"{yptags[i]}_unc"].data
                yvals_unc = np.concatenate((yvals_unc, ydata_unc))
            else:
                ydata_unc = None

            if args.nouncs & (len(cdata[xptags[i]]) > 50):
                xdata_unc = None
                ydata_unc = None

            # for covariances
            ebv = np.concatenate((ebv, cdata["EBV"].data))
            ebv_unc = np.concatenate((ebv_unc, cdata["EBV_unc"].data))

            colstr = ptype[0]
            symstr = ptype[1]

            if (nrows == 1) | (ncols == 1):
                tax = ax[i]
            else:
                px, py = divmod(pi[i], ncols)
                tax = ax[px, py]

            tax.errorbar(
                xdata,
                ydata,
                xerr=xdata_unc,
                yerr=ydata_unc,
                color=colstr,
                marker=symstr,
                linestyle="",
                label=plabel,
                alpha=palpha,
            )

            # special code to fit for the gas-to-dust ratio
            if (xptags[i] in ["EBV", "AV"]) & (yptags[i] == "NHI") & (show_gd is not None):

                if show_gd[i] & (len(xdata) > 10):
                    xlim = tax.get_xlim()
                    ylim = tax.get_ylim()
                    x = np.arange(0.0, xlim[1], 0.01)
                    gdratio = summarystats[cname][f"NHI_{xptags[i]}"][0]
                    line_orig = models.Linear1D(slope=gdratio, intercept=0.0)
                    tax.plot(x, line_orig(x), linestyle=":", label=f"$N(HI)$/{xplabels[i]} = {gdratio:.2f}: {clabel}", color=colstr)
                    tax.legend(fontsize=0.7*fontsize, loc="upper right")
                if show_gd[i] & (cname == allnames[-1]):
                    tax.set_ylim(ylim)
                    # temp hack for SMC UV paper
                    if xptags[i] == "AV":
                        tax.set_xlim(xlim[0], 3.5)
                        tax.set_ylim(ylim[0], 20.)

        tax.set_xlabel(xplabels[i], fontsize=1.3 * fontsize)
        tax.set_ylabel(yplabels[i], fontsize=1.3 * fontsize)
        if yptags[i] in ["NHI_EBV", "NHI_AV"]:
            tax.set_yscale("log")

        # now fit a line to the data
        npts = len(xvals)
        covs = np.zeros((npts, 2, 2))

        # needed for plotting and fitting
        covtags = ["C1", "C2", "B3", "C4"]
        for k in range(npts):
            #if (xptags[i] in covtags) & (yptags[i] in covtags):
            #    # approximation following Gordon et al. (2023)
            #    cov_xy = xvals[k] * yvals[k] * ((ebv_unc[k] / ebv[k]) ** 2)
            #    corr_xy = np.min([cov_xy / (xvals_unc[k] * yvals_unc[k]), 0.99])
            #    cov_xy *= xvals_unc[k] * yvals_unc[k]
            #else:
            #    cov_xy = 0.0
            cov_xy = 0.0
            covs[k, 0, 0] = xvals_unc[k] ** 2
            covs[k, 0, 1] = cov_xy
            covs[k, 1, 0] = cov_xy
            covs[k, 1, 1] = yvals_unc[k] ** 2

            if not np.all(np.linalg.eigvals(covs[k, :, :]) > 0):
                print("eigvals")
                print(xptags[i], yptags[i])
                print(k, np.all(np.linalg.eigvals(covs[k, :, :]) > 0))
                print(covs[k, :, :])

            if np.linalg.cond(covs[k, :, :]) > 1/sys.float_info.epsilon:
                print("cond")
                print(xptags[i], yptags[i])
                print(k, np.all(np.linalg.cond(covs[k, :, :])))
                print(covs[k, :, :])   

        #draw_ellipses(
        #    tax, xvals, yvals, covs, color="black", alpha=0.1
        #)

        if fitlines[i] & args.fit:
      
            xlim = tax.get_xlim()
            ylim = tax.get_ylim()
            dxlim = xlim[1] - xlim[0]
            intinfo = [xlim[0] - dxlim, xlim[1] + dxlim, dxlim/20]

            def nll(*args):
                return -lnlike_correlated(*args)

            x = np.arange(xlim[0], xlim[1], 0.01)
            line_orig = models.Linear1D(slope=-10., intercept=20.)
            fit = fitting.LinearLSQFitter()
            or_fit = fitting.FittingWithOutlierRemoval(fit, sigma_clip, niter=3, sigma=3.0)
            fitted_line= fit(line_orig, xvals, yvals, weights=1/yvals_unc)
            nparams = fitted_line.parameters
            #tax.plot(x, fitted_line(x), "k:", label=f"Fit: {nparams[1]:.2f} + {nparams[0]:.2f}x")

            #fitted_line, mask = or_fit(fitted_line, xvals, yvals, weights=1/yvals_unc)
            # print(fitted_line.parameters)

            #masked_data = np.ma.masked_array(yvals, mask=~mask)
            #print(xvals[mask])
            #tax.plot(xvals, masked_data, "ko", fillstyle="none", ms=10, label="Not used in fit")

            result = op.minimize(nll, fitted_line.parameters, args=(yvals, fitted_line, covs, intinfo, xvals))
            nparams = result["x"]
            # print(nparams)
            fitted_line = models.Linear1D(slope=nparams[0], intercept=nparams[1])
            tax.plot(x, fitted_line(x), "k--", label=f"Fit: {nparams[1]:.2f} + {nparams[0]:.2f}x")

            tax.set_ylim(ylim)

            tax.legend(fontsize=0.7*fontsize)

        elif i == 0: 
            tax.legend(fontsize=0.7*fontsize)

    fig.tight_layout()

    fname = f"ensemble_{ostr}_vs_fm90_params"
    if args.png:
        fig.savefig(f"{fname}.png")
    elif args.pdf:
        fig.savefig(f"{fname}.pdf")
    else:
        plt.show()
