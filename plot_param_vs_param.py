import argparse
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import QTable

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
    parser.add_argument("--rvebv", help="sample properties", action="store_true")
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
    for cset in args.datasets:
        fname = f"data/{cset}_ensemble_params.dat"
        allnames.append(cset)
        tdata = QTable.read(fname, format="ascii.ipac")

        # now add data if missing and derivable from expected columns
        if "B3" not in tdata.colnames:
            tdata["B3"] = tdata["C3"] / (tdata["gamma"] ** 2)
            if "C3_unc" in tdata.colnames:
                tdata["B3_unc"] = np.absolute(tdata["B3"]) * np.sqrt(tdata["C3_unc"] ** 2 +  2.0 * (tdata["gamma_unc"].value ** 2))
 
        if "IRV" not in tdata.colnames:
            tdata["IRV"] = 1. / tdata["RV"]
            tdata["IRV_unc"] = tdata["IRV"] * tdata["RV_unc"] / tdata["RV"]

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

        if args.showstats:
            print(f"Summary statistics for {cset}")
            for ccol in tdata.colnames:
                if ("name" not in ccol.lower()) & ("unc" not in ccol.lower()):
                    ave = np.average(tdata[ccol].data)
                    std = np.std(tdata[ccol].data)
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
    yplabels = ["$C_1$", "$C_2$", "$B_3 = C_3/\gamma^2$", "$C_4$", "$x_o$", r"$\gamma$"]
    yptags = ["C1", "C2", "B3", "C4", "x0", "gamma"]
    if args.sprops:
        ostr = "sprops"
        fsize = (12, 10)
        nrows = 2
        ncols = 2
        pi = [0, 1, 2, 3]
        xplabels = ["$E(B-V)$", "$A(V)$", "$E(B-V)$", "$A(V)$"]
        xptags = ["EBV", "AV", "EBV", "AV"]
        yplabels = ["$R(V)$", "$R(V)$", "$N(HI)/E(B-V)$", "$N(HI)/A(V)$"]
        yptags = ["RV", "RV", "NHI_EBV", "NHI_AV"]
    elif args.rvebv:
        ostr = "rvebv"
        fsize = (12, 6)
        nrows = 1
        ncols = 2
        pi = [0, 1]
        xplabels = ["$E(B-V)$", "$E(B-V)$"]
        xptags = ["EBV", "EBV"]
        yplabels = ["$R(V)$", "$N(HI)/E(B-V)$"]
        yptags = ["RV", "NHI_EBV"]
    elif args.gdprops:
        ostr = "gdprops"
        fsize = (12, 10)
        nrows = 2
        ncols = 2
        pi = [0, 1, 2, 3]
        xplabels = ["$A(V)$", "$C_2$", "$B_3$", "$C_4$"]
        xptags = ["AV", "C2", "B3", "C4"]
        yplabels = ["$N(HI)$", "$N(HI)/A(V)$", "$N(HI)/A(V)$", "$N(HI)/A(V)$"]
        yptags = ["NHI", "NHI_AV", "NHI_AV", "NHI_AV"]
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
        xplabels = ["$C_2$", "$C_2$", "$C_2$", "$C_4$"]
        xptags = ["C2", "C2", "C2", "C4"]
        yplabels = ["$C_1$", "$B_3 = C_3/\gamma^2$", "$C_4$", "$B_3 = C_3/\gamma^2$"]
        yptags = ["C1", "B3", "C4", "B3"]
    else:  # plot fm90 vs fm90
        ostr = "fm90"
        xplabels = ["$C_2$", "$C_2$", "$C_2$", "$C_2$", "$B_3$", "$x_0$"] #, "$C_4$"]
        xptags = ["C2", "C2", "C2", "C2", "B3", "x0"]   # , "C4"]
        yplabels = ["$C_1$", "$B_3 = C_3/\gamma^2$", "$C_4$", "$x_0$", r"$\gamma$", "$\gamma$"]   # , "$B_3 = C_3/\gamma^2$"]
        yptags = ["C1", "B3", "C4", "x0", "gamma", "gamma"]  # , "B3"]

    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=fsize)

    for cname, cdata in zip(allnames, alldata):
        ptype, palpha, plabel = ptypes[cname]
        for i in range(nrows * ncols):
            # check if uncertainties are included
            if f"{xptags[i]}_unc" in cdata.colnames:
                xdata_unc = cdata[f"{xptags[i]}_unc"]
            else:
                xdata_unc = None
            if f"{yptags[i]}_unc" in cdata.colnames:
                ydata_unc = cdata[f"{yptags[i]}_unc"]
            else:
                ydata_unc = None

            if args.nouncs & (len(cdata[xptags[i]]) > 50):
                xdata_unc = None
                ydata_unc = None

            colstr = ptype[0]
            symstr = ptype[1]

            px, py = divmod(pi[i], ncols)
            ax[px, py].errorbar(
                cdata[xptags[i]],
                cdata[yptags[i]],
                xerr=xdata_unc,
                yerr=ydata_unc,
                color=colstr,
                marker=symstr,
                linestyle="",
                label=plabel,
                alpha=palpha,
            )

    for i in range(nrows * ncols):
        px, py = divmod(pi[i], ncols)
        ax[px, py].set_xlabel(xplabels[i], fontsize=1.3 * fontsize)
        ax[px, py].set_ylabel(yplabels[i], fontsize=1.3 * fontsize)
        if yptags[i] in ["NHI_EBV", "NHI_AV"]:
            ax[px, py].set_yscale("log")

    ax[0, 0].legend(fontsize=0.7*fontsize)

    fig.tight_layout()

    fname = f"ensemble_{ostr}_vs_fm90_params"
    if args.png:
        fig.savefig(f"{fname}.png")
    elif args.pdf:
        fig.savefig(f"{fname}.pdf")
    else:
        plt.show()
