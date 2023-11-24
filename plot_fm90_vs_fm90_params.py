import argparse
import matplotlib.pyplot as plt
from astropy.table import QTable

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--datasets",
        help="give the datasets to plot",
        nargs="+",
        default=["val04", "gor03_smc", "gor03_lmc"],
        choices=["val04", "gor03_smc", "gor03_lmc", "gor24_smc", "gor24_smc_forecor"],
    )
    parser.add_argument("--av", help="plot versus A(V)", action="store_true")
    parser.add_argument("--rv", help="plot versus R(V)", action="store_true")
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    # get the data to plot
    allnames = []
    alldata = []
    for cset in args.datasets:
        fname = f"data/{cset}_ensemble_params.dat"
        allnames.append(cset)
        tdata = QTable.read(fname, format="ascii.ipac")
        if "B3" not in tdata.colnames:
            tdata["B3"] = tdata["C3"] / (tdata["gamma"] ** 2)
        alldata.append(tdata)

    # make the plots
    fontsize = 14
    font = {"size": fontsize}
    plt.rc("font", **font)
    plt.rc("lines", linewidth=2)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("ytick.major", width=2)

    fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(18, 10))

    xplabels = ["$C_2$", "$C_2$", "$C_2$", "$C_2$", "$C_2$", "$C_4$"]
    xptags = ["C2", "C2", "C2", "C2", "C2", "C4"]
    yplabels = ["$C_1$", "$B_3$", "$C_4$", "$x_o$", r"$\gamma$", "$B_3$"]
    yptags = ["C1", "B3", "C4", "x0", "gamma", "B3"]
    pi = [0, 1, 2, 3, 4, 5]

    # plot types, colors and alphas
    ptypes = {
        "val04": ("ko", 0.25),
        "gor03_smc": ("bv", 0.5),
        "gor03_lmc": ("c^", 0.5),
        "gor24_smc": ("r>", 0.2),
        "gor24_smc_forecor": ("g<", 0.5),
    }

    for cname, cdata in zip(allnames, alldata):
        ptype, palpha = ptypes[cname]
        for i in range(6):
            # check if uncertainties are included
            if f"{xptags[i]}_unc" in cdata.colnames:
                xdata_unc = cdata[f"{xptags[i]}_unc"]
            else:
                xdata_unc = None
            if f"{yptags[i]}_unc" in cdata.colnames:
                ydata_unc = cdata[f"{yptags[i]}_unc"]
            else:
                ydata_unc = None

            px, py = divmod(pi[i], 3)
            ax[px, py].errorbar(
                cdata[xptags[i]],
                cdata[yptags[i]],
                xerr=xdata_unc,
                yerr=ydata_unc,
                fmt=ptype,
                label=cname,
                alpha=palpha,
            )

    for i in range(6):
        px, py = divmod(pi[i], 3)
        ax[px, py].set_xlabel(xplabels[i], fontsize=1.3 * fontsize)
        ax[px, py].set_ylabel(yplabels[i], fontsize=1.3 * fontsize)

    ax[0, 0].legend()

    fig.tight_layout()

    fname = "ensemble_fm90_vs_fm90_params"
    if args.png:
        fig.savefig(f"{fname}.png")
    elif args.pdf:
        fig.savefig(f"{fname}.pdf")
    else:
        plt.show()
