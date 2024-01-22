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
        choices=["val04", "gor03_smc", "gor03_lmc", "gor24_smc", "gor24_smc_forecor", "m31", "m33", "cla15"],
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
        alldata.append(QTable.read(fname, format="ascii.ipac"))

    # make the plots
    fontsize = 14
    font = {"size": fontsize}
    plt.rc("font", **font)
    plt.rc("lines", linewidth=2)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("ytick.major", width=2)

    fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(18, 10))

    plabels = ["$C_1$", "$C_2$", "$C_3$", "$C_4$", "$x_o$", r"$\gamma$"]
    yptags = ["C1", "C3", "C4", "x0", "gamma", "C3"]
    xptags = ["C2", "C2", "C2", "C2", "C2", "C4"]
    pi = [0, 1, 2, 3, 4, 5]

    # plot types, colors and alphas
    ptypes = {
        "val04": ("ko", 0.25),
        "gor03_smc": ("bv", 0.5),
        "gor03_lmc": ("c^", 0.5),
        "gor24_smc": ("r>", 0.2),
        "gor24_smc_forecor": ("g<", 0.5),
        "m31": ("yD", 0.8),
        "m33": ("mo", 0.9),
        "cla15": ("rd", 0.8)
    }

    for cname, cdata in zip(allnames, alldata):

        ptype, palpha = ptypes[cname]
        for i in range(6):
            px, py = divmod(pi[i], 3)

            if "C3" in yptags[i]:
                print("c3", cdata[yptags[i]])
                print("g**2", cdata[yptags[4]] ** 2)
                c3g2 = cdata[yptags[i]] / (cdata[yptags[4]] ** 2)
                ax[px, py].plot(
                    cdata[xptags[i]], c3g2, ptype, label=cname, alpha=palpha
                )
            else:
                ax[px, py].plot(
                    cdata[xptags[i]], cdata[yptags[i]], ptype, label=cname, alpha=palpha
                )

    for i in range(6):
        px, py = divmod(pi[i], 3)

        ax[px, py].set_xlabel(xptags[i], fontsize=1.3 * fontsize)
        ax[px, py].set_ylabel(yptags[i], fontsize=1.3 * fontsize)

        if "C3" in xptags[i]:
            ax[px, py].set_xlabel("C3/$\gamma^2$", fontsize=1.3 * fontsize)
        if "C3" in yptags[i]:
            ax[px, py].set_ylabel("C3/$\gamma^2$", fontsize=1.3 * fontsize)

    ax[0, 0].legend()

    fig.tight_layout()

    fname = "ensemble_fm90_vs_fm90_params"
    if args.png:
        fig.savefig(f"{fname}.png")
    elif args.pdf:
        fig.savefig(f"{fname}.pdf")
    else:
        plt.show()
