import numpy as np
from astropy.table import QTable
import astropy.units as u

# read in the Gordon et al. 2009 and Van De Putte 2024 
# both needed for full file for MW FUSE sample
if __name__ == "__main__":
    tab1 = QTable.read("Gordon09/apj308925t2_ascii.txt", format="ascii.commented_header")
    tab2 = QTable.read("Gordon09/apj308925t3_ascii.txt", format="ascii.commented_header")
    # tab1 and tab2 verified to be the same length and have the same order

    tab3 = QTable.read("VanDePutte23/apjac9902t2_mrt.txt", format="cds")

    otab = QTable()

    # column info
    otab["Name"] = tab1["Name"]
    colnames = ["AV", "EBV", "RV"]
    tabnames = ["A(V)", "E(B-V)", "R(V)"]
    for ccol, ctab in zip(colnames, tabnames):
        otab[ccol] = np.array(tab1[ctab].data)
        otab[f"{ccol}_unc"] = np.sqrt(tab1[f"{ctab}_unc1"].data**2
                                      + tab1[f"{ctab}_unc2"].data**2)
    otab["AV"] = otab["AV"] * u.mag
    otab["AV_unc"] = otab["AV_unc"] * u.mag
    otab["EBV"] = otab["EBV"] * u.mag
    otab["EBV_unc"] = otab["EBV_unc"] * u.mag

    # FM90 parameters
    rv = otab["RV"].data
    colnames = ["C1", "C2", "C3", "C4", "x0", "gamma"]
    tabnames = ["c1", "c2", "c3", "c4", "xo", "gamma"]
    for ccol, ctab in zip(colnames, tabnames):
        otab[ccol] = np.array(tab2[ctab].data)
        otab[f"{ccol}_unc"] = np.sqrt(tab2[f"{ctab}_unc1"].data**2
                                      + tab2[f"{ctab}_unc2"].data**2)
        # convert from FM90 values in A(lambda)/A(V) to E(lambda-V)/E(B-V)
        if ctab[0:1] == "c":
            if ctab == "c1":
                otab[ccol] -= 1.0
            otab[ccol] *= rv
            otab[f"{ccol}_unc"] *= rv

    # go through van23 and get the HI and H2 columns
    npts = len(otab["Name"].data)
    nhi = np.full(npts, 0.0)
    nhi_unc = np.full(npts, 0.0)
    nh2 = np.full(npts, 0.0)
    nh2_unc = np.full(npts, 0.0)
    nh = np.full(npts, 0.0)
    nh_unc = np.full(npts, 0.0)
    for k, cname in enumerate(otab["Name"].data):
        mindx, = np.where(cname == tab3["Star"].data)
        nhi[k] = 10**tab3["logNHI"].data[mindx[0]]
        nhi_unc[k] = 0.5 * (10**(tab3["logNHI"].data[mindx[0]] + tab3["e_logNHI"].data[mindx[0]])
                            - 10**(tab3["logNHI"].data[mindx[0]] - tab3["e_logNHI"].data[mindx[0]]))
        nh2[k] = 10**tab3["logNH2"].data[mindx[0]]
        nh2_unc[k] = 0.5 * (10**(tab3["logNH2"].data[mindx[0]] + tab3["e_logNH2"].data[mindx[0]])
                            - 10**(tab3["logNH2"].data[mindx[0]] - tab3["e_logNH2"].data[mindx[0]]))
        nh[k] = 10**tab3["logNHI"].data[mindx[0]]
        nh_unc[k] = 0.5 * (10**(tab3["logNH"].data[mindx[0]] + tab3["e_logNH"].data[mindx[0]])
                            - 10**(tab3["logNH"].data[mindx[0]] - tab3["e_logNH"].data[mindx[0]]))
    otab["NHI"] = nhi
    otab["NHI_unc"] = nhi_unc
    otab["NH2"] = nh2
    otab["NH2_unc"] = nh2_unc
    otab["NH"] = nh
    otab["NH_unc"] = nh_unc

    otab.write("gor09_ensemble_params.dat", format="ascii.ipac", overwrite=True)