import numpy as np
from astropy.table import QTable
import astropy.units as u

# read in the Valencic et al. 2004 data and write the common format table
if __name__ == "__main__":
    # get the Gordon03 results
    gor03_smc = QTable.read(
        "Gordon03/gordon03_smc.dat", format="ascii.commented_header", header_start=0
    )
    gor03_lmc = QTable.read(
        "Gordon03/gordon03_lmc.dat", format="ascii.commented_header", header_start=0
    )

    # SMC
    gor03 = gor03_smc
    otab = QTable()
    ebv = np.array(gor03["ebv"]) * u.mag
    otab["EBV"] = ebv
    rv = np.array(gor03["rv"])
    otab["RV"] = rv
    otab["AV"] = ebv * rv

    otab["C1"] = np.array(gor03["c1"])
    otab["C2"] = np.array(gor03["c2"])
    otab["C3"] = np.array(gor03["c3"])
    otab["C4"] = np.array(gor03["c4"])
    otab["x0"] = np.array(gor03["x0"]) / u.micron
    otab["gamma"] = np.array(gor03["gamma"]) / u.micron

    otab.write("gor03_smc_ensemble_params.dat", format="ascii.ipac", overwrite=True)

    # LMC
    gor03 = gor03_lmc
    otab = QTable()
    ebv = np.array(gor03["ebv"]) * u.mag
    otab["EBV"] = ebv
    rv = np.array(gor03["rv"])
    otab["RV"] = rv
    otab["AV"] = ebv * rv

    otab["C1"] = np.array(gor03["c1"])
    otab["C2"] = np.array(gor03["c2"])
    otab["C3"] = np.array(gor03["c3"])
    otab["C4"] = np.array(gor03["c4"])
    otab["x0"] = np.array(gor03["x0"]) / u.micron
    otab["gamma"] = np.array(gor03["gamma"]) / u.micron

    otab.write("gor03_lmc_ensemble_params.dat", format="ascii.ipac", overwrite=True)
