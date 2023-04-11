import numpy as np
from astropy.table import QTable
import astropy.units as u

# read in the Valencic et al. 2004 data and write the common format table
if __name__ == "__main__":
    # get the Valencic04 results
    val04_dust = QTable.read("Valencic04/valencic04_datafile4.txt", format="ascii.cds")
    val04_fm = QTable.read("Valencic04/valencic04_datafile5.txt", format="ascii.cds")

    otab = QTable()
    otab["AV"] = np.array(val04_dust["A(V)"]) * u.mag
    otab["EBV"] = np.array(val04_dust["E(B-V)"]) * u.mag
    rv = np.array(val04_dust["R(V)"])
    otab["RV"] = rv

    otab["C1"] = (np.array(val04_fm["c1/R(V)+1.0"]) - 1.0) * rv
    otab["C2"] = np.array(val04_fm["c2/R(V)"]) * rv
    otab["C3"] = np.array(val04_fm["c3/R(V)"]) * rv
    otab["C4"] = np.array(val04_fm["c4/R(V)"]) * rv
    otab["x0"] = np.array(val04_fm["x0"]) / u.micron
    otab["gamma"] = np.array(val04_fm["gamma"]) / u.micron

    otab.write("val04_ensemble_params.dat", format="ascii.ipac", overwrite=True)
