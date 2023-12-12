import numpy as np

from astropy.table import QTable
import astropy.units as u

if __name__ == "__main__":
    #table = ascii.read("table3.dat", readme="ReadMe.txt")

    tab3 = QTable.read("Fitzpatrick07/table3.dat",
                       readme="Fitzpatrick07/ReadMe.txt",
                       format="ascii.cds")
    tab4 = QTable.read("Fitzpatrick07/table4.dat",
                       readme="Fitzpatrick07/ReadMe.txt",
                       format="ascii.cds")

    # build the table in the common format
    otab = QTable()
    otab["Name"] = np.array(tab3["Name"])
    otab["EBV"] = np.array(tab3["E(B-V)"]) * u.mag
    otab["EBV_unc"] = np.array(tab3["e_E(B-V)"]) * u.mag
    otab["RV"] = np.array(tab4["R(V)"])
    otab["RV_unc"] = np.array(tab4["e_R(V)"])
    
    otab["AV"] = otab["EBV"] * otab["RV"]
    otab["AV_unc"] = otab["AV"] * np.sqrt(otab["EBV_unc"].value**2 + otab["RV_unc"]**2)

    otab["C1"] = np.array(tab4["c1"])
    otab["C1_unc"] = np.array(tab4["e_c1"])
    otab["C2"] = np.array(tab4["c2"])
    otab["C2_unc"] = np.array(tab4["e_c2"])
    otab["C3"] = np.array(tab4["c3"])
    otab["C3_unc"] = np.array(tab4["e_c3"])
    otab["C4"] = np.array(tab4["c4"])
    otab["C4_unc"] = np.array(tab4["e_c4"])
    otab["x0"] = np.array(tab4["x0"]) / u.micron
    otab["x0_unc"] = np.array(tab4["e_x0"]) / u.micron
    otab["gamma"] = np.array(tab4["gamma"]) / u.micron
    otab["gamma_unc"] = np.array(tab4["e_gamma"]) / u.micron

    otab.write("fit07_ensemble_params.dat", format="ascii.ipac", overwrite=True)