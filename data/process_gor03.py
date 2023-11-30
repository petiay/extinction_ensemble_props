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

    for gor03, oname in zip([gor03_smc, gor03_lmc], ["smc", "lmc"]):
        otab = QTable()
        otab["EBV"] = gor03["ebv"] * u.mag
        otab["EBV_unc"] = gor03["ebv_u"] * u.mag
        otab["RV"] = gor03["rv"]
        otab["RV_unc"] = gor03["rv_u"]
        otab["AV"] = otab["EBV"] * otab["RV"]
        otab["AV_unc"] = otab["AV"] * np.sqrt(otab["EBV_unc"].value**2 + otab["RV_unc"]**2)

        otab["C1"] = np.array(gor03["c1"])
        otab["C1_unc"] = np.array(gor03["c1_u"])
        otab["C2"] = np.array(gor03["c2"])
        otab["C2_unc"] = np.array(gor03["c2_u"])
        otab["C3"] = np.array(gor03["c3"])
        otab["C3_unc"] = np.array(gor03["c3_u"])
        otab["C4"] = np.array(gor03["c4"])
        otab["C4_unc"] = np.array(gor03["c4_u"])
        otab["x0"] = np.array(gor03["x0"]) / u.micron
        otab["x0_unc"] = np.array(gor03["x0_u"]) / u.micron
        otab["gamma"] = np.array(gor03["gamma"]) / u.micron
        otab["gamma_unc"] = np.array(gor03["gamma_u"]) / u.micron

        otab.write(f"gor03_{oname}_ensemble_params.dat", format="ascii.ipac", overwrite=True)
