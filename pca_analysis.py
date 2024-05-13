import numpy as np
from astropy.table import QTable

files = ["gor09", "gor03_lmc", "gor24_smc_nolowebv"]
# files = ["gor09"]

c1 = []
c2 = []
b3 = []
gamma = []
c4 = []
rv = []
nhiebv = []
c1_unc = []
c2_unc = []
b3_unc = []
gamma_unc = []
c4_unc = []
rv_unc = []
nhiebv_unc = []
for cfile in files:

    itab = QTable.read(f"data/{cfile}_ensemble_params.dat", format="ascii.ipac")
    npts = len(itab)

    if "B3" not in itab.colnames:
        itab["B3"] = itab["C3"] / (itab["gamma"] ** 2)
        if "C3_unc" in itab.colnames:
            itab["B3_unc"] = np.absolute(itab["B3"]) * np.sqrt(itab["C3_unc"] ** 2 +  2.0 * (itab["gamma_unc"].value ** 2))

    c1 = np.concatenate((c1, itab["C1"].data))
    c1_unc = np.concatenate((c1_unc, itab["C1_unc"].data))
    c2 = np.concatenate((c2, itab["C2"].data))
    c2_unc = np.concatenate((c2_unc, itab["C2_unc"].data))
    b3 = np.concatenate((b3, itab["B3"].data))
    b3_unc = np.concatenate((b3_unc, itab["B3_unc"].data))
    gamma = np.concatenate((gamma, itab["gamma"].data))
    gamma_unc = np.concatenate((gamma_unc, itab["gamma_unc"].data))
    c4 = np.concatenate((c4, itab["C4"].data))
    c4_unc = np.concatenate((c4_unc, itab["C4_unc"].data))
    rv = np.concatenate((rv, itab["RV"].data))
    rv_unc = np.concatenate((rv_unc, itab["RV_unc"].data))
    tval = itab["NHI"].data / itab["EBV"].data
    nhiebv = np.concatenate((nhiebv, tval))
    tunc = tval * np.sqrt((itab["NHI_unc"].data / itab["NHI"].data)**2 
                          + (np.array(itab["EBV_unc"].data) / itab["EBV"].data)**2)
    nhiebv_unc = np.concatenate((nhiebv_unc, itab["NHI"].data / itab["EBV"].data))

# for use in this program
npts = len(c2)
# X is your data table, where the features (bump strength, C2, B3, etc) are columns and 
# individual sightlines are rows
X = np.zeros((npts, 14))
X[:, 0] = c1
X[:, 1] = c2
X[:, 2] = b3
X[:, 3] = gamma
X[:, 4] = c4
X[:, 5] = rv
X[:, 6] = nhiebv
X[:, 7] = c1_unc
X[:, 8] = c2_unc
X[:, 9] = b3_unc
X[:, 10] = gamma_unc
X[:, 11] = c4_unc
X[:, 12] = rv_unc
X[:, 13] = nhiebv_unc

np.save("pca_input.dat", X)

X = X[:, 1:7]

N_sightlines, N_features = X.shape

# normalize your data, i.e., give all features zero mean and unit variance
X_norm = (X - X.mean(axis=0)) / X.std(axis=0)

# compute covariance matrix between features; note that "@" is matrix multiply
# cov_X.shape should equal (N_features, N_features)
cov_X = np.cov(X_norm.transpose())

print(cov_X)

# decompose covariance matrix into eigenvalues/eigenvectors
# eigenvalues_X.shape should equal (N_features,)
# eigenvectors_X.shape should equal (N_features, N_features)
eigenvalues_X, eigenvectors_X = np.linalg.eig(cov_X)

# sort them in order from highest to lowest eigenvalues; now your highest eigenvalue  
# matches the eigenvector that explains most of the variance, aka your first principal component
# second highest eigenvalue is your second principal component, and etc
order = np.argsort(np.abs(eigenvalues_X))[::-1]

# reorder the eigenvalues and eigenvectors
eigenvalues_X = eigenvalues_X[order]
principal_components = eigenvectors_X[:, order]

# explained variance is simply the eigenvalue divided by sum of all eigenvalues
explained_variance = eigenvalues_X / eigenvalues_X.sum()

for k in range(len(explained_variance)):
    print(f"Component #{k+1} explains {np.round(100*explained_variance[k], 3)} % of the variance, and is in the direction", principal_components[:, k])
