import numpy as np
from astropy.table import QTable

files = ["gor09", "gor03_lmc", "gor24_smc_nolowebv"]
# files = ["gor09"]

c2 = []
b3 = []
c4 = []
rv = []
nhiebv = []
for cfile in files:

    itab = QTable.read(f"data/{cfile}_ensemble_params.dat", format="ascii.ipac")
    npts = len(itab)

    if "B3" not in itab.colnames:
        itab["B3"] = itab["C3"] / (itab["gamma"] ** 2)
        if "C3_unc" in itab.colnames:
            itab["B3_unc"] = np.absolute(itab["B3"]) * np.sqrt(itab["C3_unc"] ** 2 +  2.0 * (itab["gamma_unc"].value ** 2))

    c2 = np.concatenate((c2, itab["C2"].data))
    b3 = np.concatenate((b3, itab["B3"].data))
    c4 = np.concatenate((c4, itab["C4"].data))
    rv = np.concatenate((rv, itab["RV"].data))
    nhiebv = np.concatenate((nhiebv, itab["NHI"].data / itab["EBV"].data))


npts = len(c2)
# X is your data table, where the features (bump strength, C2, B3, etc) are columns and 
# individual sightlines are rows
X = np.zeros((npts, 5))
#X[:, 0] = itab["C1"]
X[:, 0] = c2
X[:, 1] = b3
X[:, 2] = c4
X[:, 3] = rv
X[:, 4] = nhiebv

np.save("pca_input.dat", X)

N_sightlines, N_features = X.shape

# normalize your data, i.e., give all features zero mean and unit variance
X_norm = (X - X.mean(axis=0)) / X.std(axis=0)

# compute covariance matrix between features; note that "@" is matrix multiply
# cov_X.shape should equal (N_features, N_features)
cov_X = np.cov(X_norm.transpose())

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
