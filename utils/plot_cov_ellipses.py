import numpy as np
from math import sqrt, cos, sin
from matplotlib.patches import Polygon
from scipy.linalg import eigh

# from Dries' dust_fuse_h2 repository
def cov_ellipse(x, y, cov, num_sigma=1, **kwargs):
    """
    Create an ellipse at the coordinates (x,y), that represents the
    covariance. The style of the ellipse can be adjusted using the
    kwargs.

    Returns
    -------
    ellipse: matplotlib.patches.Ellipse
    """

    position = [x, y]

    if cov[0, 1] != 0:
        # length^2 and orientation of ellipse axes is determined by
        # eigenvalues and vectors, respectively. Eigh is more stable for
        # symmetric / hermitian matrices.
        values, vectors = eigh(cov)
        width, height = np.sqrt(np.abs(values)) * num_sigma * 2
    else:
        width = sqrt(cov[0, 0]) * 2
        height = sqrt(cov[1, 1]) * 2
        vectors = np.array([[1, 0], [0, 1]])

    # I ended up using a Polygon just like Karl's plotting code. The
    # ellipse is buggy when the difference in axes is extreme (1e22). I
    # think it is because even a slight rotation will make the ellipse
    # look extremely strechted, as the extremely long axis (~1e22)
    # rotates into the short coordinate (~1).

    # two vectors representing the axes of the ellipse
    vw = vectors[:, 0] * width / 2
    vh = vectors[:, 1] * height / 2

    # generate corners
    num_corners = 64
    angles = np.linspace(0, 2 * np.pi, num_corners, endpoint=False)
    corners = np.row_stack([position + vw * cos(a) + vh * sin(a) for a in angles])

    return Polygon(corners, **kwargs)


def draw_ellipses(ax, xs, ys, covs, num_sigma=1, sigmas=None, **kwargs):
    for k, (x, y, cov) in enumerate(zip(xs, ys, covs)):
        # if sigmas is not None:
        #     color = cm.viridis(sigmas[k] / 3.0)[0]
        # ax.add_patch(cov_ellipse(x, y, cov, num_sigma, color=color, **kwargs))
        ax.add_patch(cov_ellipse(x, y, cov, num_sigma, **kwargs))