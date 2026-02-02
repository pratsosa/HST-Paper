"""
Shared plotting utilities for figure generation.

This module contains common plotting functions used across multiple figures.
"""

import numpy as np
from scipy import stats
from scipy.ndimage import gaussian_filter
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt


def plot_contour(xdata, ydata, c="k", mark=".", nlevels=3, ax=None,
                 linewidths=0.5, s=3, alpha=1, label=None):
    """
    Plot density contours with scatter points for low-density regions.
    
    Uses Gaussian KDE to estimate density and plots contours at specified levels.
    Points in low-density regions are highlighted as scatter points. This is an
    optimized version using Gaussian filtering and interpolation for performance.
    
    Parameters:
    -----------
    xdata, ydata : array-like
        Data coordinates
    c : str or color
        Color for contours and points
    mark : str
        Marker style for scatter points
    nlevels : int or list
        Number or list of density levels for contours
    ax : matplotlib axis
        Axis to plot on (default: current axis)
    linewidths : float
        Width of contour lines
    s : float
        Size of scatter points
    alpha : float
        Transparency of scatter points
    label : str
        Label for legend
    """
    if ax is None:
        ax = plt.gca()

    dx = 0.1 * (xdata.max() - xdata.min())
    dy = 0.1 * (ydata.max() - ydata.min())

    xmin, xmax = xdata.min() - dx, xdata.max() + dx
    ymin, ymax = ydata.min() - dy, ydata.max() + dy

    X, Y = np.mgrid[xmin:xmax:60j, ymin:ymax:60j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([xdata, ydata])

    kernel = stats.gaussian_kde(values, bw_method=0.2)
    Z = kernel(positions).reshape(X.shape)

    Z = gaussian_filter(Z, sigma=1.0)  # enforce smooth topology
    Z /= Z.max()                       # normalize safely

    cs = ax.contour(X, Y, Z, levels=nlevels,
                    linewidths=linewidths, colors=[c])
    levels = cs.levels

    interp = RegularGridInterpolator(
        (X[:, 0], Y[0, :]), Z,
        bounds_error=False, fill_value=Z.min()
    )

    z = interp(np.column_stack((xdata, ydata)))
    mask = z > levels[0]

    ax.scatter(xdata[~mask], ydata[~mask],
               s=s, color=c, marker=mark, alpha=alpha, label=label)
