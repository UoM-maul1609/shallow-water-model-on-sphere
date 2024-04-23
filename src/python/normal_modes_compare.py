"""
Module Name: normal_modes_compare.py

Description:
Compares the normal modes for a range of u_jets.
"""

import sys
import os
import getpass
import matplotlib

import matplotlib.pyplot as plt
import numpy as np
import scipy as scy
import scipy.interpolate as sci

matplotlib.use("Agg")

username = getpass.getuser()

# TODO: FIX AND TEST THIS MODULE pylint: disable=W0511

def normal_modes_compare(u_jets, flag): # pylint: disable=R0914 disable=R0915
    """
    normal mode analysis of a gaussian jet - see Mak, Atmospheric Dynamics, page 229
    """
    cmap_lev = 64
    basemap = plt.get_cmap(lut=cmap_lev, name="viridis")
    u1 = np.linspace(u_jets[0], u_jets[-1], cmap_lev)
    int1 = sci.interp1d(u1, np.mgrid[1 : cmap_lev + 1], kind="nearest")
    # L = 1000e3  # measure distance in units of 1000 km
    # U = 10  # measure speed in units of 10 m/s
    # L_on_U = L / U
    # UL = U * L
    time_period = 10.55 * 3600  # time-period of rotation

    # set up domain.
    ip = 200  # finite difference grid
    jp = 100
    lat_high = 85  # highest latitude
    lat_low = 65  # lowest
    re = 5.4155760e7  # radius of saturn in this region (due to squashed spheriod)
    h_jet = 1.0  # standard deviation of jet
    lat_jet = 77  # latitude of the jet
    n_ks = 36  # calculate the growth factor for this many k-values

    for i, u_jet in enumerate(u_jets):
        # y-grid
        y = np.linspace(re * lat_low * np.pi / 180.0, re * lat_high * np.pi / 180, jp)
        # distance around planet at this latitude
        x_len = 2.0 * np.pi * np.cos((lat_jet) * np.pi / 180.0) * re
        # x_len2 = 2.0 * np.pi * np.cos((lat_jet) * np.pi / 180.0) * re * 4 / 7
        # x-grid
        x = np.linspace(0, x_len, ip)
        np.meshgrid(x, y)
        dy = y[1] - y[0]
        # the jet:
        a = np.sqrt(2) * h_jet * np.pi / 180.0 * re
        b = lat_jet * np.pi / 180.0 * re
        u = u_jet * np.exp(-(((y - b) / a) ** 2))
        # derivative of vorticity wrt y (i.e. d/dy(-du/dy):
        zeta_y = 2.0 * u / a**2.0 * (1.0 - 2.0 * (y - b) ** 2.0 / a**2)
        # beta (df/dy)
        beta1 = 2.0 * 2.0 * np.pi / time_period * np.cos(y / re) / re
        sigma_max = np.zeros(n_ks)
        sigmas_max = np.zeros((n_ks, 3))
        sigmas_max_i = np.zeros((n_ks, 3))
        sigmas_min = np.zeros((n_ks, 3))
        sigmas_min_i = np.zeros((n_ks, 3))
        # now set-up matrix problem
        for n in range(1, n_ks + 1):
            k = 2 * np.pi * n / x_len  # *(x_len/x_len2);

            # A matrix:
            # top diag elements + diagonal elements + bottom diag elements
            a_matrix = (
                np.diag(1j * u[0:-1] * k / dy**2, 1)
                + np.diag(1j * k * ((zeta_y + beta1) - 2.0 * u / dy**2 - u * k**2))
                + np.diag(1j * u[1:] * k / dy**2, -1)
            )

            # B matrix:
            # top diag elements + diagonal elements + bottom diag elements
            b_matrix = (
                np.diag(-1.0 / dy**2 * np.ones(len(u) - 1), 1)
                + np.diag((k**2 + 2.0 / dy**2) * np.ones(len(u)))
                + np.diag(-1.0 / dy**2 * np.ones(len(u) - 1), -1)
            )
            # find eigenvalues:
            eigenvalues = scy.linalg.eig(a_matrix, b_matrix, left=True)
            # E,D = scy.linalg.eig(A,B)
            sigma1 = eigenvalues[0]
            # sigma1=np.diag(D)

            sigma_max[n - 1] = np.max(np.real(sigma1))
            a = np.flip(np.sort(np.real(sigma1)))
            ii = np.flip(np.argsort(np.real(sigma1)))
            sigmas_max[n - 1, 0:3] = np.real(sigma1[ii[0:3]])
            sigmas_max_i[n - 1, 0:3] = np.imag(sigma1[ii[0:3]])
            a = np.sort(np.real(sigma1))
            ii = np.argsort(np.real(sigma1))
            sigmas_min[n - 1, 0:3] = np.real(sigma1[ii[0:3]])
            sigmas_min_i[n - 1, 0:3] = np.imag(sigma1[ii[0:3]])
            # if growth rate small, set to NaN
            if sigmas_max[n - 1, 0] < 1.0e-9:
                sigmas_max_i[n - 1, 0] = np.nan

            if sigmas_max[n - 1, 1] < 1.0e-9:
                sigmas_max_i[n - 1, 1] = np.nan

            if sigmas_max[n - 1, 2] < 1.0e-9:
                sigmas_max_i[n - 1, 2] = np.nan

        h = []
        if flag == 1:
            # sigma / k is wave speed in m/s
            # k is 2*pi/lambda and lambda is x_len/n
            # so this is rotations per year
            h.append(
                plt.plot(
                    np.mgrid[1 : n_ks + 1],
                    (
                        -sigmas_max_i[:, 0]
                        / (2.0 * np.pi * np.mgrid[1 : n_ks + 1] / x_len)
                    )
                    * 86400.0
                    * 365.25
                    / x_len,
                    "-x",
                )
            )
            h.append(
                plt.plot(
                    np.mgrid[1 : n_ks + 1],
                    (
                        -sigmas_max_i[:, 1]
                        / (2.0 * np.pi * np.mgrid[1 : n_ks + 1] / x_len)
                    )
                    * 86400.0
                    * 365.25
                    / x_len,
                    "--x",
                )
            )
            h.append(
                plt.plot(
                    np.mgrid[1 : n_ks + 1],
                    (
                        -sigmas_max_i[:, 2]
                        / (2.0 * np.pi * np.mgrid[1 : n_ks + 1] / x_len)
                    )
                    * 86400.0
                    * 365.25
                    / x_len,
                    ":x",
                )
            )
            plt.ylabel("speed of wave (rotations per year)")
            plt.xlabel("number of peaks")
        elif flag == 2:
            # largest growth rate
            h.append(plt.plot(np.mgrid[1 : n_ks + 1], sigmas_max[:, 0]))
            # second largest
            h.append(plt.plot(np.mgrid[1 : n_ks + 1], sigmas_max[:, 1], "--"))
            # third largest
            h.append(plt.plot(np.mgrid[1 : n_ks + 1], sigmas_max[:, 2], ":"))
            # sum of all
            h.append(
                plt.plot(np.mgrid[1 : n_ks + 1], np.sum(sigmas_max[:, :], axis=1), lw=3)
            )
            plt.grid("on")
            plt.xlabel("number of peaks")
            plt.ylabel("Growth rate")
            plt.legend(["Largest", "2nd largest", "3rd largest", "sum"])
        elif flag == 3:
            h.append(
                plt.plot(
                    np.mgrid[1 : n_ks + 1],
                    (
                        -sigmas_max_i[:, 0]
                        / (2.0 * np.pi * np.mgrid[1 : n_ks + 1] / x_len)
                    ),
                )
            )
            h.append(
                plt.plot(
                    np.mgrid[1 : n_ks + 1],
                    (
                        -sigmas_max_i[:, 1]
                        / (2.0 * np.pi * np.mgrid[1 : n_ks + 1] / x_len)
                    ),
                    "--",
                )
            )
            plt.ylabel("speed of wave (rotations per year)")
            plt.xlabel("number of peaks")

        if len(u_jets) > 1:
            row = int(int1(u_jets[i]))
        else:
            row = 1

        for k, h_val in enumerate(h):
            h_val[0].set_color(basemap(row))

        if flag == 2:
            plt.legend(["Largest", "2nd largest", "3rd largest", "sum"])


if __name__ == "__main__":
    plt.ion()
    fig = plt.figure()
    U_JETS = [50.0, 100.0, 150.0]
    # u_jets=[100.]

    if len(sys.argv) > 1:
        FLAG = int(sys.argv[1])
    else:
        FLAG = 1

    # 1, plot speed vs wave number
    # 2, plot growth rate vs wave number
    # 3, same as 1, but only first two
    normal_modes_compare(U_JETS, FLAG)

    if FLAG in (2,3):
        sm = plt.cm.ScalarMappable(
            cmap="viridis", norm=plt.Normalize(vmin=np.min(U_JETS), vmax=np.max(U_JETS))
        )
        h_fig = fig.colorbar(sm, ax=fig.gca())
        h_fig.set_label("$U_{jet}$ (m s$^{-1}$)")

    if not os.path.exists("../../output"):
        os.mkdir("../../output")

    plt.savefig("../../output/fourier_wave_number.png", format="png", dpi=300)
