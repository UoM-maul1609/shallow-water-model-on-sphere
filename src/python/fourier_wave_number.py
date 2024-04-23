"""
Module Name: fourier_wave_number.py

Description:
This module contains code to analyse the fourier
wave number of the output.nc netCDF file.
"""

import os
import getpass
import matplotlib

import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as sci
import scipy.signal as scs

# pylint: disable=E0611
# (pylint can't find Dataset in the netCDF4 package for some reason)
from netCDF4 import Dataset as NetCDFFile

matplotlib.use("Agg")

username = getpass.getuser()

# TODO: FIX AND TEST THIS MODULE pylint: disable=W0511


def fourier_wave_number(file_name):  # pylint: disable=R0914 disable=R0915
    """
    uses findpeaks and fourier analysis to track the current wave number and
    the rotation speed (i.e. using phase information from the fft).
    """

    n_files = len(file_name)
    nc = NetCDFFile(file_name[0])
    np1 = len(nc["time"][:])
    nc.close()

    wave_number = np.zeros((n_files, np1))
    rotation_rate = np.zeros((n_files, np1 - 1))
    phase = np.zeros((n_files, np1))
    phaseold = np.zeros((n_files, np1))
    for k in range(n_files):
        nc = NetCDFFile(file_name[k])
        dt_sec = nc["time"][1] - nc["time"][0]
        _, c = np.shape(nc["vort"][0, :, :])
        fs = c
        # time_period = 1.0 / fs
        length = c
        # t = np.mgrid[0:L] * time_period
        f = fs * np.mgrid[0 : (length / 2) + 1] / length
        for l in range(np1):  # noqa: E741
            x = np.mean(nc["v"][l, 50 : 60 + 1, :], axis=0)
            y = np.fft.fft(x)
            p2 = np.abs(y / length)
            p1 = p2[0 : int(length / 2) + 1]
            p1[1:-1] = 2 * p1[1:-1]
            locs1 = scs.find_peaks(np.concatenate([x, x[0:3]]), prominence=1)
            locs = locs1[0]
            if len(locs) > 0:
                if locs[0] == locs[-1] - c:
                    locs = np.delete(locs, len(locs) - 1)

            ind = len(locs)
            phs = np.unwrap(np.angle(y)) * 180.0 / np.pi  # phase angle in degrees

            print("wave number is : " + str(f[ind]) + "; phase: " + str(phs[ind]))
            phase[k, l] = phs[ind]
            wave_number[k, l] = f[ind]

        nc.close()
        phaseold[k, :] = phase[k, :]
        # phase is the phase of the wave through it's cycle, but
        # we want the movement of the wave train around the planet,
        # so divide by number of waves.
        phase[k, :] = (
            np.unwrap(phaseold[k, :] * np.pi / 180.0)
            * 180.0
            / np.pi
            / (wave_number[k, :])
        )
        rotation_rate[k, :] = (
            -np.diff(phase[k, :], n=1, axis=0) / dt_sec * 86400.0 * 365.25 / 360
        )  # full rotations per year

    return (phase, rotation_rate, wave_number)


def do_analysis01(file_names, u_jets):
    """
    Does analysis

    Args:
        file_names: names of files.
        u_jets: beginning speed of the jets.
    """
    n_files = len(file_names)
    cmap_lev = 64
    # basemap = plt.get_cmap(lut=cmap_lev, name="ocean")
    plt.get_cmap(lut=cmap_lev, name="ocean")
    u1 = np.linspace(u_jets[0], u_jets[-1], cmap_lev)
    # int1 = sci.interp1d(u1, np.mgrid[1 : cmap_lev + 1], kind="nearest")
    sci.interp1d(u1, np.mgrid[1 : cmap_lev + 1], kind="nearest")
    mean_rot = np.zeros((n_files, 9))
    std_rot = np.zeros((n_files, 9))
    for k in range(n_files):
        (_, rotation_rate, wave_number) = fourier_wave_number([file_names[k]])
        for l in range(0, 9):  # noqa: E741
            print(k)
            (ind1,) = np.where(np.diff(wave_number[0, :]) == 0)
            (ind,) = np.where(wave_number[0, ind1] == (l + 1))
            ind = ind1[ind]
            if len(ind) > 2:
                mean_rot[k, l] = np.mean(rotation_rate[0, ind[1:-2]])
                std_rot[k, l] = np.std(rotation_rate[0, ind[1:-2]])
            else:
                mean_rot[k, l] = np.nan
                std_rot[k, l] = np.nan

        h = plt.scatter(np.mgrid[1:10], mean_rot[k, :], s=40, c=u_jets[k] * np.ones(9))
        plt.clim((u_jets[0], u_jets[-1]))
    plt.xlabel("wave number (per full rotation)")
    plt.ylabel("mean rotation rate (rotations per year)")
    h = plt.colorbar()
    h.set_label("Jet speed (m/s)")


if __name__ == "__main__":
    FILE_NAME = ["../../tests/output.nc"]
    N_FILES = len(FILE_NAME)
    U_JETS = [50.0, 100.0, 150.0]
    THIS_ONE = 2
    CMAP_LEV = 64

    (phase_main, rotation_rate_main, wave_number_main) = fourier_wave_number(FILE_NAME)

    basemap = plt.get_cmap(lut=CMAP_LEV, name="ocean")
    u1_main = np.linspace(U_JETS[0], U_JETS[-1], CMAP_LEV)
    int1_main = sci.interp1d(u1_main, np.mgrid[1 : CMAP_LEV + 1], kind="nearest")
    mean_rot_main = np.zeros((N_FILES, 9))
    std_rot_main = np.zeros((N_FILES, 9))
    for j in range(N_FILES):
        for i in range(0, 9):
            (ind1_main,) = np.where(np.diff(wave_number_main[j, :]) == 0)
            (ind_main,) = np.where(wave_number_main[j, ind1_main] == (i + 1))
            ind_main = ind1_main[ind_main]
            if len(ind_main) > 2:
                mean_rot_main[j, i] = np.mean(rotation_rate_main[j, ind_main[1:-2]])
                #                 print(rotation_rate[j,ind[1:-2]])
                std_rot_main[j, i] = np.std(rotation_rate_main[j, ind_main[1:-2]])
            else:
                mean_rot_main[j, i] = np.nan
                std_rot_main[j, i] = np.nan

    h_plot = plt.scatter(
        np.mgrid[1:10], mean_rot_main[j, :], s=40, c=U_JETS[THIS_ONE] * np.ones(9)
    )
    plt.xlabel("wave number (per full rotation)")
    plt.ylabel("mean rotation rate (rotations per year)")
    plt.clim((np.min(U_JETS), np.max(U_JETS)))
    if THIS_ONE == 0:
        h_plot = plt.colorbar()
        h_plot.set_label("Jet speed (m/s)")

    if not os.path.exists("../../output/images"):
        os.mkdir("../../output/images")

    plt.savefig("../../output/fourier_wave_number.png", format="png", dpi=300)
