"""
Module Name: analyse_output.py

Description:
This module analyses the data from the output.nc
netCDF file.
"""

import os
import getpass
import matplotlib

import matplotlib.pyplot as plt
import numpy as np
import fourier_wave_number
import normal_modes_compare

matplotlib.use("Agg")

username = getpass.getuser()

""" do the fourier analysis
"""

u_jet = [50.0, 100.0, 150]
c_vis = [0.0, 0.1, 0.2, 1.0]

spy = int(np.floor(np.sqrt(len(c_vis))))
spx = int(np.ceil(len(c_vis) / spy))

plt.ion()
fig = plt.figure(figsize=(6 * spx, 5 * spy))
for j, value in enumerate(c_vis):
    fileNames = [
        "../../tests/output_" + str(i) + "_" + str(j) + ".nc" for i in range(len(u_jet))
    ]

    ax = plt.subplot(spy, spx, j + 1)

    fourier_wave_number.do_analysis01(fileNames, u_jet)
    # same for each value of c_cvis
    normal_modes_compare.normal_modes_compare(u_jet, 1)
    plt.title("$C_{vis}$=" + str(value))
    plt.text(0.15, 0.9, "(" + chr(97 + j) + ")", transform=ax.transAxes)

if not os.path.exists("../../output"):
    os.mkdir("../../output")

plt.savefig("../../output/full_analysis.png")
