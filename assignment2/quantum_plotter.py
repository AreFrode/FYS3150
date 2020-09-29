#!/usr/bin/env python3

import matplotlib
import matplotlib.pyplot as plt
import sys

matplotlib.rcParams['text.usetex'] = True

outfile = sys.argv[1]
ifile = sys.argv[2:]

for file in ifile:
    rho = [0]
    u = [0]
    with open(file, "r") as infile:
        infile.readline()
        lines = infile.readlines()
        for line in lines:
            vals = line.split()
            rho.append(float(vals[0]))
            u.append(float(vals[1]))

    rho.append(5.)
    u.append(0)

    name = file.split("_")
    plt.plot(rho, u, label=r"$\omega_r$" + ": " + name[-1][0:4])

plt.xlabel(r"$\rho$")
plt.ylabel(r"$\psi(\rho)$")
plt.legend()
plt.savefig(outfile)
