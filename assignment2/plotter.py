#!/usr/bin/env python3

import matplotlib
import matplotlib.pyplot as plt
import sys

matplotlib.rcParams['text.usetex'] = True

outfile = sys.argv[1]
ifile = sys.argv[2]

rho = [0]
u = [0]

with open(ifile, "r") as infile:
    infile.readline()
    lines = infile.readlines()
    for line in lines:
        vals = line.split()
        rho.append(float(vals[0]))
        u.append(float(vals[1]))

rho.append(1)
u.append(0)

plt.plot(rho, u, label=ifile.split("/")[2])
plt.xlabel(r"$\rho$")
plt.ylabel(r"u($\rho$)")
plt.legend()
plt.savefig(outfile)
