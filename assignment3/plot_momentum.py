#!/usr/bin/env python3

import matplotlib
import matplotlib.pyplot as plt
import sys

matplotlib.rcParams['text.usetex'] = True

if (len(sys.argv) < 3):
    print("Not enough kwargs given")
    sys.exit(1)

outfile: str = sys.argv[1]
ifile: str = sys.argv[2:]

for file in ifile:
    momx = []
    momy = []
    momz = []
    t = []
    with open(file, "r") as infile:
        lines = infile.readlines()
        for line in lines:
            line: str = line.split()
            if int(line[1]) == 0:
                t.append(float(line[0]))
                momx.append(float(line[-3]))
                momy.append(float(line[-2]))
                momz.append(float(line[-1]))

plt.plot(t, momx, color='m', label="Angular momentum x-direction Earth")
plt.plot(t, momy, color='c', label="Angular momentum y-direction Earth")
plt.plot(t, momz, color='r', label="Angular momentum z-direction Earth")
plt.xlabel("t [Yrs]")
plt.ylabel(r"L [$\frac{kgm^2}{s}$]")
plt.legend()
plt.tight_layout()
plt.savefig(outfile)
