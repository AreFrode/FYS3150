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
    kin = []
    pot = []
    tot = []
    t = []
    read: bool = True
    with open(file, "r") as infile:
        lines = infile.readlines()
        for line in lines:
            if read:
                line: str = line.split()
                t.append(float(line[0]))
                kin.append(float(line[-3]))
                pot.append(float(line[-2]))
                tot.append(float(line[-1]))
                read = False
            else:
                read = True


plt.plot(t, kin, label="Kinetic Energy Earth")
plt.plot(t, pot, label="Potential Energy Earth")
plt.plot(t, tot, label="Total Energy Earth")
plt.xlabel("t [Yrs]")
plt.ylabel("E [J]")
plt.legend()
plt.tight_layout()
plt.savefig(outfile)
