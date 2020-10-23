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
    x = []
    y = []
    with open(file, "r") as infile:
        lines = infile.readlines()
        for line in lines:
            line: str = line.split()
            x.append(float(line[3]))
            y.append(float(line[4]))


plt.plot(x, y, '--', label="Earth")
plt.plot(0.0, 0.0, 'ro', label="Sun")
plt.xlabel("x [AU]")
plt.ylabel("y [AU]")
plt.legend()
plt.axis('equal')
plt.savefig(outfile)
