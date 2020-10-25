#!/usr/bin/env python3

import matplotlib
import matplotlib.pyplot as plt
import sys

matplotlib.rcParams['text.usetex'] = True

if (len(sys.argv) < 3):
    print("Not enough kwargs given")
    sys.exit(1)

outfile: str = sys.argv[1]
file: str = sys.argv[2]

x = [[] for _ in range(3)]
y = [[] for _ in range(3)]

with open(file, "r") as infile:
    lines = infile.readlines()
    for line in lines:
        line: str = line.split()
        idx = int(line[1]) - 1
        x[idx].append(float(line[3]))
        y[idx].append(float(line[4]))

for i in range(3):
    plt.plot(x[i], y[i], '--', label=f"Planet {i}")

plt.xlabel("x [AU]")
plt.ylabel("y [AU]")
plt.legend()
plt.axis('equal')
plt.savefig(outfile)
