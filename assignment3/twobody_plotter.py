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
names = ["Earth", "Jupiter"]

x = [[] for _ in range(2)]
y = [[] for _ in range(2)]

with open(file, "r") as infile:
    lines = infile.readlines()
    for line in lines:
        line: str = line.split()
        idx = int(line[1]) - 1
        x[idx].append(float(line[3]))
        y[idx].append(float(line[4]))

for i in range(2):
    plt.plot(x[i], y[i], '--', label=f"{names[i]}")

plt.plot(0.0, 0.0, 'ro', label="Sun")
plt.xlabel("x [AU]")
plt.ylabel("y [AU]")
plt.legend()
plt.axis('equal')
plt.savefig(outfile)
