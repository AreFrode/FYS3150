#!/usr/bin/env python3

from typing import List
import numpy as np
import matplotlib.pyplot as plt
import sys


def plotter(infile: str, exponent: int, read_u: bool = False) -> None:
    """
    Creates a plot of the output from the c++ function write_file.cpp

    Args:
        infile: /path/to/file
        exponent: the number of gridpoints the equation is solved for
        read_u (optional): plots the analytical solution over the given gridpoints

    Returns:
        As is stanard with matplotlib functions, nothing gets returned
    """

    n: int = 10**exponent
    counter: int = 1
    x = np.empty(n + 2)
    v = np.empty(n + 2)
    u = np.empty(n + 2)

    x[-1] = 1.
    v[-1] = 0.

    lines: List[str] = infile.readlines()
    for line in lines:
        line = line.split()
        x[counter] = line[0]
        v[counter] = line[1]
        if read_u:
            u[counter] = line[2]

        counter += 1

    plt.plot(x, v, label=f"numerical, n = {n}")
    if read_u:
        plt.plot(x, u, label=f"analytical, n = {n}")

if __name__ == "__main__":
    try:
        plotter(open(sys.argv[1], "r"), int(sys.argv[4]))
        plotter(open(sys.argv[2], "r"), int(sys.argv[4]) + 1)
        plotter(open(sys.argv[3], "r"), int(sys.argv[4]) + 2, read_u=True)

        plt.legend()
        plt.xlabel("x")
        plt.ylabel("v(x)")
        plt.savefig("../comparison.png")

    except IndexError as e:
        print("Wrong usage of program, no file was given as command line argument\n\
            usage: ./plotter infile1.ext infile2.ext infile3.ext exponent")
