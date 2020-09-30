#!/usr/bin/env python3

import os
import sys

N = input("Number of mesh points: ")
algo = input(
    "Choose a program \n* jacobi\n* arma\n* analytical\
    \n* jacobi-arma\n* one-electron\n* two-electrons\n: ")

filename_plot = "_".join([algo, "solution", N]) + ".pdf"
filename_data = "_".join([algo, "N", N]) + ".txt"
plot_path = "/".join([".", "plots", algo])
data_path = "/".join([".", "results", algo])

cpp_scripts = "./cpp_scripts/*.cpp"
compiler_flags = "-O3 -larmadillo"

os.system("echo compiling...")
os.system("g++ -o main.out" + " " + cpp_scripts + " " + compiler_flags)

os.system("echo executing...")
os.system("./main.out " + N + " " + algo)

if not os.path.exists(data_path):
    os.makedirs(data_path)

if algo == "one-electron":
    for i in range(4):
        os.system(
            "mv" + " " + filename_data[:-4] + "_lambda_" + str(i) + ".txt" +
            " " + data_path)

elif algo == "two-electrons":
    os.system("mv" + " " + filename_data[:-4] + "*" + " " + data_path)
else:
    os.system("mv" + " " + filename_data + " " + data_path)


if not os.path.exists(plot_path):
    os.makedirs(plot_path)

if algo == "one-electron":
    os.system("echo skipping plots...")

elif algo == "two-electrons":
    os.system("echo generating plots...")
    os.system("python3 quantum_plotter.py" + " " + filename_plot +
              " " + "/".join([data_path, filename_data[:-4] + "*"]))

    os.system(" ".join(["mv", filename_plot, plot_path]))
else:
    os.system("echo generating plots...")
    os.system("python3 plotter.py" + " " + filename_plot +
              " " + "/".join([data_path, filename_data]))

    os.system(" ".join(["mv", filename_plot, plot_path]))


os.system("echo done.")
