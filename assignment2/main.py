#!/usr/bin/env python3

import os
import sys

N = input("Number of mesh points: ")
algo = input("Choose a program - jacobi/arma/analytical/jacobi-arma/quantum: ")

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
os.system("mv" + " " + filename_data + " " + data_path)


os.system("echo generating plots...")
os.system("python3 plotter.py" + " " + filename_plot +
          " " + "/".join([data_path, filename_data]))

if not os.path.exists(plot_path):
    os.makedirs(plot_path)

os.system(" ".join(["mv", filename_plot, plot_path]))
os.system("echo done.")
