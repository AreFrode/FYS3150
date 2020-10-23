#!/usr/bin/env python3

import os
import sys

N: str = input("Number of mesh points: ")
t: str = input("Final calculation time: ")
prog: str = input("Choose a program \n* earth-fe\n* earth-vv\n* ellipse\n: ")

filename_plot: str = "_".join([prog, "solution", N, t]) + ".pdf"
filename_plot_en: str = "_".join([prog, "solution", N, t]) + "_energy.pdf"
filename_plot_mom: str = "_".join([prog, "solution", N, t]) + "_moment.pdf"
filename_data_pos: str = "_".join([prog, "N", N, "t", t]) + ".txt"
filename_data_en: str = "_".join([prog, "N", N, "t", t]) + "_energy.txt"
filename_data_mom: str = "_".join([prog, "N", N, "t", t]) + "_moment.txt"

plot_path: str = "/".join([".", "plots", prog])
data_path: str = "/".join([".", "results", prog])

cpp_scripts: str = "./cpp_scripts/*.cpp"
compiler_flags: str = "-O3"

os.system("echo compiling...")
os.system("g++ -o main.out" + " " + cpp_scripts + " " + compiler_flags)

os.system("echo executing...")
os.system("./main.out " + N + " " + t + " " + prog)

if not os.path.exists(data_path):
    os.makedirs(data_path)

os.system("mv" + " " + filename_data_pos + " " + data_path)
os.system("mv" + " " + filename_data_en + " " + data_path)
os.system("mv" + " " + filename_data_mom + " " + data_path)

if not os.path.exists(plot_path):
    os.makedirs(plot_path)

os.system("echo generating plots...")
os.system("python3.8 plotter.py" + " " + filename_plot +
          " " + "/".join([data_path, filename_data_pos]))
os.system("python3.8 plot_energy.py" + " " + filename_plot_en +
          " " + "/".join([data_path, filename_data_en]))

os.system(" ".join(["mv", filename_plot, plot_path]))
os.system(" ".join(["mv", filename_plot_en, plot_path]))

os.system("echo done.")
