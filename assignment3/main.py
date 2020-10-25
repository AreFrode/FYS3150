#!/usr/bin/env python3

import os
import sys

N: str = input("Number of mesh points: ")
t: str = input("Final calculation time: ")
prog: str = input(
    "Choose a program \n* earth-fe\n* earth-vv\n* ellipse\n* beta-3.5-circ\n* beta-3.9-circ\n* beta-4.0-circ\n* beta-4.0-ellipse\n* escape\n* earth-jupiter\n* earth-10jupiter\n* earth-1000jupiter\n* three-body\n* full-system\n: ")

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

if not os.path.exists(plot_path):
  os.makedirs(plot_path)

os.system("echo generating plots...")

if (prog == "earth-jupiter" or prog == "earth-10jupiter" or prog == "earth-1000jupiter"):
  os.system("python3.8 twobody_plotter.py" + " " + filename_plot +
            " " + "/".join([data_path, filename_data_pos]))

elif (prog == "three-body"):
  os.system("python3.8 threebody_plotter.py" + " " + filename_plot +
            " " + "/".join([data_path, filename_data_pos]))

elif (prog == "full-system"):
  os.system("python3.8 manybody_plotter.py" + " " + filename_plot +
            " " + "/".join([data_path, filename_data_pos]))

else:
  os.system("python3.8 plotter.py" + " " + filename_plot +
            " " + "/".join([data_path, filename_data_pos]))

os.system("python3.8 plot_energy.py" + " " + filename_plot_en +
          " " + "/".join([data_path, filename_data_en]))

if not (prog == "earth-fe"):
  os.system("mv" + " " + filename_data_mom + " " + data_path)
  os.system("python3.8 plot_momentum.py" + " " + filename_plot_mom +
            " " + "/".join([data_path, filename_data_mom]))
  os.system(" ".join(["mv", filename_plot_mom, plot_path]))

os.system(" ".join(["mv", filename_plot, plot_path]))
os.system(" ".join(["mv", filename_plot_en, plot_path]))

os.system("rm main.out")
os.system("echo done.")
