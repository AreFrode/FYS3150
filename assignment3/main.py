#!/usr/bin/env python3

import os
import sys

N = input("Number of mesh points: ")
t = input("Final calculation time: ")
prog = input("Choose a program \n* binary\n: ")

filename_data = "_".join([prog, "N", N, "t", t]) + ".txt"
data_path = "/".join([".", "results", prog])

cpp_scripts = "./cpp_scripts/*.cpp"

os.system("echo compiling...")
os.system("g++ -o main.out" + " " + cpp_scripts)

os.system("echo executing...")
os.system("./main.out " + N + " " + t + " " + prog)

if not os.path.exist(data_path):
    os.makedirs(data_path)

os.system("echo done.")
