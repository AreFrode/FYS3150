#!/usr/bin/env python3

import os
import sys

cpp_scripts = "./cpp_scripts/*.cpp"
compiler_flags = "-O3 -larmadillo"

os.system("echo compiling...")
os.system("g++ -o main.out" + " " + cpp_scripts + " " + compiler_flags)

os.system("echo executing...")
os.system("./main.out")

os.system("echo done.")
