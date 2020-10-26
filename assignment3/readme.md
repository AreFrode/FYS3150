# Repository containing source code for Project3 - FYS3150

All source code have been written and compiled on a desktop running ubunut 20.01.1, g++ 9.3.0. The Python file main.py serves as the main file to run the program, like a "Makefile".
The file Makefile in the tests directory is used to generate unit tests.

All code is compiled and executed on the laptop specfied above, no guarantees that it will work elsewhere.

NOTE that the Python "Makefile" main.py performs commands directly to the operating system using the "os" package to store results and plots. Make sure the Python file "Main.py" have sufficient privileges.

Specific for the plotting is that it is called specifically calling python3.8 which invokes Python 3.8.5 on the current system. This is due to the Python 3.9 installation on the current system being a bit problematic regarding Matplotlib. Thus, for this to run, it is required that python3.8 is installed

## Contents
* main.py
* plotter.py
* plot_energy.py
* plot_momentum.py
* twobody_plotter.py
* threebody_plotter.py
* manybody_plotter.py
* positions_and_vel.txt
* masses.txt
* cpp_scripts/main.cpp
* cpp_scripts/planet.cpp
* cpp_scripts/planet.hpp
* cpp_scripts/solarsystem.cpp
* cpp_scripts/solarsystem.hpp
* tests/test.cpp
* tests/Makefile
* report/report.pdf

