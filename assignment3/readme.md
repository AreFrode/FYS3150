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

## Compilation and Execution
supplied with the repository is a Python script main.py. Running it with the command "Python main.py" results in a promot to enter number of mesh points for the simulation, then a final time. (MAKE SURE the final time is enetered without a period, i.e. as an int).
The program will then present 13 options:

* earth-fe            - Solves Earth-Sun using forward Euler
* earth-vv            - Solved Earth-Sun using velocity Verlet
* ellipse             - Solved Earth-Sun with IC resulting in an ellipse
* beta-3.5-circ       - Modifies beta = 2.5 for circular orbit
* beta-3.9-circ       - Modifies beta = 2.9 for circular orbit
* beta-4.0-circ       - Modifies beta = 3.0 for circular orbit
* beta-4.0-ellipse    - Modifies beta = 3.0 for elliptical orbit
* escape              - Solved Eart-Sun with IC resulting in Earth's escape velocity
* earth-jupiter       - Solves Earth-Jupiter around a static Sun
* earth-10jupiter     - Solves Earth-Jupiter around a static Sun. Jupter has 10 times more mass
* earth-1000jupiter   - Solved Earth-Jupiter around a static Sun. Jupiter has 1000 times more mass
* three-body          - Solves Earth-Sun-Jupiter system
* full-system         - Solved system with all planets in the solar system, and Pluto

## Tests
Supplied is a directory tests, it contains a Makefile with the options
* all   - compiles all programs to runnable .out.files
* test  - compiles unit tests
* clean - removes all .out files

## Bugs
There is an unknown bug in solarsystem.cpp, in the velocity Verlet solver. This results in multiplanetary systems with several moving bodies not simulating correctly.
