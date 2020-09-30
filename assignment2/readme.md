# Repository contaning source code for Project2 - FYS3150

The code has been compiled and written on a laptop running
ubuntu 20.04.1, g++ 9.3.0. The Python file main.py "Makefile" and the
Makefile for the tests is developed to work on the laptop
as just described, with no guarantees of it working elsewhere.

NOTE that the Python "Makefile" generates results and plots and
manipulates directory structures for storing them. Make sure
that the Python file "main.py" have the sufficient privileges.

## Contents
    main.py
    plotter.py
    quantum_plotter.py
    cpp_scripts/armadillosolver.cpp
    cpp_scripts/jacobisolver.cpp
    cpp_scripts/main.cpp
    cpp_scripts/tridiagonalmatrix.cpp
    cpp_scripts/tridiagonalmatrix.hpp
    tests/print_5x5.cpp
    tests/test.cpp

## Compilation and Execution
    Supplied with the repository is a Python script main.py.
    Running it with the commmand "Python main.py" results in
    6 different programs after a prompt asking for the number
    of grid points.
    The programs are:
    jacobi          - Solves the problem and times the jacobi solver
    arma            - Solves the problem and times the Armadillo solver
    analytical      - Solves the analytical buckling beam eigenvectors
    jacobi-arma     - Takes the difference between jacobi and Armadillo solver
    one-electron    - Solves the one-electron quantum case and prints out first 4 eigvals
    two-electrons   - Solves the two-electron quantum case and prints out first 4 eigvals

## Tests
    Supplied is a directory tests it conatins a Makefile with four options
    all         - compiles all programs to runnable .out files
    test        - compiles only the unit tests
    results     - compiles only the runtime examples
    clean       - removes all .out files

