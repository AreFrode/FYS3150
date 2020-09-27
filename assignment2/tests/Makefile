CC := g++
FLAGS = -larmadillo

PROGRAMS = test results
DEP = ../cpp_scripts/jacobisolver.cpp ../cpp_scripts/tridiagonalmatrix.cpp

.PHONY : all clean

all : $(PROGRAMS)

test:
	$(CC) -o test.out test.cpp $(DEP) $(FLAGS)

results:
	$(CC) -o results.out print_5x5.cpp $(DEP) $(FLAGS)

clean:
	rm -f *.out $(PROGRAMS)
