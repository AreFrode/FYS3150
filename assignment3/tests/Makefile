CC := g++

PROGRAMS = test
DEP = ../cpp_scripts/planet.cpp ../cpp_scripts/solarsystem.cpp

.PHONY : all clean

all : $(PROGRAMS)

test:
		$(CC) -o test.out test.cpp $(DEP)

clean:
		rm -f *.out $(PROGRAMS)
