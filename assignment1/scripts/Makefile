CC := g++
FLAGS = -O3

PROGRAMS := write error table_gen_spec lu_decomposition
DEP = functions.hpp lib.hpp

.PHONY : all clean

all : $(PROGRAMS)

create_png : write topng

write:
	$(CC) -o write_to_file.out write_file.cpp functions.cpp $(FLAGS)

topng:
	./write_to_file.out ../results/solution10to 1
	./write_to_file.out ../results/solution10to 2
	./write_to_file.out ../results/solution10to 3 && chmod a+x plotter.py && ./plotter.py ../results/solution10to1.txt ../results/solution10to2.txt ../results/solution10to3.txt 1

error:
	$(CC) -o relative_error.out error.cpp functions.cpp $(FLAGS)

table_gen_spec:
	$(CC) -o gen_spec_table.out runtime.cpp functions.cpp $(FLAGS)

lu_decomposition:
	$(CC) -o lu.out lu_compare.cpp lib.cpp functions.cpp $(FLAGS)

test:
	$(CC) -o test.out tests/test_functions.cpp functions.cpp
	./test.out
	rm -f test.out

clean:
	rm -f *.out $(PROGRAMS)
