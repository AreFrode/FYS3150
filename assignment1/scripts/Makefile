CC := g++
FLAGS = -O3

create_png : write tofile

.PHONY : all clean

write:
	$(CC) -o write_to_file.out write_file.cpp functions.cpp $(FLAGS)

tofile:
	./write_to_file.out ../solution10to 1
	./write_to_file.out ../solution10to 2
	./write_to_file.out ../solution10to 3 && chmod a+x plotter.py && ./plotter.py ../solution10to1.txt ../solution10to2.txt ../solution10to3.txt 1

table_gen_spec:
	$(CC) -o gen_spec_table.out runtime.cpp functions.cpp $(FLAGS)

