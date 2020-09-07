CC := g++

write_to_file : write tofile

.PHONY : all clean

write:
	$(CC) -o write_to_file.out write_file.cpp functions.cpp -O3

tofile:
	./write_to_file.out ../solution10to 1
	./write_to_file.out ../solution10to 2
	./write_to_file.out ../solution10to 3 && chmod a+x plotter.py && ./plotter.py ../solution10to1.txt ../solution10to2.txt ../solution10to3.txt 1

thomas:
	$(CC) -o thomas.out thomas.cpp functions.cpp
