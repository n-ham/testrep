CC=g++
CXXFLAGS=-O3 -std=c++11 -Wall -Wextra -pedantic

default: egalgo1 egalgo2 egalgo3

egalgo1: egalgo1.cc algorithms.h
	$(CC) $(CXXFLAGS) -o egalgo1 egalgo1.cc

egalgo2: egalgo2.cc algorithms.h
	$(CC) $(CXXFLAGS) -o egalgo2 egalgo2.cc

egalgo3: egalgo3.cc algorithms.h
	$(CC) $(CXXFLAGS) -o egalgo3 egalgo3.cc

clean:
	rm -f egalgo1 egalgo2 egalgo3

test:
	./egalgo1
	./egalgo2
	./egalgo3
	
