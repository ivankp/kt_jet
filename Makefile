CC := g++
CFLAGS := -Wall -O3

.PHONY: all clean

EXE := test

all: $(EXE)

test: %: %.cc
	$(CC) $(CFLAGS) $(filter %.cc,$^) $(filter %.o,$^) -o $@

fjcore.o: %.o: fjcore.cc fjcore.hh
	$(CC) $(CFLAGS) -c $(filter %.cc,$^) -o $@

test: cluster.hh Timer.h fjcore.o

clean:
	rm -f $(EXE) $(wildcard *.o)
