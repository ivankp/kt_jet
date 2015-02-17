CC := g++
CFLAGS := -Wall -g

ROOT_CFLAGS := $(shell root-config --cflags)
ROOT_LIBS := $(shell root-config --libs)

.PHONY: all clean

EXE := test

all: $(EXE)

test: %: %.cc
	$(CC) $(CFLAGS) $(ROOT_CFLAGS) $(filter %.cc,$^) -o $@ $(ROOT_LIBS)

test: cluster.hh

clean:
	rm -f $(EXE) $(wildcard *.o)
