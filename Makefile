CC := g++
CCFLAGS := -Wall -O2

.PHONY: all clean

EXE := test

all: $(EXE)

$(EXE): %: %.cc
	$(CC) $(CCFLAGS) $^ -o $@
	
%.o: %.cc %.hh
	$(CC) $(CCFLAGS) -c $(filter %.cc,$^) -o $@
	
test: cluster.o

clean:
	rm -f $(EXE) *.o

