CC=g++
OS := $(shell uname)
ifeq ($(OS),Darwin)
  # Run MacOS commands
  CFLAGS=  -Wall  -Ofast -std=c++11  -flto -pipe -funit-at-a-time  -Wfatal-errors
  LDFLAGS=-flto -lpthread
else
  # check for Linux and add use of OpenMP
  CFLAGS=  -Wall  -Ofast -std=c++11  -flto -pipe -funit-at-a-time  -Wfatal-errors -fopenmp
  LDFLAGS=-flto -lpthread -fopenmp
endif


ifeq ($(gprof),1)
CFLAGS=-std=c++0x -pg -O4   -march=native
LDFLAGS=-pg
endif

ifeq ($(valgrind),1)
CFLAGS=-std=c++0x -O4 -g
LDFLAGS=-g
endif



EXEC=ES_simulation

all: $(EXEC)



ES_simulation.o: ES_simulation.cpp
	$(CC) -o $@ -c $< $(CFLAGS)


ES_simulation: ES_simulation.o
	$(CC) -o $@ $^ $(LDFLAGS)


clean:
	rm -rf *.o
	rm -rf $(EXEC)


rebuild: clean $(EXEC)
