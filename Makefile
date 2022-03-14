CC=mpicxx
CXXFLAGS=-Wall -O3 -pedantic -g
LDLIBS=-lboost_program_options
TARGET=main
NP=6
ARGS1=--dt 0.001 --T 100 --Nx 101 --Ny 101 --a 0.75 --b 0.06   --eps 50.0 --mu1 5.0 --mu2 0.0
ARGS2=--dt 0.001 --T 100 --Nx 251 --Ny 251 --a 0.75 --b 0.06   --eps 13.0 --mu1 5.0 --mu2 0.0
ARGS3=--dt 0.001 --T 100 --Nx 101 --Ny 101 --a 0.50 --b 0.10   --eps 50.0 --mu1 5.0 --mu2 0.0
ARGS4=--dt 0.001 --T 100 --Nx 151 --Ny 81  --a 0.75 --b 0.0001 --eps 12.5 --mu1 1.0 --mu2 0.01

default: $(TARGET)
all: $(TARGET)

$(TARGET): $(TARGET).o ReactionDiffusion.o

%.o: %.cpp
	$(CC) $(CXXFLAGS) -o $@ -c $< $(LDLIBS)

.PHONY: clean

test1: $(TARGET)
	mpiexec -np $(NP) ./$(TARGET) $(ARGS1)

test2: $(TARGET)
	mpiexec -np $(NP) ./$(TARGET) $(ARGS2)

test3: $(TARGET)
	mpiexec -np $(NP) ./$(TARGET) $(ARGS3)

test4: $(TARGET)
	mpiexec -np $(NP) ./$(TARGET) $(ARGS4)

clean:
	rm -rf $(TARGET) *.o output.txt
