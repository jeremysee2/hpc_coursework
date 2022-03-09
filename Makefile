CC=g++
CXXFLAGS=-Wall -O3 -pedantic
LDLIBS=-llapack -lblas
TARGET=main
ARGS1=--dt 0.001 --T 100 --Nx 101 --Ny 101 --a 0.75 --b 0.06   --mu1 50.0 --mu2 5.0 --eps 0.0
ARGS2=--dt 0.001 --T 100 --Nx 251 --Ny 251 --a 0.75 --b 0.06   --mu1 13.0 --mu2 5.0 --eps 0.0
ARGS3=--dt 0.001 --T 100 --Nx 101 --Ny 101 --a 0.50 --b 0.10   --mu1 50.0 --mu2 5.0 --eps 0.0
ARGS4=--dt 0.001 --T 100 --Nx 151 --Ny 81  --a 0.75 --b 0.0001 --mu1 12.5 --mu2 1.0 --eps 0.01

default: $(TARGET)
all: $(TARGET)

$(TARGET): $(TARGET).o

%.o: %.cpp
	$(CC) $(CXXFLAGS) -o $@ -c $< $(LDLIBS)

.PHONY: clean

test1: $(TARGET)
	./$(TARGET) $(ARGS1)

test2: $(TARGET)
	./$(TARGET) $(ARGS2)

test3: $(TARGET)
	./$(TARGET) $(ARGS3)

test4: $(TARGET)
	./$(TARGET) $(ARGS4)

clean:
	rm -rf $(TARGET) *.o
