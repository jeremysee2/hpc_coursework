CC=g++
CXXFLAGS=-Wall -O3 -pedantic
LDLIBS=-llapack -lblas
TARGET=main

default: $(TARGET)
all: $(TARGET)

$(TARGET): $(TARGET).o

%.o: %.cpp
	$(CC) $(CXXFLAGS) -o $@ -c $< $(LDLIBS)

.PHONY: clean

clean:
	rm -rf $(TARGET) *.o
