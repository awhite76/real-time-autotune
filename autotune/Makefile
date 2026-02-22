CXX      ?= g++
TARGET   ?= main

SRC      := main.cpp pitch.cpp time_stretch.cpp
OBJ      := $(SRC:.cpp=.o)

CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -Wpedantic
LDLIBS   ?= -lasound

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJ)
