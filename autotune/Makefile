CXX      ?= g++
TARGET   ?= main

SRCS     := main.cpp pv.cpp pitch.cpp time_stretch.cpp util.cpp 
OBJS     := $(SRCS:.cpp=.o)

CXXFLAGS ?= -g -O0 -std=c++17 -Wall -Wextra -Wpedantic
LDLIBS   ?= -lasound -lfftw3f -lm -lspeexdsp

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJS)
