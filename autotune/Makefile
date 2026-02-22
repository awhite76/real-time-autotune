CXX      ?= g++
TARGET   ?= main

SRCS     := main.cpp pv.cpp
OBJS     := $(SRCS:.cpp=.o)

CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -Wpedantic
LDLIBS   ?= -lasound -lfftw3f -lm

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJS)