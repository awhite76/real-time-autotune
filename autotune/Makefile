CXX      ?= g++
TARGET   ?= main

SRCS     := main.cpp pv.cpp pitch.cpp time_stretch.cpp util.cpp 
OBJS     := $(SRCS:.cpp=.o)

TEST_RESAMPLE_SRCS := test_resample.cpp time_stretch.cpp util.cpp
TEST_RESAMPLE_OBJS := $(TEST_RESAMPLE_SRCS:.cpp=.o)

TEST_PV_SRCS := test_pv.cpp pv.cpp util.cpp
TEST_PV_OBJS := $(TEST_PV_SRCS:.cpp=.o)

CXXFLAGS ?= -g -O0 -std=c++17 -Wall -Wextra -Wpedantic
LDLIBS   ?= -lasound -lfftw3f -lm -lspeexdsp

.PHONY: all clean test

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

test_resample: $(TEST_RESAMPLE_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

test_pv: $(TEST_PV_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJS)
