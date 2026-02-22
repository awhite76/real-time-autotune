CXX      ?= g++
TARGET   ?= main
SRC      ?= main.cpp

CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -Wpedantic
LDLIBS   ?= -lasound -lfftw3f -lm

# CXXFLAGS += $(shell pkg-config --cflags alsa)
# LDLIBS   += $(shell pkg-config --libs alsa)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

clean:
	rm -f $(TARGET)
