# Makefile

CXX      := g++
CXXFLAGS := -O2 -std=c++17 -Wall -Wextra -Wno-unused-parameter
LDFLAGS  := 

# Files
HDRS := wsbm_io_bin.hpp WSBMGenerator.hpp
SRCS := WSBMGenerator.cpp wsbm_write.cpp community_assignment.cpp
OBJS := WSBMGenerator.o

# Binaries
WRITER := wsbm_write
READER := main

.PHONY: all clean run_write run_main

all: $(WRITER) $(READER)

$(WRITER): wsbm_write.cpp $(OBJS) $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ wsbm_write.cpp $(OBJS) $(LDFLAGS)

$(READER): community_assignment.cpp $(OBJS) $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ community_assignment.cpp $(OBJS) $(LDFLAGS)

WSBMGenerator.o: WSBMGenerator.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -c WSBMGenerator.cpp -o $@

clean:
	rm -f $(WRITER) $(READER) $(OBJS) wsbm_data.bin

N    ?= 100
K    ?= 4
PIN  ?= 0.6
POUT ?= 0.1
OUT  ?= wsbm_data.bin
IN   ?= $(OUT)

run_write: $(WRITER)
	./$(WRITER) $(N) $(K) $(PIN) $(POUT) $(OUT)

run_main: $(READER)
	./$(READER) $(IN) $(N) $(K)
