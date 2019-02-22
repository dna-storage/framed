.PHONY: clean help .depend

EXE=system_sim
SOURCES=$(wildcard ../dnastorage/system_sim/src/*.cpp)
OBJECTS=$(SOURCES:.cpp=.o)
INCLUDES=../dnastorage/system_sim/headers

CFLAGS=-O3 -std=c++11 -I$(INCLUDES) -I../dnastorage/system_sim/utils

all: $(EXE)

$(EXE):$(OBJECTS)
	$(CXX) $(CFLAGS) $^ -o $@

%.o: %.cpp
	$(CXX) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(EXE)
