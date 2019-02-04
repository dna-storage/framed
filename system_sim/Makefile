.PHONY: clean help .depend

EXE=system_sim
SOURCES=$(wildcard ../dnastorage/system_sim/src/*.cpp)
OBJECTS=$(SOURCES:.cpp=.o)
INCLUDES=../dnastorage/system_sim/headers
INCLUDES+=../dnastorage/system_sim/utils

CFLAGS=-O3 -I$(INCLUDES)

all: $(EXE)

$(EXE):$(OBJECTS)
	$(CXX) $(CFLAGS) $^ -o $@

%.o: %.cpp
	$(CXX) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(EXE)
