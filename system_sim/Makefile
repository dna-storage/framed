.PHONY: clean help .depend
BUILD_DIR=../dnastorage/system_sim/build
EXE_DIR=../dnastorage/system_sim/bin
SOURCE_DIR=../dnastorage/system_sim/src


SOURCES=$(wildcard ../dnastorage/system_sim/src/*.cpp)
OBJECTS=$(patsubst $(SOURCE_DIR)/%,$(BUILD_DIR)/%,$(SOURCES:.cpp=.o))




INCLUDE_DIR=../dnastorage/system_sim/include
UTILS_DIR=../dnastorage/system_sim/utils
CFLAGS=-O3 -std=c++11 -I$(INCLUDE_DIR) -I$(UTILS_DIR)

all: $(EXE_DIR)/system_sim

$(EXE_DIR)/system_sim:$(OBJECTS)
	mkdir -p $(EXE_DIR)
	$(CXX) $(CFLAGS) $^ -o $@

$(BUILD_DIR)/%.o: $(SOURCE_DIR)/%.cpp
	mkdir -p $(BUILD_DIR)
	$(CXX) $(CFLAGS) -c $< -o $@

clean:
	rm -r $(EXE_DIR) $(BUILD_DIR)  
