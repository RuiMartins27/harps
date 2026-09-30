# Compiler and flags
CXX = mpicxx
CXXFLAGS = -Wall -std=c++17

# Source and output directories
SRC_DIR = src
BIN_DIR = bin
ANALYSIS_DIR = analysis/analysis_output
OUTPUTS_DIR = Outputs
MUSCLE3_HOME ?= /home/martins/muscle3

# PETSc location
PETSC_DIR := $(shell spack location -i petsc)

# Include and library directories
BASE_INCLUDE_DIRS = -I/usr/include -I./src -I$(PETSC_DIR)/include
MUSCLE_INCLUDE_DIRS = -I$(MUSCLE3_HOME)/include

BASE_LIBRARY_DIRS = -L$(PETSC_DIR)/lib
MUSCLE_LIBRARY_DIRS = -L$(MUSCLE3_HOME)/lib

BASE_LIBS = -lpetsc
MUSCLE_LIBS = -lmuscle -lymmsl -pthread

# Targets
TARGET = harps.exe
MUSCLE_TARGET = harps_muscle.exe
POLAR_TARGET = harps_polar.exe

COMMON_SOURCES = $(wildcard $(SRC_DIR)/*.cpp)
COMMON_OBJECTS = $(patsubst $(SRC_DIR)/%.cpp, $(BIN_DIR)/%.o, $(COMMON_SOURCES))

.PHONY: all muscle polar clean clean-png clean-txt clean-all

all: $(TARGET)

muscle: $(MUSCLE_TARGET)

polar: $(POLAR_TARGET)

# Ensure bin directory exists
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# Main builds
$(TARGET): $(BIN_DIR) $(COMMON_OBJECTS) $(BIN_DIR)/main_fields.o
	$(CXX) $(CXXFLAGS) $(COMMON_OBJECTS) $(BIN_DIR)/main_fields.o -o $@ \
	$(BASE_LIBRARY_DIRS) $(BASE_LIBS)

$(MUSCLE_TARGET): $(BIN_DIR) $(COMMON_OBJECTS) $(BIN_DIR)/main_muscle_harps.o
	$(CXX) $(CXXFLAGS) $(COMMON_OBJECTS) $(BIN_DIR)/main_muscle_harps.o -o $@ \
	$(BASE_LIBRARY_DIRS) $(MUSCLE_LIBRARY_DIRS) $(BASE_LIBS) $(MUSCLE_LIBS)

$(POLAR_TARGET): $(BIN_DIR) $(COMMON_OBJECTS) $(BIN_DIR)/main_muscle_harps_polar.o
	$(CXX) $(CXXFLAGS) $(COMMON_OBJECTS) $(BIN_DIR)/main_muscle_harps_polar.o -o $@ \
	$(BASE_LIBRARY_DIRS) $(MUSCLE_LIBRARY_DIRS) $(BASE_LIBS) $(MUSCLE_LIBS)


# Compile common sources
$(BIN_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(BASE_INCLUDE_DIRS) -c $< -o $@

# Compile main files
$(BIN_DIR)/main_fields.o: main_fields.cpp
	$(CXX) $(CXXFLAGS) $(BASE_INCLUDE_DIRS) -c $< -o $@

$(BIN_DIR)/main_muscle_harps.o: main_muscle_harps.cpp
	$(CXX) $(CXXFLAGS) $(BASE_INCLUDE_DIRS) $(MUSCLE_INCLUDE_DIRS) -c $< -o $@

$(BIN_DIR)/main_muscle_harps_polar.o: main_muscle_harps_polar.cpp
	$(CXX) $(CXXFLAGS) $(BASE_INCLUDE_DIRS) $(MUSCLE_INCLUDE_DIRS) -c $< -o $@


# Cleaning rules
clean:
	rm -rf $(BIN_DIR) $(TARGET) $(MUSCLE_TARGET) $(POLAR_TARGET)

clean-png:
	rm -f $(ANALYSIS_DIR)/*.{png,svg,pdf}

clean-txt:
	rm -f $(OUTPUTS_DIR)/*.txt

clean-all: clean clean-png clean-txt
