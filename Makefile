# ------------------------------------------------------------------
# Compiler and flags
# ------------------------------------------------------------------
CXX := g++

# Base language standard
CXXFLAGS := -std=c++17

# Choose debug vs release flags.
# - Xcode sets CONFIGURATION=Debug or CONFIGURATION=Release.
# - From the command line, "make debug" sets MAKECMDGOALS=debug.
ifneq (,$(filter Debug debug,$(CONFIGURATION) $(MAKECMDGOALS)))
    # Debug build: no optimisation, with debug info
    CXXFLAGS += -g
else
    # Release build: optimised
    CXXFLAGS += -O3 -march=native -Wno-deprecated-declarations
endif



# Include directories
INC_FLAGS := -Isrc/Phylib -Isrc/DecompositionTree -Isrc/StandardLikelihood -Isrc/Eigen
CXXFLAGS += $(INC_FLAGS)


# Source files
SRCS := \
	src/TidyTree.cpp \
	src/TreeHeights.cpp \
	src/RunSimulation.cpp \
	src/LvDLikelihoods.cpp \
	src/DecompositionTree/decompositionTree.cpp \
	src/DecompositionTree/decompositionLikelihood.cpp \
	src/StandardLikelihood/standardLikelihood.cpp \
	src/Phylib/sequences/phylip_seq.cpp

# Object files
OBJS := \
	build/TidyTree.o \
	build/RunSimulation.o \
	build/TreeHeights.o \
	build/LvDLikelihoods.o \
	build/DecompositionTree/decompositionTree.o \
	build/DecompositionTree/decompositionLikelihood.o \
	build/StandardLikelihood/standardLikelihood.o \
	build/Phylib/sequences/phylip_seq.o

# Object files excluding main files
OBJS_NO_MAIN := \
	build/DecompositionTree/decompositionTree.o \
	build/DecompositionTree/decompositionLikelihood.o \
	build/StandardLikelihood/standardLikelihood.o \
	build/Phylib/sequences/phylip_seq.o


# Executable names
TARGETS := RunSimulation TreeHeights LvDLikelihoods TidyTree


# Default rule: Compile in release mode
all: $(TARGETS)

# Xcode uses "build" as ACTION → map that to "all"
build: all

# Debug build (command line)
debug: clean
debug: $(TARGETS)

# Build TidyTree executable
TidyTree: build/TidyTree.o $(OBJS_NO_MAIN)
	$(CXX) $(CXXFLAGS)  $^ -o $@

# Build TreeHeights executable
TreeHeights: build/TreeHeights.o $(OBJS_NO_MAIN)
	$(CXX) $(CXXFLAGS)  $^ -o $@
	
RunSimulation: build/RunSimulation.o $(OBJS_NO_MAIN)
	$(CXX) $(CXXFLAGS)  $^ -o $@
	
LvDLikelihoods: build/LvDLikelihoods.o $(OBJS_NO_MAIN)
	$(CXX) $(CXXFLAGS)  $^ -o $@

build/TidyTree.o: src/TidyTree.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS)  -c $< -o $@
	
build/TreeHeights.o: src/TreeHeights.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS)  -c $< -o $@
	
build/RunSimulation.o: src/RunSimulation.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS)  -c $< -o $@
	
build/LvDLikelihoods.o: src/LvDLikelihoods.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS)  -c $< -o $@

build/DecompositionTree/decompositionTree.o: src/DecompositionTree/decompositionTree.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS)  -c $< -o $@
	
build/DecompositionTree/decompositionLikelihood.o: src/DecompositionTree/decompositionLikelihood.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS)  -c $< -o $@

build/StandardLikelihood/standardLikelihood.o: src/StandardLikelihood/standardLikelihood.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS)  -c $< -o $@

build/Phylib/sequences/phylip_seq.o: src/Phylib/sequences/phylip_seq.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS)  -c $< -o $@


# Clean rule: Remove object files and executables
clean:
	rm -rf build $(TARGETS)

# Phony targets
.PHONY: all clean
