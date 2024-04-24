PATHS = src/
PATHB = build/
PATHD = build/depends/
PATHO = build/objs/

CXX = hipcc
EXEC = $(PATHB)run.out
SRCS = tools.cpp loaders.cpp hamiltonian.cpp basis.cpp basic_solver.cpp generate_indices.cpp hamiltonian_bitset_representation.cpp hamiltonian_device.cpp diagnostics.cpp hip_wrappers.cpp
OBJS = $(addprefix $(PATHO), $(SRCS:.cpp=.o))
DEPS = $(addprefix $(PATHD), $(SRCS:.cpp=.d))
CXXFLAGS = -std=c++20 -MMD -MP -Wall -fopenmp -Ofast# -Wno-unused-result

TEST_SRC = ../tests/hash_tests.cpp
TEST_OBJ = $(TEST_SRC:.cpp=.o)
TEST_EXEC = hash_test.out

PLAYGROUND_SRC = playground.cpp tools.cpp
PLAYGROUND_OBJ = $(PLAYGROUND_SRC:.cpp=.o)
PLAYGROUND_EXEC = playground.out

# Build all
all: setup $(EXEC)

setup:
	mkdir -p $(PATHB) $(PATHD) $(PATHO)

# Link
$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Rule to generate object files from cpp files
$(PATHO)%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ -MF $(PATHD)$*.d

# Rule to compile test source file
$(TEST_OBJ): $(TEST_SRC)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Rule to create the test executable
test: $(TEST_OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $(TEST_EXEC)

play:
	hipcc $(PLAYGROUND_SRC) -o $(PLAYGROUND_EXEC) -fopenmp -Wno-unused-result

# Include the dependency files
-include $(DEPS)

clean:
	rm -f $(OBJS) $(EXEC) $(DEPS) $(TEST_OBJ) $(TEST_EXEC) $(PLAYGROUND_OBJ) $(PLAYGROUND_EXEC)

.PHONY: all clean test play setup
