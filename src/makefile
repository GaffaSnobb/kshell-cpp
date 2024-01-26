# CXX = clang++
# CXX = /usr/lib/llvm-14/bin/clang++
CXX = hipcc
EXEC = run.out
SRCS = tools.cpp loaders.cpp hamiltonian.cpp basis.cpp basic_solver.cpp
OBJS = $(SRCS:.cpp=.o)
DEPS = $(OBJS:.o=.d)
CXXFLAGS = -std=c++20 -MMD -Wall -fopenmp -Ofast

TEST_SRC = ../tests/hash_tests.cpp
TEST_OBJ = $(TEST_SRC:.cpp=.o)
TEST_EXEC = hash_test.out

PLAYGROUND_SRC = playground.cpp
PLAYGROUND_OBJ = $(PLAYGROUND_SRC:.cpp=.o)
PLAYGROUND_EXEC = playground.out

# Build all
all: $(EXEC)

# Link
$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(EXEC) $(OBJS)

# Rule to generate object files from cpp files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

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

.PHONY: all clean
