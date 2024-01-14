CXX = clang++
TARGET = run.out
SRCS = tools.cpp loaders.cpp hamiltonian.cpp basis.cpp basic_solver.cpp
OBJS = $(SRCS:.cpp=.o)
DEPS = $(OBJS:.o=.d)
CXXFLAGS = -std=c++20 -MMD

TEST_SRC = ../tests/hash_tests.cpp
TEST_OBJ = $(TEST_SRC:.cpp=.o)
TEST_EXEC = hash_test.out

# Build all
all: $(TARGET)

# Link
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

# Rule to generate object files from cpp files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to compile test source file
$(TEST_OBJ): $(TEST_SRC)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Rule to create the test executable
test: $(TEST_OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $(TEST_EXEC)

# Include the dependency files
-include $(DEPS)

# Clean rule
clean:
	rm -f $(OBJS) $(TARGET) $(DEPS) $(TEST_OBJ) $(TEST_EXEC)

.PHONY: all clean
