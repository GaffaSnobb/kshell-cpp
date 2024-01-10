CXX = clang++
TARGET = run.out
SRCS = tools.cpp loaders.cpp hamiltonian.cpp basis.cpp basic_solver.cpp
OBJS = $(SRCS:.cpp=.o)
DEPS = $(OBJS:.o=.d)
CXXFLAGS = -std=c++20 -MMD# -I../external/cppitertools

# Build all
all: $(TARGET)

# Link
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

# Rule to generate object files from cpp files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Include the dependency files
-include $(DEPS)

# Clean rule
clean:
	rm -f $(OBJS) $(TARGET) $(DEPS)

.PHONY: all clean
