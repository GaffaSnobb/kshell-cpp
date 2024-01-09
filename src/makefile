CXX = clang++
TARGET = run.out
OBJS = loaders.o tools.o basic_solver.o
CXXFLAGS = -std=c++20

# Build all
all: $(TARGET)

# Link
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

tools.o: tools.cpp tools.hpp
	$(CXX) $(CXXFLAGS) -c tools.cpp

loaders.o: loaders.cpp loaders.hpp tools.hpp data_structures.hpp
	$(CXX) $(CXXFLAGS) -c loaders.cpp

basic_solver.o: basic_solver.cpp loaders.hpp tools.hpp data_structures.hpp
	$(CXX) $(CXXFLAGS) -c basic_solver.cpp

clean:
	rm -f $(OBJS) $(TARGET)
