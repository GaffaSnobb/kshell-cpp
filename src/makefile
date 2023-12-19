CXX = clang++
TARGET = run
OBJS = loaders.o tools.o
CXXFLAGS = -std=c++20

# Build all
all: $(TARGET)

# Link
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

tools.o: tools.cpp tools.h
	$(CXX) $(CXXFLAGS) -c tools.cpp

loaders.o: loaders.cpp tools.h
	$(CXX) $(CXXFLAGS) -c loaders.cpp

clean:
	rm -f $(OBJS) $(TARGET)
