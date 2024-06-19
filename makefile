PATHS = src/
PATHB = build/
PATHD = build/depends/
PATHO = build/objs/

PATHS_TEST = tests/
PATHB_TEST = $(PATHB)
PATHD_TEST = $(PATHD)
PATHO_TEST = $(PATHO)

CXX = hipcc
CXXFLAGS = -std=c++20 -MMD -MP -Wall -fopenmp -Ofast# -Wno-unused-result
CXXFLAGS_TEST = -std=c++20 -MMD -MP -Wall -fopenmp -O0
# CXXFLAGS = -std=c++20 -MMD -MP -Wall -fopenmp -O0 -g	# Debug flags.

SRCS = $(addprefix $(PATHS), tools.cpp loaders.cpp hamiltonian.cpp basis.cpp basic_solver.cpp generate_indices.cpp hamiltonian_bitset_representation.cpp hamiltonian_device.cpp diagnostics.cpp hip_wrappers.cpp lanczos.cpp)
OBJS = $(addprefix $(PATHO), $(notdir $(SRCS:.cpp=.o)))
DEPS = $(addprefix $(PATHD), $(notdir $(SRCS:.cpp=.d)))
EXEC = $(PATHB)run.out

SRCS_TEST = $(addprefix $(PATHS_TEST), hash_tests.cpp)
OBJS_TEST = $(addprefix $(PATHO), $(notdir $(SRCS_TEST:.cpp=.o)))
EXEC_TEST = $(PATHB)tests.out

# Build all
all: setup $(EXEC)

test: setup $(EXEC_TEST)

setup:
	mkdir -p $(PATHB) $(PATHD) $(PATHO)

# Link
$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Link (tests)
$(EXEC_TEST): $(OBJS_TEST)
	$(CXX) $(CXXFLAGS_TEST) -o $@ $^

# Rule to generate object files from cpp files
$(PATHO)%.o: $(PATHS)%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ -MF $(PATHD)$(notdir $*).d

# Rule to generate object files from cpp files (tests)
$(PATHO_TEST)%.o: $(PATHS_TEST)%.cpp
	$(CXX) $(CXXFLAGS_TEST) -c $< -o $@ -MF $(PATHD_TEST)$(notdir $*).d

-include $(DEPS)

clean:
	rm -rf $(PATHB)

.PHONY: all clean test play setup
