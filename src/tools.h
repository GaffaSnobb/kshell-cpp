#ifndef TOOLS_H
#define TOOLS_H

#include <vector>
#include "data_structures.h"

template <typename T>
void print_vector(const std::vector<T>& vec);
template <typename T>
void print_vector(std::string name, const std::vector<T>& vec);
template <typename T>
void print(std::string name, T value);

template <typename T>
std::vector<T> range(int start, int stop, int step);

void print(std::vector<OrbitalParameters> orbitals);

#include "tools.tpp"

#endif // TOOLS_H
