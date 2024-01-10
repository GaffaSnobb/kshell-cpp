#ifndef TOOLS_HPP
#define TOOLS_HPP

#include <vector>
#include "data_structures.hpp"

template <typename T>
void print_vector(const std::vector<T>& vec);
template <typename T>
void print_vector(const std::vector<std::vector<T>>& nested_vector);
template <typename T>
void print_vector(std::string name, const std::vector<T>& vec);
template <typename T>
void print(std::string name, T value);

template <typename T0, typename T1, typename T2, typename T3>
std::vector<T0> range(T1 start, T2 stop, T3 step);

void print(std::vector<OrbitalParameters> orbitals);

#include "tools.tpp"

#endif // TOOLS_HPP