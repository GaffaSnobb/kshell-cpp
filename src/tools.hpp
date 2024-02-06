#ifndef TOOLS_HPP
#define TOOLS_HPP

#include <vector>
#include <chrono>
#include <bitset>
#include "parameters.hpp"
#include "data_structures.hpp"
#include "../external/eigen-3.4.0/Eigen/Dense"

// void print_bit_representation(unsigned short state);
void complete_hermitian_matrix(Eigen::MatrixXd& matrix);
template <typename T>
void print_bit_representation(const T& value);
void print_vector(const std::vector<std::bitset<n_bits_bitset>>& vec);
void print(const Key6& key);
void print_vector(const std::vector<Key5>& vec);
template <typename T>
void print_vector(const std::vector<T>& vec);
template <typename T>
void print_vector(const std::vector<std::vector<T>>& nested_vector);
template <typename T>
void print_vector(std::string name, const std::vector<T>& vec);
template <typename T>
void print(std::string name, T value);
void print(std::vector<OrbitalParameters> orbitals);
template <typename T0, typename T1, typename T2, typename T3>
std::vector<T0> range(T1 start, T2 stop, T3 step);
template <typename T>
double mean(std::vector<T> vec);
short index(const std::vector<unsigned short>& vec, const unsigned short value);
short check_existence_and_bisect(const std::vector<unsigned short>& vec, const unsigned short value);
std::chrono::milliseconds timer(std::chrono::time_point<std::chrono::high_resolution_clock> start, std::string name);
std::chrono::time_point<std::chrono::high_resolution_clock> timer();
long long timer(std::chrono::time_point<std::chrono::high_resolution_clock> start);
template <typename T1, typename T2>
bool compare_arrays(T1* arr1, T1* arr2, T2 size);
template <typename T1, typename T2>
void print_flattened_2d_array(T1* arr, T2 size);
template <typename T0, typename T1, typename T2>
void print_loop_timer(std::vector<T0>& loop_timings, T1 idx, T2 n_iterations);

#include "tools.tpp"
#endif // TOOLS_HPP