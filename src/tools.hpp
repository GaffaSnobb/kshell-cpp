#ifndef TOOLS_HPP
#define TOOLS_HPP

#include <vector>
#include <chrono>
#include <bitset>
#include <filesystem>
#include <stdint.h>

#include "parameters.hpp"
#include "data_structures.hpp"
#include "../external/eigen-3.4.0/Eigen/Dense"

using std::vector;
using std::chrono::time_point;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;
using std::string;

// void print_bit_representation(uint16_t state);
bool compare_with_ref_matrix(double* arr, string ref_name, size_t size);
void complete_symmetric_matrix(Eigen::MatrixXd& matrix);
vector<double> read_symmetric_matrix(const string& filename);
void print(vector<OrbitalParameters> orbitals);
void print_vector(const vector<std::bitset<n_bits_bitset>>& vec);
void print_vector(const vector<Key5>& vec);
void print(const Key6& key);
int16_t index(const vector<uint16_t>& vec, const uint16_t value);
int16_t check_existence_and_bisect(const vector<uint16_t>& vec, const uint16_t value);
time_point<high_resolution_clock> timer();
milliseconds timer(time_point<high_resolution_clock> start, string name);
int64_t timer(time_point<high_resolution_clock> start);
void check_if_file_exists(const std::filesystem::path& path);

template <typename T1, typename T2>
void print_flattened_2d_array(const T1* arr, const T2 size);

template <typename T1, typename T2>
void write_flattened_2d_array_to_file(const T1* arr, const T2 size, const string& filename);

template <typename T1, typename T2, typename T3>
bool compare_arrays(T1* arr1, T2* arr2, T3 size);

template <typename T1, typename T2, typename T3>
bool compare_arrays_upper_triangle(T1* arr1, vector<T2> arr2, T3 size);

template <typename T>
void print_bit_representation(const T& value);

template <typename T>
void print_vector(const vector<T>& vec);

template <typename T>
void print_vector(const vector<vector<T>>& nested_vector);

template <typename T>
void print_vector(string name, const vector<T>& vec);

template <typename T>
void print(string name, T value);

template <typename T0, typename T1, typename T2, typename T3>
vector<T0> range(T1 start, T2 stop, T3 step);

template <typename T>
double mean(vector<T> vec);

template <typename T0, typename T1, typename T2>
void print_loop_timer(vector<T0>& loop_timings, T1 idx, T2 n_iterations);

#include "tools.tpp"
#endif // TOOLS_HPP