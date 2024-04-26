#include <iostream>
#include <vector>
#include <chrono>
#include <bitset>
#include <filesystem>
#include <stdexcept>
#include <stdint.h>

#include "parameters.hpp"
#include "tools.hpp"
#include "data_structures.hpp"
#include "../external/eigen-3.4.0/Eigen/Dense"

using std::cout;
using std::endl;
using std::vector;
using std::chrono::time_point;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;
using std::string;

// void complete_hermitian_matrix(Eigen::MatrixXcd& matrix) // For complex values.
void complete_hermitian_matrix(Eigen::MatrixXd& matrix)
{
    /*
    Copies the values from the upper triangle to the lower triangle.
    Does not copy the diagonal.

    This function will probably be removed soon since its not really
    necessary to explicitly store the other triangle when we know that
    the matrix is Hermitian.
    */
    auto start = timer();
    int32_t rows = matrix.rows();
    int32_t cols = matrix.cols();

    for (size_t row_idx = 0; row_idx < rows; row_idx++)
    {
        for (size_t col_idx = row_idx + 1; col_idx < cols; col_idx++)
        {
            // matrix(col_idx, row_idx) = std::conj(matrix(row_idx, col_idx));  // For complex values.
            matrix(col_idx, row_idx) = matrix(row_idx, col_idx);
        }
    }
    timer(start, "complete_hermitian_matrix");
}

void print(vector<OrbitalParameters> orbitals)
{
    for (const OrbitalParameters& orb : orbitals)
    {
        print("n", orb.n);
        print("l", orb.l);
        print("j", orb.j);
        print("tz", orb.tz);
        print_vector("jz", orb.jz);
        cout << endl;
    }
}

void print_vector(const vector<std::bitset<n_bits_bitset>>& vec)
{
    for (const std::bitset<n_bits_bitset> val : vec)
    {
        bool first_value = true;
        cout << val << ": [";
        for (size_t i = 0; i < n_bits_bitset; i++)
        {
            if (val.test(i))
            {
                if (not first_value) cout << ", ";
                cout << i;
                first_value = false;
            }
        }
        cout << "]" << endl;
    }
}

void print_vector(const vector<Key5>& vec)
{
    for (Key5 elem : vec)
    {
        cout << "{" << elem.a << ", " << elem.b << ", " << elem.c << ", " << elem.d << ", " << elem.e << "}," << endl;
    }
}

void print(const Key6& key)
{
    cout << "{" << key.j1 << ", " << key.m1 << ", " << key.j2 << ", " << key.m2 << ", " << key.J << ", " << key.M << "}," << endl;
}

int16_t index(const vector<uint16_t>& vec, const uint16_t value)
{
    /*
    Returns the index of the first occurence of `value` in `vec`. If
    `value` is not found, -1 is returned. Assumes that `vec` is short
    enough to be indexed by a int16_t.

    Using std::find and std::distance has practically identical
    computation time at -Ofast.
    
    Parameters
    ----------
    vec : vector<uint16_t>
        The vector to search in.
    
    value : uint16_t
        The value to search for.

    Returns
    -------
    int16_t
        The index of the first occurence of `value` in `vec`. If
        `value` is not found, -1 is returned.
    */
    // auto it = std::find(vec.begin(), vec.end(), value);

    // if (it != vec.end())
    // {
    //     int16_t idx = std::distance(vec.begin(), it);
    //     return idx;
    // }
    // else return -1;
    int16_t size = static_cast<int16_t>(vec.size());
    for (int16_t res = 0; res < size; res++)
    {
        if (vec[res] == value) return res;
    }
    return -1;
}

int16_t check_existence_and_bisect(const vector<uint16_t>& vec, const uint16_t value)
{
    int16_t size = static_cast<int16_t>(vec.size());
    for (int16_t res = 0; res < size; res++)
    {
        if (vec[res] > value) return res;
        else if (vec[res] == value) return -1;
    }
    return size;
}

time_point<high_resolution_clock> timer()
{
    return high_resolution_clock::now();
}

milliseconds timer(time_point<high_resolution_clock> start, string name)
{
    time_point<high_resolution_clock> stop = high_resolution_clock::now();
    milliseconds duration = std::chrono::duration_cast<milliseconds>(stop - start);
    cout << name << ": " << duration.count()/1000.0 << " s" << endl;
    return duration;
}

int64_t timer(time_point<high_resolution_clock> start)
{
    time_point<high_resolution_clock> stop = high_resolution_clock::now();
    milliseconds duration = std::chrono::duration_cast<milliseconds>(stop - start);
    return duration.count();
}

void check_if_file_exists(const std::filesystem::path& path)
{
    if (!std::filesystem::exists(path))
    {
        throw std::runtime_error("File does not exist: " + path.string());
    }
    if (!std::filesystem::is_regular_file(path))
    {
        throw std::runtime_error("The path is not a file: " + path.string());
    }
}