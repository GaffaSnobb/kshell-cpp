#include <iostream>
#include <vector>
#include <chrono>
#include <bitset>
#include "parameters.hpp"
#include "tools.hpp"
#include "data_structures.hpp"
#include "../external/eigen-3.4.0/Eigen/Dense"

using std::cout;
using std::endl;
// using std::chrono::high_resolution_clock;
// using std::chrono::duration_cast;
// using std::chrono::milliseconds;

// void print_bit_representation(unsigned short state)
// {
//     cout << "value: " << state << endl;
//     std::bitset<16> bits(state);
//     cout << "bits: ";
//     for (int i = 15; i >= 0; --i) cout << bits[i] << " ";
//     cout << endl;
// }

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
    int rows = matrix.rows();
    int cols = matrix.cols();

    for (int row_idx = 0; row_idx < rows; row_idx++)
    {
        for (int col_idx = row_idx + 1; col_idx < cols; col_idx++)
        {
            // matrix(col_idx, row_idx) = std::conj(matrix(row_idx, col_idx));  // For complex values.
            matrix(col_idx, row_idx) = matrix(row_idx, col_idx);
        }
    }
    timer(start, "complete_hermitian_matrix");
}

void print(std::vector<OrbitalParameters> orbitals)
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

void print_vector(const std::vector<std::bitset<n_bits_bitset>>& vec)
{
    for (const std::bitset<n_bits_bitset> val : vec)
    {
        bool first_value = true;
        cout << val << ": [";
        for (int i = 0; i < n_bits_bitset; i++)
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

void print_vector(const std::vector<Key5>& vec)
{
    for (Key5 elem : vec)
    {
        cout << "{" << elem.a << ", " << elem.b << ", " << elem.c << ", " << elem.d << ", " << elem.e << "}," << endl;
    }
}

short index(const std::vector<unsigned short>& vec, const unsigned short value)
{
    /*
    Returns the index of the first occurence of `value` in `vec`. If
    `value` is not found, -1 is returned. Assumes that `vec` is short
    enough to be indexed by a short.

    Using std::find and std::distance has practically identical
    computation time at -Ofast.
    
    Parameters
    ----------
    vec : std::vector<unsigned short>
        The vector to search in.
    
    value : unsigned short
        The value to search for.

    Returns
    -------
    short
        The index of the first occurence of `value` in `vec`. If
        `value` is not found, -1 is returned.
    */
    // auto it = std::find(vec.begin(), vec.end(), value);

    // if (it != vec.end())
    // {
    //     short idx = std::distance(vec.begin(), it);
    //     return idx;
    // }
    // else return -1;
    short size = static_cast<short>(vec.size());
    for (short res = 0; res < size; res++)
    {
        if (vec[res] == value) return res;
    }
    return -1;
}

short check_existence_and_bisect(const std::vector<unsigned short>& vec, const unsigned short value)
{
    short size = static_cast<short>(vec.size());
    for (short res = 0; res < size; res++)
    {
        if (vec[res] > value) return res;
        else if (vec[res] == value) return -1;
    }
    return size;
}

std::chrono::time_point<std::chrono::high_resolution_clock> timer()
{
    return std::chrono::high_resolution_clock::now();
}

std::chrono::milliseconds timer(std::chrono::time_point<std::chrono::high_resolution_clock> start, std::string name)
{
    std::chrono::time_point<std::chrono::high_resolution_clock> stop = std::chrono::high_resolution_clock::now();
    std::chrono::milliseconds duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    cout << name << ": " << duration.count()/1000.0 << " s" << endl;
    return duration;
}

long long timer(std::chrono::time_point<std::chrono::high_resolution_clock> start)
{
    std::chrono::time_point<std::chrono::high_resolution_clock> stop = std::chrono::high_resolution_clock::now();
    std::chrono::milliseconds duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    return duration.count();
}