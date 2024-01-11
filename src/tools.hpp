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
void print(std::vector<OrbitalParameters> orbitals);

template <typename T0, typename T1, typename T2, typename T3>
std::vector<T0> range(T1 start, T2 stop, T3 step);

short index(const std::vector<unsigned short>& vec, const unsigned short value);
inline short negative_one_pow(short exponent)   // Dont really know if the compiler cares about my inline... : (
{
    /*
    Calculates -1 raised to the power of `exponent` efficiently using
    bitwise operations. This function  uses the property that -1 to an
    even power is 1, and to an odd power is -1. It checks the least
    significant bit of `exponent` to determine if `exponent` is odd or
    even.

    Parameters
    ----------
    exponent : unsigned short
        The exponent to raise -1 to. Should be a small, non-negative
        integer.

    Returns
    -------
    int
        The result of raising -1 to the power of `exponent`. Returns 1
        if `exponent` is even, and -1 if `exponent` is odd.
    */
    return (exponent & 1) ? -1 : 1;
}

short check_existence_and_bisect(const std::vector<unsigned short>& vec, const unsigned short value);

#include "tools.tpp"

#endif // TOOLS_HPP