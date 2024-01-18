#ifndef TOOLS_HPP
#define TOOLS_HPP

#include <vector>
#include <chrono>
#include <bitset>
#include "parameters.hpp"
#include "data_structures.hpp"

void print_vector(const std::vector<std::bitset<n_bits_bitset>>& vec);
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

short index(const std::vector<unsigned short>& vec, const unsigned short value);
inline short negative_one_pow(short exponent)   // Dont really know if the compiler cares about my inline... : (
{
    /*
    Calculates -1 raised to the power of `exponent` efficiently using
    bitwise operations. This function  uses the property that -1 to an
    even power is 1, and to an odd power is -1. It checks the least
    significant bit of `exponent` to determine if `exponent` is odd or
    even.

    NOTE: std::pow is practically identically fast as the bit-check and
    modulo solutions with -O0 and -O1. There is however a slight
    advantage at -O2 and higher:

    np.mean(pow_O0)/1000, np.mean(bit_O0)/1000, np.mean(modulo_O0)/1000
    (19.54375, 19.70725, 19.9485)

    np.mean(pow_O1)/1000, np.mean(bit_O1)/1000, np.mean(modulo_O1)/1000
    (1.49025, 1.48375, 1.4715)

    np.mean(pow_O2)/1000, np.mean(bit_O2)/1000, np.mean(modulo_O2)/1000
    (1.50675, 1.461, 1.4745)

    np.mean(pow_O3)/1000, np.mean(bit_O3)/1000, np.mean(modulo_O3)/1000
    (1.50575, 1.4795, 1.4725)

    np.mean(pow_Ofast)/1000, np.mean(bit_Ofast)/1000, np.mean(modulo_Ofast)/1000
    (1.542, 1.4765, 1.475)

    Parameters
    ----------
    exponent : unsigned short
        The exponent to raise -1 to. Should be a small, non-negative
        integer.

    Returns
    -------
    short
        The result of raising -1 to the power of `exponent`. Returns 1
        if `exponent` is even, and -1 if `exponent` is odd.
    */
    // return std::pow(-1, exponent);
    return (exponent & 1) ? -1 : 1;
    // return (exponent % 2 == 0) ? 1 : -1;
}

short check_existence_and_bisect(const std::vector<unsigned short>& vec, const unsigned short value);

inline unsigned short set_bit_and_count_swaps(std::bitset<n_bits_bitset>& state, unsigned short bit_to_set)
{
    /*
    Set bit number `bit_to_set` and count how many bits before
    `bit_to_set` are set / count how many operator swaps must be
    performed to place the creation operator correctly.

    We need to know how many bits before the bit we want
    to annihilate are set. This corresponds to how many
    times we must swap the positions of operators and
    thus if we get a phase of +1 or -1. This is because
    we have to make sure that the annihilation operator
    is placed next to the creation operator it tries to
    annihilate:

        c_0 | (0, 3) > = c_0 c_0^\dagger c_3^\dagger | core >
                        = c_3^\dagger | core >
                        = | (3) >

        c_0 | (3, 0) > = c_0 c_3^\dagger c_0^\dagger | core >
                        = - c_0 c_0^\dagger c_3^\dagger | core >
                        = - c_3^\dagger | core >
                        = - | (3) >
    */
    unsigned short count = 0;
    // for (unsigned short i = 0; i < bit_to_set; i++) if (state.test(i)) count++;
    for (unsigned short i = 0; i < bit_to_set; i++) if (state[i]) count++;
    state.set(bit_to_set);
    return count;
}

inline unsigned short reset_bit_and_count_swaps(std::bitset<n_bits_bitset>& state, unsigned short bit_to_reset)
{
    /*
    Reset bit number `bit_to_reset` and count how many bits before
    `bit_to_reset` are set / count how many operator swaps must be
    performed to place the annihilation operator correctly.
    */
    unsigned short count = 0;
    // for (unsigned short i = 0; i < bit_to_reset; i++) if (state.test(i)) count++;
    for (unsigned short i = 0; i < bit_to_reset; i++) if (state[i]) count++;
    state.reset(bit_to_reset);
    return count;
}

std::chrono::milliseconds timer(std::chrono::time_point<std::chrono::high_resolution_clock> start, std::string name);
std::chrono::time_point<std::chrono::high_resolution_clock> timer();

#include "tools.tpp"

#endif // TOOLS_HPP