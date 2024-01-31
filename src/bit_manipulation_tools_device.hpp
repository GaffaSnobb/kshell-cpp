#ifndef BIT_MANIPULATION_TOOLS_DEVICE_HPP
#define BIT_MANIPULATION_TOOLS_DEVICE_HPP

#include <hip/hip_runtime.h>

namespace bittools_device
{
__device__ inline unsigned short reset_bit_and_count_swaps(unsigned long long& state, const unsigned short bit_to_reset)
{
    const unsigned long long mask = (1ULL << bit_to_reset) - 1; // Create a mask for all bits before the bit_to_reset.
    const unsigned int count = __popcll(state & mask);          // Count set bits before bit_to_set using popcll.
    state &= ~(1ULL << bit_to_reset);                           // Reset the bit at bit_to_reset.
    return count;
}

__device__ inline bool is_bit_set(unsigned long long state, const unsigned short bit_to_check)
{
    /*
    Check if a specific bit is set.
    */
    return (state & (1ULL << bit_to_check)) != 0;
}

__device__ inline unsigned short set_bit_and_count_swaps(unsigned long long& state, const unsigned short bit_to_set)
{
    const unsigned long long mask = (1ULL << bit_to_set) - 1; // Create a mask for bits before bit_to_set.
    const unsigned int count = __popcll(state & mask);        // Count set bits before bit_to_set using popcll.
    state |= (1ULL << bit_to_set);                            // Set the bit at bit_to_set.
    return count;
}

__device__ inline short negative_one_pow(const short exponent)
{
    /*
    Calculates -1 raised to the power of `exponent` efficiently using
    bitwise operations. This function  uses the property that -1 to an
    even power is 1, and to an odd power is -1. It checks the least
    significant bit of `exponent` to determine if `exponent` is odd or
    even.
    */
    return (exponent & 1) ? -1 : 1;
    // return (exponent % 2 == 0) ? 1 : -1;
}
}

#endif