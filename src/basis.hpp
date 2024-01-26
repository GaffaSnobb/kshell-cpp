#ifndef BASIS_HPP
#define BASIS_HPP

#include <vector>
#include <bitset>
#include "parameters.hpp"
#include "data_structures.hpp"

namespace basis
{
std::vector<std::vector<unsigned short>> calculate_m_basis_states_vector_representation(const Interaction& interaction);
std::vector<std::bitset<n_bits_bitset>> calculate_m_basis_states_bitset_representation(const Interaction& interaction);
std::vector<unsigned long long> calculate_m_basis_states_primitive_bit_representation(const Interaction& interaction);
}
#endif