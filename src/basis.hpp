#ifndef BASIS_HPP
#define BASIS_HPP

#include <vector>
#include <bitset>
#include "parameters.hpp"
#include "data_structures.hpp"

std::vector<std::vector<unsigned short>> calculate_m_basis_states(const Interaction& interaction);
std::vector<std::bitset<n_bits_bitset>> calculate_m_basis_states_bit_representation(const Interaction& interaction);

#endif