#ifndef BASIS_HPP
#define BASIS_HPP

#include <vector>
#include <bitset>
#include <stdint.h>

#include "parameters.hpp"
#include "data_structures.hpp"

namespace basis
{
std::vector<std::vector<uint16_t>> calculate_m_basis_states_vector_representation(const Interaction& interaction);
std::vector<std::bitset<n_bits_bitset>> calculate_m_basis_states_bitset_representation(const Interaction& interaction);
const std::vector<uint64_t> calculate_m_basis_states_primitive_bit_representation(
    const uint16_t n_proton_m_substates,
    const uint16_t n_neutron_m_substates,
    const uint16_t n_valence_protons,
    const uint16_t n_valence_neutrons,
    const std::vector<int16_t>& all_jz_values
);
}
#endif