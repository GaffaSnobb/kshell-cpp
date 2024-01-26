#ifndef HAMILTONIAN_BITSET_REPRESENTATION_HPP
#define HAMILTONIAN_BITSET_REPRESENTATION_HPP

#include "data_structures.hpp"
#include "../external/eigen-3.4.0/Eigen/Dense"

namespace hamiltonian_bitset
{
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> create_hamiltonian_bit_representation(const Interaction& interaction);
}

#endif