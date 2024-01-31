#ifndef HAMILTONIAN_HPP
#define HAMILTONIAN_HPP

#include "data_structures.hpp"
#include "../external/eigen-3.4.0/Eigen/Dense"

namespace hamiltonian
{
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> create_hamiltonian_primitive_bit_representation(const Interaction& interaction);
void create_hamiltonian_device_dispatcher(const Interaction& interaction);
}
#endif