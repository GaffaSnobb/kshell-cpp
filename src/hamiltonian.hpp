#ifndef HAMILTONIAN_HPP
#define HAMILTONIAN_HPP

#include "data_structures.hpp"
#include "../external/eigen-3.4.0/Eigen/Dense"

void create_hamiltonian(const Interaction& interaction);
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> create_hamiltonian_bit_representation(const Interaction& interaction);
#endif