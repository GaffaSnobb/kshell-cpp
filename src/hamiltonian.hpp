#ifndef HAMILTONIAN_HPP
#define HAMILTONIAN_HPP

#include "data_structures.hpp"

namespace hamiltonian
{
void create_hamiltonian_primitive_bit_representation(const Interaction& interaction, const Indices& indices, double* H);
void create_hamiltonian_primitive_bit_representation_tmp(const Interaction& interaction, const Indices& indices, double* H);
}
#endif