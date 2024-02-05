#ifndef HAMILTONIAN_HPP
#define HAMILTONIAN_HPP

#include "data_structures.hpp"

namespace hamiltonian
{
void create_hamiltonian_primitive_bit_representation_new(const Interaction& interaction, const Indices& indices, double* H);
void create_hamiltonian_primitive_bit_representation_reference(const Interaction& interaction, const Indices& indices, double* H);
}
#endif