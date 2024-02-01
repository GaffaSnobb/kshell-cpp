#ifndef HAMILTONIAN_DEVICE_HPP
#define HAMILTONIAN_DEVICE_HPP

#include "data_structures.hpp"

namespace hamiltonian_device
{
void create_hamiltonian_device_dispatcher(const Interaction& interaction, const Indices& indices, double* H);
}

#endif