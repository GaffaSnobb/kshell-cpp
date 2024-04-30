#ifndef DIAGNOSTICS_HPP
#define DIAGNOSTICS_HPP

#include "data_structures.hpp"

namespace diagnostics
{
void print_hamiltonian_info(const Interaction& interaction);
void print_dtype_sizes();
void print_gpu_mem_usage(const Interaction& interaction, const Indices& indices);
}

#endif