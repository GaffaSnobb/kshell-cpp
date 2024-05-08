#ifndef DIAGNOSTICS_HPP
#define DIAGNOSTICS_HPP

#include "data_structures.hpp"

namespace diagnostics
{
const string DIAG_STR_START = "\n[DIAGNOSTICS]--------";
const string DIAG_STR_END   = "[DIAGNOSTICS END]----\n";
void print_hamiltonian_info(const Interaction& interaction);
void print_dtype_sizes();
void print_gpu_diagnostics(const Interaction& interaction, const Indices& indices);
}

#endif