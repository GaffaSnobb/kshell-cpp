#include <iostream>
#include "data_structures.hpp"
#include "tools.hpp"

using std::cout;
using std::endl;

namespace diagnostics
{
void print_hamiltonian_info(const Interaction& interaction)
{
    const unsigned int m_dim = interaction.basis_states.size();
    cout << "---------------------" << endl;
    print("n_valence_protons", interaction.model_space_protons.n_valence_nucleons);
    print("n_valence_neutrons", interaction.model_space_neutrons.n_valence_nucleons);
    print("m_dim", m_dim);
    print("m_dim**2", m_dim*m_dim);
    print("H size (MB): ", m_dim*m_dim*sizeof(double)/1000./1000.);
    print("H diag size (MB): ", m_dim*sizeof(double)/1000./1000.);
    cout << "---------------------" << endl;
}
}