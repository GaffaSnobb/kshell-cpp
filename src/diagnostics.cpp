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
    cout << "\n[DIAGNOSTICS]--------" << endl;
    print("n_valence_protons", interaction.model_space_protons.n_valence_nucleons);
    print("n_valence_neutrons", interaction.model_space_neutrons.n_valence_nucleons);
    print("m_dim", m_dim);
    print("m_dim**2", m_dim*m_dim);
    print("H size (MB): ", m_dim*m_dim*sizeof(double)/1000./1000.);
    print("H diag size (MB): ", m_dim*sizeof(double)/1000./1000.);
    cout << "[DIAGNOSTICS END]----\n" << endl;
}

void print_dtype_sizes()
{
    cout << "\n[DIAGNOSTICS]--------" << endl;
    cout << "bool              : " << sizeof(bool)               << "B (" << sizeof(bool)*8               << "b)" << endl;
    cout << "char              : " << sizeof(char)               << "B (" << sizeof(char)*8               << "b)" << endl;
    cout << "unsigned short    : " << sizeof(unsigned short)     << "B (" << sizeof(unsigned short)*8     << "b)" << " [0, " << (2ULL << (8*sizeof(unsigned short) - 1)) - 1       << "]"  << endl;
    cout << "short             : " << sizeof(short)              << "B (" << sizeof(short)*8              << "b)" << " [-"   << (2ULL << (8*sizeof(short) - 2))                    << ", " << (2ULL << (8*sizeof(short) - 2)) - 1     << "]" << endl;
    cout << "unsigned int      : " << sizeof(unsigned int)       << "B (" << sizeof(unsigned int)*8       << "b)" << " [0, " << (2ULL << (8*sizeof(unsigned int) - 1)) - 1         << "]"  << endl;
    cout << "int               : " << sizeof(int)                << "B (" << sizeof(int)*8                << "b)" << " [-"   << (2ULL << (8*sizeof(int) - 2))                      << ", " << (2ULL << (8*sizeof(int) - 2)) - 1       << "]" << endl;
    cout << "unsigned long     : " << sizeof(unsigned long)      << "B (" << sizeof(unsigned long)*8      << "b)" << " [0, " << (2ULL << (8*sizeof(unsigned long) - 1)) - 1        << "]"  << endl;
    cout << "long              : " << sizeof(long)               << "B (" << sizeof(long)*8               << "b)" << " [-"   << (2ULL << (8*sizeof(long) - 2))                     << ", " << (2ULL << (8*sizeof(long) - 2)) - 1      << "]" << endl;
    cout << "unsigned long long: " << sizeof(unsigned long long) << "B (" << sizeof(unsigned long long)*8 << "b)" << " [0, " << (2ULL << (8*sizeof(unsigned long long) - 1)) - 1   << "]"  << endl;
    cout << "long long         : " << sizeof(long long)          << "B (" << sizeof(long long)*8          << "b)" << " [-"   << (2ULL << (8*sizeof(long long) - 2))                << ", " << (2ULL << (8*sizeof(long long) - 2)) - 1 << "]" << endl;
    cout << "float             : " << sizeof(float)              << "B (" << sizeof(float)*8              << "b)" << endl;
    cout << "double            : " << sizeof(double)             << "B (" << sizeof(double)*8             << "b)" << endl;
    cout << "long double       : " << sizeof(long double)        << "B (" << sizeof(long double)*8        << "b)" << endl;
    cout << "[DIAGNOSTICS END]----\n" << endl;
}
}