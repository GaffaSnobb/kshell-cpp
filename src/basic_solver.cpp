#include <iostream>
#include <stdint.h>
#include "generate_indices.hpp"
#include "data_structures.hpp"
#include "parameters.hpp"
#include "loaders.hpp"
#include "tools.hpp"
#include "hamiltonian.hpp"
#include "hamiltonian_device.hpp"
#include "hamiltonian_bitset_representation.hpp"
#include "diagnostics.hpp"
#include "../external/eigen-3.4.0/Eigen/Dense"
#include "../external/eigen-3.4.0/Eigen/Eigenvalues"
// #include "../external/boost_1_84_0/boost/container_hash/hash.hpp"

using std::cout;
using std::endl;
using std::string;

int main(int argc, char* argv[])
{
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
    exit(0);

    const uint16_t n_valence_protons = std::stoi(argv[1]);
    const uint16_t n_valence_neutrons = std::stoi(argv[2]);
    const string interaction_filename = "./snt/w.snt";
    
    const Interaction interaction = loaders::load_interaction(interaction_filename, n_valence_protons, n_valence_neutrons);
    const Indices indices = indices::generate_indices(interaction);
    diagnostics::print_hamiltonian_info(interaction);
    const int m_dim = interaction.basis_states.size();

    double* H_host_reference = new double[m_dim*m_dim];
    for (size_t idx = 0; idx < m_dim*m_dim; idx++) H_host_reference[idx] = 0; // Not needed for the calculation, using this so that printing H is less cluttered.

    double* H_host = new double[m_dim*m_dim];
    for (size_t idx = 0; idx < m_dim*m_dim; idx++) H_host[idx] = 0; // Not needed for the calculation, using this so that printing H is less cluttered.

    std::vector<int64_t> timing;
    for (size_t i = 0; i < 1; i++)
    {
        auto start = timer();
        hamiltonian::create_hamiltonian_primitive_bit_representation_reference(interaction, indices, H_host_reference);
        cout << endl;
        // hamiltonian::create_hamiltonian_primitive_bit_representation_new(interaction, indices, H_host);
        hamiltonian_device::create_hamiltonian_device_dispatcher(interaction, indices, H_host);
        timing.push_back(timer(start));
    }
    cout << endl;
    print_vector(timing);
    cout << mean(timing) << endl;
    cout << "H_host_reference == H_host: " << compare_arrays(H_host_reference, H_host, m_dim*m_dim) << endl;

    // print_flattened_2d_array(H_host, m_dim);
    // print_flattened_2d_array(H_host_reference, m_dim);
    // Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> H_host_mapped(H_host, m_dim, m_dim);
    // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(H_host_mapped.transpose());
    // std::cout << "Eigenvalues:\n" << solver.eigenvalues() << std::endl;
    
    delete[] H_host_reference;
    delete[] H_host;
    return 0;
}