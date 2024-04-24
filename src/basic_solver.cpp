#include <iostream>
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

int main(int argc, char* argv[])
{
    unsigned short n_valence_protons = std::stoi(argv[1]);
    unsigned short n_valence_neutrons = std::stoi(argv[2]);
    std::string interaction_filename = "../snt/w.snt";
    
    const Interaction interaction = loaders::load_interaction(interaction_filename, n_valence_protons, n_valence_neutrons);
    const Indices indices = indices::generate_indices(interaction);
    diagnostics::print_hamiltonian_info(interaction);
    int m_dim = interaction.basis_states.size();

    double* H_host_reference = new double[m_dim*m_dim];
    for (int idx = 0; idx < m_dim*m_dim; idx++) H_host_reference[idx] = 0; // Not needed for the calculation, using this so that printing H is less cluttered.

    double* H_host = new double[m_dim*m_dim];
    for (int idx = 0; idx < m_dim*m_dim; idx++) H_host[idx] = 0; // Not needed for the calculation, using this so that printing H is less cluttered.

    std::vector<int> timing;
    for (int i = 0; i < 1; i++)
    {
        auto start = timer();
        hamiltonian::create_hamiltonian_primitive_bit_representation_reference(interaction, indices, H_host_reference);
        cout << endl;
        // hamiltonian::create_hamiltonian_primitive_bit_representation_new(interaction, indices, H_host);
        // hamiltonian_device::create_hamiltonian_device_dispatcher(interaction, indices, H_host);
        timing.push_back(timer(start));
    }
    cout << endl;
    print_vector(timing);
    cout << mean(timing) << endl;
    cout << "H_host_reference == H_host: " << compare_arrays(H_host_reference, H_host, m_dim*m_dim) << endl;

    // print_flattened_2d_array(H_host, m_dim);
    // Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> H_host_mapped(H_host, m_dim, m_dim);
    // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(H_host_mapped.transpose());
    // std::cout << "Eigenvalues:\n" << solver.eigenvalues() << std::endl;
    
    delete[] H_host_reference;
    delete[] H_host;
    return 0;
}