#include <iostream>
#include <stdint.h>
#include <hip/hip_runtime.h>

#include "generate_indices.hpp"
#include "data_structures.hpp"
#include "parameters.hpp"
#include "loaders.hpp"
#include "tools.hpp"
#include "hip_wrappers.hpp"
#include "hamiltonian.hpp"
#include "hamiltonian_device.hpp"
#include "hamiltonian_bitset_representation.hpp"
#include "diagnostics.hpp"
#include "macros.hpp"
#include "../external/eigen-3.4.0/Eigen/Dense"
#include "../external/eigen-3.4.0/Eigen/Eigenvalues"
// #include "../external/boost_1_84_0/boost/container_hash/hash.hpp"

using std::cout;
using std::endl;
using std::string;
using std::vector;

int main(int argc, char* argv[])
{
    diagnostics::print_dtype_sizes();
    const uint16_t n_valence_protons = std::stoi(argv[1]);
    const uint16_t n_valence_neutrons = std::stoi(argv[2]);
    const string interaction_filename = "./resources/interaction_files/w.snt";
    // const string interaction_filename = "./resources/interaction_files/usda.snt";
    // const string interaction_filename = "./resources/interaction_files/gxpf1a.snt";
    
    string ref_matrix_filename = "p";
    ref_matrix_filename.append(argv[1]);
    ref_matrix_filename.append("n");
    ref_matrix_filename.append(argv[2]);

    const Interaction interaction = loaders::load_interaction(interaction_filename, n_valence_protons, n_valence_neutrons);
    const Indices indices = indices::generate_indices(interaction);
    const size_t m_dim = interaction.basis_states.size();
    diagnostics::print_hamiltonian_info(interaction);
    diagnostics::print_gpu_diagnostics(interaction, indices);

    // double* H_host_reference = new double[m_dim*m_dim];
    // for (size_t idx = 0; idx < m_dim*m_dim; idx++) H_host_reference[idx] = 0; // Not needed for the calculation, using this so that printing H is less cluttered.

    double* H_from_device = new double[m_dim*m_dim];
    for (size_t idx = 0; idx < m_dim*m_dim; idx++) H_from_device[idx] = 0; // Not needed for the calculation, using this so that printing H is less cluttered.

    std::vector<int64_t> timing;
    for (size_t i = 0; i < 3; i++)
    {
        auto start = timer();
        // hamiltonian::create_hamiltonian_primitive_bit_representation_reference(interaction, indices, H_host_reference); cout << endl;
        // hamiltonian::create_hamiltonian_primitive_bit_representation_new(interaction, indices, H_host);
        hamiltonian_device::create_hamiltonian_device_dispatcher(interaction, indices, H_from_device);
        timing.push_back(timer(start));
    }
    cout << endl;
    print_vector(timing);
    cout << "Average over " << timing.size() << " runs: " << mean(timing) << " milliseconds" << endl;

    // cout << "compare_arrays(H_host_reference, H_from_device, m_dim*m_dim): " << compare_arrays_upper_triangle(H_host_reference, H_from_device, m_dim) << endl;
    // cout << "H_host_reference == reference: " << compare_with_ref_matrix(H_host_reference, ref_matrix_filename, m_dim) << endl;
    cout << "H_from_device == reference: " << compare_with_ref_matrix(H_from_device, ref_matrix_filename, m_dim) << endl;

    // print_flattened_2d_array(H_host_reference, m_dim);
    // print_flattened_2d_array(H_from_device, m_dim);
    // print_flattened_2d_array(H_host_reference, m_dim);
    // Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> H_host_mapped(H_host, m_dim, m_dim);
    // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(H_host_mapped.transpose());
    // std::cout << "Eigenvalues:\n" << solver.eigenvalues() << std::endl;
    
    // delete[] H_host_reference;
    delete[] H_from_device;
    return 0;
}