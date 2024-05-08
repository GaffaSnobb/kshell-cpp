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

// namespace device_arrays
// {
// __constant__ uint16_t orbital_idx_to_composite_m_idx_map_flattened_indices_const[CONST_MEM_ARR_LEN_N_ORBITALS];
// __constant__ uint16_t creation_orb_indices_0_const[CONST_MEM_ARR_LEN_INDICES];
// __constant__ uint16_t creation_orb_indices_1_const[CONST_MEM_ARR_LEN_INDICES];
// __constant__ uint16_t annihilation_orb_indices_0_const[CONST_MEM_ARR_LEN_INDICES];
// __constant__ uint16_t annihilation_orb_indices_1_const[CONST_MEM_ARR_LEN_INDICES];
// __constant__ uint16_t j_coupled_const[CONST_MEM_ARR_LEN_INDICES];
// __constant__ int16_t m_coupled_const[CONST_MEM_ARR_LEN_INDICES];
// __constant__ double tbme_const[CONST_MEM_ARR_LEN_INDICES];
// } // namespace device_arrays

// __global__ void tmpkern()
// {
//     for (size_t i = 0; i < 6; i++)
//     {
//         printf("%d, ", device_arrays::orbital_idx_to_composite_m_idx_map_flattened_indices_const[i]);
//     }
// }

// void gpu_init(const Interaction& interaction, const Indices& indices)
// {
//     /*
//     Add shit to constant memory.
//     */
//     hipDeviceProp_t prop;
//     const size_t device_id = 0;
//     HIP_ASSERT(hipGetDeviceProperties(&prop, device_id));

//     const size_t coi_0 = indices.creation_orb_indices_0.size()*sizeof(uint16_t);
//     const size_t coi_1 = indices.creation_orb_indices_1.size()*sizeof(uint16_t);
//     const size_t aoi_0 = indices.annihilation_orb_indices_0.size()*sizeof(uint16_t);
//     const size_t aoi_1 = indices.annihilation_orb_indices_1.size()*sizeof(uint16_t);
//     const size_t jc = indices.j_coupled.size()*sizeof(uint16_t);
//     const size_t mc = indices.m_coupled.size()*sizeof(int16_t);
//     const size_t tbme = indices.tbme.size()*sizeof(double);
//     const size_t oitcmimf = interaction.model_space.orbitals.size()*sizeof(uint16_t);
//     const size_t total = coi_0 + coi_1 + aoi_0 + aoi_1 + jc + mc + tbme + oitcmimf;
    
//     hip_wrappers::hipMemcpyToSymbol(device_arrays::creation_orb_indices_0_const, indices.creation_orb_indices_0);
//     hip_wrappers::hipMemcpyToSymbol(device_arrays::creation_orb_indices_1_const, indices.creation_orb_indices_1);
//     hip_wrappers::hipMemcpyToSymbol(device_arrays::annihilation_orb_indices_0_const, indices.annihilation_orb_indices_0);
//     hip_wrappers::hipMemcpyToSymbol(device_arrays::annihilation_orb_indices_1_const, indices.annihilation_orb_indices_1);
//     hip_wrappers::hipMemcpyToSymbol(device_arrays::j_coupled_const, indices.j_coupled);
//     hip_wrappers::hipMemcpyToSymbol(device_arrays::m_coupled_const, indices.m_coupled);
//     hip_wrappers::hipMemcpyToSymbol(device_arrays::orbital_idx_to_composite_m_idx_map_flattened_indices_const, indices.orbital_idx_to_composite_m_idx_map_flattened_indices, oitcmimf);

//     for (size_t i = 0; i < interaction.model_space.orbitals.size(); i++)
//     {
//         printf("%d, ", indices.orbital_idx_to_composite_m_idx_map_flattened_indices[i]);
//     }
//     printf("\n");
//     // tmpkern<<<1, 1>>>();
    
//     cout << diagnostics::DIAG_STR_START << endl;
//     cout << coi_0/1e3 << " kB" << " (" << indices.creation_orb_indices_0.size() << " elements)" << endl;
//     cout << coi_1/1e3 << " kB" << " (" << indices.creation_orb_indices_1.size() << " elements)" << endl;
//     cout << aoi_0/1e3 << " kB" << " (" << indices.annihilation_orb_indices_0.size() << " elements)" << endl;
//     cout << aoi_1/1e3 << " kB" << " (" << indices.annihilation_orb_indices_1.size() << " elements)" << endl;
//     cout << jc/1e3 << " kB" << " (" << indices.j_coupled.size() << " elements)" << endl;
//     cout << mc/1e3 << " kB" << " (" << indices.m_coupled.size() << " elements)" << endl;
//     cout << tbme/1e3 << " kB" << " (" << indices.tbme.size() << " elements)" << endl;
//     cout << oitcmimf/1e3 << " kB" << " (" << interaction.model_space.orbitals.size() << " elements)" << endl;
//     cout << total/1e6 << " MB __constant__ mem used of " << prop.totalConstMem/1e6 << " MB total (" << total*100.0/prop.totalConstMem << "%)" << endl;
//     cout << diagnostics::DIAG_STR_END << endl;

//     // exit(0);

//     diagnostics::print_gpu_diagnostics(interaction, indices);
// }

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

    double* H_host_reference = new double[m_dim*m_dim];
    for (size_t idx = 0; idx < m_dim*m_dim; idx++) H_host_reference[idx] = 0; // Not needed for the calculation, using this so that printing H is less cluttered.

    double* H_host = new double[m_dim*m_dim];
    for (size_t idx = 0; idx < m_dim*m_dim; idx++) H_host[idx] = 0; // Not needed for the calculation, using this so that printing H is less cluttered.

    std::vector<int64_t> timing;
    for (size_t i = 0; i < 1; i++)
    {
        auto start = timer();
        hamiltonian::create_hamiltonian_primitive_bit_representation_reference(interaction, indices, H_host_reference); cout << endl;
        // hamiltonian::create_hamiltonian_primitive_bit_representation_new(interaction, indices, H_host);
        hamiltonian_device::create_hamiltonian_device_dispatcher(interaction, indices, H_host);
        timing.push_back(timer(start));
    }
    cout << endl;
    print_vector(timing);
    cout << mean(timing) << endl;

    cout << "compare_arrays(H_host_reference, H_host, m_dim*m_dim): " << compare_arrays(H_host_reference, H_host, m_dim*m_dim) << endl;
    // cout << "H_host_reference == reference: " << compare_with_ref_matrix(H_host_reference, ref_matrix_filename, m_dim) << endl;
    // cout << "H_host == reference: " << compare_with_ref_matrix(H_host, ref_matrix_filename, m_dim) << endl;

    // print_flattened_2d_array(H_host, m_dim);
    // print_flattened_2d_array(H_host_reference, m_dim);
    // Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> H_host_mapped(H_host, m_dim, m_dim);
    // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(H_host_mapped.transpose());
    // std::cout << "Eigenvalues:\n" << solver.eigenvalues() << std::endl;
    
    delete[] H_host_reference;
    delete[] H_host;
    return 0;
}