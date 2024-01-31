#include <vector>
#include <iostream>
#include "bit_manipulation_tools_device.hpp"
#include "data_structures.hpp"
#include "generate_indices.hpp"
#include "basis.hpp"
#include "tools.hpp"

using std::cout;
using std::endl;

namespace hamiltonian_device
{
__device__ double calculate_onebody_matrix_element_primitive_bit_representation_device(
    const unsigned short n_orbitals,
    const double* spe,
    const unsigned short* orbital_idx_to_composite_m_idx_map_flattened_indices,
    const unsigned long long& left_state,
    const unsigned long long& right_state
)
{
    double onebody_res = 0;
    unsigned short creation_start_m_idx = 0;
    unsigned short annihilation_start_m_idx = 0;
    
    for (unsigned short creation_and_annihilation_orb_idx = 0; creation_and_annihilation_orb_idx < n_orbitals; creation_and_annihilation_orb_idx++)
    {
        const unsigned short creation_end_m_idx = orbital_idx_to_composite_m_idx_map_flattened_indices[creation_and_annihilation_orb_idx];
        const unsigned short annihilation_end_m_idx = orbital_idx_to_composite_m_idx_map_flattened_indices[creation_and_annihilation_orb_idx];

        for (unsigned short creation_comp_m_idx = creation_start_m_idx; creation_comp_m_idx < creation_end_m_idx; creation_comp_m_idx++)
        {
            for (unsigned short annihilation_comp_m_idx = annihilation_start_m_idx; annihilation_comp_m_idx < annihilation_end_m_idx; annihilation_comp_m_idx++)
            {
                // double switch_ = 1; // To eliminate if-statements.
                unsigned long long new_right_state = right_state;   // The contents of right_state is copied, not referenced.

                // switch_ = switch_*bittools_device::is_bit_set(new_right_state, annihilation_comp_m_idx);
                if (not bittools_device::is_bit_set(new_right_state, annihilation_comp_m_idx)) continue;
                const unsigned short n_operator_swaps_annihilation = bittools_device::reset_bit_and_count_swaps(new_right_state, annihilation_comp_m_idx);
                const short annihilation_sign = bittools_device::negative_one_pow(n_operator_swaps_annihilation);

                // switch_ = switch_*(not bittools_device::is_bit_set(new_right_state, creation_comp_m_idx));
                if (bittools_device::is_bit_set(new_right_state, creation_comp_m_idx)) continue;
                const unsigned short n_operator_swaps_creation = bittools_device::set_bit_and_count_swaps(new_right_state, creation_comp_m_idx);
                const short creation_sign = bittools_device::negative_one_pow(n_operator_swaps_creation);

                // switch_ = switch_*(left_state == new_right_state);
                if (left_state != new_right_state) continue;
                onebody_res += annihilation_sign*creation_sign*spe[creation_and_annihilation_orb_idx];//*switch_;   // Or annihilation_orb_idx, they are the same.
            }
        }
        annihilation_start_m_idx = annihilation_end_m_idx; // Update annihilation_start_m_idx to the beginning of the next section of the map.
        creation_start_m_idx = creation_end_m_idx; // Update creation_start_m_idx to the beginning of the next section of the map.
    }
    return onebody_res;
}

__global__ void matrix_element_dispatcher(
    double* H,
    const unsigned long long* basis_states,
    const double* spe,
    const unsigned short* orbital_idx_to_composite_m_idx_map_flattened_indices,
    const unsigned int m_dim,
    const unsigned short n_orbitals
)
{
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    int row_idx = idx/m_dim;
    int col_idx = idx%m_dim;

    if ((row_idx < m_dim) and (col_idx < m_dim))
    {
        const unsigned long long left_state = basis_states[row_idx];
        const unsigned long long right_state = basis_states[col_idx];

        H[idx] = calculate_onebody_matrix_element_primitive_bit_representation_device(
            n_orbitals,
            spe,
            orbital_idx_to_composite_m_idx_map_flattened_indices,
            left_state,
            right_state
        );
    }
}

void create_hamiltonian_device_dispatcher(const Interaction& interaction)
{
    const Indices indices = indices::generate_indices(interaction);
    const std::vector<unsigned long long> basis_states = basis::calculate_m_basis_states_primitive_bit_representation(interaction);
    const unsigned int m_dim = basis_states.size();
    cout << "---------------------" << endl;
    print("n_valence_protons", interaction.model_space_protons.n_valence_nucleons);
    print("n_valence_neutrons", interaction.model_space_neutrons.n_valence_nucleons);
    print("m_dim", m_dim);
    print("m_dim**2", m_dim*m_dim);
    print("H size (MB): ", m_dim*m_dim*sizeof(double)/1000./1000.);
    cout << "---------------------" << endl;

    unsigned long long* basis_states_array = new unsigned long long[m_dim];
    double* H = new double[m_dim*m_dim];
    const unsigned int n_cols = m_dim;  // Yes they are the same, I just want to explicitly see what the flattened 2D -> 1D indexing looks like.
    // const unsigned int n_rows = m_dim;

    for (int row_idx = 0; row_idx < m_dim; row_idx++)
    {
        basis_states_array[row_idx] = basis_states[row_idx];
        // for (int col_idx = 0; col_idx < m_dim; col_idx++)
        // {
        //     H[row_idx*n_cols + col_idx] = 0;    // I dont think this is needed.
        // }
    }
    const int threads_per_block = 128;
    const int blocks_per_grid = (m_dim*m_dim + threads_per_block - 1) / threads_per_block;

    double* H_device = nullptr;
    unsigned long long* basis_states_device = nullptr;
    double* spe_array_device = nullptr;
    unsigned short* orbital_idx_to_composite_m_idx_map_flattened_indices_device = nullptr;
    hipMalloc(&H_device, m_dim*m_dim*sizeof(double));
    hipMalloc(&basis_states_device, m_dim*sizeof(unsigned long long));
    hipMalloc(&spe_array_device, interaction.model_space.n_orbitals*sizeof(double));
    hipMalloc(&orbital_idx_to_composite_m_idx_map_flattened_indices_device, interaction.model_space.n_orbitals*sizeof(unsigned short));
    
    // hipMemcpy(H_device, H, m_dim*m_dim*sizeof(double), hipMemcpyHostToDevice);
    hipMemcpy(basis_states_device, basis_states_array, m_dim*sizeof(unsigned long long), hipMemcpyHostToDevice);
    hipMemcpy(spe_array_device, interaction.spe_array, interaction.model_space.n_orbitals*sizeof(double), hipMemcpyHostToDevice);
    hipMemcpy(orbital_idx_to_composite_m_idx_map_flattened_indices_device, indices.orbital_idx_to_composite_m_idx_map_flattened_indices, interaction.model_space.n_orbitals*sizeof(unsigned short), hipMemcpyHostToDevice);

    hipLaunchKernelGGL(
        matrix_element_dispatcher, dim3(blocks_per_grid), dim3(threads_per_block), 0, 0,
        H_device,
        basis_states_device,
        spe_array_device,
        orbital_idx_to_composite_m_idx_map_flattened_indices_device,
        m_dim,
        interaction.model_space.n_orbitals
    );
    hipDeviceSynchronize();
    cout << "\n" << endl;

    hipMemcpy(H, H_device, m_dim*m_dim*sizeof(double), hipMemcpyDeviceToHost);
    hipFree(H_device);
    hipFree(basis_states_device);
    hipFree(spe_array_device);
    hipFree(orbital_idx_to_composite_m_idx_map_flattened_indices_device);

    // for (int row_idx = 0; row_idx < m_dim; row_idx++)
    // {
    //     for (int col_idx = 0; col_idx < m_dim; col_idx++)
    //     {
    //         cout << H[row_idx*n_cols + col_idx] << ", ";
    //     }
    //     cout << endl;
    // }
    // cout << endl;

    
    delete[] basis_states_array;
    delete[] H;
}
}