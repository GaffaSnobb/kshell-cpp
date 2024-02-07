#include <vector>
#include <iostream>
#include <hip/hip_runtime.h>
#include "hip_wrappers.hpp"
#include "bit_manipulation_tools_device.hpp"
#include "data_structures.hpp"
#include "generate_indices.hpp"
#include "basis.hpp"
#include "tools.hpp"

using std::cout;
using std::endl;

namespace hamiltonian_device
{
__device__ double calculate_onebody_matrix_element(
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

__global__ void onebody_matrix_element_dispatcher(
    double* H_diag,
    const unsigned long long* basis_states,
    const double* spe,
    const unsigned short* orbital_idx_to_composite_m_idx_map_flattened_indices,
    const unsigned int m_dim,
    const unsigned short n_orbitals
)
{
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    // int row_idx = idx/m_dim;
    // int col_idx = idx%m_dim;

    // if ((row_idx < m_dim) and (col_idx < m_dim) and (row_idx == col_idx))
    if (idx < m_dim)
    {
        const unsigned long long left_state = basis_states[idx];
        const unsigned long long right_state = basis_states[idx];
        
        H_diag[idx] = calculate_onebody_matrix_element(
            n_orbitals,
            spe,
            orbital_idx_to_composite_m_idx_map_flattened_indices,
            left_state,
            right_state
        );
    }
}

__global__ void twobody_matrix_element_dispatcher(
    double* H_diag,
    const unsigned long long* basis_states,
    const double* spe,
    const unsigned short* orbital_idx_to_composite_m_idx_map_flattened_indices,
    const unsigned int m_dim,
    const unsigned short n_orbitals
)
{
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    // int row_idx = idx/m_dim;
    // int col_idx = idx%m_dim;

    // if ((row_idx < m_dim) and (col_idx < m_dim) and (row_idx == col_idx))
    if (idx < m_dim)
    {
        const unsigned long long left_state = basis_states[idx];
        const unsigned long long right_state = basis_states[idx];
        
        H_diag[idx] = calculate_onebody_matrix_element(
            n_orbitals,
            spe,
            orbital_idx_to_composite_m_idx_map_flattened_indices,
            left_state,
            right_state
        );
    }
}

void create_hamiltonian_device_dispatcher(const Interaction& interaction, const Indices& indices, double* H)
{
    const unsigned int m_dim = interaction.basis_states.size();
    const unsigned short n_orbitals = interaction.model_space.n_orbitals;

    const int threads_per_block = 128;
    const int blocks_per_grid = (m_dim + threads_per_block - 1) / threads_per_block;
    double* H_diag_device = nullptr;
    unsigned long long* basis_states_device = nullptr;
    double* spe_array_device = nullptr;
    unsigned short* orbital_idx_to_composite_m_idx_map_flattened_indices_device = nullptr;
    
    auto start = timer();
        double* H_diag_tmp = new double[m_dim];
        hip_wrappers::hipMalloc(&H_diag_device, m_dim*sizeof(double));
        hip_wrappers::hipMalloc(&basis_states_device, m_dim*sizeof(unsigned long long));
        hip_wrappers::hipMalloc(&spe_array_device, n_orbitals*sizeof(double));
        hip_wrappers::hipMalloc(&orbital_idx_to_composite_m_idx_map_flattened_indices_device, n_orbitals*sizeof(unsigned short));
        
        hip_wrappers::hipMemcpy(basis_states_device, interaction.basis_states.data(), m_dim*sizeof(unsigned long long), hipMemcpyHostToDevice);
        hip_wrappers::hipMemcpy(spe_array_device, interaction.spe_array, n_orbitals*sizeof(double), hipMemcpyHostToDevice);
        hip_wrappers::hipMemcpy(orbital_idx_to_composite_m_idx_map_flattened_indices_device, indices.orbital_idx_to_composite_m_idx_map_flattened_indices, n_orbitals*sizeof(unsigned short), hipMemcpyHostToDevice);

        hipLaunchKernelGGL(
            onebody_matrix_element_dispatcher, dim3(blocks_per_grid), dim3(threads_per_block), 0, 0,
            H_diag_device,
            basis_states_device,
            spe_array_device,
            orbital_idx_to_composite_m_idx_map_flattened_indices_device,
            m_dim,
            n_orbitals
        );
        hip_wrappers::hipDeviceSynchronize();

        hip_wrappers::hipMemcpy(H_diag_tmp, H_diag_device, m_dim*sizeof(double), hipMemcpyDeviceToHost);
        hip_wrappers::hipFree(H_diag_device);
        hip_wrappers::hipFree(basis_states_device);
        hip_wrappers::hipFree(spe_array_device);
        hip_wrappers::hipFree(orbital_idx_to_composite_m_idx_map_flattened_indices_device);
    timer(start, "[DEVICE] one-body calc, alloc, copy, and free time");

    for (int diag_idx = 0; diag_idx < m_dim; diag_idx++) H[diag_idx*m_dim + diag_idx] = H_diag_tmp[diag_idx]; // Copy data to the diagonal of H.
    delete[] H_diag_tmp;
}
}