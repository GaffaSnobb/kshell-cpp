#include <vector>
#include <iostream>
#include <stdint.h>
#include <hip/hip_runtime.h>

#include "hip_wrappers.hpp"
#include "bit_manipulation_tools_device.hpp"
#include "data_structures.hpp"
#include "generate_indices.hpp"
#include "basis.hpp"
#include "tools.hpp"
#include "parameters.hpp"
#include "macros.hpp"
#include "diagnostics.hpp"

using std::cout;
using std::endl;

/*
I cannot for the life of me figure out how to put these device constant
memory arrays in a header so that I can use them in several translation
units. `extern` just makes the entire compilation crash, and if I
declare the arrays in a header, values I put into them do not persist
across translation units.
*/
__constant__ uint16_t dev_const_creation_orb_indices_0[CONST_MEM_ARR_LEN_INDICES];
__constant__ uint16_t dev_const_creation_orb_indices_1[CONST_MEM_ARR_LEN_INDICES];
__constant__ uint16_t dev_const_annihilation_orb_indices_0[CONST_MEM_ARR_LEN_INDICES];
__constant__ uint16_t dev_const_annihilation_orb_indices_1[CONST_MEM_ARR_LEN_INDICES];
__constant__ uint16_t dev_const_j_coupled[CONST_MEM_ARR_LEN_INDICES];
__constant__ int16_t dev_const_m_coupled[CONST_MEM_ARR_LEN_INDICES];
__constant__ uint16_t dev_const_orbital_idx_to_composite_m_idx_map_flattened_indices[CONST_MEM_ARR_LEN_N_ORBITALS];
__constant__ double dev_const_tbme[CONST_MEM_ARR_LEN_INDICES];
__constant__ double dev_const_spe[CONST_MEM_ARR_LEN_N_ORBITALS];

__host__ void dev_init(
    const Interaction& interaction,
    const Indices& indices,
    uint64_t*& dev_basis_states
)
{
    /*
    Add shit to constant memory.
    */
    hipDeviceProp_t prop;
    const size_t device_id = 0;
    HIP_ASSERT(hipGetDeviceProperties(&prop, device_id));

    const size_t m_dim = interaction.basis_states.size();

    /*
    Byte sizes. These are the amounts of bytes which have actual data in
    them. The constant arrays are likely reserving more space than this,
    as it is difficult to change the size of the constant arrays when
    the size has to be known at compile-time. These sizes should however
    never be larger than `CONST_MEM_ARR_LEN_INDICES` and
    `CONST_MEM_ARR_LEN_N_ORBITALS`.
    */
    const size_t coi_0 = indices.creation_orb_indices_0.size()*sizeof(uint16_t);
    const size_t coi_1 = indices.creation_orb_indices_1.size()*sizeof(uint16_t);
    const size_t aoi_0 = indices.annihilation_orb_indices_0.size()*sizeof(uint16_t);
    const size_t aoi_1 = indices.annihilation_orb_indices_1.size()*sizeof(uint16_t);
    const size_t jc = indices.j_coupled.size()*sizeof(uint16_t);
    const size_t mc = indices.m_coupled.size()*sizeof(int16_t);
    const size_t oitcmimf = interaction.model_space.orbitals.size()*sizeof(uint16_t);
    const size_t tbme = indices.tbme.size()*sizeof(double);
    const size_t spe = interaction.spe.size()*sizeof(double);
    const size_t total = coi_0 + coi_1 + aoi_0 + aoi_1 + jc + mc + tbme + oitcmimf + spe;

    assert(coi_0 <= (CONST_MEM_ARR_LEN_INDICES*sizeof(uint16_t)));
    assert(coi_1 <= (CONST_MEM_ARR_LEN_INDICES*sizeof(uint16_t)));
    assert(aoi_0 <= (CONST_MEM_ARR_LEN_INDICES*sizeof(uint16_t)));
    assert(aoi_1 <= (CONST_MEM_ARR_LEN_INDICES*sizeof(uint16_t)));
    assert(jc <= (CONST_MEM_ARR_LEN_INDICES*sizeof(uint16_t)));
    assert(mc <= (CONST_MEM_ARR_LEN_INDICES*sizeof(int16_t)));
    assert(oitcmimf <= (CONST_MEM_ARR_LEN_N_ORBITALS*sizeof(uint16_t)));
    assert(tbme <= (CONST_MEM_ARR_LEN_INDICES*sizeof(double)));
    assert(spe <= (CONST_MEM_ARR_LEN_N_ORBITALS*sizeof(double)));
    
    // Constant device arrays.
    hip_wrappers::hipMemcpyToSymbol(dev_const_creation_orb_indices_0, indices.creation_orb_indices_0);
    hip_wrappers::hipMemcpyToSymbol(dev_const_creation_orb_indices_1, indices.creation_orb_indices_1);
    hip_wrappers::hipMemcpyToSymbol(dev_const_annihilation_orb_indices_0, indices.annihilation_orb_indices_0);
    hip_wrappers::hipMemcpyToSymbol(dev_const_annihilation_orb_indices_1, indices.annihilation_orb_indices_1);
    hip_wrappers::hipMemcpyToSymbol(dev_const_j_coupled, indices.j_coupled);
    hip_wrappers::hipMemcpyToSymbol(dev_const_m_coupled, indices.m_coupled);
    hip_wrappers::hipMemcpyToSymbol(dev_const_orbital_idx_to_composite_m_idx_map_flattened_indices, indices.orbital_idx_to_composite_m_idx_map_flattened_indices, oitcmimf);
    hip_wrappers::hipMemcpyToSymbol(dev_const_tbme, indices.tbme);
    hip_wrappers::hipMemcpyToSymbol(dev_const_spe, interaction.spe);
    
    cout << diagnostics::DIAG_STR_START << endl;
    cout << coi_0/1e3 << " kB" << " (" << indices.creation_orb_indices_0.size() << " elements)" << endl;
    cout << coi_1/1e3 << " kB" << " (" << indices.creation_orb_indices_1.size() << " elements)" << endl;
    cout << aoi_0/1e3 << " kB" << " (" << indices.annihilation_orb_indices_0.size() << " elements)" << endl;
    cout << aoi_1/1e3 << " kB" << " (" << indices.annihilation_orb_indices_1.size() << " elements)" << endl;
    cout << jc/1e3 << " kB" << " (" << indices.j_coupled.size() << " elements)" << endl;
    cout << mc/1e3 << " kB" << " (" << indices.m_coupled.size() << " elements)" << endl;
    cout << oitcmimf/1e3 << " kB" << " (" << interaction.model_space.orbitals.size() << " elements)" << endl;
    cout << tbme/1e3 << " kB" << " (" << indices.tbme.size() << " elements)" << endl;
    cout << spe/1e3 << " kB" << " (" << interaction.spe.size() << " elements)" << endl;
    cout << total/1e6 << " MB __constant__ mem used of " << prop.totalConstMem/1e6 << " MB total (" << total*100.0/prop.totalConstMem << "%)" << endl;
    cout << diagnostics::DIAG_STR_END << endl;

    diagnostics::print_gpu_diagnostics(interaction, indices);
    
    // Device arrays.
    hip_wrappers::hipMalloc(&dev_basis_states, m_dim*sizeof(uint64_t));
    hip_wrappers::hipMemcpy(dev_basis_states, interaction.basis_states.data(), m_dim*sizeof(uint64_t), hipMemcpyHostToDevice);
}

__host__ void dev_uninit(uint64_t*& dev_basis_states)
{
    hip_wrappers::hipFree(dev_basis_states);
}

namespace hamiltonian_device
{
__device__ double calculate_onebody_matrix_element(
    const size_t n_orbitals,
    const uint64_t& left_state,
    const uint64_t& right_state
)
{
    double onebody_res = 0;
    uint16_t creation_start_m_idx = 0;
    uint16_t annihilation_start_m_idx = 0;
    
    for (uint16_t creation_and_annihilation_orb_idx = 0; creation_and_annihilation_orb_idx < n_orbitals; creation_and_annihilation_orb_idx++)
    {
        const uint16_t creation_end_m_idx =
            dev_const_orbital_idx_to_composite_m_idx_map_flattened_indices[creation_and_annihilation_orb_idx];
        const uint16_t annihilation_end_m_idx = creation_end_m_idx;

        for (uint16_t creation_comp_m_idx = creation_start_m_idx; creation_comp_m_idx < creation_end_m_idx; creation_comp_m_idx++)
        {
            for (uint16_t annihilation_comp_m_idx = annihilation_start_m_idx; annihilation_comp_m_idx < annihilation_end_m_idx; annihilation_comp_m_idx++)
            {
                // double switch_ = 1; // To eliminate if-statements.
                uint64_t new_right_state = right_state;   // The contents of right_state is copied, not referenced.

                // switch_ = switch_*bittools_device::is_bit_set(new_right_state, annihilation_comp_m_idx);
                if (not bittools_device::is_bit_set(new_right_state, annihilation_comp_m_idx)) continue;
                const uint16_t n_operator_swaps_annihilation = bittools_device::reset_bit_and_count_swaps(new_right_state, annihilation_comp_m_idx);
                const int16_t annihilation_sign = bittools_device::negative_one_pow(n_operator_swaps_annihilation);

                // switch_ = switch_*(not bittools_device::is_bit_set(new_right_state, creation_comp_m_idx));
                if (bittools_device::is_bit_set(new_right_state, creation_comp_m_idx)) continue;
                const uint16_t n_operator_swaps_creation = bittools_device::set_bit_and_count_swaps(new_right_state, creation_comp_m_idx);
                const int16_t creation_sign = bittools_device::negative_one_pow(n_operator_swaps_creation);

                // switch_ = switch_*(left_state == new_right_state);
                if (left_state != new_right_state) continue;
                onebody_res += annihilation_sign*creation_sign*dev_const_spe[creation_and_annihilation_orb_idx];//*switch_;   // Or annihilation_orb_idx, they are the same.
            }
        }
        annihilation_start_m_idx = annihilation_end_m_idx; // Update annihilation_start_m_idx to the beginning of the next section of the map.
        creation_start_m_idx = creation_end_m_idx; // Update creation_start_m_idx to the beginning of the next section of the map.
    }
    return onebody_res;
}

__global__ void onebody_matrix_element_dispatcher(
    double* H_diag,
    const uint64_t* dev_basis_states,
    const uint32_t m_dim,
    const size_t n_orbitals
)
{
    const size_t idx = blockIdx.x*blockDim.x + threadIdx.x;

    if (idx < m_dim)
    {
        const uint64_t left_state = dev_basis_states[idx];
        const uint64_t right_state = dev_basis_states[idx];
        
        H_diag[idx] = calculate_onebody_matrix_element(
            n_orbitals,
            left_state,
            right_state
        );
    }
}

__device__ double calculate_twobody_matrix_element(
    const size_t n_orbitals,
    const uint64_t& left_state,
    const uint64_t& right_state
)
{
    uint16_t j1_idx, m1_idx, j2_idx, m2_idx, j3_idx;  // Indices for the Clebsch-Gordan 5D array...
    uint16_t flat_cg_idx;                             // ...translated to a flattened 1D array.

    double twobody_res = 0;
    // double creation_res_switch = 0; // For removing if statements in the inner creation loop. Multiply the result by 0 instead of using if (condition) continue;.
    const uint32_t n_indices = indices.creation_orb_indices_0.size();
    for (size_t i = 0; i < n_indices; i++)
    {
        const uint16_t creation_orb_idx_0 = indices.creation_orb_indices_0[i];
        const uint16_t creation_orb_idx_1 = indices.creation_orb_indices_1[i];
        const uint16_t annihilation_orb_idx_0 = indices.annihilation_orb_indices_0[i];
        const uint16_t annihilation_orb_idx_1 = indices.annihilation_orb_indices_1[i];
        const uint16_t j_coupled = indices.j_coupled[i];
        const int16_t m_coupled = indices.m_coupled[i];
        const double tbme = indices.tbme[i];

        const double creation_norm = inverse_sqrt_2[creation_orb_idx_0 == creation_orb_idx_1];
        const double annihilation_norm = inverse_sqrt_2[annihilation_orb_idx_0 == annihilation_orb_idx_1];

        // Annihilation terms
        for (uint16_t annihilation_comp_m_idx_0 : indices.orbital_idx_to_composite_m_idx_map[annihilation_orb_idx_0])
        {
            for (uint16_t annihilation_comp_m_idx_1 : indices.orbital_idx_to_composite_m_idx_map[annihilation_orb_idx_1])
            {
                if (
                    (indices.composite_m_idx_to_m_map[annihilation_comp_m_idx_0] +  // m1
                    indices.composite_m_idx_to_m_map[annihilation_comp_m_idx_1]) != // m2
                    m_coupled                                                       // M
                )
                {
                    /*
                    The Clebsch-Gordan coefficient is always zero when
                    m1 + m2 != M. This is because of conservation of the
                    z-component of angular momentum in the coupling
                    process.
                    */
                    continue;
                }
                
                // if (not right_state.test(annihilation_comp_m_idx_0))
                if (not bittools_device::is_bit_set(right_state, annihilation_comp_m_idx_0)) continue;

                uint64_t new_right_state_annihilation = right_state;
                const uint16_t n_operator_swaps_annihilation_0 = bittools_device::reset_bit_and_count_swaps(new_right_state_annihilation, annihilation_comp_m_idx_0);
                int16_t annihilation_sign = bittools_device::negative_one_pow(n_operator_swaps_annihilation_0);

                // if (not new_right_state_annihilation.test(annihilation_comp_m_idx_1)) continue;
                if (not bittools_device::is_bit_set(new_right_state_annihilation, annihilation_comp_m_idx_1)) continue;
                const uint16_t n_operator_swaps_annihilation_1 = bittools_device::reset_bit_and_count_swaps(new_right_state_annihilation, annihilation_comp_m_idx_1);
                annihilation_sign *= bittools_device::negative_one_pow(n_operator_swaps_annihilation_1);
                
                j1_idx = (indices.orbital_idx_to_j_map[annihilation_orb_idx_0] + 1)/2 - 1;
                m1_idx = (indices.composite_m_idx_to_m_map[annihilation_comp_m_idx_0] + 5)/2;
                j2_idx = (indices.orbital_idx_to_j_map[annihilation_orb_idx_1] + 1)/2 - 1;
                m2_idx = (indices.composite_m_idx_to_m_map[annihilation_comp_m_idx_1] + 5)/2;
                j3_idx = j_coupled/2;
                flat_cg_idx = ((((j1_idx*6 + m1_idx)*3 + j2_idx)*6 + m2_idx)*6 + j3_idx);    // Indexing a 1D array as if it was a 5D array.

                const double cg_annihilation = clebsch_gordan_array[flat_cg_idx];

                // Creation terms
                double creation_res = 0.0;

                for (uint16_t creation_comp_m_idx_0 : indices.orbital_idx_to_composite_m_idx_map[creation_orb_idx_0])
                {
                    for (uint16_t creation_comp_m_idx_1 : indices.orbital_idx_to_composite_m_idx_map[creation_orb_idx_1])
                    {   
                        if (
                            (indices.composite_m_idx_to_m_map[creation_comp_m_idx_0] +  // m1
                            indices.composite_m_idx_to_m_map[creation_comp_m_idx_1]) != // m2
                            m_coupled                                                   // M
                        )
                        {
                            /*
                            The Clebsch-Gordan coefficient is always zero when
                            m1 + m2 != M. This is because of conservation of the
                            z-component of angular momentum in the coupling
                            process.
                            */
                            continue;
                        }
                        if (bittools_device::is_bit_set(new_right_state_annihilation, creation_comp_m_idx_1)) continue;
                        // creation_res_switch = (not new_right_state_annihilation[creation_comp_m_idx_1]);
                        
                        uint64_t new_right_state_creation = new_right_state_annihilation;
                        const uint16_t n_operator_swaps_creation_1 = bittools_device::set_bit_and_count_swaps(new_right_state_creation, creation_comp_m_idx_1);
                        int16_t creation_sign = bittools_device::negative_one_pow(n_operator_swaps_creation_1);

                        if (bittools_device::is_bit_set(new_right_state_annihilation, creation_comp_m_idx_0)) continue;
                        // creation_res_switch *= (not new_right_state_annihilation[creation_comp_m_idx_0]);
                        
                        const uint16_t n_operator_swaps_creation_0 = bittools_device::set_bit_and_count_swaps(new_right_state_creation, creation_comp_m_idx_0);
                        creation_sign *= bittools_device::negative_one_pow(n_operator_swaps_creation_0);

                        if (left_state != new_right_state_creation) continue;
                        // creation_res_switch *= (left_state == new_right_state_creation);

                        j1_idx = (indices.orbital_idx_to_j_map[creation_orb_idx_0] + 1)/2 - 1;
                        m1_idx = (indices.composite_m_idx_to_m_map[creation_comp_m_idx_0] + 5)/2;
                        j2_idx = (indices.orbital_idx_to_j_map[creation_orb_idx_1] + 1)/2 - 1;
                        m2_idx = (indices.composite_m_idx_to_m_map[creation_comp_m_idx_1] + 5)/2;
                        flat_cg_idx = ((((j1_idx*6 + m1_idx)*3 + j2_idx)*6 + m2_idx)*6 + j3_idx);    // Indexing a 1D array as if it was a 5D array.

                        const double cg_creation = clebsch_gordan_array[flat_cg_idx];

                        // if (cg_creation == 0) continue;  // Might be faster to just multiply with 0 instead of checking.
                        creation_res += creation_sign*cg_creation;//*creation_res_switch;
                    }
                }
                twobody_res += annihilation_norm*creation_norm*tbme*creation_res*annihilation_sign*cg_annihilation;
            }
        }
    }
    return twobody_res;
}

__global__ void twobody_matrix_element_dispatcher(
    double* H,
    const uint64_t* dev_basis_states,
    const uint32_t m_dim,
    const size_t n_orbitals
)
{
    const size_t col_idx = blockIdx.x*blockDim.x + threadIdx.x;
    const size_t row_idx = blockIdx.y*blockDim.y + threadIdx.y;
    const size_t idx = row_idx*m_dim + col_idx;
    // printf("(r, c): (%zu, %zu)\n", row_idx, col_idx);

    if ((col_idx < m_dim) and (row_idx < m_dim))
    {
        const uint64_t left_state = dev_basis_states[row_idx];
        const uint64_t right_state = dev_basis_states[col_idx];
        
        H[idx] = calculate_twobody_matrix_element(
            n_orbitals,
            left_state,
            right_state
        );
    }
}

void create_hamiltonian_device_dispatcher(const Interaction& interaction, const Indices& indices, double* H)
{
    static __device__ uint64_t* dev_basis_states = nullptr;
    dev_init(interaction, indices, dev_basis_states);
    const size_t m_dim = interaction.basis_states.size();
    const size_t n_orbitals = interaction.model_space.n_orbitals;
    static __device__ double* dev_H_diag = nullptr;
    double* H_diag_tmp = new double[m_dim];
    
    auto start = timer();
        const size_t threads_per_block_onebody = 2;
        const size_t blocks_per_grid_onebody = (m_dim + threads_per_block_onebody - 1)/threads_per_block_onebody;
        
        hip_wrappers::hipMalloc(&dev_H_diag, m_dim*sizeof(double));
        
        hipLaunchKernelGGL(
            onebody_matrix_element_dispatcher, dim3(blocks_per_grid_onebody), dim3(threads_per_block_onebody), 0, 0,
            dev_H_diag,
            dev_basis_states,
            m_dim,
            n_orbitals
        );
        hip_wrappers::hipDeviceSynchronize();

        const size_t block_dim = 16;
        const dim3 threads_per_block_twobody(block_dim, block_dim);
        const dim3 blocks_per_grid_twobody(ceil(m_dim/((double)block_dim)), ceil(m_dim/((double)block_dim)));

        hipLaunchKernelGGL(
            twobody_matrix_element_dispatcher, blocks_per_grid_twobody, threads_per_block_twobody, 0, 0,
            dev_H_diag,
            dev_basis_states,
            m_dim,
            n_orbitals
        );
        hip_wrappers::hipDeviceSynchronize();

        hip_wrappers::hipMemcpy(H_diag_tmp, dev_H_diag, m_dim*sizeof(double), hipMemcpyDeviceToHost);
        hip_wrappers::hipFree(dev_H_diag);
    timer(start, "[DEVICE] one-body calc, alloc, copy, and free time");

    for (size_t diag_idx = 0; diag_idx < m_dim; diag_idx++) H[diag_idx*m_dim + diag_idx] = H_diag_tmp[diag_idx]; // Copy data to the diagonal of H.
    delete[] H_diag_tmp;
    dev_uninit(dev_basis_states);
}
}