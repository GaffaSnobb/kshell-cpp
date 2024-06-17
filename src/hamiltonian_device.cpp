#include <vector>
#include <iostream>
#include <cstddef> // stdints
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

bool is_show_dev_const_mem_usage = true;    // For making sure that the info is only displayed once per program execution.

__constant__ uint16_t dev_const_creation_orb_indices_0[CONST_MEM_ARR_LEN_INDICES];
__constant__ uint16_t dev_const_creation_orb_indices_1[CONST_MEM_ARR_LEN_INDICES];
__constant__ uint16_t dev_const_annihilation_orb_indices_0[CONST_MEM_ARR_LEN_INDICES];
__constant__ uint16_t dev_const_annihilation_orb_indices_1[CONST_MEM_ARR_LEN_INDICES];
__constant__ uint16_t dev_const_j_coupled[CONST_MEM_ARR_LEN_INDICES];
__constant__ int16_t dev_const_m_coupled[CONST_MEM_ARR_LEN_INDICES];
__constant__ uint16_t dev_const_orbital_idx_to_composite_m_idx_map_flattened_indices[CONST_MEM_ARR_LEN_N_ORBITALS];
__constant__ double dev_const_tbme[CONST_MEM_ARR_LEN_INDICES];
__constant__ double dev_const_spe[CONST_MEM_ARR_LEN_N_ORBITALS];

__constant__ uint16_t dev_const_annihilation_comp_m_start_idx_0[CONST_MEM_ARR_LEN_INDICES];
__constant__ uint16_t dev_const_annihilation_comp_m_end_idx_0[CONST_MEM_ARR_LEN_INDICES];
__constant__ uint16_t dev_const_annihilation_comp_m_start_idx_1[CONST_MEM_ARR_LEN_INDICES];
__constant__ uint16_t dev_const_annihilation_comp_m_end_idx_1[CONST_MEM_ARR_LEN_INDICES];
__constant__ uint16_t dev_const_creation_comp_m_start_idx_0[CONST_MEM_ARR_LEN_INDICES];
__constant__ uint16_t dev_const_creation_comp_m_end_idx_0[CONST_MEM_ARR_LEN_INDICES];
__constant__ uint16_t dev_const_creation_comp_m_start_idx_1[CONST_MEM_ARR_LEN_INDICES];
__constant__ uint16_t dev_const_creation_comp_m_end_idx_1[CONST_MEM_ARR_LEN_INDICES];
__constant__ int16_t dev_const_composite_m_idx_to_m_map[CONST_MEM_ARR_LEN_MODEL_SPACE_DEGENERACY];
__constant__ uint16_t dev_const_orbital_idx_to_j_map[CONST_MEM_ARR_LEN_N_ORBITALS];

__host__ void dev_init(
    const Interaction& interaction,
    const Indices& indices,
    uint64_t*& dev_basis_states,
    double*& dev_H_diag,
    double*& dev_H_upper_triangle
)
{
    /*
    Add shit to device memory.
    */
    auto start = timer();
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

    const size_t acmsi_0 = indices.annihilation_comp_m_start_idx_0.size()*sizeof(uint16_t);
    const size_t acmei_0 = indices.annihilation_comp_m_end_idx_0.size()*sizeof(uint16_t);
    const size_t acmsi_1 = indices.annihilation_comp_m_start_idx_1.size()*sizeof(uint16_t);
    const size_t acmei_1 = indices.annihilation_comp_m_end_idx_1.size()*sizeof(uint16_t);
    const size_t ccmsi_0 = indices.creation_comp_m_start_idx_0.size()*sizeof(uint16_t);
    const size_t ccmei_0 = indices.creation_comp_m_end_idx_0.size()*sizeof(uint16_t);
    const size_t ccmsi_1 = indices.creation_comp_m_start_idx_1.size()*sizeof(uint16_t);
    const size_t ccmei_1 = indices.creation_comp_m_end_idx_1.size()*sizeof(uint16_t);
    const size_t cmitmm = indices.composite_m_idx_to_m_map.size()*sizeof(int16_t);
    const size_t oitjm = indices.orbital_idx_to_j_map.size()*sizeof(uint16_t);

    const size_t total = coi_0 + coi_1 + aoi_0 + aoi_1 + jc + mc + tbme + oitcmimf + spe + acmsi_0 + acmei_0 + acmsi_1 + acmei_1 + ccmsi_0 + ccmei_0 + ccmsi_1 + ccmei_1 + cmitmm + oitjm;

    assert(coi_0 <= (CONST_MEM_ARR_LEN_INDICES*sizeof(uint16_t)));
    assert(coi_1 <= (CONST_MEM_ARR_LEN_INDICES*sizeof(uint16_t)));
    assert(aoi_0 <= (CONST_MEM_ARR_LEN_INDICES*sizeof(uint16_t)));
    assert(aoi_1 <= (CONST_MEM_ARR_LEN_INDICES*sizeof(uint16_t)));
    assert(jc <= (CONST_MEM_ARR_LEN_INDICES*sizeof(uint16_t)));
    assert(mc <= (CONST_MEM_ARR_LEN_INDICES*sizeof(int16_t)));
    assert(oitcmimf <= (CONST_MEM_ARR_LEN_N_ORBITALS*sizeof(uint16_t)));
    assert(tbme <= (CONST_MEM_ARR_LEN_INDICES*sizeof(double)));
    assert(spe <= (CONST_MEM_ARR_LEN_N_ORBITALS*sizeof(double)));

    assert(acmsi_0 <= (CONST_MEM_ARR_LEN_INDICES*sizeof(uint16_t)));
    assert(acmei_0 <= (CONST_MEM_ARR_LEN_INDICES*sizeof(uint16_t)));
    assert(acmsi_1 <= (CONST_MEM_ARR_LEN_INDICES*sizeof(uint16_t)));
    assert(acmei_1 <= (CONST_MEM_ARR_LEN_INDICES*sizeof(uint16_t)));
    assert(ccmsi_0 <= (CONST_MEM_ARR_LEN_INDICES*sizeof(uint16_t)));
    assert(ccmei_0 <= (CONST_MEM_ARR_LEN_INDICES*sizeof(uint16_t)));
    assert(ccmsi_1 <= (CONST_MEM_ARR_LEN_INDICES*sizeof(uint16_t)));
    assert(ccmei_1 <= (CONST_MEM_ARR_LEN_INDICES*sizeof(uint16_t)));
    assert(cmitmm <= (CONST_MEM_ARR_LEN_MODEL_SPACE_DEGENERACY*sizeof(int16_t)));
    assert(oitjm <= (CONST_MEM_ARR_LEN_N_ORBITALS*sizeof(uint16_t)));
    
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

    hip_wrappers::hipMemcpyToSymbol(dev_const_annihilation_comp_m_start_idx_0, indices.annihilation_comp_m_start_idx_0);
    hip_wrappers::hipMemcpyToSymbol(dev_const_annihilation_comp_m_end_idx_0, indices.annihilation_comp_m_end_idx_0);
    hip_wrappers::hipMemcpyToSymbol(dev_const_annihilation_comp_m_start_idx_1, indices.annihilation_comp_m_start_idx_1);
    hip_wrappers::hipMemcpyToSymbol(dev_const_annihilation_comp_m_end_idx_1, indices.annihilation_comp_m_end_idx_1);
    hip_wrappers::hipMemcpyToSymbol(dev_const_creation_comp_m_start_idx_0, indices.creation_comp_m_start_idx_0);
    hip_wrappers::hipMemcpyToSymbol(dev_const_creation_comp_m_end_idx_0, indices.creation_comp_m_end_idx_0);
    hip_wrappers::hipMemcpyToSymbol(dev_const_creation_comp_m_start_idx_1, indices.creation_comp_m_start_idx_1);
    hip_wrappers::hipMemcpyToSymbol(dev_const_creation_comp_m_end_idx_1, indices.creation_comp_m_end_idx_1);
    hip_wrappers::hipMemcpyToSymbol(dev_const_composite_m_idx_to_m_map, indices.composite_m_idx_to_m_map);
    hip_wrappers::hipMemcpyToSymbol(dev_const_orbital_idx_to_j_map, indices.orbital_idx_to_j_map);
    
    if (is_show_dev_const_mem_usage)
    {
        is_show_dev_const_mem_usage = false;
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

        cout << acmsi_0/1e3 << " kB" << " (" << indices.annihilation_comp_m_start_idx_0.size() << " elements)" << endl;
        cout << acmei_0/1e3 << " kB" << " (" << indices.annihilation_comp_m_end_idx_0.size() << " elements)" << endl;
        cout << acmsi_1/1e3 << " kB" << " (" << indices.annihilation_comp_m_start_idx_1.size() << " elements)" << endl;
        cout << acmei_1/1e3 << " kB" << " (" << indices.annihilation_comp_m_end_idx_1.size() << " elements)" << endl;
        cout << ccmsi_0/1e3 << " kB" << " (" << indices.creation_comp_m_start_idx_0.size() << " elements)" << endl;
        cout << ccmei_0/1e3 << " kB" << " (" << indices.creation_comp_m_end_idx_0.size() << " elements)" << endl;
        cout << ccmsi_1/1e3 << " kB" << " (" << indices.creation_comp_m_start_idx_1.size() << " elements)" << endl;
        cout << ccmei_1/1e3 << " kB" << " (" << indices.creation_comp_m_end_idx_1.size() << " elements)" << endl;
        cout << cmitmm/1e3 << " kB" << " (" << indices.composite_m_idx_to_m_map.size() << " elements)" << endl;
        cout << oitjm/1e3 << " kB" << " (" << indices.orbital_idx_to_j_map.size() << " elements)" << endl;
        cout << diagnostics::DIAG_STR_END << endl;
    }

    
    // Device arrays.
    hip_wrappers::hipMalloc(&dev_basis_states, m_dim*sizeof(uint64_t));
    hip_wrappers::hipMemcpy(dev_basis_states, interaction.basis_states.data(), m_dim*sizeof(uint64_t), hipMemcpyHostToDevice);

    hip_wrappers::hipMalloc(&dev_H_diag, m_dim*sizeof(double));
    hip_wrappers::hipMemset(dev_H_diag, 0, m_dim*sizeof(double));

    const size_t n_elements_upper_triangle = m_dim*(m_dim + 1)/2;

    hip_wrappers::hipMalloc(&dev_H_upper_triangle, n_elements_upper_triangle*sizeof(double));
    hip_wrappers::hipMemset(dev_H_upper_triangle, 0, n_elements_upper_triangle*sizeof(double));

    timer(start, "[DEVICE_INFO] dev_init time");
}

__host__ void dev_uninit(uint64_t*& dev_basis_states, double*& dev_H_diag, double*& dev_H_upper_triangle)
{
    hip_wrappers::hipFree(dev_basis_states);
    hip_wrappers::hipFree(dev_H_diag);
    hip_wrappers::hipFree(dev_H_upper_triangle);
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
    const size_t n_indices,
    const uint64_t& left_state,
    const uint64_t& right_state
)
{
    uint16_t j1_idx, m1_idx, j2_idx, m2_idx, j3_idx;  // Indices for the Clebsch-Gordan 5D array...
    uint16_t flat_cg_idx;                             // ...translated to a flattened 1D array.

    double twobody_res = 0;
    // double creation_res_switch = 0; // For removing if statements in the inner creation loop. Multiply the result by 0 instead of using if (condition) continue;.
    for (size_t i = 0; i < n_indices; i++)
    {
        const uint16_t creation_orb_idx_0 = dev_const_creation_orb_indices_0[i];
        const uint16_t creation_orb_idx_1 = dev_const_creation_orb_indices_1[i];
        const uint16_t annihilation_orb_idx_0 = dev_const_annihilation_orb_indices_0[i];
        const uint16_t annihilation_orb_idx_1 = dev_const_annihilation_orb_indices_1[i];
        const uint16_t j_coupled = dev_const_j_coupled[i];
        const int16_t m_coupled = dev_const_m_coupled[i];
        const double tbme = dev_const_tbme[i];

        const uint16_t annihilation_comp_m_start_idx_0 = dev_const_annihilation_comp_m_start_idx_0[i];
        const uint16_t annihilation_comp_m_end_idx_0 = dev_const_annihilation_comp_m_end_idx_0[i];
        const uint16_t annihilation_comp_m_start_idx_1 = dev_const_annihilation_comp_m_start_idx_1[i];
        const uint16_t annihilation_comp_m_end_idx_1 = dev_const_annihilation_comp_m_end_idx_1[i];
        const uint16_t creation_comp_m_start_idx_0 = dev_const_creation_comp_m_start_idx_0[i];
        const uint16_t creation_comp_m_end_idx_0 = dev_const_creation_comp_m_end_idx_0[i];
        const uint16_t creation_comp_m_start_idx_1 = dev_const_creation_comp_m_start_idx_1[i];
        const uint16_t creation_comp_m_end_idx_1 = dev_const_creation_comp_m_end_idx_1[i];

        const double creation_norm = dev_const_inverse_sqrt_2[creation_orb_idx_0 == creation_orb_idx_1];
        const double annihilation_norm = dev_const_inverse_sqrt_2[annihilation_orb_idx_0 == annihilation_orb_idx_1];

        // Annihilation terms
        // for (uint16_t annihilation_comp_m_idx_0 : indices.orbital_idx_to_composite_m_idx_map[annihilation_orb_idx_0])
        for (uint16_t annihilation_comp_m_idx_0 = annihilation_comp_m_start_idx_0; annihilation_comp_m_idx_0 < annihilation_comp_m_end_idx_0; annihilation_comp_m_idx_0++)
        {
            // for (uint16_t annihilation_comp_m_idx_1 : indices.orbital_idx_to_composite_m_idx_map[annihilation_orb_idx_1])
            for (uint16_t annihilation_comp_m_idx_1 = annihilation_comp_m_start_idx_1; annihilation_comp_m_idx_1 < annihilation_comp_m_end_idx_1; annihilation_comp_m_idx_1++)
            {
                if (
                    (dev_const_composite_m_idx_to_m_map[annihilation_comp_m_idx_0] +  // m1
                    dev_const_composite_m_idx_to_m_map[annihilation_comp_m_idx_1]) != // m2
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
                
                j1_idx = (dev_const_orbital_idx_to_j_map[annihilation_orb_idx_0] + 1)/2 - 1;
                m1_idx = (dev_const_composite_m_idx_to_m_map[annihilation_comp_m_idx_0] + 5)/2;
                j2_idx = (dev_const_orbital_idx_to_j_map[annihilation_orb_idx_1] + 1)/2 - 1;
                m2_idx = (dev_const_composite_m_idx_to_m_map[annihilation_comp_m_idx_1] + 5)/2;
                j3_idx = j_coupled/2;
                flat_cg_idx = ((((j1_idx*6 + m1_idx)*3 + j2_idx)*6 + m2_idx)*6 + j3_idx);    // Indexing a 1D array as if it was a 5D array.

                const double cg_annihilation = dev_const_clebsch_gordan_array[flat_cg_idx];

                // Creation terms
                double creation_res = 0.0;

                // for (uint16_t creation_comp_m_idx_0 : indices.orbital_idx_to_composite_m_idx_map[creation_orb_idx_0])
                for (uint16_t creation_comp_m_idx_0 = creation_comp_m_start_idx_0; creation_comp_m_idx_0 < creation_comp_m_end_idx_0; creation_comp_m_idx_0++)
                {
                    // for (uint16_t creation_comp_m_idx_1 : indices.orbital_idx_to_composite_m_idx_map[creation_orb_idx_1])
                    for (uint16_t creation_comp_m_idx_1 = creation_comp_m_start_idx_1; creation_comp_m_idx_1 < creation_comp_m_end_idx_1; creation_comp_m_idx_1++)
                    {   
                        if (
                            (dev_const_composite_m_idx_to_m_map[creation_comp_m_idx_0] +  // m1
                            dev_const_composite_m_idx_to_m_map[creation_comp_m_idx_1]) != // m2
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

                        j1_idx = (dev_const_orbital_idx_to_j_map[creation_orb_idx_0] + 1)/2 - 1;
                        m1_idx = (dev_const_composite_m_idx_to_m_map[creation_comp_m_idx_0] + 5)/2;
                        j2_idx = (dev_const_orbital_idx_to_j_map[creation_orb_idx_1] + 1)/2 - 1;
                        m2_idx = (dev_const_composite_m_idx_to_m_map[creation_comp_m_idx_1] + 5)/2;
                        flat_cg_idx = ((((j1_idx*6 + m1_idx)*3 + j2_idx)*6 + m2_idx)*6 + j3_idx);    // Indexing a 1D array as if it was a 5D array.

                        const double cg_creation = dev_const_clebsch_gordan_array[flat_cg_idx];

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

__device__ inline void flat_tri_idx_to_2d_indices_ascending(const size_t n_rows, const size_t flat_tri_idx, size_t &row_idx, size_t &col_idx)
{   /*
    Find the row, col indices of the original matrix from the flat index
    of the upper triangle (with diagonal). Search for row in ascending
    order. Heres the Python code for the descending version:

def flat_tri_idx_to_2d_indices_descending(n, flat_tri_idx):
    """
    Find the row, col indices of the original matrix from the flat index
    of the upper triangle (with diagonal). Search for row in descending
    order.
    """
    cum_n_triangular_elements = n*(n + 1)//2
    for row_idx in reversed(range(n)):
        
        n_triangular_elements_row = n - row_idx
        cum_n_triangular_elements -= n_triangular_elements_row    # Cumulative number of triangular elements at row `row_idx`.

        if flat_tri_idx >= cum_n_triangular_elements:
            cum_n_triangular_elements += n_triangular_elements_row
            break

    col_idx = n - (cum_n_triangular_elements - flat_tri_idx)
    return row_idx, col_idx
    */
    size_t cum_n_triangular_elements = 0;
    for (row_idx = 0; row_idx < n_rows; row_idx++)
    {
        const size_t n_triangular_elements_row = n_rows - row_idx;
        cum_n_triangular_elements += n_triangular_elements_row;  // Cumulative number of triangular elements at row `row_idx`.

        if (flat_tri_idx < cum_n_triangular_elements) break;

    }
    col_idx = n_rows - (cum_n_triangular_elements - flat_tri_idx);
}

__global__ void twobody_matrix_element_dispatcher(
    double* H_diag,
    double* H_upper_triangle,
    const uint64_t* dev_basis_states,
    const uint32_t m_dim,
    const size_t n_orbitals,
    const size_t n_indices
)
{
    const size_t n_elements_upper_triangle = m_dim*(m_dim + 1)/2;
    const size_t idx = blockIdx.x*blockDim.x + threadIdx.x;

    if (idx < n_elements_upper_triangle)
    {
        size_t row_idx;
        size_t col_idx;

        flat_tri_idx_to_2d_indices_ascending(m_dim, idx, row_idx, col_idx);

        const uint64_t left_state = dev_basis_states[row_idx];
        const uint64_t right_state = dev_basis_states[col_idx];
        
        H_upper_triangle[idx] = calculate_twobody_matrix_element(
            n_orbitals,
            n_indices,
            left_state,
            right_state
        );

        if (col_idx == row_idx)
        {   /*
            Add one-body results. They are located only on the diagonal.
            */
            H_upper_triangle[idx] += H_diag[col_idx];  // Or row_idx.
        }
    }
}

void create_hamiltonian_device_dispatcher(const Interaction& interaction, const Indices& indices, double* H)
{
    static __device__ uint64_t* dev_basis_states = nullptr;
    static __device__ double* dev_H_upper_triangle = nullptr;
    static __device__ double* dev_H_diag = nullptr;
    dev_init(interaction, indices, dev_basis_states, dev_H_diag, dev_H_upper_triangle);
    const size_t m_dim = interaction.basis_states.size();
    const size_t n_orbitals = interaction.model_space.n_orbitals;
    const size_t n_indices = indices.creation_orb_indices_0.size();
    const size_t n_elements_upper_triangle = m_dim*(m_dim + 1)/2;

    double* H_upper_triangle = new double[n_elements_upper_triangle];
    
    auto start = timer();
        const size_t threads_per_block_onebody = 256;
        const size_t blocks_per_grid_onebody = (m_dim + threads_per_block_onebody - 1)/threads_per_block_onebody;
        cout << "[DEVICE_INFO] threads_per_block_onebody: " << threads_per_block_onebody << endl;
        cout << "[DEVICE_INFO] blocks_per_grid_onebody: " << blocks_per_grid_onebody << endl;
        
        hipLaunchKernelGGL(
            onebody_matrix_element_dispatcher, dim3(blocks_per_grid_onebody), dim3(threads_per_block_onebody), 0, 0,
            dev_H_diag,
            dev_basis_states,
            m_dim,
            n_orbitals
        );
        hip_wrappers::hipGetLastError();
        hip_wrappers::hipDeviceSynchronize();

        const size_t threads_per_block_twobody = 64;
        const size_t blocks_per_grid_twobody = (n_elements_upper_triangle + threads_per_block_twobody - 1)/threads_per_block_twobody;
        cout << "\n[DEVICE_INFO] threads_per_block_twobody: " << threads_per_block_twobody << endl;
        cout << "[DEVICE_INFO] blocks_per_grid_twobody: " << blocks_per_grid_twobody << endl;
        cout << "[DEVICE_INFO] total threads: " << threads_per_block_twobody*blocks_per_grid_twobody << endl;
        cout << "[DEVICE_INFO] n_elements_upper_triangle: " << n_elements_upper_triangle << endl;

        hipLaunchKernelGGL(
            twobody_matrix_element_dispatcher, blocks_per_grid_twobody, threads_per_block_twobody, 0, 0,
            dev_H_diag,
            dev_H_upper_triangle,
            dev_basis_states,
            m_dim,
            n_orbitals,
            n_indices
        );
        hip_wrappers::hipGetLastError();
        hip_wrappers::hipDeviceSynchronize();
        hip_wrappers::hipMemcpy(H_upper_triangle, dev_H_upper_triangle, n_elements_upper_triangle*sizeof(double), hipMemcpyDeviceToHost);

        dev_uninit(dev_basis_states, dev_H_diag, dev_H_upper_triangle);
    timer(start, "[DEVICE_INFO] H calc, alloc, copy, and free time");

    size_t lim = m_dim;
    size_t counter = 0;
    size_t prev_offset = 0;
    size_t row_idx = 0;

    for (size_t flat_tri_idx = 0; flat_tri_idx < n_elements_upper_triangle; flat_tri_idx++)
    {   /*
        Convert an index from the flat upper triangle array into an
        index of the flat H array. If the H array is

            H = [[ 0  1  2  3]
                 [ 4  5  6  7]
                 [ 8  9 10 11]
                 [12 13 14 15]]

        the flat triangle array will be

            tri = [0, 1, 2, 3, 5, 6, 7, 10, 11, 15]

        tri[5] with value 6 corresponds to element H[1][2].
        */
        if (counter >= lim)
        {   /*
            Row 0 has lim number of upper triangular elements, row 1 has
            lim - 1 elements, row 2 lim - 2, ...

            The first triangular element of row 0 is offset by 0
            The first triangular element of row 1 is offset by 0 + 1
            The first triangular element of row 2 is offset by 0 + 1 + 2
            ...
            */
            lim--;
            counter = 0;
            prev_offset = row_idx + prev_offset;
        }

        row_idx = m_dim - lim;
        size_t flat_idx = flat_tri_idx + row_idx + prev_offset;
        counter++;

        H[flat_idx] = H_upper_triangle[flat_tri_idx];
    }

    delete[] H_upper_triangle;
}
}