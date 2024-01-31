#include <bitset>
#include <iomanip>
#include "omp.h"
#include "data_structures.hpp"
#include "parameters.hpp"
#include "basis.hpp"
#include "generate_indices.hpp"
#include "tools.hpp"
#include "bit_manipulation_tools.hpp"
#include "../external/eigen-3.4.0/Eigen/Dense"

namespace hamiltonian_bitset
{
double calculate_onebody_matrix_element_bitset_representation(
    const Interaction& interaction,
    const Indices& indices,
    const std::bitset<n_bits_bitset>& left_state,
    const std::bitset<n_bits_bitset>& right_state
)
{
    double onebody_res = 0;

    for (unsigned short creation_orb_idx = 0; creation_orb_idx < interaction.model_space.n_orbitals; creation_orb_idx++)
    {
        /*
        More generally, there should also be a loop over the same values
        for annihilation_idx, but the SPEs for most, if not all of the
        interaction files, are only defined for when the two indices are
        the same.
        */
        const unsigned short annihilation_orb_idx = creation_orb_idx;
    
        for (unsigned short creation_comp_m_idx : indices.orbital_idx_to_composite_m_idx_map[creation_orb_idx])
        {
            for (unsigned short annihilation_comp_m_idx : indices.orbital_idx_to_composite_m_idx_map[annihilation_orb_idx])
            {
                /*
                Index scheme for the sd model space:

                              22    23        
                            -  O  -  O  -             neutron s1/2: 5
                             -1/2   1/2

                  16    17    18    19    20    21    
                -  O  -  O  -  O  -  O  -  O  -  O    neutron d5/2: 4
                 -5/2  -3/2  -1/2   1/2   3/2   5/2

                        12    13    14    15
                      -  O  -  O  -  O  -  O  -       neutron d3/2: 3
                       -3/2  -1/2   1/2   3/2

                              10    11        
                            -  O  -  O  -             proton s1/2: 2
                             -1/2   1/2

                   4     5     6     7     8     9    
                -  O  -  O  -  O  -  O  -  O  -  O    proton d5/2: 1
                 -5/2  -3/2  -1/2   1/2   3/2   5/2

                         0     1     2     3
                      -  O  -  O  -  O  -  O  -       proton d3/2: 0
                       -3/2  -1/2   1/2   3/2

                orbital_idx_to_composite_m_idx_map translates the orbital
                indices to the composite m indices of the magnetic
                substates. For example, proton d3/2 has orbital index 0
                and composite m substate indices 0, 1, 2, 3.
                */
                std::bitset<n_bits_bitset> new_right_state = right_state;   // The contents of right_state is copied, not referenced.

                // if (not new_right_state.test(annihilation_comp_m_idx))
                if (not new_right_state[annihilation_comp_m_idx])
                {
                    /*
                    If the index cannot be found, then it does not exist
                    in the list. Aka. we are trying to annihilate a
                    substate which is un-occupied and the result is
                    zero.

                    We need to know how many bits before the bit we want
                    to annihilate are set. This corresponds to how many
                    times we must swap the positions of operators and
                    thus if we get a phase of +1 or -1. This is because
                    we have to make sure that the annihilation operator
                    is placed next to the creation operator it tries to
                    annihilate:

                        c_0 | (0, 3) > = c_0 c_0^\dagger c_3^\dagger | core >
                                       = c_3^\dagger | core >
                                       = | (3) >

                        c_0 | (3, 0) > = c_0 c_3^\dagger c_0^\dagger | core >
                                       = - c_0 c_0^\dagger c_3^\dagger | core >
                                       = - c_3^\dagger | core >
                                       = - | (3) >
                    */
                    continue;
                }

                const unsigned short n_operator_swaps_annihilation = bittools::reset_bit_and_count_swaps(new_right_state, annihilation_comp_m_idx);
                const short annihilation_sign = bittools::negative_one_pow(n_operator_swaps_annihilation);
                
                // if (new_right_state.test(creation_comp_m_idx))
                if (new_right_state[creation_comp_m_idx])
                {
                    /*
                    `creation_comp_m_idx` already exists in
                    `new_right_state`. Creating a state which already
                    exists yields 0.
                    */
                    continue;
                }
                const unsigned short n_operator_swaps_creation = bittools::set_bit_and_count_swaps(new_right_state, creation_comp_m_idx);
                const short creation_sign = bittools::negative_one_pow(n_operator_swaps_creation);

                if (left_state != new_right_state)
                {
                    /*
                    If the states are not equal, the inner product will
                    be zero:

                        < i | j > = delta_ij
                    */
                    continue;
                }
                onebody_res += annihilation_sign*creation_sign*interaction.spe[creation_orb_idx];   // Or annihilation_orb_idx, they are the same.
            }
        }
    }
    return onebody_res;
}

double calculate_twobody_matrix_element_bitset_representation(
    const Interaction& interaction,
    const Indices& indices,
    const std::bitset<n_bits_bitset>& left_state,
    const std::bitset<n_bits_bitset>& right_state
)
{
    double twobody_res = 0;
    // double creation_res_switch = 0; // For removing if statements in the inner creation loop. Multiply the result by 0 instead of using if (condition) continue;.
    const unsigned int n_indices = indices.creation_orb_indices_0.size();
    for (unsigned int i = 0; i < n_indices; i++)
    {
        /*
        The values in the following vectors corresponds to using these
        nested loops:

        for creation_orb_idx_0 in range(n_orbitals):
            for creation_orb_idx_1 in range(creation_orb_idx_0, n_orbitals):
                
                for annihilation_orb_idx_0 in range(n_orbitals):
                    for annihilation_orb_idx_1 in range(annihilation_orb_idx_0, n_orbitals):
                    
                        for j_coupled in range(j_min, j_max + 2, 2):
                            for m_coupled in range(-j_coupled, j_coupled + 2, 2):
        
        It gives good reduction in program run time by using
        pre-calculated indices instead of four nested loops.
        */
        const unsigned short creation_orb_idx_0 = indices.creation_orb_indices_0[i];
        const unsigned short creation_orb_idx_1 = indices.creation_orb_indices_1[i];
        const unsigned short annihilation_orb_idx_0 = indices.annihilation_orb_indices_0[i];
        const unsigned short annihilation_orb_idx_1 = indices.annihilation_orb_indices_1[i];
        const unsigned short j_coupled = indices.j_coupled[i];
        const short m_coupled = indices.m_coupled[i];
        const double tbme = indices.tbme[i];

        const double creation_norm = 1/std::sqrt(1 + (creation_orb_idx_0 == creation_orb_idx_1));
        const double annihilation_norm = 1/std::sqrt(1 + (annihilation_orb_idx_0 == annihilation_orb_idx_1));

        // Annihilation terms
        for (unsigned short annihilation_comp_m_idx_0 : indices.orbital_idx_to_composite_m_idx_map[annihilation_orb_idx_0])
        {
            for (unsigned short annihilation_comp_m_idx_1 : indices.orbital_idx_to_composite_m_idx_map[annihilation_orb_idx_1])
            {
                // if (not right_state.test(annihilation_comp_m_idx_0))
                if (not right_state[annihilation_comp_m_idx_0])
                {
                    /*
                    If the index cannot be found, then it does not exist
                    in the list. Aka. we are trying to annihilate a
                    substate which is un-occupied and the result is
                    zero.

                    We need to know how many bits before the bit we want
                    to annihilate are set. This corresponds to how many
                    times we must swap the positions of operators and
                    thus if we get a phase of +1 or -1. This is because
                    we have to make sure that the annihilation operator
                    is placed next to the creation operator it tries to
                    annihilate:

                        c_0 | (0, 3) > = c_0 c_0^\dagger c_3^\dagger | core >
                                       = c_3^\dagger | core >
                                       = | (3) >

                        c_0 | (3, 0) > = c_0 c_3^\dagger c_0^\dagger | core >
                                       = - c_0 c_0^\dagger c_3^\dagger | core >
                                       = - c_3^\dagger | core >
                                       = - | (3) >
                    */
                    continue;
                }

                std::bitset<n_bits_bitset> new_right_state_annihilation = right_state;
                const unsigned short n_operator_swaps_annihilation_0 = bittools::reset_bit_and_count_swaps(new_right_state_annihilation, annihilation_comp_m_idx_0);
                short annihilation_sign = bittools::negative_one_pow(n_operator_swaps_annihilation_0);

                // if (not new_right_state_annihilation.test(annihilation_comp_m_idx_1)) continue;
                if (not new_right_state_annihilation[annihilation_comp_m_idx_1]) continue;
                const unsigned short n_operator_swaps_annihilation_1 = bittools::reset_bit_and_count_swaps(new_right_state_annihilation, annihilation_comp_m_idx_1);
                annihilation_sign *= bittools::negative_one_pow(n_operator_swaps_annihilation_1);

                const Key6 annihilation_key = {
                    indices.orbital_idx_to_j_map[annihilation_orb_idx_0],
                    indices.composite_m_idx_to_m_map[annihilation_comp_m_idx_0],
                    indices.orbital_idx_to_j_map[annihilation_orb_idx_1],
                    indices.composite_m_idx_to_m_map[annihilation_comp_m_idx_1],
                    j_coupled,
                    m_coupled
                };

                const double cg_annihilation = clebsch_gordan.at(annihilation_key);
                if (cg_annihilation == 0) continue;

                // Creation terms
                double creation_res = 0.0;

                for (unsigned short creation_comp_m_idx_0 : indices.orbital_idx_to_composite_m_idx_map[creation_orb_idx_0])
                {
                    for (unsigned short creation_comp_m_idx_1 : indices.orbital_idx_to_composite_m_idx_map[creation_orb_idx_1])
                    {
                        /*
                        Performance notes
                        -----------------
                        Using a "switch" variable instead of if-statements
                        gives worse performance, but that might be due to
                        increased number of CG dictionary lookups.

                        creation_res_switch = (not new_right_state_annihilation[creation_comp_m_idx_1]);
                        creation_res_switch *= (not new_right_state_annihilation[creation_comp_m_idx_0]);
                        creation_res_switch *= (left_state == new_right_state_creation);

                        The multiplications themselves might be better than ifs
                        for potential GPU acceleration.
                        */
                        
                        // if (new_right_state_annihilation.test(creation_comp_m_idx_1)) continue;
                        if (new_right_state_annihilation[creation_comp_m_idx_1]) continue;
                        // creation_res_switch = (not new_right_state_annihilation[creation_comp_m_idx_1]);
                        
                        std::bitset<n_bits_bitset> new_right_state_creation = new_right_state_annihilation;
                        const unsigned short n_operator_swaps_creation_1 = bittools::set_bit_and_count_swaps(new_right_state_creation, creation_comp_m_idx_1);
                        short creation_sign = bittools::negative_one_pow(n_operator_swaps_creation_1);

                        // if (new_right_state_annihilation.test(creation_comp_m_idx_0)) continue;
                        if (new_right_state_annihilation[creation_comp_m_idx_0]) continue;
                        // creation_res_switch *= (not new_right_state_annihilation[creation_comp_m_idx_0]);
                        
                        const unsigned short n_operator_swaps_creation_0 = bittools::set_bit_and_count_swaps(new_right_state_creation, creation_comp_m_idx_0);
                        creation_sign *= bittools::negative_one_pow(n_operator_swaps_creation_0);

                        if (left_state != new_right_state_creation) continue;
                        // creation_res_switch *= (left_state == new_right_state_creation);

                        const Key6 creation_key = {
                            indices.orbital_idx_to_j_map[creation_orb_idx_0],
                            indices.composite_m_idx_to_m_map[creation_comp_m_idx_0],
                            indices.orbital_idx_to_j_map[creation_orb_idx_1],
                            indices.composite_m_idx_to_m_map[creation_comp_m_idx_1],
                            j_coupled,
                            m_coupled
                        };

                        const double cg_creation = clebsch_gordan.at(creation_key);
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

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> create_hamiltonian_bitset_representation(const Interaction& interaction)
{
    const Indices indices = indices::generate_indices(interaction);
    const std::vector<std::bitset<n_bits_bitset>> basis_states = basis::calculate_m_basis_states_bitset_representation(interaction);
    const unsigned int m_dim = basis_states.size();

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> H;
    H.resize(m_dim, m_dim);
    H.setZero();

    auto start = timer();
    for (unsigned int left_idx = 0; left_idx < m_dim; left_idx++)
    {
        for (unsigned int right_idx = left_idx; right_idx < m_dim; right_idx++)
        {
            /*
            Calculate only the upper triangle of the Hamiltonian matrix.
            H is hermitian so we dont have to calculate both triangles.
            */
            H(left_idx, right_idx) += calculate_onebody_matrix_element_bitset_representation(
                interaction,
                indices,
                basis_states[left_idx],
                basis_states[right_idx]
            );
        }
    }
    timer(start, "calculate_onebody_matrix_element_bit_representation");
    start = timer();
    int thread_id = omp_get_thread_num();
    std::vector<long long> loop_timings;
    auto loop_timer = timer();

    for (unsigned int left_idx = 0; left_idx < m_dim; left_idx++)
    {   
        if (thread_id == 0) loop_timer = timer();
        
        #pragma omp parallel for
        for (int right_idx = left_idx; right_idx < m_dim; right_idx++)
        {    
            H(left_idx, right_idx) += calculate_twobody_matrix_element_bitset_representation(
                interaction,
                indices,
                basis_states[left_idx],
                basis_states[right_idx]
            );
        }
        if (thread_id == 0)
        {
            loop_timings.push_back(timer(loop_timer));
            double mean_time = mean(loop_timings)/1000;
            int num_threads = omp_get_num_threads();
            
            cout << "\r[" << left_idx << " of ≈ " << (double)m_dim/num_threads << "]";
            cout << " [loop time: ";
            cout << std::setfill(' ') << std::setw(5) << loop_timings.back()/1000. << " s";
            cout << " - mean loop time: ";
            cout << std::setfill(' ') << std::setw(10) << mean_time << " s]";
            cout << " [est. time left: ";
            cout << std::setfill(' ') << std::setw(10) << (m_dim - (left_idx + 1))*mean_time << " s]" << std::flush;
        }
    }
    cout << endl;   // For the progress bar.
    timer(start, "calculate_twobody_matrix_element_bit_representation");
    complete_hermitian_matrix(H);    
    return H;
}
}