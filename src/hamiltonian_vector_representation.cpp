#include <vector>
#include <cmath>
#include "data_structures.hpp"
#include "tools.hpp"
#include "basis.hpp"
#include "../external/eigen-3.4.0/Eigen/Dense"
#include "../external/tqdm-cpp/tqdm.hpp"

double calculate_onebody_matrix_element(
    const Interaction& interaction,
    const Indices& indices,
    const std::vector<unsigned short>& left_state,
    const std::vector<unsigned short>& right_state
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
        unsigned short annihilation_orb_idx = creation_orb_idx;
    
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
                std::vector<unsigned short> new_right_state = right_state;
                short annihilated_substate_idx = index(new_right_state, annihilation_comp_m_idx);

                if (annihilated_substate_idx == -1)
                {
                    /*
                    If the index cannot be found, then it does not exist
                    in the list. Aka. we are trying to annihilate a
                    substate which is un-occupied and the result is
                    zero.

                    We need the index of the substate which will be
                    annihilated because there might be a phase of -1.
                    This is because we have to make sure that the
                    annihilation operator is placed next to the creation
                    operator it tries to annihilate:

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
                
                short annihilation_sign = negative_one_pow(annihilated_substate_idx);
                new_right_state.erase(new_right_state.begin() + annihilated_substate_idx);

                short created_substate_idx = check_existence_and_bisect(new_right_state, creation_comp_m_idx);
                if (created_substate_idx == -1)
                {
                    /*
                    -1 means that `creation_comp_m_idx` already exists
                    in `new_right_state`. Creating a state which already
                    exists yields 0.
                    */
                    continue;
                }

                new_right_state.insert(new_right_state.begin() + created_substate_idx, creation_comp_m_idx);
                short creation_sign = negative_one_pow(created_substate_idx);

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

double calculate_twobody_matrix_element(
    const Interaction& interaction,
    const Indices& indices,
    const std::vector<unsigned short>& left_state,
    const std::vector<unsigned short>& right_state
)
{
    std::vector<unsigned short> new_right_state_creation;
    double twobody_res = 0;
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
                const short annihilated_substate_idx_0 = index(right_state, annihilation_comp_m_idx_0);

                if (annihilated_substate_idx_0 == -1)
                {
                    /*
                    If the index cannot be found, then it does not exist
                    in the list. Aka. we are trying to annihilate a
                    substate which is un-occupied and the result is
                    zero.

                    We need the index of the substate which will be
                    annihilated because there might be a phase of -1.
                    This is because we have to make sure that the
                    annihilation operator is placed next to the creation
                    operator it tries to annihilate:

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
                
                std::vector<unsigned short> new_right_state_annihilation = right_state;
                short annihilation_sign = negative_one_pow(annihilated_substate_idx_0);
                new_right_state_annihilation.erase(new_right_state_annihilation.begin() + annihilated_substate_idx_0);

                short annihilated_substate_idx_1 = index(new_right_state_annihilation, annihilation_comp_m_idx_1);
                if (annihilated_substate_idx_1 == -1) continue;
                annihilation_sign *= negative_one_pow(annihilated_substate_idx_1);
                new_right_state_annihilation.erase(new_right_state_annihilation.begin() + annihilated_substate_idx_1);

                Key6 annihilation_key = {
                    indices.orbital_idx_to_j_map[annihilation_orb_idx_0],
                    indices.composite_m_idx_to_m_map[annihilation_comp_m_idx_0],
                    indices.orbital_idx_to_j_map[annihilation_orb_idx_1],
                    indices.composite_m_idx_to_m_map[annihilation_comp_m_idx_1],
                    j_coupled,
                    m_coupled
                };

                double cg_annihilation = clebsch_gordan.at(annihilation_key);
                // double cg_annihilation = 1;
                if (cg_annihilation == 0) continue;

                // Creation terms
                double creation_res = 0.0;

                for (unsigned short creation_comp_m_idx_0 : indices.orbital_idx_to_composite_m_idx_map[creation_orb_idx_0])
                {
                    for (unsigned short creation_comp_m_idx_1 : indices.orbital_idx_to_composite_m_idx_map[creation_orb_idx_1])
                    {
                        std::vector<unsigned short> new_right_state_creation = new_right_state_annihilation;
                        
                        short created_substate_idx_1 = check_existence_and_bisect(new_right_state_creation, creation_comp_m_idx_1);
                        if (created_substate_idx_1 == -1) continue;
                        short creation_sign = negative_one_pow(created_substate_idx_1);
                        new_right_state_creation.insert(new_right_state_creation.begin() + created_substate_idx_1, creation_comp_m_idx_1);

                        short created_substate_idx_0 = check_existence_and_bisect(new_right_state_creation, creation_comp_m_idx_0);
                        if (created_substate_idx_0 == -1) continue;
                        creation_sign *= negative_one_pow(created_substate_idx_0);
                        new_right_state_creation.insert(new_right_state_creation.begin() + created_substate_idx_0, creation_comp_m_idx_0);

                        if (left_state != new_right_state_creation) continue;

                        Key6 creation_key = {
                            indices.orbital_idx_to_j_map[creation_orb_idx_0],
                            indices.composite_m_idx_to_m_map[creation_comp_m_idx_0],
                            indices.orbital_idx_to_j_map[creation_orb_idx_1],
                            indices.composite_m_idx_to_m_map[creation_comp_m_idx_1],
                            j_coupled,
                            m_coupled
                        };

                        double cg_creation = clebsch_gordan.at(creation_key);
                        // if (cg_creation == 0) continue;  // Might be faster to just multiply with 0 instead of checking.

                        creation_res += creation_sign*cg_creation;
                    }
                }
                twobody_res += annihilation_norm*creation_norm*tbme*creation_res*annihilation_sign*cg_annihilation;
            }
        }
    }
    return twobody_res;
}

void create_hamiltonian(const Interaction& interaction)
{
    const Indices indices = generate_indices(interaction);
    const std::vector<std::vector<unsigned short>> basis_states = calculate_m_basis_states(interaction);
    const unsigned int m_dim = basis_states.size();

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> H;
    H.resize(m_dim, m_dim);
    H.setZero();

    print_vector(basis_states);

    auto start = timer();
    for (unsigned int left_idx = 0; left_idx < m_dim; left_idx++)
    {
        for (unsigned int right_idx = left_idx; right_idx < m_dim; right_idx++)
        {
            /*
            Calculate only the upper triangle of the Hamiltonian matrix.
            H is hermitian so we dont have to calculate both triangles.
            */
            H(left_idx, right_idx) += calculate_onebody_matrix_element(
                interaction,
                indices,
                basis_states[left_idx],
                basis_states[right_idx]
            );
        }
    }
    timer(start, "calculate_onebody_matrix_element");
    start = timer();
    // for (int left_idx = 0; left_idx < m_dim; left_idx++)
    for (int left_idx : tq::trange(m_dim))
    {
        for (int right_idx = left_idx; right_idx < m_dim; right_idx++)
        {
            H(left_idx, right_idx) += calculate_twobody_matrix_element(
                interaction,
                indices,
                basis_states[left_idx],
                basis_states[right_idx]
            );
        }
    }
    cout << endl;
    timer(start, "calculate_twobody_matrix_element");
    complete_hermitian_matrix(H);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
    es.compute(H);

    // // for (auto val : es.eigenvalues())
    // // {
    // //     cout << val << endl;
    // // }
    cout << "The eigenvalues of A are: " << es.eigenvalues().transpose() << endl;
    print("m_dim", m_dim);
    return;
}