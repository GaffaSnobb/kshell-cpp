#include <vector>
#include <iostream>
#include <algorithm>
#include "data_structures.hpp"
#include "tools.hpp"
#include "basis.hpp"
#include "../external/eigen-3.4.0/Eigen/Dense"

using std::cout;
using std::endl;

const Indices generate_indices(const Interaction& interaction)
{
    /*
    Generate values for all the attributes of the Indices data
    structure.
    */
    std::vector<unsigned short> orbital_idx_to_j_map;
    std::vector<std::vector<unsigned short>> orbital_idx_to_composite_m_idx_map;
    unsigned short previous_degeneracy = 0;

    for (int i = 0; i < interaction.model_space.orbitals.size(); i++)
    {   
        unsigned short j = interaction.model_space.orbitals[i].j;
        unsigned short current_degeneracy = interaction.model_space.orbitals[i].degeneracy;

        orbital_idx_to_j_map.push_back(j);
        orbital_idx_to_composite_m_idx_map.push_back(
            range<unsigned short>(
                previous_degeneracy,
                previous_degeneracy + current_degeneracy,
                1
            )
        );
        previous_degeneracy = current_degeneracy + previous_degeneracy;
    }
    std::vector<unsigned short> creation_orb_indices_0;
    std::vector<unsigned short> creation_orb_indices_1;
    std::vector<unsigned short> annihilation_orb_indices_0;
    std::vector<unsigned short> annihilation_orb_indices_1;
    std::vector<unsigned short> j_coupled_list;
    std::vector<short> m_coupled_list;
    std::vector<double> tbme_list;

    for (unsigned short creation_orb_idx_0 = 0; creation_orb_idx_0 < interaction.model_space.orbitals.size(); creation_orb_idx_0++)
    {
        for (unsigned short creation_orb_idx_1 = creation_orb_idx_0; creation_orb_idx_1 < interaction.model_space.orbitals.size(); creation_orb_idx_1++)
        {
            for (unsigned short annihilation_orb_idx_0 = 0; annihilation_orb_idx_0 < interaction.model_space.orbitals.size(); annihilation_orb_idx_0++)
            {
                for (unsigned short annihilation_orb_idx_1 = annihilation_orb_idx_0; annihilation_orb_idx_1 < interaction.model_space.orbitals.size(); annihilation_orb_idx_1++)
                {
                    unsigned short j_min = std::max(
                        std::abs(orbital_idx_to_j_map[creation_orb_idx_0] - orbital_idx_to_j_map[creation_orb_idx_1]),
                        std::abs(orbital_idx_to_j_map[annihilation_orb_idx_0] - orbital_idx_to_j_map[annihilation_orb_idx_1])
                    );
                    unsigned short j_max = std::min(
                        orbital_idx_to_j_map[creation_orb_idx_0] + orbital_idx_to_j_map[creation_orb_idx_1],
                        orbital_idx_to_j_map[annihilation_orb_idx_0] + orbital_idx_to_j_map[annihilation_orb_idx_1]
                    );
                    for (unsigned short j_coupled = j_min; j_coupled <= j_max; j_coupled += 2)
                    {
                        /*
                        j_coupled is the total angular momentum to which 
                        (creation_orb_idx_0, creation_orb_idx_1) and
                        (annihilation_orb_idx_0, annihilation_orb_idx_1)
                        couple. Follows the standard angular momentum
                        coupling rules:

                            j = | j1 - j2 |, | j1 - j2 | + 1, ..., j1 + j2

                        but we have to respect the allowed range for
                        both pairs of total angular momentum values so
                        that j_coupled is contained in both ranges.

                        Step length of 2 because all angular momentum
                        values are multiplied by 2 to avoid fractions.
                        + 2 so that the end point is included.
                        */
                        double tbme_tmp;
                        try
                        {   
                            Key key = {
                                creation_orb_idx_0,
                                creation_orb_idx_1,
                                annihilation_orb_idx_0,
                                annihilation_orb_idx_1,
                                j_coupled
                            };
                            tbme_tmp = interaction.tbme_map.at(key);
                        }
                        catch (const std::out_of_range& e)
                        {
                            /*
                            If the current interaction file is not
                            defining any two-body matrix elements
                            for this choice of orbitals and coupled
                            j then the result is 0.
                            */
                            continue;
                        }
                        for (short m_coupled = -j_coupled; m_coupled <= j_coupled; m_coupled += 2)
                        {
                            /*
                            m_coupled is simply the z component of the
                            coupled total angular momentum, j_coupled.
                            */
                            creation_orb_indices_0.push_back(creation_orb_idx_0);
                            creation_orb_indices_1.push_back(creation_orb_idx_1);
                            annihilation_orb_indices_0.push_back(annihilation_orb_idx_0);
                            annihilation_orb_indices_1.push_back(annihilation_orb_idx_1);
                            j_coupled_list.push_back(j_coupled);
                            m_coupled_list.push_back(m_coupled);
                            tbme_list.push_back(tbme_tmp);
                        }
                    }
                }
            }
        }
    }
    const Indices indices(
        interaction.model_space.all_jz_values,
        orbital_idx_to_j_map,
        orbital_idx_to_composite_m_idx_map,
        creation_orb_indices_0,
        creation_orb_indices_1,
        annihilation_orb_indices_0,
        annihilation_orb_indices_1,
        j_coupled_list,
        m_coupled_list,
        tbme_list
    );
    return indices;
}

void create_hamiltonian(const Interaction& interaction)
{
    const Indices indices = generate_indices(interaction);
    const std::vector<std::vector<unsigned short>> basis_states = calculate_m_basis_states(interaction);
    const unsigned int m_dim = basis_states.size();

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> H;
    H.resize(m_dim, m_dim);

    for (int left_idx = 0; left_idx < m_dim; left_idx++)
    {
        for (int right_idx = left_idx; right_idx < m_dim; right_idx++)
        {
            
        }
    }

    return;
}