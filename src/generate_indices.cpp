#include "data_structures.hpp"
#include "tools.hpp"

namespace indices
{
const Indices generate_indices(const Interaction& interaction)
{
    /*
    Generate values for all the attributes of the Indices data
    structure.

    Pre-generated indices for the two-body Hamiltonian matrix element
    calculations are used to reduce the computational load.
    */
    auto start = timer();
    std::vector<uint16_t> orbital_idx_to_j_map;
    std::vector<std::vector<uint16_t>> orbital_idx_to_composite_m_idx_map;
    uint16_t* orbital_idx_to_composite_m_idx_map_flattened_indices = new uint16_t[interaction.model_space.orbitals.size()]; // End index of each map section.
    uint16_t previous_degeneracy = 0;
    uint16_t previous_size = 0;
    
    for (int32_t i = 0; i < interaction.model_space.orbitals.size(); i++)
    {   
        uint16_t j = interaction.model_space.orbitals[i].j;
        uint16_t current_degeneracy = interaction.model_space.orbitals[i].degeneracy;

        orbital_idx_to_j_map.push_back(j);
        orbital_idx_to_composite_m_idx_map.push_back(
            range<uint16_t>(
                previous_degeneracy,
                previous_degeneracy + current_degeneracy,
                1
            )
        );
        previous_degeneracy = current_degeneracy + previous_degeneracy;

        /*
        Make a flattened orbital index -> m-substate index array which
        stores only the end index of the current map. Since the
        m-substate indices are in increasing order and incremeted by 1,
        there is no need to explicitly store all the numbers. The
        flattened array is used for passing the index mapping to a
        kernel.
        */
        uint16_t current_size = orbital_idx_to_composite_m_idx_map[i].size();
        orbital_idx_to_composite_m_idx_map_flattened_indices[i] = current_size + previous_size;

        previous_size = previous_size + current_size;
    }

    std::vector<uint16_t> creation_orb_indices_0;
    std::vector<uint16_t> creation_orb_indices_1;
    std::vector<uint16_t> annihilation_orb_indices_0;
    std::vector<uint16_t> annihilation_orb_indices_1;
    std::vector<uint16_t> j_coupled_list;
    std::vector<int16_t> m_coupled_list;
    std::vector<double> tbme_list;

    for (uint16_t creation_orb_idx_0 = 0; creation_orb_idx_0 < interaction.model_space.orbitals.size(); creation_orb_idx_0++)
    {
        for (uint16_t creation_orb_idx_1 = creation_orb_idx_0; creation_orb_idx_1 < interaction.model_space.orbitals.size(); creation_orb_idx_1++)
        {
            for (uint16_t annihilation_orb_idx_0 = 0; annihilation_orb_idx_0 < interaction.model_space.orbitals.size(); annihilation_orb_idx_0++)
            {
                for (uint16_t annihilation_orb_idx_1 = annihilation_orb_idx_0; annihilation_orb_idx_1 < interaction.model_space.orbitals.size(); annihilation_orb_idx_1++)
                {
                    uint16_t j_min = std::max(
                        std::abs(orbital_idx_to_j_map[creation_orb_idx_0] - orbital_idx_to_j_map[creation_orb_idx_1]),
                        std::abs(orbital_idx_to_j_map[annihilation_orb_idx_0] - orbital_idx_to_j_map[annihilation_orb_idx_1])
                    );
                    uint16_t j_max = std::min(
                        orbital_idx_to_j_map[creation_orb_idx_0] + orbital_idx_to_j_map[creation_orb_idx_1],
                        orbital_idx_to_j_map[annihilation_orb_idx_0] + orbital_idx_to_j_map[annihilation_orb_idx_1]
                    );
                    for (uint16_t j_coupled = j_min; j_coupled <= j_max; j_coupled += 2)
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
                            Key5 key = {
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
                        for (int16_t m_coupled = -j_coupled; m_coupled <= j_coupled; m_coupled += 2)
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
        orbital_idx_to_composite_m_idx_map_flattened_indices,
        creation_orb_indices_0,
        creation_orb_indices_1,
        annihilation_orb_indices_0,
        annihilation_orb_indices_1,
        j_coupled_list,
        m_coupled_list,
        tbme_list
    );
    timer(start, "generate_indices");
    return indices;
}
}