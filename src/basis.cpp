#include <vector>
#include <iostream>
#include "data_structures.hpp"
#include "tools.hpp"
#include "../external/cppitertools/combinations.hpp"

using std::cout;
using std::endl;

std::vector<std::vector<unsigned short>> calculate_m_basis_states(const Interaction& interaction)
{
    /*
    Calculate the M-scheme basis states.

    Given a model space and a number of valence nucleons, calculate all
    possible combinations of nucleons whose magnetic substates sum up
    to `M_target`. Consider the following figure which shows the sd
    model space. Note that the ordering of the orbitals is the same
    order as they appear in usda.snt.

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

    This function will iterate over all the possible configurations of
    the valence particles and save the configurations whose m values sum
    up to M_target.
    */
    unsigned short n_proton_m_substates = interaction.model_space_protons.all_jz_values.size();
    unsigned short n_neutron_m_substates = interaction.model_space_neutrons.all_jz_values.size();

    std::vector<unsigned short> proton_indices = range<unsigned short>( // All possible proton m substate indices.
        0, n_proton_m_substates, 1
    );

    auto proton_index_combinations = iter::combinations(
        proton_indices, // Iterable.
        interaction.model_space_protons.n_valence_nucleons  // Number of elements to choose from the iterable.
    );

    std::vector<unsigned short> neutron_indices = range<unsigned short>(    // All possible neutron m substate indices.
        n_proton_m_substates, n_proton_m_substates + n_neutron_m_substates, 1
    );

    auto neutron_index_combinations = iter::combinations(
        neutron_indices,
        interaction.model_space_neutrons.n_valence_nucleons
    );

    unsigned short M_target;
    if (interaction.model_space.n_valence_nucleons%2 == 0)
    {
        /*
        For an even number of valence nucleons, the M = 0 basis states
        are enough to describe all angular momenta.
        */
        M_target = 0;
    }
    else
    {
        /*
        For odd-numbered M = 1/2 is enough, but remember, all angular
        momenta in this code are multiplied by 2.
        */
        M_target = 1;
    }

    std::vector<std::vector<unsigned short>> basis_states;

    for (auto&& proton_combination : proton_index_combinations)
    {
        for (auto&& neutron_combination : neutron_index_combinations)
        {   
            std::vector<unsigned short> combined_index_combinations;
            short M = 0;    // For summing the m value of each nucleon.
            for (unsigned short p : proton_combination)
            {
                combined_index_combinations.push_back(p);
                M += interaction.model_space.all_jz_values[p];
            }
            for (unsigned short n : neutron_combination)
            {
                combined_index_combinations.push_back(n);
                M += interaction.model_space.all_jz_values[n];
            }
            if (M == M_target)
            {
                /*
                Keep the combination only if the sum of the m values of
                all the nucleons in the combination is equal to the
                target M value.
                */
                basis_states.push_back(combined_index_combinations);
            }
        }
    }
    return basis_states;
}