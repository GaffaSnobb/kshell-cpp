#include <vector>
#include <iostream>
#include <bitset>
#include <stdint.h>

#include "data_structures.hpp"
#include "parameters.hpp"
#include "tools.hpp"
#include "bit_manipulation_tools.hpp"
#include "../external/cppitertools/combinations.hpp"

using std::cout;
using std::endl;

namespace basis
{
std::vector<std::vector<uint16_t>> calculate_m_basis_states_vector_representation(const Interaction& interaction)
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
    uint16_t n_proton_m_substates = interaction.model_space_protons.all_jz_values.size();
    uint16_t n_neutron_m_substates = interaction.model_space_neutrons.all_jz_values.size();

    std::vector<uint16_t> proton_indices = range<uint16_t>( // All possible proton m substate indices.
        0, n_proton_m_substates, 1
    );

    auto proton_index_combinations = iter::combinations(
        proton_indices, // Iterable.
        interaction.model_space_protons.n_valence_nucleons  // Number of elements to choose from the iterable.
    );

    std::vector<uint16_t> neutron_indices = range<uint16_t>(    // All possible neutron m substate indices.
        n_proton_m_substates, n_proton_m_substates + n_neutron_m_substates, 1
    );

    auto neutron_index_combinations = iter::combinations(
        neutron_indices,
        interaction.model_space_neutrons.n_valence_nucleons
    );

    uint16_t M_target;
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

    std::vector<std::vector<uint16_t>> basis_states;

    for (auto&& proton_combination : proton_index_combinations)
    {
        for (auto&& neutron_combination : neutron_index_combinations)
        {   
            std::vector<uint16_t> combined_index_combinations;
            int16_t M = 0;    // For summing the m value of each nucleon.
            for (uint16_t p : proton_combination)
            {
                combined_index_combinations.push_back(p);
                M += interaction.model_space.all_jz_values[p];
            }
            for (uint16_t n : neutron_combination)
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

std::vector<std::bitset<n_bits_bitset>> calculate_m_basis_states_bitset_representation(const Interaction& interaction)
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
    uint16_t n_proton_m_substates = interaction.model_space_protons.all_jz_values.size();
    uint16_t n_neutron_m_substates = interaction.model_space_neutrons.all_jz_values.size();

    std::vector<uint16_t> proton_indices = range<uint16_t>( // All possible proton m substate indices.
        0, n_proton_m_substates, 1
    );

    auto proton_index_combinations = iter::combinations(
        proton_indices, // Iterable.
        interaction.model_space_protons.n_valence_nucleons  // Number of elements to choose from the iterable.
    );

    std::vector<uint16_t> neutron_indices = range<uint16_t>(    // All possible neutron m substate indices.
        n_proton_m_substates, n_proton_m_substates + n_neutron_m_substates, 1
    );

    auto neutron_index_combinations = iter::combinations(
        neutron_indices,
        interaction.model_space_neutrons.n_valence_nucleons
    );

    uint16_t M_target;
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

    std::vector<std::bitset<n_bits_bitset>> basis_states;

    for (auto&& proton_combination : proton_index_combinations)
    {
        for (auto&& neutron_combination : neutron_index_combinations)
        {   
            std::bitset<n_bits_bitset> basis_state;
            int16_t M = 0;    // For summing the m value of each nucleon.
            for (uint16_t p : proton_combination)
            {
                basis_state.set(p);
                M += interaction.model_space.all_jz_values[p];
            }
            for (uint16_t n : neutron_combination)
            {
                basis_state.set(n);
                M += interaction.model_space.all_jz_values[n];
            }
            if (M == M_target)
            {
                /*
                Keep the combination only if the sum of the m values of
                all the nucleons in the combination is equal to the
                target M value.
                */
                basis_states.push_back(basis_state);
            }
        }
    }
    return basis_states;
}

const std::vector<uint64_t> calculate_m_basis_states_primitive_bit_representation(
    const uint16_t n_proton_m_substates,
    const uint16_t n_neutron_m_substates,
    const uint16_t n_valence_protons,
    const uint16_t n_valence_neutrons,
    const std::vector<int16_t>& all_jz_values
)
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
    const uint16_t n_valence_nucleons = n_valence_protons + n_valence_neutrons;

    std::vector<uint16_t> proton_indices = range<uint16_t>( // All possible proton m substate indices.
        0, n_proton_m_substates, 1
    );

    auto proton_index_combinations = iter::combinations(
        proton_indices, // Iterable.
        n_valence_protons  // Number of elements to choose from the iterable.
    );

    std::vector<uint16_t> neutron_indices = range<uint16_t>(    // All possible neutron m substate indices.
        n_proton_m_substates, n_proton_m_substates + n_neutron_m_substates, 1
    );

    auto neutron_index_combinations = iter::combinations(
        neutron_indices,
        n_valence_neutrons
    );

    uint16_t M_target;
    if (n_valence_nucleons%2 == 0)
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

    std::vector<uint64_t> basis_states;

    for (auto&& proton_combination : proton_index_combinations)
    {
        for (auto&& neutron_combination : neutron_index_combinations)
        {   
            uint64_t basis_state = 0;
            int16_t M = 0;    // For summing the m value of each nucleon.
            for (uint16_t p : proton_combination)
            {
                bittools::set_bit(basis_state, p);
                M += all_jz_values[p];
            }
            for (uint16_t n : neutron_combination)
            {
                bittools::set_bit(basis_state, n);
                M += all_jz_values[n];
            }
            if (M == M_target)
            {
                /*
                Keep the combination only if the sum of the m values of
                all the nucleons in the combination is equal to the
                target M value.
                */
                basis_states.push_back(basis_state);
            }
        }
    }
    return basis_states;
}
}