#include <vector>
#include <stdint.h>

#include "../src/data_structures.hpp"
#include "../src/lanczos.hpp"

void lanczos_test()
{   
    constexpr size_t m_dim = 5;
    std::vector<int16_t> dummy_jz;
    OrbitalParameters dummy_orbital_parameter(
        0, // uint16_t n_,
        0, // uint16_t l_,
        0, // uint16_t j_,
        0, // uint16_t degeneracy_,
        0, // int16_t tz_,
        dummy_jz // std::vector<int16_t> jz_
    );

    std::vector<OrbitalParameters> dummy_orbitals;
    dummy_orbitals.push_back(dummy_orbital_parameter);

    std::vector<int16_t> all_jz_values;
    ModelSpace dummy_model_space(
        dummy_orbitals, // std::vector<OrbitalParameters> orbitals_,
        all_jz_values, // std::vector<int16_t> all_jz_values_,
        0, // uint16_t n_valence_nucleons_,
        0 // uint16_t n_orbitals_
    );

    std::vector<double> dummy_spe;
    std::unordered_map<Key5, double> dummy_tbme_map;
    std::vector<Key5> dummy_tbme_keys;
    std::vector<uint64_t> dummy_basis_states;
    for (size_t i = 0; i < m_dim; i++) dummy_basis_states.push_back(i);    // Doesnt matter what the contents are, only the length is important.
    Interaction dummy_interaction(
        0, //uint16_t tbme_mass_dependence_method_,
        0, //uint16_t n_core_protons_,
        0, //uint16_t n_core_neutrons_,
        0, //uint16_t n_core_nucleons_,
        0, //double tbme_mass_dependence_exponent_,
        0, //double tbme_mass_dependence_denominator_,
        dummy_spe, //std::vector<double> spe_,
        nullptr, //double* spe_array_,
        dummy_model_space, // ModelSpace model_space_,
        dummy_model_space, // ModelSpace model_space_protons_,
        dummy_model_space, // ModelSpace model_space_neutrons_,
        dummy_tbme_map, // std::unordered_map<Key5, double> tbme_map_,
        dummy_tbme_keys, // std::vector<Key5> tbme_keys_,
        dummy_basis_states // std::vector<uint64_t> basis_states_
    );
    
    double H[5*5] = {
        -2.65526241, -1.73658919,  1.05043732, -1.35836282, -0.60596862,
        -1.73658919, -1.04257302, -0.38122495,  0.67562902, -0.56439201,
        1.05043732, -0.38122495,  3.95116467, -0.66926132,  0.58965748,
        -1.35836282,  0.67562902, -0.66926132, -0.28581319, -0.37952717,
        -0.60596862, -0.56439201,  0.58965748, -0.37952717, -1.22605036
    };

    lanczos::lanczos(dummy_interaction, H);
}