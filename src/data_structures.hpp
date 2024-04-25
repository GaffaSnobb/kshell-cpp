#ifndef DATA_STRUCTURES_HPP
#define DATA_STRUCTURES_HPP

#include <vector>
#include <unordered_map>
#include <stdint.h>

struct Key5
{
    /*
    Custom struct to be used as a key in a std::unordered_map.
    */
    uint16_t a, b, c, d, e;

    bool operator==(const Key5& other) const
    {
        return a == other.a && b == other.b && c == other.c && d == other.d && e == other.e;
    }
};

struct Key6
{
    /*
    Custom struct to be used as a key in a std::unordered_map for the
    Clebsh-Gordan coefficients.
    */
    uint16_t j1;
    int16_t m1;
    uint16_t j2;
    int16_t m2;
    uint16_t J;
    int16_t M;

    bool operator==(const Key6& other) const
    {
        return j1 == other.j1 && m1 == other.m1 && j2 == other.j2 && m2 == other.m2 && J == other.J && M == other.M;
    }
};

namespace std
{
    template <>
    struct hash<Key5>
    {
        std::size_t operator()(const Key5& k) const
        {   // Compute individual hash values for two data members and combine them using XOR and bit shifting
            return ((std::hash<uint16_t>()(k.a)
                ^ (std::hash<uint16_t>()(k.b) << 1)) >> 1)
                ^ (std::hash<uint16_t>()(k.c)
                ^ (std::hash<uint16_t>()(k.d) << 1))
                ^ (std::hash<uint16_t>()(k.e) << 1);
        }
    };
    // template<>
    // struct hash<Key5>
    // {
    //     std::size_t operator()(const Key5& k) const
    //     {
    //         std::size_t h1 = std::hash<uint16_t>()(k.a);
    //         std::size_t h2 = std::hash<uint16_t>()(k.b);
    //         std::size_t h3 = std::hash<uint16_t>()(k.c);
    //         std::size_t h4 = std::hash<uint16_t>()(k.d);
    //         std::size_t h5 = std::hash<uint16_t>()(k.e);

    //         // Combine hash values
    //         return ((((h1 ^ (h2 << 1)) >> 1) ^ (h3 << 1)) >> 1) ^ ((h4 ^ (h5 << 1)) << 1);
    //     }
    // };

    // template<>
    // struct hash<Key5>
    // {
    //     std::size_t operator()(const Key5& k) const
    //     {
    //         std::size_t seed = 0;
    //         auto hash_combine = [&seed](std::size_t h)
    //         {
    //             h ^= seed + 0x9e3779b9 + (h << 6) + (h >> 2);
    //             seed ^= h;
    //         };

    //         hash_combine(std::hash<uint16_t>()(k.a));
    //         hash_combine(std::hash<uint16_t>()(k.b));
    //         hash_combine(std::hash<uint16_t>()(k.c));
    //         hash_combine(std::hash<uint16_t>()(k.d));
    //         hash_combine(std::hash<uint16_t>()(k.e));

    //         return seed;
    //     }
    // };

    template <>
    struct hash<Key6>
    {
        size_t operator()(const Key6& k) const
        {   // Hash individual members and combine them
            size_t hash_j1 = std::hash<uint16_t>()(k.j1);
            size_t hash_m1 = std::hash<int16_t>()(k.m1);
            size_t hash_j2 = std::hash<uint16_t>()(k.j2);
            size_t hash_m2 = std::hash<int16_t>()(k.m2);
            size_t hash_J = std::hash<uint16_t>()(k.J);
            size_t hash_M = std::hash<int16_t>()(k.M);

            // Combine these hashes together
            return ((((hash_j1 ^ (hash_m1 << 1)) >> 1) ^ (hash_j2 << 1)) >> 1) ^
                   ((((hash_m2 ^ (hash_J << 1)) >> 1) ^ (hash_M << 1)) >> 1);
        }
    };
}

struct Indices
{
    /*
    Attributes
    ----------
    orbital_idx_to_m_idx_map : list[tuple[int, ...]]
        Map the index of an orbital in a model space to the magnetic
        substate indices of that orbital. How the orbitals are ordered
        is defined by the interaction file (.snt). For example, in w.snt
        the indices are as follows:

            0 = p 0d_3/2
            1 = p 0d_5/2
            2 = p 1s_1/2
            3 = n 0d_3/2
            4 = n 0d_5/2
            5 = n 1s_1/2

        Please note that I have changed the indices to start at 0 while
        in the actual interaction file the indices start at 1. This
        means that orbital 0 (d3/2) has four magnetic substates and thus

            orbital_idx_to_m_idx_map[0] = (0, 1, 2, 3)
            orbital_idx_to_m_idx_map[1] = (0, 1, 2, 3, 4, 5)
            orbital_idx_to_m_idx_map[2] = (0, 1)
            ...

    orbital_m_pair_to_composite_m_idx_map : dict[tuple[int, int], int]
        Map a pair of (orbital index, m substate index) to a composite
        m substate index. Please consider the following scheme sd model
        space scheme:

                          10    11        
                        -  O  -  O  -             s1/2: 2
                         -1/2   1/2

               4     5     6     7     8     9    
            -  O  -  O  -  O  -  O  -  O  -  O    d5/2: 1
             -5/2  -3/2  -1/2   1/2   3/2   5/2

                     0     1     2     3
                  -  O  -  O  -  O  -  O  -       d3/2: 0
                   -3/2  -1/2   1/2   3/2

        The orbitals are indexed 0, 1, 2. The m substates are indexed by
        a composite / cumulative / whatever you wanna call it, m
        substate index, 0, 1, ..., 11. This map translates the orbital
        index and the orbital's 'local' m substate index into the
        composite m substate index. The 'local' m indices are for d3/2,
        d5/2, s1/2 respectively,

            [0, 1, 2, 3],
            [0, 1, 2, 3, 4, 5],
            [0, 1].

        orbital_idx_to_j_map : list[int]
            Map the index of an orbital to its 2*j value. Considering
            the above scheme,

                orbital_idx_to_j_map = [3, 5, 1]

        m_composite_idx_to_m_map : list[int]
            Map a composite m index to the 2*m value it corresponds to.
            Considering the above scheme,

                m_composite_idx_to_m_map = [
                    -3, -1, +1, +3, -5, -3, -1, +1, +3, +5, -1, +1,
                ]

        orbital_idx_to_comp_m_idx_map : list[tuple[int, ...]]
            Translate the orbital indices to the composite m indices of
            the magnetic substates. For example, proton d3/2 has orbital
            index 0 and composite m substate indices 0, 1, 2, 3.

    */
    const std::vector<int16_t> composite_m_idx_to_m_map;
    const std::vector<uint16_t> orbital_idx_to_j_map;
    const std::vector<std::vector<uint16_t>> orbital_idx_to_composite_m_idx_map;
    const uint16_t* orbital_idx_to_composite_m_idx_map_flattened_indices; // End index of each section.
    const std::vector<uint16_t> creation_orb_indices_0;
    const std::vector<uint16_t> creation_orb_indices_1;
    const std::vector<uint16_t> annihilation_orb_indices_0;
    const std::vector<uint16_t> annihilation_orb_indices_1;
    const std::vector<uint16_t> j_coupled;
    const std::vector<int16_t> m_coupled;
    const std::vector<double> tbme;

    Indices(
        std::vector<int16_t> composite_m_idx_to_m_map_,
        std::vector<uint16_t> orbital_idx_to_j_map_,
        std::vector<std::vector<uint16_t>> orbital_idx_to_composite_m_idx_map_,
        uint16_t* orbital_idx_to_composite_m_idx_map_flattened_indices_,
        std::vector<uint16_t> creation_orb_indices_0_,
        std::vector<uint16_t> creation_orb_indices_1_,
        std::vector<uint16_t> annihilation_orb_indices_0_,
        std::vector<uint16_t> annihilation_orb_indices_1_,
        std::vector<uint16_t> j_coupled_,
        std::vector<int16_t> m_coupled_,
        std::vector<double> tbme_
    ) :
    composite_m_idx_to_m_map(composite_m_idx_to_m_map_),
    orbital_idx_to_j_map(orbital_idx_to_j_map_),
    orbital_idx_to_composite_m_idx_map(orbital_idx_to_composite_m_idx_map_),
    orbital_idx_to_composite_m_idx_map_flattened_indices(orbital_idx_to_composite_m_idx_map_flattened_indices_),
    creation_orb_indices_0(creation_orb_indices_0_),
    creation_orb_indices_1(creation_orb_indices_1_),
    annihilation_orb_indices_0(annihilation_orb_indices_0_),
    annihilation_orb_indices_1(annihilation_orb_indices_1_),
    j_coupled(j_coupled_),
    m_coupled(m_coupled_),
    tbme(tbme_) {}

    ~Indices()
    {
        delete[] orbital_idx_to_composite_m_idx_map_flattened_indices;
    }
};  // Indices

struct OrbitalParameters
{
    const uint16_t n;   // The "principal quantum number".
    const uint16_t l;   // Orbital angular momentum.
    const uint16_t j;   // Total angular momentum.
    const uint16_t degeneracy;
    const int16_t tz;           // Isospin.
    const std::vector<int16_t> jz;  // All possible z projections of the total angular momentum vector.

    OrbitalParameters(
        uint16_t n_,
        uint16_t l_,
        uint16_t j_,
        uint16_t degeneracy_,
        int16_t tz_,
        std::vector<int16_t> jz_
    ) : n(n_), l(l_), j(j_), degeneracy(degeneracy_), tz(tz_), jz(jz_) {}
};

struct ModelSpace
{
    const std::vector<OrbitalParameters> orbitals;
    const std::vector<int16_t> all_jz_values;
    const uint16_t n_valence_nucleons;
    const uint16_t n_orbitals;

    ModelSpace(
        std::vector<OrbitalParameters> orbitals_,
        std::vector<int16_t> all_jz_values_,
        uint16_t n_valence_nucleons_,
        uint16_t n_orbitals_
    ) : orbitals(orbitals_), all_jz_values(all_jz_values_), n_valence_nucleons(n_valence_nucleons_), n_orbitals(n_orbitals_) {}
};

struct Interaction
{
    const uint16_t tbme_mass_dependence_method;
    const uint16_t n_core_protons;
    const uint16_t n_core_neutrons;
    const uint16_t n_core_nucleons;
    const double tbme_mass_dependence_exponent;
    const double tbme_mass_dependence_denominator;
    const std::vector<double> spe;    // Single-particle energies. The orbital to which they belong is the same as the index of the SPE.
    const double* spe_array;
    const ModelSpace model_space;
    const ModelSpace model_space_protons;
    const ModelSpace model_space_neutrons;
    const std::unordered_map<Key5, double> tbme_map;
    const std::vector<Key5> tbme_keys;
    const std::vector<uint64_t> basis_states;

    Interaction(
        uint16_t tbme_mass_dependence_method_,
        uint16_t n_core_protons_,
        uint16_t n_core_neutrons_,
        uint16_t n_core_nucleons_,
        double tbme_mass_dependence_exponent_,
        double tbme_mass_dependence_denominator_,
        std::vector<double> spe_,
        double* spe_array_,
        ModelSpace model_space_,
        ModelSpace model_space_protons_,
        ModelSpace model_space_neutrons_,
        std::unordered_map<Key5, double> tbme_map_,
        std::vector<Key5> tbme_keys_,
        std::vector<uint64_t> basis_states_
    ) :
    tbme_mass_dependence_method(tbme_mass_dependence_method_),
    n_core_protons(n_core_protons_),
    n_core_neutrons(n_core_neutrons_),
    n_core_nucleons(n_core_nucleons_),
    tbme_mass_dependence_exponent(tbme_mass_dependence_exponent_),
    tbme_mass_dependence_denominator(tbme_mass_dependence_denominator_),
    spe(spe_),
    spe_array(spe_array_),
    model_space(model_space_),
    model_space_protons(model_space_protons_),
    model_space_neutrons(model_space_neutrons_),
    tbme_map(tbme_map_),
    tbme_keys(tbme_keys_),
    basis_states(basis_states_) {}

    ~Interaction()
    {
        delete[] spe_array;
    }
};

#endif // DATA_STRUCTURES_HPP