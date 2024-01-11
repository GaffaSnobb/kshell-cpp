#ifndef DATA_STRUCTURES_HPP
#define DATA_STRUCTURES_HPP

#include <vector>
#include <unordered_map>

struct Key5
{
    /*
    Custom struct to be used as a key in a std::unordered_map.
    */
    unsigned short a, b, c, d, e;

    bool operator==(const Key5& other) const
    {
        return a == other.a && b == other.b && c == other.c && d == other.d && e == other.e;
    }
};

struct Key6
{
    /*
    Custom struct to be used as a key in a std::unordered_map.
    */
    short a, b, c, d, e, f;

    bool operator==(const Key6& other) const
    {
        return a == other.a && b == other.b && c == other.c && d == other.d && e == other.e && f == other.f;
    }
};

namespace std
{
    template <>
    struct hash<Key5>
    {
        std::size_t operator()(const Key5& k) const
        {   // Compute individual hash values for two data members and combine them using XOR and bit shifting
            return ((std::hash<unsigned short>()(k.a)
                ^ (std::hash<unsigned short>()(k.b) << 1)) >> 1)
                ^ (std::hash<unsigned short>()(k.c)
                ^ (std::hash<unsigned short>()(k.d) << 1))
                ^ (std::hash<unsigned short>()(k.e) << 1);
        }
    };

    template <>
    struct hash<Key6>
    {
        std::size_t operator()(const Key6& k) const
        {
            return ((std::hash<short>()(k.a)
                ^ (std::hash<short>()(k.b) << 1)) >> 1)
                ^ (std::hash<short>()(k.c)
                ^ (std::hash<short>()(k.d) << 1))
                ^ (std::hash<short>()(k.e)
                ^ (std::hash<short>()(k.f) << 1));
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
    const std::vector<short> composite_m_idx_to_m_map;
    const std::vector<unsigned short> orbital_idx_to_j_map;
    const std::vector<std::vector<unsigned short>> orbital_idx_to_composite_m_idx_map;
    const std::vector<unsigned short> creation_orb_indices_0;
    const std::vector<unsigned short> creation_orb_indices_1;
    const std::vector<unsigned short> annihilation_orb_indices_0;
    const std::vector<unsigned short> annihilation_orb_indices_1;
    const std::vector<unsigned short> j_coupled;
    const std::vector<short> m_coupled;
    const std::vector<double> tbme;

    Indices(
        std::vector<short> composite_m_idx_to_m_map_,
        std::vector<unsigned short> orbital_idx_to_j_map_,
        std::vector<std::vector<unsigned short>> orbital_idx_to_composite_m_idx_map_,
        std::vector<unsigned short> creation_orb_indices_0_,
        std::vector<unsigned short> creation_orb_indices_1_,
        std::vector<unsigned short> annihilation_orb_indices_0_,
        std::vector<unsigned short> annihilation_orb_indices_1_,
        std::vector<unsigned short> j_coupled_,
        std::vector<short> m_coupled_,
        std::vector<double> tbme_
    ) :
    composite_m_idx_to_m_map(composite_m_idx_to_m_map_),
    orbital_idx_to_j_map(orbital_idx_to_j_map_),
    orbital_idx_to_composite_m_idx_map(orbital_idx_to_composite_m_idx_map_),
    creation_orb_indices_0(creation_orb_indices_0_),
    creation_orb_indices_1(creation_orb_indices_1_),
    annihilation_orb_indices_0(annihilation_orb_indices_0_),
    annihilation_orb_indices_1(annihilation_orb_indices_1_),
    j_coupled(j_coupled_),
    m_coupled(m_coupled_),
    tbme(tbme_) {}
};

struct OrbitalParameters
{
    const unsigned short n;   // The "principal quantum number".
    const unsigned short l;   // Orbital angular momentum.
    const unsigned short j;   // Total angular momentum.
    const unsigned short degeneracy;
    const short tz;           // Isospin.
    const std::vector<short> jz;  // All possible z projections of the total angular momentum vector.

    OrbitalParameters(
        unsigned short n_,
        unsigned short l_,
        unsigned short j_,
        unsigned short degeneracy_,
        short tz_,
        std::vector<short> jz_
    ) : n(n_), l(l_), j(j_), degeneracy(degeneracy_), tz(tz_), jz(jz_) {}
};

struct ModelSpace
{
    const std::vector<OrbitalParameters> orbitals;
    const std::vector<short> all_jz_values;
    const unsigned short n_valence_nucleons;
    const unsigned short n_orbitals;

    ModelSpace(
        std::vector<OrbitalParameters> orbitals_,
        std::vector<short> all_jz_values_,
        unsigned short n_valence_nucleons_,
        unsigned short n_orbitals_
    ) : orbitals(orbitals_), all_jz_values(all_jz_values_), n_valence_nucleons(n_valence_nucleons_), n_orbitals(n_orbitals_) {}
};

struct Interaction
{
    const unsigned short tbme_mass_dependence_method;
    const unsigned short n_core_protons;
    const unsigned short n_core_neutrons;
    const unsigned short n_core_nucleons;
    const double tbme_mass_dependence_exponent;
    const double tbme_mass_dependence_denominator;
    const std::vector<double> spe;    // Single-particle energies. The orbital to which they belong is the same as the index of the SPE.
    const ModelSpace model_space;
    const ModelSpace model_space_protons;
    const ModelSpace model_space_neutrons;
    const std::unordered_map<Key5, double> tbme_map;
    const std::vector<Key5> tbme_keys;

    Interaction(
        unsigned short tbme_mass_dependence_method_,
        unsigned short n_core_protons_,
        unsigned short n_core_neutrons_,
        unsigned short n_core_nucleons_,
        double tbme_mass_dependence_exponent_,
        double tbme_mass_dependence_denominator_,
        std::vector<double> spe_,
        ModelSpace model_space_,
        ModelSpace model_space_protons_,
        ModelSpace model_space_neutrons_,
        std::unordered_map<Key5, double> tbme_map_,
        std::vector<Key5> tbme_keys_
    ) :
    tbme_mass_dependence_method(tbme_mass_dependence_method_),
    n_core_protons(n_core_protons_),
    n_core_neutrons(n_core_neutrons_),
    n_core_nucleons(n_core_nucleons_),
    tbme_mass_dependence_exponent(tbme_mass_dependence_exponent_),
    tbme_mass_dependence_denominator(tbme_mass_dependence_denominator_),
    spe(spe_),
    model_space(model_space_),
    model_space_protons(model_space_protons_),
    model_space_neutrons(model_space_neutrons_),
    tbme_map(tbme_map_),
    tbme_keys(tbme_keys_) {}
};

#endif // DATA_STRUCTURES_HPP