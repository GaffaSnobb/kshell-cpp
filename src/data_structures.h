#ifndef DATA_STRUCTURES_H
#define DATA_STRUCTURES_H

struct Key
{
    /*
    Custom struct to be used as a key in a std::unordered_map.
    */
    unsigned short a, b, c, d, e;

    bool operator==(const Key& other) const
    {
        return a == other.a && b == other.b && c == other.c && d == other.d && e == other.e;
    }
};

namespace std
{
    template <>
    struct hash<Key>
    {
        std::size_t operator()(const Key& k) const
        {   // Compute individual hash values for two data members and combine them using XOR and bit shifting
            return ((std::hash<unsigned short>()(k.a) 
                 ^ (std::hash<unsigned short>()(k.b) << 1)) >> 1)
                 ^ (std::hash<unsigned short>()(k.c) 
                 ^ (std::hash<unsigned short>()(k.d) << 1))
                 ^ (std::hash<unsigned short>()(k.e) << 1);
        }
    };
}

struct OrbitalParameters
{
    const unsigned short n;   // The "principal quantum number".
    const unsigned short l;   // Orbital angular momentum.
    const unsigned short j;   // Total angular momentum.
    const short tz;           // Isospin.
    const std::vector<short> jz;  // All possible z projections of the total angular momentum vector.

    OrbitalParameters(
        unsigned short n_,
        unsigned short l_,
        unsigned short j_,
        short tz_,
        std::vector<short> jz_
    ) : n(n_), l(l_), j(j_), tz(tz_), jz(jz_) {}
};

struct ModelSpace
{
    const std::vector<OrbitalParameters> orbitals;
    // const unsigned short n_valence_nucleons;

    ModelSpace(
        std::vector<OrbitalParameters> orbitals_
        // unsigned short n_valence_nucleons_
    ) : orbitals(orbitals_) {}

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
    const std::unordered_map<Key, double> tbme_map;
    const std::vector<Key> tbme_keys;

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
        std::unordered_map<Key, double> tbme_map_,
        std::vector<Key> tbme_keys_
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

#endif // DATA_STRUCTURES_H