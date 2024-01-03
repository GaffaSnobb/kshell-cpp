#include <iostream>
#include "data_structures.hpp"
#include "loaders.hpp"
#include "tools.hpp"

using std::cout;
using std::endl;

int main()
{
    std::string interaction_filename = "../snt/w.snt";
    const Interaction interaction = load_interaction(interaction_filename);
    cout << interaction.model_space_protons.orbitals[0].j << endl;

    print(interaction.model_space.orbitals);
    cout << interaction.tbme_mass_dependence_denominator << endl;

    const Key key = {0, 0, 0, 0, 0};
    cout << interaction.tbme_map.at(key) << endl;

    // Key key = {i0, i1, i2, i3, j};
    // tbme_keys.push_back(key);
    // tbme_map[key] = tbme;

    // print("n_proton_orbitals", interaction.model_space_protons.n_orbitals);
    // print("n_neutron_orbitals", interaction.model_space_neutrons.n_orbitals);
    // print("n_core_protons", interaction.n_core_protons);
    // print("n_core_neutrons", interaction.n_core_neutrons);

    // print(interaction.model_space_protons.orbitals);
    // print(interaction.model_space_neutrons.orbitals);
    // print_vector(j_couple);
    return 0;
}