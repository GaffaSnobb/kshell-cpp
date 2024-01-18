#include <iostream>
#include <algorithm>
#include <chrono>
#include <bitset>
#include "data_structures.hpp"
#include "parameters.hpp"
#include "loaders.hpp"
#include "tools.hpp"
#include "hamiltonian.hpp"
#include "../external/cppitertools/combinations.hpp"
#include "../external/eigen-3.4.0/Eigen/Dense"
#include "../external/eigen-3.4.0/Eigen/Eigenvalues"
// #include "../external/boost_1_84_0/boost/container_hash/hash.hpp"

using std::cout;
using std::endl;
// using std::chrono::high_resolution_clock;
// using std::chrono::duration_cast;
// using std::chrono::milliseconds;

int main()
{

    // std::bitset<16> state_0;
    // state_0.set(14);
    // cout << (not 0) << endl;
    // cout << "state_0: " << state_0 << endl;
    // state_0.reset(14);
    // cout << "state_0: " << state_0 << endl;

    // double a = 98.5;
    // cout << a << endl;
    // cout << a*state_0.test(14) << endl;

    // cout << res << endl;
    // state.set(3);
    // if (state.test(4)) cout << "Sub-state 5 is occupied." << endl;
    // state.reset(4);
    // state.flip(9);

    // std::hash<std::bitset<16>> hasher;
    // size_t hash_value = hasher(state_0);

    // cout << hash_value << endl;

    // return 0;

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> H;

    std::vector<int> timing;
    for (int i = 0; i < 10; i++)
    {
        auto start = timer();
        unsigned short n_valence_protons = 2;
        unsigned short n_valence_neutrons = 1;
        std::string interaction_filename = "../snt/w.snt";
        const Interaction interaction = load_interaction(interaction_filename, n_valence_protons, n_valence_neutrons);
        // create_hamiltonian(interaction);
        H = create_hamiltonian_bit_representation(interaction);
        timing.push_back(timer(start, "main").count());
    }
    print_vector(timing);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
    es.compute(H);
    cout << "The eigenvalues of A are: " << es.eigenvalues().transpose() << endl;
    
    print("m_dim", H.rows());
    // print("n_proton_orbitals", interaction.model_space_protons.orbitals.size());
    // print("n_neutron_orbitals", interaction.model_space_neutrons.orbitals.size());
    // print("n_core_protons", interaction.n_core_protons);
    // print("n_core_neutrons", interaction.n_core_neutrons);
    return 0;
}