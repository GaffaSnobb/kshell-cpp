#include <fstream>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <bitset>
#include "data_structures.hpp"
#include "parameters.hpp"
#include "loaders.hpp"
#include "tools.hpp"
#include "hamiltonian.hpp"
#include "hamiltonian_device.hpp"
#include "hamiltonian_bitset_representation.hpp"
#include "../external/cppitertools/combinations.hpp"
#include "../external/eigen-3.4.0/Eigen/Dense"
#include "../external/eigen-3.4.0/Eigen/Eigenvalues"
// #include "../external/boost_1_84_0/boost/container_hash/hash.hpp"

using std::cout;
using std::endl;
// using std::chrono::high_resolution_clock;
// using std::chrono::duration_cast;
// using std::chrono::milliseconds;

int main(int argc, char* argv[])
{
    // unsigned long long state = 0;
    // print_bit_representation(state);
    // state |= (1ULL << 63);
    // print_bit_representation(state);
    // return 0;

    // unsigned short n_valence_protons = 1;
    // unsigned short n_valence_neutrons = 2;
    unsigned short n_valence_protons = std::stoi(argv[1]);
    unsigned short n_valence_neutrons = std::stoi(argv[2]);

    // Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> H_bitset;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> H;

    std::vector<int> timing;
    for (int i = 0; i < 5; i++)
    {
        auto start = timer();
        std::string interaction_filename = "../snt/w.snt";
        const Interaction interaction = loaders::load_interaction(interaction_filename, n_valence_protons, n_valence_neutrons);
        // H_bitset = hamiltonian_bitset::create_hamiltonian_bitset_representation(interaction);
        // H = hamiltonian::create_hamiltonian_primitive_bit_representation(interaction);
        hamiltonian_device::create_hamiltonian_device_dispatcher(interaction);
        timing.push_back(timer(start, "main").count());
    }


    print_vector(timing);
    cout << mean(timing) << endl;
    // cout << (H == H_bitset) << endl;

    // auto start = timer();
    // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
    // es.compute(H);
    // es.eigenvalues();
    // timer(start, "eigensolver");

    // std::ofstream outfile("matrix.txt");
    // if (outfile.is_open())
    // {
    //     outfile << H << std::endl;
    //     outfile.close();
    // }

    // cout << "The eigenvalues of A are: " << es.eigenvalues().transpose() << endl;
    
    return 0;
}