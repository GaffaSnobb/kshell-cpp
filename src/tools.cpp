#include <iostream>
#include <vector>
#include "tools.hpp"
#include "data_structures.hpp"

using std::cout;
using std::endl;

void print(std::vector<OrbitalParameters> orbitals)
{
    for (const OrbitalParameters& orb : orbitals)
    {
        print("n", orb.n);
        print("l", orb.l);
        print("j", orb.j);
        print("tz", orb.tz);
        print_vector("jz", orb.jz);
        cout << endl;
    }
}