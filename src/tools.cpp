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

short index(const std::vector<unsigned short>& vec, const unsigned short value)
{
    /*
    Returns the index of the first occurence of `value` in `vec`. If
    `value` is not found, -1 is returned. Assumes that `vec` is short
    enough to be indexed by a short.
    
    Parameters
    ----------
    vec : std::vector<unsigned short>
        The vector to search in.
    
    value : unsigned short
        The value to search for.

    Returns
    -------
    short
        The index of the first occurence of `value` in `vec`. If
        `value` is not found, -1 is returned.
    */
    for (short res = 0; res < static_cast<short>(vec.size()); res++)
    {
        if (vec[res] == value) return res;
    }
    return -1;
}

short check_existence_and_bisect(const std::vector<unsigned short>& vec, const unsigned short value)
{
    short size = static_cast<short>(vec.size());
    for (short res = 0; res < size; res++)
    {
        if (vec[res] == value) return -1;
        else if (vec[res] > value) return res;
    }
    return size;
}