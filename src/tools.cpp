#include <iostream>
#include <vector>
#include <chrono>
#include "tools.hpp"
#include "data_structures.hpp"

using std::cout;
using std::endl;
// using std::chrono::high_resolution_clock;
// using std::chrono::duration_cast;
// using std::chrono::milliseconds;

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

void print_vector(const std::vector<Key5>& vec)
{
    for (Key5 elem : vec)
    {
        cout << "{" << elem.a << ", " << elem.b << ", " << elem.c << ", " << elem.d << ", " << elem.e << "}," << endl;
    }
}

short index(const std::vector<unsigned short>& vec, const unsigned short value)
{
    /*
    Returns the index of the first occurence of `value` in `vec`. If
    `value` is not found, -1 is returned. Assumes that `vec` is short
    enough to be indexed by a short.

    Using std::find and std::distance has practically identical
    computation time at -Ofast.
    
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
    // auto it = std::find(vec.begin(), vec.end(), value);

    // if (it != vec.end())
    // {
    //     short idx = std::distance(vec.begin(), it);
    //     return idx;
    // }
    // else return -1;
    short size = static_cast<short>(vec.size());
    for (short res = 0; res < size; res++)
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
        if (vec[res] > value) return res;
        else if (vec[res] == value) return -1;
    }
    return size;
}

std::chrono::time_point<std::chrono::high_resolution_clock> timer()
{
    return std::chrono::high_resolution_clock::now();
}

std::chrono::milliseconds timer(std::chrono::time_point<std::chrono::high_resolution_clock> start, std::string name)
{
    std::chrono::time_point<std::chrono::high_resolution_clock> stop = std::chrono::high_resolution_clock::now();
    std::chrono::milliseconds duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    cout << name << ": " << duration.count()/1000.0 << " s" << endl;
    return duration;
}