#include <fstream>
#include <vector>
#include <sstream>

#include "tools.hpp"
#include "data_structures.hpp"

using std::cout;
using std::endl;

void make_tbme_map(
    std::unordered_map<Key, double>& tbme_map,
    std::vector<Key>& tbme_keys,
    std::vector<unsigned short>& orb_0,
    std::vector<unsigned short>& orb_1,
    std::vector<unsigned short>& orb_2,
    std::vector<unsigned short>& orb_3,
    std::vector<unsigned short>& j_couple,
    std::vector<double>& tbme_list,
    std::vector<OrbitalParameters> model_space_orbitals
)
{
    /*
    Organise the TBMEs into an unordered map which takes a tuple of five
    unsigned shorts as a key and a double as the value. The first four
    numbers in each key represent the orbitals which are interacting and
    the fifth number is the total angular momentum (j) they couple to.
    The value of the unordered map is the corresponding TBME.
    */    
    for (int i = 0; i < tbme_list.size(); i++)
    {   
        unsigned short i0 = orb_0[i];
        unsigned short i1 = orb_1[i];
        unsigned short i2 = orb_2[i];
        unsigned short i3 = orb_3[i];
        unsigned short j = j_couple[i];
        double tbme = tbme_list[i];

        Key key = {i0, i1, i2, i3, j};
        tbme_keys.push_back(key);
        tbme_map[key] = tbme;

        short sign_01 = std::pow(-1, (model_space_orbitals[i0].j + model_space_orbitals[i1].j)/2 - j + 1);
        short sign_23 = std::pow(-1, (model_space_orbitals[i2].j + model_space_orbitals[i3].j)/2 - j + 1);

        /*
        The following network of if statements is taken from KSHELL's
        espe.py.
        */
        if (i0 != i1)
        {
            Key key = {i1, i0, i2, i3, j};
            tbme_keys.push_back(key);
            tbme_map[key] = tbme*sign_01;
        }
        if (i2 != i3)
        {
            Key key = {i0, i1, i3, i2, j};
            tbme_keys.push_back(key);
            tbme_map[key] = tbme*sign_23;
        }
        if ((i0 != i1) and (i2 != i3))
        {
            Key key = {i1, i0, i3, i2, j};
            tbme_keys.push_back(key);
            tbme_map[key] = tbme*sign_01*sign_23;
        }
        // if (i0, i1) != (i2, i3)
        if ((i0 != i2) or (i1 != i3))
        {
            Key key = {i2, i3, i0, i1, j};
            tbme_keys.push_back(key);
            tbme_map[key] = tbme;
            if (i0 != i1)
            {
                Key key = {i2, i3, i1, i0, j};
                tbme_keys.push_back(key);
                tbme_map[key] = tbme*sign_01;
            }
            if (i2 != i3)
            {
                Key key = {i3, i2, i0, i1, j};
                tbme_keys.push_back(key);
                tbme_map[key] = tbme*sign_23;
            }
            if ((i0 != i1) and (i2 != i3))
            {
                Key key = {i3, i2, i1, i0, j};
                tbme_keys.push_back(key);
                tbme_map[key] = tbme*sign_01*sign_23;
            }
        }
    }
    return;
}

const Interaction load_interaction(const std::string& interaction_filename)
{
    /*
    Load raw interaction data from the interaction file.
    */
    std::vector<unsigned short> orb_0, orb_1, orb_2, orb_3, j_couple;
    std::vector<double> tbme, spe;

    std::vector<OrbitalParameters> model_space_protons_orbitals;
    std::vector<OrbitalParameters> model_space_neutrons_orbitals;
    std::vector<OrbitalParameters> model_space_orbitals;

    std::ifstream infile(interaction_filename);
    std::string line;
    unsigned int n_tbme;    // gs8.snt has over 10k TBMEs. short is enough for this, but other interactions might have even more.
    unsigned short tbme_mass_dependence_method;
    unsigned short n_proton_orbitals;
    unsigned short n_neutron_orbitals;
    unsigned short n_spe;
    unsigned short n_core_neutrons;
    unsigned short n_core_protons;
    double tbme_mass_dependence_exponent;
    double tbme_mass_dependence_denominator;
    int _;  // Garbage.

    while (std::getline(infile, line))
    {
        std::string target_string = "! model space";
        if (!line.empty() && (line == target_string))
        {
            /*
            Read model space info. Example from w.snt:
            ...
            ! model space
            3   3     8   8  <-- This line!
            1       0   2   3  -1    !  1 = p 0d_3/2
            2       0   2   5  -1    !  2 = p 0d_5/2
            3       1   0   1  -1    !  3 = p 1s_1/2
            4       0   2   3   1    !  4 = n 0d_3/2
            5       0   2   5   1    !  5 = n 0d_5/2
            6       1   0   1   1    !  6 = n 1s_1/2
            ...
            */
            std::getline(infile, line);
            std::istringstream iss(line);
            iss >> n_proton_orbitals >> n_neutron_orbitals >> n_core_protons >> n_core_neutrons;
            break;
        }
    }

    for (int i = 0; i < n_proton_orbitals; i++)
    {
        /*
        ! model space
        3   3     8   8
        1       0   2   3  -1    !  1 = p 0d_3/2  <-- These lines!
        2       0   2   5  -1    !  2 = p 0d_5/2  <-- These lines!
        3       1   0   1  -1    !  3 = p 1s_1/2  <-- These lines!
        4       0   2   3   1    !  4 = n 0d_3/2
        5       0   2   5   1    !  5 = n 0d_5/2
        6       1   0   1   1    !  6 = n 1s_1/2
        ...
        */
        std::getline(infile, line);
        std::istringstream iss(line);
        unsigned short n_tmp, l_tmp, j_tmp;
        short tz_tmp;
        iss >> _ >> n_tmp >> l_tmp >> j_tmp >> tz_tmp;
        std::vector<short> jz_tmp = range<short>(-j_tmp, j_tmp + 1, 2);
        model_space_protons_orbitals.emplace_back(n_tmp, l_tmp, j_tmp, tz_tmp, jz_tmp);
        model_space_orbitals.emplace_back(n_tmp, l_tmp, j_tmp, tz_tmp, jz_tmp);
    }

    for (int i = 0; i < n_neutron_orbitals; i++)
    {
        /*
        ! model space
        3   3     8   8
        1       0   2   3  -1    !  1 = p 0d_3/2
        2       0   2   5  -1    !  2 = p 0d_5/2
        3       1   0   1  -1    !  3 = p 1s_1/2
        4       0   2   3   1    !  4 = n 0d_3/2  <-- These lines!
        5       0   2   5   1    !  5 = n 0d_5/2  <-- These lines!
        6       1   0   1   1    !  6 = n 1s_1/2  <-- These lines!
        ...
        */
        std::getline(infile, line);
        std::istringstream iss(line);
        unsigned short n_tmp, l_tmp, j_tmp;
        short tz_tmp;
        iss >> _ >> n_tmp >> l_tmp >> j_tmp >> tz_tmp;
        std::vector<short> jz_tmp = range<short>(-j_tmp, j_tmp + 1, 2);
        model_space_neutrons_orbitals.emplace_back(n_tmp, l_tmp, j_tmp, tz_tmp, jz_tmp);
        model_space_orbitals.emplace_back(n_tmp, l_tmp, j_tmp, tz_tmp, jz_tmp);
    }

    while (std::getline(infile, line))
    {
        std::string target_string = "! interaction";
        if (!line.empty() && (line == target_string))
        {
            /*
            Read the number of single-particle energies. Example from w.snt:
            ...
            ! interaction
                6   0   <--- This line!
            1   1      1.64658
            2   2     -3.94780
            3   3     -3.16354
            4   4      1.64658
            5   5     -3.94780
            6   6     -3.16354
            ...
            */
            std::getline(infile, line);
            std::istringstream iss(line);
            iss >> n_spe >> _;
            break;
        }
    }

    for (int i = 0; i < n_spe; i++)
    {
        /*
        Read the single-particle energies. Example:
        ...
        ! interaction
            6   0
        1   1      1.64658  <--- These lines!
        2   2     -3.94780  <--- These lines!
        3   3     -3.16354  <--- These lines!
        4   4      1.64658  <--- These lines!
        5   5     -3.94780  <--- These lines!
        6   6     -3.16354  <--- These lines!
        ...
        */
        std::getline(infile, line);
        std::istringstream iss(line);
        double spe_tmp;
        iss >> _ >> _ >> spe_tmp;
        spe.push_back(spe_tmp);
    }
    
    std::getline(infile, line);
    std::istringstream iss(line);
    iss >> n_tbme >> tbme_mass_dependence_method >> tbme_mass_dependence_denominator >> tbme_mass_dependence_exponent;

    for (int i = 0; i < n_tbme; i++)
    {
        /*
        Read the two-body matrix elements. Example:
        ...
        158   1  18  -0.30000
        1   1   1   1     0     -2.18450  <--- These lines!
        1   1   1   1     2     -0.06650  <--- These lines!
        1   1   1   2     2      0.61490  <--- These lines!
        1   1   1   3     2      0.51540  <--- These lines!
        1   1   2   2     0     -3.18560  <--- These lines etc!
        ...
        */
        std::getline(infile, line);
        std::istringstream iss(line);
        unsigned short orb_0_tmp, orb_1_tmp, orb_2_tmp, orb_3_tmp, j_couple_tmp;
        double tbme_tmp;
        
        iss >> orb_0_tmp >> orb_1_tmp >> orb_2_tmp >> orb_3_tmp >> j_couple_tmp >> tbme_tmp;
        orb_0_tmp--;    // Start indices from 0.
        orb_1_tmp--;
        orb_2_tmp--;
        orb_3_tmp--;
        j_couple_tmp = j_couple_tmp*2;
        orb_0.push_back(orb_0_tmp);
        orb_1.push_back(orb_1_tmp);
        orb_2.push_back(orb_2_tmp);
        orb_3.push_back(orb_3_tmp);
        j_couple.push_back(j_couple_tmp);
        tbme.push_back(tbme_tmp);
    }
    infile.close();

    std::unordered_map<Key, double> tbme_map;
    std::vector<Key> tbme_keys;
    make_tbme_map(
        tbme_map,
        tbme_keys,
        orb_0,
        orb_1,
        orb_2,
        orb_3,
        j_couple,
        tbme,
        model_space_orbitals
    );

    const ModelSpace model_space_protons(model_space_protons_orbitals);
    const ModelSpace model_space_neutrons(model_space_neutrons_orbitals);
    const ModelSpace model_space(model_space_orbitals);

    const Interaction interaction(
        tbme_mass_dependence_method,
        n_core_protons,
        n_core_neutrons,
        n_core_protons + n_core_neutrons,
        tbme_mass_dependence_exponent,
        tbme_mass_dependence_denominator,
        spe,
        model_space,
        model_space_protons,
        model_space_neutrons,
        tbme_map,
        tbme_keys      
    );
    return interaction;
}