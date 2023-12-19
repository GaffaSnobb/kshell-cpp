#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

#include "tools.h"

using std::cout;
using std::endl;

void load_interaction()
{
    std::ifstream file("../snt/w.snt");
    std::vector<unsigned short> orb_0, orb_1, orb_2, orb_3, j_couple;
    std::vector<double> tbme, spe;
    std::string line;
    int n_spe, n_tbme;
    int _;

    while (std::getline(file, line))
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
            std::getline(file, line);
            std::istringstream iss(line);
            iss >> n_spe >> _;
            break;
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
            std::getline(file, line);
            std::istringstream iss(line);
            double spe_tmp;
            iss >> spe_tmp;
            spe.push_back(spe_tmp);
        }
        
        // std::istringstream iss(line);
        // double val1, val2, val3, val4, val5, val6;
        // if (iss >> val1 >> val2 >> val3 >> val4 >> val5 >> val6) {
        //     column1.push_back(val1);
        //     column2.push_back(val2);
        //     column3.push_back(val3);
        //     column4.push_back(val4);
        //     column5.push_back(val5);
        //     column6.push_back(val6);
        // }
    }
    cout << n_spe << endl;
}

int main() {

    load_interaction();
    return 0;
}