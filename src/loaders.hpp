#ifndef LOADERS_HPP
#define LOADERS_HPP

#include <iostream>
#include "data_structures.hpp"

Interaction load_interaction(const std::string& interaction_filename, const unsigned short n_valence_protons, const unsigned short n_valence_neutrons);

#endif