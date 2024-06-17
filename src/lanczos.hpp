#ifndef LANCZOS_HPP
#define LANCZOS_HPP

#include <stdint.h>

#include "data_structures.hpp"

namespace lanczos
{
void lanczos(const Interaction &interaction, const double *H);
} // namespace lanczos


#endif