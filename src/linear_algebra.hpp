#ifndef LINEAR_ALGEBRA_HPP
#define LINEAR_ALGEBRA_HPP

#include <stdint.h>
#include <cstddef>

namespace lalg
{
void mat_vec_mul(const double* mat, const double* vec, const size_t m_dim, double* res);
double vec_vec_dot(const double *vec_0, const double *vec_1, const size_t size);
void normalise_vector(double *vec, const size_t size);
} // namespace lalg


#endif