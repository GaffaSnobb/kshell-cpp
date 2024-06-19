#include <stdint.h>
#include <cstddef>
#include <cmath>
#include <stdexcept>

namespace lalg
{
void mat_vec_mul(const double* mat, const double* vec, const size_t m_dim, double* res)
{
    for (size_t row_idx = 0; row_idx < m_dim; row_idx++)
    {
        res[row_idx] = 0;
        for (size_t col_idx = 0; col_idx < m_dim; col_idx++)
        {
            res[row_idx] += mat[row_idx*m_dim + col_idx]*vec[col_idx];
        }
    }
}

double vec_vec_dot(const double *vec_0, const double *vec_1, const size_t size)
{
    double sum = 0;
    for (size_t i = 0; i < size; i++) sum += vec_0[i]*vec_1[i];
    return sum;
}

void normalise_vector(double *vec, const size_t size)
{
    const double norm = std::sqrt(vec_vec_dot(vec, vec, size));
    if (norm <= 0) throw std::runtime_error("Vector norm has to be larger than 0!");
    for (size_t i = 0; i < size; i++) vec[i] /= norm;
}
} // namespace lalg