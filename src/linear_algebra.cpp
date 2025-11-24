#include <stdint.h>
#include <cstddef>
#include <cmath>
#include <stdexcept>
#include <lapacke.h>
#include <spdlog/spdlog.h>

#include "tools.hpp"

namespace lalg
{
void mat_vec_mul(const double* mat, const double* vec, const size_t m_dim, double* res)
{   /*
    Left multiply a (m_dim x m_dim) matrix on a vector, save the
    resulting vector in `res`:

        res = mat@vec
    */
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
{   /*
    Take the dot product of `vec_0` and `vec_1`.
    */
    double sum = 0;
    for (size_t i = 0; i < size; i++) sum += vec_0[i]*vec_1[i];
    return sum;
}

void normalise_vector(double *vec, const size_t size)
{   /*
    Normalise `vec`, aka. make its length / norm / magnitude equal to 1.
    */
    const double norm = std::sqrt(vec_vec_dot(vec, vec, size));
    if (norm <= 0) throw std::runtime_error("Vector norm has to be larger than 0!");
    for (size_t i = 0; i < size; i++) vec[i] /= norm;
}

void vec_div(const double *vec, const double factor, const size_t size, double *res)
{   /*
    Divide every element of `vec` by `factor`. Store the resulting
    vector in `res`.
    */
    for (size_t i = 0; i < size; i++) res[i] = vec[i]/factor;
}

void diagonalise(const size_t n_rows_or_cols, double *diagonal_then_eigenvalues, double *off_diagonal_then_garbage)
{   /*
    DSTERF computes all eigenvalues of a symmetric tridiagonal matrix
    using the Pal-Walker-Kahan variant of the QL or QR algorithm.

    https://netlib.org/lapack/explore-html/d4/d9d/group__sterf_gad293bb81da1c7785b42796d1e197f08c.html
    
    lapack_int LAPACKE_dsterf( lapack_int n, double* d, double* e );

    Consider turning off nan check: #ifndef LAPACK_DISABLE_NAN_CHECK
    */

    lapack_int info = LAPACKE_dsterf(
        n_rows_or_cols,             // rows or cols, the matrix is square.
        diagonal_then_eigenvalues,  // Contains diagonal now, eigenvalues after.
        off_diagonal_then_garbage   // Contains off-diagonal now, garbage after.
    );

    if (info != 0)
    {
        spdlog::critical("dsterf failed with info={}", info);
        throw std::runtime_error("dsterf failed with info = " + std::to_string(info));
    }
}

} // namespace lalg