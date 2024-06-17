#include <iostream>
#include <cstddef> // stdints
#include <random>
#include <cmath>
#include <stdexcept>

#include "data_structures.hpp"

using std::cout;
using std::endl;

namespace lanczos
{
void mat_vec_mul(
    const double* mat,
    const double* vec,
    double* res,
    const size_t m_dim,
    const size_t n_steps
)
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
    for (size_t i = 0; i < size; i++) sum += vec_0[i] + vec_1[i];
    return sum;
}

double vector_norm(const double *vec, const size_t size)
{
    double norm = 0;
    for (size_t i = 0; i < size; i++) norm += vec[i]*vec[i];
    norm = std::sqrt(norm);

    return norm;
}

void normalise_vector(double *vec, const size_t size)
{
    const double norm = vector_norm(vec, size);
    if (norm <= 0) throw std::runtime_error("Vector norm has to be larger than 0!");
    for (size_t i = 0; i < size; i++) vec[i] /= norm;
}

// void lanczos(const Interaction interaction, double *H_from_device)
void lanczos(const Interaction &interaction, const double *H)
{   /*
    m_dim = A.shape[0]
    Q = np.zeros((m_dim, num_steps + 1))
    alpha = np.zeros(num_steps)
    beta = np.zeros(num_steps + 1)
    q = np.random.normal(size=m_dim)
    q /= np.linalg.norm(q)
    Q[:, 0] = q

    for k in range(num_steps):
        v = A @ Q[:, k]
        alpha[k] = Q[:, k].T @ v
        v -= alpha[k] * Q[:, k] + beta[k] * Q[:, k-1]
        beta[k+1] = np.linalg.norm(v)

        if beta[k+1] < 1e-10:  # A convergence tolerance
            print("Early convergence at step:", k)
            return Q[:, :k+1], alpha[:k], beta[:k+1]

        Q[:, k+1] = v / beta[k+1]

    return Q[:, :num_steps], alpha, beta[:num_steps]

    https://stackoverflow.com/questions/68745223/eigenvectors-of-lanczos-algorithm-differ-from-numpy-linal-eig-for-complex-matr
    */
    const size_t m_dim = interaction.basis_states.size();
    const size_t n_steps = 10;
    
    double H_krylov[m_dim*n_steps];     // Flat 2D array.
    double lanc_vecs[m_dim*n_steps];    // Flat 2D array.
    for (size_t i = 0; i < m_dim*n_steps; i++)
    {   // Just to be sure during development, might not be needed in the end.
        lanc_vecs[i] = 0;
        H_krylov[i] = 0;
    }
    double init_lanc_vec[m_dim];
    double alpha[n_steps] = {0};
    double beta[n_steps + 1] = {0};
    double w[m_dim];

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> distribution(0.0, 1.0); // mean, std

    for (size_t i = 0; i < m_dim; i++) init_lanc_vec[i] = distribution(gen);

    normalise_vector(init_lanc_vec, m_dim);

    for (size_t i = 0; i < m_dim; i++)
    // {   /*
    //     Equiv. to lanc_vecs_2d[i][0] = init_lanc_vec[i]. Broadcast init_lanc_vec to the first
    //     element of each row in lanc_vecs.
    //     */
    //     lanc_vecs[i*n_steps] = init_lanc_vec[i];
    {
        lanc_vecs[i] = init_lanc_vec[i];
    }

    // First step:
    double *current_lanc_vec = lanc_vecs + 0*m_dim;     // Placeholder for easier reading.
    mat_vec_mul(H, current_lanc_vec, w, m_dim, n_steps);
    alpha[0] = vec_vec_dot(current_lanc_vec, w, m_dim);
    for (size_t i = 0; i < m_dim; i++) w[i] = w[i] - alpha[0]*current_lanc_vec[i];
    H_krylov[0*m_dim + 0] = alpha[0];   // On the diagonal of H_krylov.

    for (size_t k = 1; k < n_steps; k++)
    {
        // mat_vec_mul(H, lanc_vecs);  // v = A @ Q[:, k]
    }

}
} // namespace lanczos
