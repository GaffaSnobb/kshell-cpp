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
void mat_vec_mul(const double** mat, const double* vec, double* res, const size_t m_dim, const size_t n_steps)
{
    for (size_t i = 0; i < m_dim; i++)
    {
        res[i] = 0;
        for (size_t j = 0; j < m_dim; j++)
        {
            res[i] += mat[i][j]*vec[j];
        }
    }
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
    */
    const size_t m_dim = interaction.basis_states.size();
    const size_t n_steps = 10;
    
    double Q[m_dim*n_steps];    // Flat 2D array.
    for (size_t i = 0; i < m_dim*n_steps; i++) Q[i] = 0;    // Just to be sure during development, might not be needed in the end.
    double q[m_dim];
    double alpha[n_steps] = {0};
    double beta[n_steps + 1] = {0};
    double v[m_dim];

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> distribution(0.0, 1.0); // mean, std

    for (size_t i = 0; i < m_dim; i++) q[i] = distribution(gen);

    normalise_vector(q, m_dim);

    for (size_t i = 0; i < m_dim; i++)
    {   /*
        Equiv. to Q_2d[i][0] = q[i]. Broadcast q to the first element of
        each row in Q.
        */
        Q[i*n_steps] = q[i];
    }

    for (size_t k = 0; k < n_steps; k++)
    {
        // mat_vec_mul(H, Q);  // v = A @ Q[:, k]
    }



}
} // namespace lanczos
