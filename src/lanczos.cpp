#include <iostream>
#include <cstddef> // stdints
#include <random>

#include "data_structures.hpp"
#include "tools.hpp"
#include "linear_algebra.hpp"
#include "../external/eigen-3.4.0/Eigen/Dense"
#include "../external/eigen-3.4.0/Eigen/Eigenvalues"

using std::cout;
using std::endl;

namespace lanczos
{
void lanczos(
    const Interaction &interaction,
    /*const*/ double *H_lol//,
    // const size_t n_lanc_steps
)
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
    // const size_t m_dim = interaction.basis_states.size();
    const size_t m_dim = 4;
    const size_t n_lanc_steps = m_dim;
    if (n_lanc_steps > m_dim) throw std::runtime_error("n_lanc_steps cannot be larger than m_dim!");

    double H[5*5] = {
        -2.65526241, -1.73658919,  1.05043732, -1.35836282, -0.60596862,
        -1.73658919, -1.04257302, -0.38122495,  0.67562902, -0.56439201,
        1.05043732, -0.38122495,  3.95116467, -0.66926132,  0.58965748,
        -1.35836282,  0.67562902, -0.66926132, -0.28581319, -0.37952717,
        -0.60596862, -0.56439201,  0.58965748, -0.37952717, -1.22605036
    };
    
    double H_krylov[m_dim*n_lanc_steps];     // Flat 2D array. This is the tridiagonal matrix, T.
    double lanc_vecs[m_dim*n_lanc_steps];    // Flat 2D array. Also called Krylov vectors. Lanczos vectors are stored as rows in this matrix.
    for (size_t i = 0; i < m_dim*n_lanc_steps; i++)
    {   // Just to be sure during development, might not be needed in the end.
        lanc_vecs[i] = 0;
        H_krylov[i] = 0;
    }
    double alpha;   // Main diagonal elements.
    double beta;    // Super- and subdiagonal elements.
    double w[m_dim];
    double *current_lanc_vec = nullptr; // tmp pointers for easier reading.
    double *prev_lanc_vec = nullptr;
    double *next_lanc_vec = nullptr;

    // Generate normalised random initial Lanczos vector.
    std::random_device rd;
    std::mt19937 gen(rd());
    gen.seed(1337);
    std::normal_distribution<> distribution(0.0, 1.0); // mean, std
    for (size_t i = 0; i < m_dim; i++) lanc_vecs[i] = distribution(gen);    // Put it into the first row.
    lalg::normalise_vector(lanc_vecs, m_dim);

    // First step:
    current_lanc_vec = lanc_vecs + 0*m_dim;
    lalg::mat_vec_mul(H, current_lanc_vec, m_dim, w);   // First Krylov subspace vector.
    alpha = lalg::vec_vec_dot(current_lanc_vec, w, m_dim);
    for (size_t i = 0; i < m_dim; i++) w[i] = w[i] - alpha*current_lanc_vec[i];
    H_krylov[0*m_dim + 0] = alpha;   // On the diagonal of H_krylov.

    for (size_t lanc_step_idx = 1; lanc_step_idx < n_lanc_steps; lanc_step_idx++)
    {   /*
        Tridiagonal matrix with `n_lanc_steps` rows and `m_dim` cols,
        where d is the main diagonal, u is the superdiagonal, and l is
        the subdiagonal:
        
            [[d0, u1,  0,  0, 0]
             [l1, d1, u2,  0, 0]
             [ 0, l2, d2, u3, 0]
             [ 0,  0, l3, d3, u4]
             [ 0,  0,  0, l4, d4]]
        */
        prev_lanc_vec = current_lanc_vec;
        current_lanc_vec = lanc_vecs + lanc_step_idx*m_dim;         // Picks the correct row in lanc_vecs.
        
        beta = std::sqrt(lalg::vec_vec_dot(w, w, m_dim));           // Norm of w.
        lalg::vec_div(w, beta, m_dim, current_lanc_vec);            // Divide w by beta, aka. normalise w and save the result in `current_lanc_vec`.
        lalg::mat_vec_mul(H, current_lanc_vec, m_dim, w);           // Calculate next Krylov subspace vector.
        alpha = lalg::vec_vec_dot(current_lanc_vec, w, m_dim);      // NOTE: Why is there no square root here?
        H_krylov[lanc_step_idx*m_dim + lanc_step_idx] = alpha;      // Main diagonal.

        // Diagonalise H_krylov here and stop lanc step iteration if eigenvalues (which?) converge.

        for (size_t i = 0; i < m_dim; i++) w[i] = w[i] - alpha*current_lanc_vec[i] - beta*prev_lanc_vec[i];

        // Re-orthogonalise w here!

        H_krylov[(lanc_step_idx - 1)*m_dim + lanc_step_idx] = beta;      // Superdiagonal, aka. the row above the main diagonal on the same col.
        H_krylov[lanc_step_idx*m_dim + (lanc_step_idx - 1)] = beta;      // Subdiagonal, aka. the col behind the main diagonal on the same row.   
    }

    printf("\n\n");
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> H_host_mapped(H_krylov, n_lanc_steps, m_dim);
    // cout << H_host_mapped << endl;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(H_host_mapped);
    std::cout << "Eigenvalues:\n" << solver.eigenvalues() << std::endl;
    printf("\n\n");
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> H_host_mapped_2(H, m_dim, m_dim);
    // cout << H_host_mapped_2 << endl;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver_2(H_host_mapped_2);
    std::cout << "Eigenvalues:\n" << solver_2.eigenvalues() << std::endl;

}
} // namespace lanczos
