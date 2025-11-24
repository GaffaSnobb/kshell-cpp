#include <iostream>
#include <cstddef> // stdints
#include <random>

#include "data_structures.hpp"
#include "tools.hpp"
#include "linear_algebra.hpp"
#include "../external/eigen-3.4.0/Eigen/Dense"
#include "../external/eigen-3.4.0/Eigen/Eigenvalues"
#include "../tests/data/symmetric_test_matrices.hpp"

using std::cout;
using std::endl;

namespace lanczos
{
void lanczos(
    const Interaction &interaction,
    /*const*/ double *Hlool//,
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
   //    const size_t m_dim = interaction.basis_states.size();
   
    double *H = (double *)testmat::seven_by_seven;
    const size_t m_dim = 7;

    // HERE!! The algorithm as it is now requires exlicit representation of the entire matrix, not just upper diag as it is now.
    // print_flattened_2d_array(H, m_dim, m_dim     );

    const size_t n_lanc_steps = 7;
    if (n_lanc_steps > m_dim) throw std::runtime_error("n_lanc_steps cannot be larger than m_dim!");
    
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
    double *current_orth_lanc_vec = nullptr;

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

    for (int64_t lanc_step_idx = 1; lanc_step_idx < n_lanc_steps; lanc_step_idx++)
    {   /*
        NOTE: The index in this loop has to be signed because of the
        comparison with a (sometimes) negative number in the
        orthogonolisation loop index.

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

        // // Re-orthogonalise w here!
        // for (int64_t lanc_orth_idx = 0; lanc_orth_idx < (lanc_step_idx - 2); lanc_orth_idx++)   // Maybe this should start at 1? 0 is the random vector.
        // {
        //     current_orth_lanc_vec = lanc_vecs + lanc_orth_idx*m_dim;     // Pick the correct row in lanc_vecs.
        //     const double tmp_dot = lalg::vec_vec_dot(current_orth_lanc_vec, w, m_dim);
        //     for (size_t i = 0; i < m_dim; i++) w[i] = w[i] - current_orth_lanc_vec[i]*tmp_dot;
        // }

        H_krylov[(lanc_step_idx - 1)*m_dim + lanc_step_idx] = beta;      // Superdiagonal, aka. the row above the main diagonal on the same col.
        H_krylov[lanc_step_idx*m_dim + (lanc_step_idx - 1)] = beta;      // Subdiagonal, aka. the col behind the main diagonal on the same row.   
        // break;
    }

    printf("\n\n");
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> H_host_mapped(H_krylov, n_lanc_steps, m_dim);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(H_host_mapped);
    
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> H_host_mapped_2(H, m_dim, m_dim);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver_2(H_host_mapped_2);

    printf("Approx, exact, diff\n");

    auto approx_vals = solver.eigenvalues();
    auto exact_vals  = solver_2.eigenvalues();

    size_t n = std::min(approx_vals.size(), exact_vals.size());

    for (size_t i = 0; i < n; ++i) {
        printf("%f,  %f,  %f\n", approx_vals[i], exact_vals[i], std::abs(approx_vals[i] - exact_vals[i]));
    }
    printf("\n");

    // printf("Approx. eigenvalues:\n");
    // for (auto elem : solver.eigenvalues()) printf("%f, ", elem);
    // printf("\n");
    // printf("Exact eigenvalues:\n");
    // for (auto elem : solver_2.eigenvalues()) printf("%f, ", elem);
    // printf("\n");

    // std::cout << "Eigenvalues:\n" << solver_2.eigenvalues() << std::endl;

}

void lanczos_from_kshell()
{
    const size_t max_lanc_vec = 10;
    const double *H = testmat::eight_by_eight;
    const size_t m_dim = 8;
    const size_t n_lanc_steps = 5;
    const size_t n_skip_diag = 2;   // Diagonalise only every n_skip_diag'th iteration.

    // double tridiagonal_matrix[max_lanc_vec*max_lanc_vec];   // The matrix whose eigenvalues approximate the eigenvalues of the original matrix. 1D representation of 2D array.
    double tridiagonal_matrix_diagonal[max_lanc_vec];
    array_tools::fill_array(tridiagonal_matrix_diagonal, max_lanc_vec, 0);
    double tridiagonal_matrix_offdiagonal[max_lanc_vec - 1];
    array_tools::fill_array(tridiagonal_matrix_offdiagonal, max_lanc_vec - 1, 0);
    double tridiagonal_matrix_offdiagonal_tmp[max_lanc_vec - 1];    // tmp array because dsterf uses the off-diag array as tmp storage for its calculations, aka. destroying the array (but I wanna keep the values).
    double tridiagonal_matrix_eigenvals[max_lanc_vec];
    array_tools::fill_array(tridiagonal_matrix_eigenvals, max_lanc_vec, 1e6);
    double tridiagonal_matrix_eigenvals_previous[max_lanc_vec];
    double all_lanczos_vectors[m_dim*max_lanc_vec];   // Each row is a lanczos vector of length `m_dim`. 1D representation of 2D array.
    double *lanczos_vec_current = nullptr; // Will be set to the start address of each row of all_lanczos_vectors in the iv loop. For the lanczos vector of the current step.
    double *lanczos_vec_next = nullptr;    // Will be set to the start address of each row of all_lanczos_vectors in the iv loop. For the lanczos vector of the next step.

    print(tridiagonal_matrix_eigenvals_previous, max_lanc_vec);

    // Generate normalised random initial Lanczos vector.
    std::random_device rd;
    std::mt19937 gen(rd());
    gen.seed(1337);
    std::normal_distribution<> distribution(0.0, 1.0); // mean, std
    for (size_t i = 0; i < m_dim; i++) all_lanczos_vectors[i] = distribution(gen);    // Put it into the first row.
    lalg::normalise_vector(all_lanczos_vectors, m_dim);

    /* 
    Here would be the "`maxiter` loop" aka. the thick-restart loop. Without it
    there are no thick restarts.

    for (size_t thick_restart_idx = 0; thick_restart_idx < max_thick_restarts; thick_restart_idx++)

    the next loop would then be inside the thick restart loop.
    */
    for (size_t iv = 0; iv < (max_lanc_vec - 1); iv++)
    {   /*
        Equivalent to the `do iv = n_iv_start, max_lv-1` loop in KSHELL's
        `lanczos.f90`.

        This is a single thick restart. The tri-diagonal matrix is constructed
        and diagonalised, then the convergence requirement is checked. If
        convergence, the loop is broken (thick-restart loop is broken too) and
        the eigenvalues are good.

        A thick restart usually keeps some Lanczos (?) vectors from the
        previous thick restart, meaning that `iv` doesn't have to start at 0.
        */
        lanczos_vec_current = all_lanczos_vectors + m_dim*iv;  // Set the `lanczos_vec_current` pointer to the start of the iv'th row.
        lanczos_vec_next = all_lanczos_vectors + m_dim*(iv + 1);

        // call matvec(vec(iv), vec(iv+1))  Multiply the Hamiltonian with vec(iv) and store the res in vec(iv+1).
        lalg::mat_vec_mul(H, lanczos_vec_current, m_dim, lanczos_vec_next); // lanczos_vec_next = H*lanczos_vec_current

        // call dotprod(vec(iv+1), vec(iv), an)
        // tridiagonal_matrix[iv, iv] = an
        const double diagonal_element = lalg::vec_vec_dot(lanczos_vec_next, lanczos_vec_current, m_dim);
        tridiagonal_matrix_diagonal[iv] = diagonal_element; // Or should it start at iv + 1?
        
        
        // te_last(:neig) = teval(:neig)   ! Keep the previous eigenvalues of the tri-diagonal matrix for checking convergence.
        std::memcpy(tridiagonal_matrix_eigenvals_previous, tridiagonal_matrix_eigenvals, max_lanc_vec*sizeof(double));

        // Diagonalise the tri-matrix (at some interval, shouldn't be every time (computationally heavy))
        if ( ((iv%n_skip_diag) == 0) or (iv == (max_lanc_vec - 2)) )
        {   /*
            Diagonalise the tri-diagonal matrix. Diagonalise only every
            n_skip_diag'th iteration or at the last iteration of the iv
            loop.

            Copy diagonal arr to eigenvalues arr because the input array
            to dsterf with diagonals will be overwritten with
            eigenvalues.

            Copy off-diag arr to a tmp arr because `dsterf` uses the
            off-diag arr for tmp storage and overwrites the data.
            */
            std::memcpy(tridiagonal_matrix_eigenvals, tridiagonal_matrix_diagonal, max_lanc_vec*sizeof(double));
            std::memcpy(tridiagonal_matrix_offdiagonal_tmp, tridiagonal_matrix_offdiagonal, (max_lanc_vec - 1)*sizeof(double));
            lalg::diagonalise(max_lanc_vec, tridiagonal_matrix_eigenvals, tridiagonal_matrix_offdiagonal_tmp);
            
        }

        // call reorth(iv+1)
        // call dotprod(vec(iv+1), vec(iv+1), bn)   ! The norm of the residual vector.

        // Check the norm bn. The norm is proportional to the error of the Lanczos approximation.
        // If it is sufficiently small, stop the Lanczos procedure because the eigenvalues of the tri-matrix are sufficiently close to the true eigenvalues.

        // Populate upper and lower diag of the tri-matrix.
        // bn = sqrt(bn)
        // tridiagonal_matrix[iv, iv+1] = bn
        // tridiagonal_matrix[iv+1, iv] = bn

        // x = 1.d0/bn
        // do mq = 1, ndim
        //     vec(iv+1)%p(mq) = x * vec(iv+1)%p(mq)
        // end do


    }
}

} // namespace lanczos
