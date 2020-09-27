/**
* jacobisolver.cpp: implementation file for JacobiSolver
*
* Author: Are Frode Kvanum
*
* Completion Date: 26.09.2020
*/

#include <iomanip>

#include "tridiagonalmatrix.hpp"

// CONSTRUCTOR

/**
* Default constructor
* wrapper for calling on parent-constructor in TridiagonalMatrix
*
* @param N Dimensionality (NxN) of matrix to be constructed
* @param diag Diagonal element of the matrix
* @param nondiag Nondiagonal element of the matrix (symmetric)
* @param rho_max (default = 1.0) maximum length of the system
*/
void JacobiSolver::init(int N, double diag, double nondiag, double rho_max) {
    initialize(N, diag, nondiag, rho_max);
}

// PUBLIC MEMBER FUNCTIONS

/**
* solve: performs the Jacobi rotation algorithm to reduce all diagonal elements
*        to eigenvalues and nondiagonal elements close to 0
*
* @return Number of iterations(rotations) performed
*/
int JacobiSolver::solve() {
    int k, l, counter = 0;
    double tol = 1.0e-8;
    double iterations = (double) m_N * (double) m_N * (double) m_N;
    double max_offdiag = maxoffdiag(&k, &l);

    while (fabs(max_offdiag) > tol && (double) counter < iterations) {
        max_offdiag = maxoffdiag(&k, &l);
        rotate(k, l);
        counter++;
    }

    return counter;
}

/**
* write_to_file: writes results to output file
*
* @param fname /path/to/file
* @param transform Number of rotations performed
* @return Nothing is returned, as this function writes to file
*/
void JacobiSolver::write_to_file(string fname, int transform) {
    int low_idx = find_lowest_eigval();
    m_ofile.open(fname);
    m_ofile << "N.O. transformations: " << transform << endl;
    for (int i = 0; i < m_N; i++) {
        m_ofile << m_x(i) << " " << m_eigmat(i, low_idx) << endl;
    }
    m_ofile.close();
}

/**
* find_lowest_eigval: Searches through the reduced tridiagonal matrix
*                     to find the lowest eigenvalue
*
* @return Index whre the lowest eigenvalue is located
*/
int JacobiSolver::find_lowest_eigval() {
    int idx = 0;
    double min = m_Toeplitz(0,0);
    for (int i = 1; i < m_N; i++) {
        if (m_Toeplitz(i,i) < min) {
            min = m_Toeplitz(i,i);
            idx = i;
        }
    }
    return idx;
}

/**
* print_eigvals: prints the diagonal elements of the tridiagonal matrix
*                USED FOR PRINTING SELECTED RESULTS
*
* @return Nothing is returned
*/
void JacobiSolver::print_eigvals() {
    vec computed = zeros<vec>(m_N);
    for (int i = 0; i < m_N; i++)
        computed(i) = m_Toeplitz(i,i);

    uvec indices = sort_index(computed);
    for (int i = 0; i < 4; i++)
        cout << setw(15) << setprecision(8) << computed(indices(i)) << endl;
}

/**
* print_toeplitz: prints the entire tridiagonal toeplitz matrix
*                 USED FOR PRINTING SELECTED RESULTS
*
* return Nothing is returned
*/
void JacobiSolver::print_toeplitz() {
    to_string(m_Toeplitz);
}

/**
* unittest_maxoffdiag: used for testing private function maxoffdiag
*                      USED FOR UNITTESTS ONLY
*
* @return computed result of maxoffdiag
*/
double JacobiSolver::unittest_maxoffdiag() {
    m_N = 5;
    int k, l;
    m_Toeplitz = zeros<mat>(m_N,m_N);
    for (int i = 0; i < m_N; i++) {
        for (int j = i; j < m_N; j++) {
            m_Toeplitz(i,j) = i*j;
        }
    }
    return maxoffdiag(&k, &l);
}

// PRIVATE MEMBER FUNCTIONS

/**
* maxoffdiag: finds the maximum nondiagonal element in a symmetric matrix
*
* @param k pointer to index k
* @param l pointer to index l
* @return maximum nondiagonal element in the matrix
*/
double JacobiSolver::maxoffdiag(int *k, int *l) {
    double max = 0.0;
    for (int i = 0; i < m_N; i++) {
        for (int j = i+1; j < m_N; j++) {
            if (fabs(m_Toeplitz(i,j)) > max) {
                max = fabs(m_Toeplitz(i,j));
                *k = i;
                *l = j;
            }
        }
    }

    return max;
}

/**
* rotate: performs a single Jacobi rotation
*
* @param k the k'th index of the highest nondiagonal element in a tridiagonal matrix
* @param l the l'th index of the highest nondiagonal element in a tridiagonal matrix
* @return Nothing is returned, m_Toeplitz is reduced and the eigenvectors are recomputed
*/
void JacobiSolver::rotate(int k, int l) {
    double s, c, a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    if (m_Toeplitz(k,l) != 0.0) {
        double t, tau;
        tau = (m_Toeplitz(l,l) - m_Toeplitz(k,k))/(2*m_Toeplitz(k,l));

        if (tau > 0) {
            t = 1.0/(tau + sqrt(1.0 + tau*tau));
        } else {
            t = -1.0/(-tau + sqrt(1.0 + tau*tau));
        }
        c = 1.0 / (sqrt(1 + t*t));
        s = t*c;
    } else {
        c = 1.0;
        s = 0.0;
    }
    a_kk = m_Toeplitz(k,k);
    a_ll = m_Toeplitz(l,l);

    m_Toeplitz(k,k) = c*c*a_kk - 2.0*c*s*m_Toeplitz(k,l) + s*s*a_ll;
    m_Toeplitz(l,l) = s*s*a_kk + 2.0*c*s*m_Toeplitz(k,l) + c*c*a_ll;
    m_Toeplitz(k,l) = 0.0;
    m_Toeplitz(l,k) = 0.0;

    for (int i = 0; i < m_N; i++) {
        if (i != k && i != l) {
            a_ik = m_Toeplitz(i,k);
            a_il = m_Toeplitz(i,l);
            m_Toeplitz(i,k) = c*a_ik - s*a_il;
            m_Toeplitz(k,i) = m_Toeplitz(i,k);
            m_Toeplitz(i,l) = c*a_il + s*a_ik;
            m_Toeplitz(l,i) = m_Toeplitz(i,l);
        }

        r_ik = m_eigmat(i,k);
        r_il = m_eigmat(i,l);
        m_eigmat(i,k) = c*r_ik - s*r_il;
        m_eigmat(i,l) = c*r_il + s*r_ik;
    }
}

