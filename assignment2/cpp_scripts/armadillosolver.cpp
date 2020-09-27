/**
* armadillosolver.cpp: implementation file for ArmadilloSolver
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
* Calls TridiagonalMatrix's constructor and sets ut vector for storing eigenvalues
*
* @param N Dimensionality (NxN) of matrix to be constructed
* @param diag Diagonal element of the matrix
* @param nondiag Nondiagonal element of the matrix (symmetric)
* @param rho_max (default = 1.0) maximum length of the system
*/
void ArmadilloSolver::init(int N, double diag, double nondiag, double rho_max) {
    initialize(N, diag, nondiag, rho_max);
    m_eigval = zeros<vec>(m_N);
}

// PUBLIC MEMBER FUNCTIONS

/**
* solve: wrapper for calling Armadillo library function eig_sym
*
* @return nothing is returned, eigenvalues are stored in m_eigval
*         and eigenvectors stored in m_eigmat
*/
void ArmadilloSolver::solve() {
    eig_sym(m_eigval, m_eigmat, m_Toeplitz);
}

/**
* print_eigvals: prints the diagonal elements of the tridiagonal matrix
*                USED FOR PRINTING SELECTED RESULTS
*
* @return Nothing is returned
*/
void ArmadilloSolver::print_eigvals() {
    for (int i = 0; i < 4; i++)
        cout << setw(15) << setprecision(8) << m_eigval(i) << endl;
}

/**
* write_to_file: writes results to output file
*
* @param fname /path/to/file
* @return Nothing is returned, as this function writes to file
*/
void ArmadilloSolver::write_to_file(string fname) {
    m_ofile.open(fname);
    m_ofile << "Armadillo, no N.O. transformations" << endl;
    for (int i = 0; i < m_N; i++) {
        m_ofile << m_x(i) << " " << m_eigmat(i, 0) << endl;
    }
    m_ofile.close();
}
