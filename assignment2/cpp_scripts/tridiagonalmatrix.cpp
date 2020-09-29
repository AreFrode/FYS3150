/**
* tridiagonalmatrix.cpp: implementation file for TridiagonalMatrix
*
* Author: Are Frode Kvanum
*
* Completion Date: 26.09.2020
*/

#include "tridiagonalmatrix.hpp"

// CONSTRUCTOR

/**
* Default constructor
* Sets up the tridiagonal matrix, defines and stores step length
* and creates a matrix for storing eigenvectors
*
* @param N Dimensionality (NxN) of matrix to be constructed
* @param diag Diagonal element of the matrix
* @param nondiag Nondiagonal element of the matrix (symmetric)
* @param rho_max (default = 1.0) maximum length of the system
*/
void TridiagonalMatrix::initialize(int N, double diag, double nondiag, double rho_max) {
    m_N = N;
    m_h = rho_max/m_N;
    m_Toeplitz = zeros<mat>(m_N,m_N);
    m_diag = diag/(m_h*m_h);
    m_nondiag = nondiag/(m_h*m_h);
    m_x = linspace(m_h, rho_max-m_h, m_N);

    m_Toeplitz(0,0) = m_diag;
    m_Toeplitz(0,1) = m_nondiag;
    for (int i = 1; i < m_N-1; i++) {
        m_Toeplitz(i,i-1) = m_nondiag;
        m_Toeplitz(i,i) = m_diag;
        m_Toeplitz(i,i+1) = m_nondiag;
    }

    m_Toeplitz(m_N-1,m_N-2) = m_nondiag;
    m_Toeplitz(m_N-1,m_N-1) = m_diag;

    m_eigmat = zeros<mat>(m_N,m_N);
    for (int i = 0; i < m_N; i++)
        for (int j = 0; j < m_N; j++)
            if (i == j)
                m_eigmat(i,j) = 1.0;
}

// PUBLIC MEMEBER FUNCTIONS

/**
* add_harmonic_potential: adds the harmonic oscillator potential to
*                         diagonal elements
*
* @return Nothing is returned, m_Toeplitz is modified
*/
void TridiagonalMatrix::add_harmonic_potential() {
    for (int i = 0; i < m_N; i++)
        m_Toeplitz(i,i) += (m_x(i)*m_x(i));
}

/**
* add_harmonic_potential_two_electrons: add the potential for a case
*                                       of two electrons to the
*                                       diagonal elements
*
* @param omega frequency parameter
* @return Nothing is returned, m_Toeplitz is modified
*/
void TridiagonalMatrix::add_harmonic_potential_two_electrons(double omega) {
    for (int i = 0; i < m_N; i++) {
        m_Toeplitz(i,i) += (omega*omega)*(m_x(i)*m_x(i)) + 1./m_x(i);
    }
}

/**
* get_x: getter function for private variable m_x
*
* @param i The index the value in m_x is loaded from
* @return The specified value at index i in m_x
*/
double TridiagonalMatrix::get_x(int i) {
    return m_x(i);
}

/**
* get_N: getter function for private variable m_N
*
* @return the integer value m_N
*/
int TridiagonalMatrix::get_N() {
    return m_N;
}

/**
* get_h: getter function for private variable m_h
*
* @return the double value m_h
*/
double TridiagonalMatrix::get_h() {
    return m_h;
}

/**
* get_eigmat: getter function for private variable m_eigmat
*
* @param i Index specifying the i'th value in m_eigmat
* @param j Index specifying the j'th value in m_eigmat
* @return The specified value at index (i,j) in m_eigmat
*/
double TridiagonalMatrix::get_eigmat(int i, int j) {
    return m_eigmat(i,j);
}

/**
* get_toeplitz: getter function for private variable m_toeplitz
*
* @param i Index specifying the i'th value in m_toeplitz
* @param j Index specifying the j'th value in m_toeplitz
* @return The specified value at index (i,j) in m_toeplitz
*/
double TridiagonalMatrix::get_toeplitz(int i, int j) {
    return m_Toeplitz(i,j);
}

// PRIVATE MEMBER FUNCTIONS

/**
* to_stiring: sends a string representation of the supplied matrix to std.out
*
* @param matrix The matrix to be printed
* @return Nothing gets returned
*/
void TridiagonalMatrix::to_string(mat matrix) {
    for (int i = 0; i < m_N; i++) {
        cout << "[ ";
        for (int j = 0; j < m_N; j++) {
            cout << matrix(i,j) << ", ";
        }
        cout << "]" << endl;
    }
}
