#include "tridiagonalmatrix.hpp"

void TridiagonalMatrix::initialize(int N, double diag, double nondiag) {
    m_N = N;
    m_h = 1.0/m_N;
    m_Toeplitz = zeros<mat>(m_N,m_N);
    m_diag = diag/(m_h*m_h);
    m_nondiag = nondiag/(m_h*m_h);

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

void TridiagonalMatrix::print(mat matrix) {
    for (int i = 0; i < m_N; i++) {
        cout << "[ ";
        for (int j = 0; j < m_N; j++) {
            cout << matrix(i,j) << ", ";
        }
        cout << "]" << endl;
    }
}
