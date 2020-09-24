#include <iomanip>

#include "tridiagonalmatrix.hpp"

void ArmadilloSolver::init(int N, double diag, double nondiag) {
    initialize(N, diag, nondiag);
    m_eigval = zeros<vec>(m_N);
}

void ArmadilloSolver::solve() {
    eig_sym(m_eigval, m_Toeplitz);
    double pi = acos(-1.0);
    cout << "RESULTS: " << endl;
    cout << setiosflags(ios::showpoint | ios::uppercase);
    cout << "Number of Eigenvalues = " << m_N << endl;
    cout << "Exact versus numerical eigenvalues: " << endl;
    for (int i = 0; i < m_N; i++) {
        double exact = m_diag + 2*m_nondiag*cos((i+1)*pi/(m_N+1));
        cout << m_eigval[i] << ", " << exact << endl;
    }
}
