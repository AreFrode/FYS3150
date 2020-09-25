#include <iomanip>

#include "tridiagonalmatrix.hpp"

void ArmadilloSolver::init(int N, double diag, double nondiag) {
    initialize(N, diag, nondiag);
    m_eigval = zeros<vec>(m_N);
}

void ArmadilloSolver::solve() {
    eig_sym(m_eigval, m_eigmat, m_Toeplitz);
    /*
    double pi = acos(-1.0);
    cout << "RESULTS: " << endl;
    cout << setiosflags(ios::showpoint | ios::uppercase);
    cout << "Number of Eigenvalues = " << m_N << endl;
    cout << "Exact versus numerical eigenvalues: " << endl;
    for (int i = 0; i < m_N; i++) {
        double exact = m_diag + 2*m_nondiag*cos((i+1)*pi/(m_N+1));
        cout << m_eigval[i] << ", " << exact << endl;
    }
    */
}

void ArmadilloSolver::write_to_file(string fname) {
    vec x = linspace(m_h, 1-m_h, m_N);
    m_ofile.open(fname);
    for (int i = 0; i < m_N; i++) {
        m_ofile << x(i) << " " << m_eigmat(i, 0) << endl;
    }
    m_ofile.close();
}
