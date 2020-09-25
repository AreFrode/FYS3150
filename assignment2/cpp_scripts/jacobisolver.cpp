#include "tridiagonalmatrix.hpp"

void JacobiSolver::init(int N, double diag, double nondiag) {
    initialize(N, diag, nondiag);
}

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

void JacobiSolver::write_to_file(string fname) {
    vec x = linspace(m_h, 1-m_h, m_N);
    int low_idx = findlowesteigval();
    m_ofile.open(fname);
    for (int i = 0; i < m_N; i++) {
        m_ofile << x(i) << " " << m_eigmat(i, low_idx) << endl;
    }
    m_ofile.close();
}

int JacobiSolver::findlowesteigval() {
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

void JacobiSolver::print_eigvals() {
    cout << "[";
    for (int i = 0; i < m_N; i++)
        cout << m_Toeplitz(i,i) << ", ";
    cout << "]" << endl;
}

void JacobiSolver::print_toeplitz() {
    print(m_Toeplitz);
}
