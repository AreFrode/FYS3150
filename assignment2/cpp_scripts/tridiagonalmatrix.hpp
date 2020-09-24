#ifndef TridiagonalMatrix_hpp
#define TridiagonalMatrix_hpp
#include <fstream>
#include <armadillo>
#include <cmath>

using namespace std;
using namespace arma;

class TridiagonalMatrix {
protected:
    int m_N;
    double m_h, m_diag, m_nondiag;
    mat m_Toeplitz;
    ofstream m_ofile;

    void print(mat matrix);

public:
    void initialize(int N, double diag, double nondiag);
};

class ArmadilloSolver : public TridiagonalMatrix {
private:
    vec m_eigval;

public:
    void init(int N, double diag, double nondiag);
    void solve();
};

class JacobiSolver : public TridiagonalMatrix {
private:
    mat m_eigmat;

    double maxoffdiag(int *k, int *l);
    void rotate(int k, int l);
    int findlowesteigval();

public:
    void init(int N, double diag, double nondiag);
    void solve();
    void print_eigvec();
    double unittest_maxoffdiag();
    void print_eigvals();
    void print_toeplitz();

};

#endif
