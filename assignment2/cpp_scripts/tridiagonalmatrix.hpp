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
    mat m_Toeplitz, m_eigmat;
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
    void write_to_file(string fname);
};

class JacobiSolver : public TridiagonalMatrix {
private:

    double maxoffdiag(int *k, int *l);
    void rotate(int k, int l);
    int findlowesteigval();

public:
    void init(int N, double diag, double nondiag);
    int solve();
    void write_to_file(string fname);
    double unittest_maxoffdiag();
    void print_eigvals();
    void print_toeplitz();

};

#endif
