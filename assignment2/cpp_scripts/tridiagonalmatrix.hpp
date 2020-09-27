/**
* tridiagonalmatrix.hpp: header file for TridiagonalMatrix, JacobiSolver and ArmadilloSolver
*
* Author: Are Frode Kvanum
*
* Completion Date: 26.09.2020
*/

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
    vec m_x;
    mat m_Toeplitz, m_eigmat;
    ofstream m_ofile;

    void to_string(mat matrix);

public:
    void initialize(int N, double diag, double nondiag, double rho_max = 1.0);
    void add_harmonic_potential();
    double get_x(int i);
    int get_N();
    double get_h();
    double get_eigmat(int i , int j);
    double get_toeplitz(int i, int j);
};

class JacobiSolver : public TridiagonalMatrix {
private:

    double maxoffdiag(int *k, int *l);
    void rotate(int k, int l);

public:
    void init(int N, double diag, double nondiag, double rho_max = 1.0);
    int solve();
    void write_to_file(string fname, int transform);
    int find_lowest_eigval();
    void print_eigvals();
    void print_toeplitz();
    double unittest_maxoffdiag();

};

class ArmadilloSolver : public TridiagonalMatrix {
private:
    vec m_eigval;

public:
    void init(int N, double diag, double nondiag, double rho_max = 1.0);
    void solve();
    void print_eigvals();
    void write_to_file(string fname);
};

#endif
