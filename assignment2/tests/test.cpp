#include "../cpp_scripts/tridiagonalmatrix.hpp"

#define TOL 1e-8

void test_maxoffdiag();
void test_eigenvalues();

int main(int argc, char *argv[]) {
    test_maxoffdiag();
    cout << "test_maxoffdiag passed" << endl;

    test_eigenvalues();
    cout << "test_eigenvalues passed" << endl;

    return 0;
}

void test_maxoffdiag() {
    double computed, expected;
    JacobiSolver my_solver;
    my_solver.init(10, 2, -1);
    computed = my_solver.unittest_maxoffdiag();
    expected = 12.;
    if (abs(expected - computed) > TOL) {
        cout << "Excpected " << expected << " != " << computed << " Computed" << endl;
        cout << "test_maxoffdiag not passed" << endl;
        exit(1);
    }
}

void test_eigenvalues() {
    double expected;
    JacobiSolver my_solver;
    my_solver.init(20, 2, -1);
    my_solver.solve();
    double step = my_solver.get_h();
    int dim = my_solver.get_N();
    vec computed = zeros<vec>(dim);
    double pi = acos(-1.);
    double diag = 2.0/(step * step);
    double nondiag = -1.0/(step * step);
    for (int i = 0; i < dim; i++)
        computed(i) = my_solver.get_toeplitz(i,i);

    uvec indices = sort_index(computed);
    for (int i = 0; i < dim; i++) {
        expected = diag + 2*nondiag*cos((i+1)*pi/(dim+1));
        if (abs(computed(indices(i)) - expected) > TOL) {
            cout << "Excpected " << expected << " != " << computed << " Computed" << endl;
            cout << "test_eigenvalues not passed" << endl;
            exit(1);
        }
    }
}
