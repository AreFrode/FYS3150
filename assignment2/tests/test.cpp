#include "../cpp_scripts/tridiagonalmatrix.hpp"

void test_maxoffdiag();

int main(int argc, char *argv[]) {
    test_maxoffdiag();

    return 0;
}

void test_maxoffdiag() {
    double computed, expected;
    JacobiSolver my_solver;
    my_solver.init(10, 2, -1);
    computed = my_solver.unittest_maxoffdiag();
    expected = 12.;
    cout << "computed: " << computed << ", expected: " << expected << endl;

    my_solver.print_toeplitz();
}
