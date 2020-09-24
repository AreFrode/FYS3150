#include "tridiagonalmatrix.hpp"

int main(int argc, char *argv[]) {
    // ArmadilloSolver my_solver;
    // my_solver.init(20, 2, -1);
    // my_solver.solve();

    JacobiSolver my_solver;
    my_solver.init(20, 2, -1);
    //my_solver.print_eigvec();
    my_solver.solve();
    //my_solver.print_toeplitz();
    my_solver.print_eigvec();

    ArmadilloSolver arma_solver;
    arma_solver.init(20, 2, -1);
    arma_solver.solve();




    return 0;
}
