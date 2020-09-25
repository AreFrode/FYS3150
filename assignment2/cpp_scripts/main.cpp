#include "tridiagonalmatrix.hpp"

int main(int argc, char *argv[]) {

    if (argc < 3) {
        cout << "Missing command line input!" << endl;
        cout << "Supply N and algorithm when prompted\nN is n.o. mesh points" << endl;
        exit(1);
    }
    // ArmadilloSolver my_solver;
    // my_solver.init(20, 2, -1);
    // my_solver.solve();

    int N = atoi(argv[1]);
    string algo = string(argv[2]);
    string fname = algo + "_N_" + to_string(N) + ".txt";

    if (algo == "jacobi") {
        JacobiSolver my_solver;
        my_solver.init(N, 2, -1);
        int no_trans = my_solver.solve();
        my_solver.write_to_file(fname);
    } else if (algo == "arma") {
        ArmadilloSolver arma_solver;
        arma_solver.init(N, 2, -1);
        arma_solver.solve();
        arma_solver.write_to_file(fname);
    }

    return 0;
}
