#include <ctime>
#include <iomanip>

#include "tridiagonalmatrix.hpp"

int main(int argc, char *argv[]) {

    if (argc < 3) {
        cout << "Missing command line input!" << endl;
        cout << "Supply N and algorithm when prompted\nN is n.o. mesh points" << endl;
        exit(1);
    }

    int N = atoi(argv[1]);
    string algo = string(argv[2]);
    string fname = algo + "_N_" + to_string(N) + ".txt";
    clock_t start, end;
    if (algo == "jacobi") {
        JacobiSolver my_solver;
        my_solver.init(N, 2, -1);
        start = clock();
        int transformations = my_solver.solve();
        end = clock();
        cout << setw(15) << setprecision(9) << "Time taken JacobiSolver: " << 1000.*((end - start)/(double)CLOCKS_PER_SEC) << "ms" << endl;
        my_solver.write_to_file(fname, transformations);
    } else if (algo == "arma") {
        ArmadilloSolver arma_solver;
        arma_solver.init(N, 2, -1);
        start = clock();
        arma_solver.solve();
        end = clock();
        cout << setw(15) << setprecision(9) << "Time taken ArmaSolver: " << 1000.*((end - start)/(double)CLOCKS_PER_SEC) << "ms" << endl;
        arma_solver.write_to_file(fname);
    } else if (algo == "analytical") {
        double h = 1./N;
        vec x = linspace(h, 1-h, N);
        ofstream ofile;
        double pi = acos(-1.0);
        ofile.open(fname);
        ofile << "Analytcial solution" << endl;
        for (int i = 1; i <= N; i++) {
            ofile << x(i-1) << " " << sin((i*1.*pi)/(N+1)) << endl;
        }
        ofile.close();
    } else if (algo == "jacobi-arma") {
        JacobiSolver my_solver;
        my_solver.init(N, 2, -1);
        int t = my_solver.solve();
        ArmadilloSolver arma_solver;
        arma_solver.init(N, 2, -1);
        arma_solver.solve();
        int low_idx = my_solver.find_lowest_eigval();
        ofstream ofile;
        ofile.open(fname);
        ofile << "abs(Jacobi - armadillo) eigenvector" << endl;
        for (int i = 0; i < N; i++) {
            ofile << my_solver.get_x(i) << " " << abs(my_solver.get_eigmat(i, low_idx) - arma_solver.get_eigmat(i, 0)) << endl;
        }
        ofile.close();
    } else if (algo == "one-electron") {
        JacobiSolver my_solver;
        my_solver.init(N, 2., -1., 5.0);
        my_solver.add_harmonic_potential();
        int rot = my_solver.solve();
        for (int i = 0; i < 4; i++) {
            fname = algo + "_N_" + to_string(N);
            my_solver.write_to_file_quantum(fname, rot, i);
        }
        my_solver.print_eigvals();
    } else if (algo == "two-electrons") {
        vec omega = {0.01, 0.5, 1., 5.};
        for (double i : omega) {
            JacobiSolver my_solver;
            fname = algo + "_N_" + to_string(N) + "_omega_" + to_string(i) + ".txt";
            my_solver.init(N, 2. , -1., 5.0);
            my_solver.add_harmonic_potential_two_electrons(i);
            int rot = my_solver.solve();
            my_solver.write_to_file(fname, rot);
        }

    }

    return 0;
}
