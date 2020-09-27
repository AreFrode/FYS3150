#include "../cpp_scripts/tridiagonalmatrix.hpp"

int main(int argc, char *argv[]) {
    JacobiSolver my_solver;
    my_solver.init(5, 2, -1);
    int rot = my_solver.solve();
    cout << "Diagonalized after " << rot << " similarity transformations" << endl;
    cout << "Tridagonal matrix after reduction" << endl;
    my_solver.print_toeplitz();
    int low = my_solver.find_lowest_eigval();
    cout << "Lowest eigenvector corresponding to eigenvalue: " << my_solver.get_toeplitz(low, low) << endl;
    cout << "[";
    for (int i = 0; i < my_solver.get_N(); i++)
        cout << my_solver.get_eigmat(i, low) << ", ";
    cout << "]" << endl;

    return 0;
}
