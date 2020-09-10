#include <iostream>
#include <cmath>

#include "../functions.hpp"

using namespace std;

int main() {
    int n = pow(10, 1);
    double h = 1./(n + 1.);
    double *c = FillArray(-1., n-1);
    double *a = FillArray(-1., n-1);
    double *b = FillArray(2., n);
    double *v1 = FillArray(0., n);
    double *v2 = FillArray(0., n);
    double *solution = FillArray(0., n);

    double xval = h;
    for (int i = 0; i < n; i++) {
        solution[i] = h*h*f(xval);
        xval += h;
    }

    GeneralTridiagSolver(a, b, c, v1, solution, n);
    SpecializedTridiagSolver(a[0], b[0], c[0], v2, solution, n);

    cout << "Solution for General Tridiag solver, n=10" << endl << "[";
    for (int i = 0; i < n; i++) {
        cout << v1[i] << ", ";
    }
    cout << "]" << endl << endl;

    cout << "Solution for Specialized Tridiag solver, n=10" << endl << "[";
    for (int i = 0; i < n; i++) {
        cout << v2[i] << ", ";
    }
    cout << "]" << endl;


    delete [] c;
    delete [] a;
    delete [] b;
    delete [] v1;
    delete [] v2;
    delete [] solution;

    return 0;
}
