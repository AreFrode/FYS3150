#include <iostream>
#include <cmath>

#include "../functions.hpp"

#define TOL 1e-10

using namespace std;

int test_f();
int test_u();
int test_FillArray();
int test_GeneralTridiagSolver();
int test_SpecializedTridiagSolver();

int main(int argc, char *argv[]) {
    if (test_f() != 0) {
        cout << "test_f() did not pass" << endl;
        exit(1);
    }
    cout << "test_f passed" << endl;

    if (test_u() != 0) {
        cout << "test_u() did not pass" << endl;
        exit(1);
    }
    cout << "test_u passed" << endl;

    if (test_FillArray() != 0) {
        cout << "test_FillArray() did not pass" << endl;
        exit(1);
    }
    cout << "test_FillArray passed" << endl;

    if (test_GeneralTridiagSolver() != 0) {
        cout << "test_GeneralTridiagSolver() did not pass" << endl;
        exit(1);
    }
    cout << "test_GeneralTridiagSolver passed" << endl;

    if (test_SpecializedTridiagSolver() != 0) {
        cout << "test_SpecializedTridiagSolver() did not pass" << endl;
        exit(1);
    }
    cout << "test_SpecializedTridiagSolver passed" << endl;

    cout << "All tests passed" << endl;
    return 0;
}

int test_f() {
    double expected = 100.;
    double computed = f(0);
    if (fabs(expected - computed) > TOL) return 1;

    return 0;
}

int test_u() {
    double expected2, expected1 = 0;
    double computed2, computed1 = u(0);
    computed2 = u(1);
    if (fabs(expected1 - computed1) > TOL || fabs(expected1 - computed2) > TOL) return 1;

    return 0;
}

int test_FillArray() {
    int n = 10;
    double *a = FillArray(2., n);
    for (int i = 1; i < n; i++)
        if (a[i] != a[0]) return 1;

    delete [] a;
    return 0;
}

int test_SpecializedTridiagSolver() {
    int n = pow(10., 2);
    double h = 1./(n + 1.);
    double a, b, c;
    a = -1; b = 2; c = -1;
    double *v = FillArray(0., n);
    double *solution = FillArray(0., n);
    double *uvec = FillArray(0., n);

    double xval = h;
    for (int i = 0; i < n; i++) {
        solution[i] = h*h*f(xval);
        uvec[i] = u(xval);
        xval += h;
    }

    SpecializedTridiagSolver(a, b, c, v, solution, n);

    for (int i = 0; i < n; i++) {
        if (fabs(v[i] - uvec[i]) > 1e-3) return 1;
    }

    delete [] v;
    delete [] solution;
    delete [] uvec;
    return 0;
}

int test_GeneralTridiagSolver() {
    int n = pow(10., 2);
    double h = 1./(n + 1.);
    double *c = FillArray(-1., n-1);
    double *a = FillArray(-1., n-1);
    double *b = FillArray(2., n);
    double *v = FillArray(0., n);
    double *solution = FillArray(0., n);
    double *uvec = FillArray(0., n);

    double xval = h;
    for (int i = 0; i < n; i++) {
        solution[i] = h*h*f(xval);
        uvec[i] = u(xval);
        xval += h;
    }

    GeneralTridiagSolver(a, b, c, v, solution, n);

    for (int i = 0; i < n; i++) {
        if (fabs(v[i] - uvec[i]) > 1e-3) return 1;
    }

    delete [] c;
    delete [] a;
    delete [] b;
    delete [] v;
    delete [] solution;
    delete [] uvec;
    return 0;
}
