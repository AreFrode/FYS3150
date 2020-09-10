#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>

#include "functions.hpp"

using namespace std;

ofstream ofile;

int main(int argc, char *argv[]) {
    string filename = "../results/error.csv";
    int exponent = 1;
    int n = pow(10., exponent);

    double a, b, c;
    a = c = -1;
    b = 2;

    ofile.open(filename);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << "Step size, max(relative error)" << endl;

    while (exponent < 8) {
        double h = 1./(n + 1.);
        double *v = FillArray(0., n);
        double *solution = FillArray(0., n);
        double *x = FillArray(0., n);
        double *u_vector = FillArray(0., n);

        double xval = h;
        for (int i = 0; i < n; i++) {
            u_vector[i] = u(xval);
            solution[i] = h*h*f(xval);
            x[i] = xval;
            xval += h;
        }

        SpecializedTridiagSolver(a, b, c, v, solution, n);
        double relative_error = MaxError(v, u_vector, n);

        ofile << setw(5) << setprecision(1) << scientific << 1./n << ",";
        ofile << setw(5) << setprecision(5) << fixed << relative_error << endl;

        exponent++;
        n = pow(10., exponent);

        delete [] v;
        delete [] solution;
        delete [] x;
        delete [] u_vector;
    }

    ofile.close();

    return 0;
}
