#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>

#include "functions.hpp"

using namespace std;

ofstream ofile;

int main(int argc, char *argv[]) {
    int exponent;
    string filename;
    if (argc <= 1) {
        cout << "Missing command line argument!" << endl;
        cout << "Usage: " << argv[0] << " filename n\nn is the power in 10^n" << endl;
        exit(1);
    } else {
        filename = argv[1];
        exponent = atoi(argv[2]);
    }

    filename.append(to_string(exponent) + ".txt");

    int n = pow(10., exponent);
    double h = 1./(n + 1.);
    double *c = InitArray(-1., n-1);
    double *a = InitArray(-1., n-1);
    double *b = InitArray(2., n);
    double *v = InitArray(0., n);
    double *solution = InitArray(0., n);
    double *x = InitArray(0., n);

    double xval = h;
    for (int i = 0; i < n; i++) {
        solution[i] = h*h*f(xval);
        x[i] = xval;
        xval += h;
    }

    GeneralSolver(a, b, c, v, solution, n);

    ofile.open(filename);
    ofile << setiosflags(ios::showpoint | ios::uppercase);

    for (int i = 0; i < n; i++) {
        ofile << setw(15) << setprecision(8) << x[i];
        ofile << setw(15) << setprecision(8) << v[i];
        ofile << setw(15) << setprecision(8) << u(x[i]) << endl;
    }
    ofile.close();

    delete [] c;
    delete [] a;
    delete [] b;
    delete [] solution;
    delete [] v;

    return 0;
}
