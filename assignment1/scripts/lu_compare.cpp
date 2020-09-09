#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <ctime>

#include "lib.h"
#include "functions.hpp"

using namespace std;

ofstream ofile;

int main(int argc, char *argv[]) {

    string filename = "../results/lu_results.csv";

    ofile.open(filename);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << "Matrix size, Runtime [ms], max(relative_error)" << endl;

    for (int n = 1e1; n < 1e4; n *= 10) {
        void **ptr = matrix(n, n, sizeof(double));
        double **matr = (double**)ptr;

        FillTridiagonalMatrix(matr, -1, 2, -1, n);
        int *index = new int[n];
        double d, h = 1./(n + 1.);
        double *u_vector = FillArray(0., n);
        double *b = FillArray(0., n);
        double xval = h;
        for (int i = 0; i < n; i++) {
            b[i] = h*h*f(xval);
            u_vector[i] = u(xval);
            xval += h;
        }

        time_t end, start = clock();
        ludcmp(matr, n, index, &d);
        lubksb(matr, n, index, b);
        end = clock();

        double error = MaxError(b, u_vector, n);

        ofile << setw(1) << setprecision(5) << scientific << n << "x" << n << ",";
        ofile << setw(5) << setprecision(3) << fixed << 1000*(double)(end - start)/CLOCKS_PER_SEC << ",";
        ofile << setw(5) << setprecision(5) << error << endl;

        delete [] b;
        delete [] u_vector;
        delete [] index;
        free_matrix(ptr);
    }

    return 0;
}
