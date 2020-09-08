#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <ctime>

#include "functions.hpp"

using namespace std;

ofstream ofile;

int main(int argc, char *argv[]) {
    int range, loop, exponent = 1;
    string filename;
    if (argc <= 3) {
        cout << "Missing command line argument!" << endl;
        cout << "Usage: " << argv[0] << "filename n loop\nn is the power in 10^n, loop is amount of runs to calculate mean from" << endl;
        exit(1);
    } else {
        filename = argv[1];
        range = atoi(argv[2]);
        loop = atoi(argv[3]);
    }

    int n = pow(10., exponent);
    filename.append(".csv");

    ofile.open(filename);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << "Grid size, General solver [ms], Specialized solver [ms]" << endl;

    while (exponent < range) {
        double h = 1./(n + 1.);
        double *c = InitArray(-1., n-1);
        double *a = InitArray(-1., n-1);
        double *b = InitArray(2., n);
        double *v1 = InitArray(0., n);
        double *v2 = InitArray(0., n);
        double *solution = InitArray(0., n);
        double *x = InitArray(0., n);

        double xval = h;
        for (int i = 0; i < n; i++) {
            solution[i] = h*h*f(xval);
            x[i] = xval;
            xval += h;
        }

        time_t start1, start2, end1, end2;
        double mean1, mean2;

        for (int i = 0; i < loop; i++) {
            start1 = clock();
            GeneralSolver(a, b, c, v1, solution, n);
            end1 = clock();
            mean1 += (double)(end1 - start1)/CLOCKS_PER_SEC;
            start2 = clock();
            SpecializedSolver(a, b, c, v2, solution, n);
            end2 = clock();
            mean2 += (double)(end2 - start2)/CLOCKS_PER_SEC;
        }
        mean1 /= loop;
        mean2 /= loop;

        ofile << setw(5) << setprecision(5) << "10e" << exponent << ",";
        ofile << setw(5) << setprecision(5) << 1000*mean1 << ",";
        ofile << setw(5) << setprecision(5) << 1000*mean2 << endl;
        //cout << "Running General solver with 10^" << exponent << "x10^" << exponent << " grid points, mean time from " << loop << " loops: " << 1000*mean1 << "ms" << endl;
        //cout << "Running Specialized solver with 10^" << exponent << "x10^" << exponent << " grid points, mean time from " << loop << " loops. " << 1000*mean2 << "ms" << endl;


        delete [] c;
        delete [] a;
        delete [] b;
        delete [] solution;
        delete [] v1;
        delete [] v2;

        exponent++;
        n = pow(10., exponent);
    }

    ofile.close();

    return 0;
}

