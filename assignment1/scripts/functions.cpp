#include <cmath>

#include "functions.hpp"

void GeneralTridiagSolver(double *a, double *b, double *c, double *v, double *sol, int n) {
    /**
    * Solves a system with a general tridiagonal matrix
    *
    *@param a Lower diagonal elements in matrix
    *@param b Diagnoal elements in matrix
    *@param c Upper diagonal elements in matrix
    *@param v l.h.s solution vector
    *@param sol r.h.s closed form solution vector
    *@param n Number of gridpoints
    *@return fills the vector v with the solution
    */
    double *diag, *temp;
    diag = new double[n];
    temp = new double[n];
    double weight;
    diag[0] = b[0];
    temp[0] = sol[0];

    for (int i = 1; i < n; i++) {
        diag[i] = b[i] - (a[i-1]*c[i-1])/diag[i-1];
        temp[i] = sol[i] - (a[i-1]*temp[i-1])/diag[i-1];
    }

    v[n-1] = temp[n-1]/diag[n-1];
    for (int i = n-2; i >= 0; i--)
        v[i] = (temp[i] - c[i]*v[i+1])/diag[i];

    delete [] diag;
    delete [] temp;
}

void SpecializedTridiagSolver(double a, double b, double c, double *v, double *sol, int n) {
    /**
    * Solves a specialized system where the diagonal elements are equal (but different bewteen diags)
    *
    *@param a Lower diagonal elements in matrix
    *@param b Diagnoal elements in matrix
    *@param c Upper diagonal elements in matrix
    *@param v l.h.s solution vector
    *@param sol r.h.s closed form solution vector
    *@param n Number of gridpoints
    *@return fills the vector v with the solution
    */
    double *diag, *temp;
    diag = new double[n];
    temp = new double[n];
    diag[0] = b;
    temp[0] = sol[0];
    double ac = a*c;

    for (int i = 1; i < n; i++) {
        diag[i] = b - (ac)/diag[i-1];
        temp[i] = sol[i] - (a*temp[i-1])/diag[i-1];
    }

    v[n-1] = temp[n-1]/diag[n-1];
    for (int i = n-2; i >= 0; i--)
        v[i] = (temp[i] - c*v[i+1])/diag[i];

    delete [] diag;
    delete [] temp;
}

double f(double x) {
    /**
    *Explicit source term
    *
    *@param x the value to be solved for
    *@return f evaluated at x
    */
    return 100.*exp(-10.*x);
}

double u(double x) {
    /**
    *Explicit closed form solution
    *
    *@param x the value to be solved for
    *@return u evaluated at x
    */
    return 1.-(1-exp(-10.))*x-exp(-10.*x);
}

double MaxError(double *v, double *u, int n) {
    /**
    *Extracts the max value of the realtive error between v and u
    *
    *@param v numerical solution vector
    *@param u closed form analytic solution
    *@param n length of v and u
    *@return the maximum error in the data set
    */
    double eps, error = log10(fabs((v[0] - u[0])/u[0]));
        for (int i = 1; i < n; i++) {
            eps = log10(fabs((v[i] - u[i])/u[i]));
            if (eps > error)
                error = eps;
        }

    return eps;
}

double *FillArray(int value, int N) {
    /**
    *Initializes an array with constant and equal elements
    *
    *@param value The value assigned to all elements
    *@param N The length of the array
    *@return An initialized array
    */
    double *outvec = new double[N];
    for (int i = 0; i < N; i++) outvec[i] = value;
    return outvec;
}

void FillTridiagonalMatrix(double **matrix, int a, int b, int c, int n) {
    /**
    *Fills a tridiagonal matrix with constant values along the diagonals
    *
    *@param matrix 2d matrix
    *@param a value along the lower diagonal
    *@param b value along the diagonal
    *@param c value along the upper diagonal
    *@param n length of array (nxn)
    *@return nothing is returned, but matrix is filled with values
    */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                matrix[i][j] = b;
            } else if (j == (i+1)) {
                matrix[i][j] = c;
            } else if (j == (i-1)) {
                matrix[i][j] = a;
            }
        }
    }
}
