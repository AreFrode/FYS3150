#include "functions.hpp"
#include <cmath>

void GeneralSolver(double *a, double *b, double *c, double *v, double *sol, int n) {
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
        weight = a[i-1] / diag[i-1];
        diag[i] = b[i] - weight*c[i-1];
        temp[i] = sol[i] - weight*temp[i-1];
    }

    v[n-1] = temp[n-1]/diag[n-1];
    for (int i = n-2; i >= 0; i--)
        v[i] = (temp[i] - c[i]*v[i+1])/diag[i];

    delete [] diag;
    delete [] temp;
}

void SpecializedSolver(double *a, double *b, double *c, double *v, double *sol, int n) {
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
    double b_const = diag[0] = b[0];
    temp[0] = sol[0];
    double ac_const = a[0]*c[0];

    for (int i = 1; i < n; i++) {
        diag[i] = b_const - (ac_const)/diag[i-1];
        temp[i] = sol[i] - (a[i-1]*temp[i-1])/diag[i-1];
    }

    v[n-1] = temp[n-1]/diag[n-1];
    for (int i = n-2; i >= 0; i--)
        v[i] = (temp[i] - c[i]*v[i+1])/diag[i];

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

double *InitArray(int value, int N) {
    /**
    *Initializes an array with all equal elements
    *
    *@param value The value assigned to all elements
    *@param N The length of the array
    *@return An initialized array
    */
    double *outvec = new double[N];
    for (int i = 0; i < N; i++) outvec[i] = value;
    return outvec;
}
