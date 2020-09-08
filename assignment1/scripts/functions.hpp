#ifndef FUNCTIONS_H
#define FUNCTIONS_H

void GeneralSolver(double* a, double* b, double* c, double* v, double* sol, int n);

void SpecializedSolver(double *a, double *b, double *c, double *v, double *sol, int n);

double f(double x);

double u(double x);

double *InitArray(int value, int N);

#endif
