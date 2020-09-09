#ifndef FUNCTIONS_H
#define FUNCTIONS_H

void GeneralSolver(double* a, double* b, double* c, double* v, double* sol, int n);

void SpecializedSolver(double a, double b, double c, double *v, double *sol, int n);

double f(double x);

double u(double x);

double MaxError(double* v, double *u, int n);

double *FillArray(int value, int N);

void FillTridiagonalMatrix(double **matrix, int a, int b, int c, int n);

#endif
