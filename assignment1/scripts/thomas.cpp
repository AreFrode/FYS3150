#include <iostream>
#include <cmath>

using namespace std;

double *initVector(int);
double **initMatrix(int);
void generalSolver(double *, double *, double *, double *, double *, int);

int main(int argc, char *argv[]) {
    int exponent;
    if (argc <= 1) {
        cout << "Missing command line argument!" << endl;
        cout << "Usage: " << argv[0] << " n\nn is the power in 10^n" << endl;
        exit(1);
    } else {
        exponent = atoi(argv[1]);
    }

    int n = pow(10, exponent);
    double h = 1.0/(n + 1);
    double *c = initVector(n);
    double *a = initVector(n);
    double *b = initVector(n);
    double *v = initVector(n+1);
    double *solution = initVector(n);

    v[0] = v[n] = 0;
    double x = h;
    for (int i = 0; i < n; i++) {
        c[i] = -1;
        a[i] = -1;
        b[i] = 2;
        solution[i] = h*h*100.0*exp(-10*x);
        x += h;
    }

    generalSolver(a, b, c, v, solution, n);

    cout << "[ ";
    for (int i = 0; i < n+1; i++)
        cout << v[i] << ", ";

    cout << "]" << endl;

    cout << "[ ";
    double tmp = h;
    for (int i = 0; i < 10; i++) {
        cout << 1 - tmp + tmp*exp(-10) - exp(-10*tmp) << ", "; // dette er skummelt...
        tmp += h;
    }

    cout << "]" << endl;

    delete [] c;
    delete [] a;
    delete [] b;
    delete [] solution;
    delete [] v;

    return 0;
}

void generalSolver(double *a, double *b, double *c, double *v, double *sol, int n) {
    double *diag, *f;
    diag = new double[n];
    f = new double[n];
    double weight;
    diag[0] = b[0];
    f[0] = sol[0];

    for (int i = 1; i < n; i++) {
        weight = a[i-1] / diag[i-1];
        diag[i] = b[i] - weight*c[i-1];
        f[i] = sol[i] - weight*f[i-1];
    }

    v[n-1] = f[n-1]/diag[n-1];
    for (int i = n-2; i >= 1; i--)
        v[i] = (f[i] - c[i]*v[i+1])/diag[i];

    delete [] diag;
    delete [] f;
}

double *initVector(int N) {
    return new double[N];
}

double **initMatrix(int N) {
    double **A;
    A = new double*[N];
    for (int i = 0; i < N; i++) {
        A[i] = new double[N];
    }
    return A;
}
