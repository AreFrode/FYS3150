#include <iostream>
#include <cmath>

using namespace std;

double *initVector(int);
double **initMatrix(int);
void destroyMatrix(double, int);

int main(int argc, char *argv[]) {
    int exp;
    if (argc <= 1) {
        cout << "Missing command line argument!" << endl;
        cout << "Usage: " << argv[0] << " n\nn is the power in 10^n" << endl;
        exit(1);
    } else {
        exp = atoi(argv[1]);
    }

    int n = pow(10, exp);
    double *upper = initVector(n);
    double *lower = initVector(n);
    double *diag = initVector(n);
    double *solution = initVector(n+1);
    double *x = initVector(n+1);

    x[0] = 0;
    x[n] = 1;
    solution[0] = solution[n] = 0;

    for (int i = 0; i < n; i++) {
        upper[i] = -1;
        lower[i] = -1;
        diag[i] = 2;
    }

    double w;
    for (int i = 1; i < n; i++) {
        w = lower[i] / diag[i-1];
        diag[i] = diag[i] - w*upper[i-1];
        solution[i] = solution[i] - w*diag[i-1];
    }

    for (int i = n-1; i >= 1; i--)
        x[i] = (solution[i] - upper[i]*x[i+1])/diag[i];

    cout << "[ ";
    for (int i = 0; i < n+1; i++)
        cout << x[i] << ", ";

    cout << "]" << endl;

    cout << "[ ";
    double tmp;
    double anal;
    for (int i = 0; i < 10; i++) {
        tmp = 0.1*i;
        cout << 100*pow(2.7, -10*tmp) << ", "; // dette er skummelt...
    }
    cout << "]" << endl;

    delete [] upper;
    delete [] lower;
    delete [] diag;
    delete [] solution;
    delete [] x;

    return 0;
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

void destroyMatrix(double **A, int N) {
    for (int i = 0; i < N; i++)
        delete [] A[i];
    delete [] A;
}
