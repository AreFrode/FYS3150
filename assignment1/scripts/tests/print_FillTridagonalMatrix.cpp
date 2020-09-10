#include <iostream>

#include "../functions.hpp"

using namespace std;

int main() {
    int n = 5;
    double **m;
    double a, b, c;
    a = -1.; b = 2.; c = -1.;
    m = new double*[n];
    for (int i = 0; i < n; i++)
        m[i] = new double[n];

    FillTridiagonalMatrix(m, a, b, c, n);

    for (int i = 0; i < n; i++) {
        cout << "[";
        for (int j = 0; j < n; j++)
                cout << m[i][j] << ", ";

        cout << "]" << endl;
    }

    for (int i = 0; i < n; i++)
        delete [] m[i];
    delete [] m;

    return 0;
}
