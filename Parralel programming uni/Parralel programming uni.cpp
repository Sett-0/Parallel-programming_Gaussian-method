#include <iostream>
#include <iomanip>
#include <random>
using namespace std;

void swap_rows(double** matrix, int row1, int row2) {
    double* tmp = matrix[row1];
    matrix[row1] = matrix[row2];
    matrix[row2] = tmp;
}

void print_matrix(double** matrix, const int& N, int limit_cols = 8, int limit_rows = 8) {
    for (int i = 0; i < N; i++) {
        if (i != limit_rows) {
            cout << "|\t";
            for (int j = 0; j < N + 1; j++) {
                if (j != limit_cols)
                    cout << setprecision(3) << matrix[i][j] << "\t";
                else {
                    cout << "\t...";
                    break;
                }
            }
            cout << "\t|\n";
        }
        else {
            cout << "|\t......\t";
            for (int j = 1; j < limit_cols + 1; j++)
                cout << "......\t";
            cout << "\t|\n";
            break;
        }
    }
    cout << endl;
}

void print_X(double* X, const int& N, int limit_cols = 8) {
    for (int i = 0; i < N - 1; i++) {
        if (i != limit_cols)
            cout << "X" << i + 1 << " = " << X[i] << ", ";
        else {
            cout << "..., ";
            break;
        }
    }
    cout << "X" << N << " = " << X[N - 1] << ".";
}

int main() {
    int N;
    double min = -100, max = 100;

    cout << "Select number of equations: ";
    cin >> N;

    // Generating random double in range
    random_device rd; // obtain a random number from hardware
    mt19937 gen(rd()); // seed the generator
    uniform_real_distribution<> distr(min, max); // define the range

    // Filling the matrix with random float numbers
    double** matrix = new double*[N];
    for (int i = 0; i < N; i++) {
        matrix[i] = new double[N + 1]; // N + 1 because of the Y vector
        for (int j = 0; j < N + 1; j++)
            matrix[i][j] = distr(gen);
    }

    // Forward elimination
    for (int k = 0; k < N - 1; k++) {
        // Searching for the maximum nonzero first multiplier
        double max_elem = 0;
        int index = k;
        for (int i = k; i < N; i++) {
            double first_elem = abs(matrix[i][k]);
            if (max_elem < first_elem) {
                max_elem = first_elem;
                index = i;
            }
        }
        if (max_elem == 0)
            return -1;

        // Making the row with the max element first
        if (index != 0)
            swap_rows(matrix, k, index);

        // Normalizing the matrix
        for (int i = k; i < N; i++) {
            for (int j = N; j >= k; j--) {
                matrix[i][j] /= matrix[i][k];
            }
        }

        // Substracting the first row from every other row
        for (int i = k + 1; i < N; i++) {
            for (int j = k; j < N + 1; j++) {
                matrix[i][j] -= matrix[k][j];
            }
        }
    }

    // Back substitution
    double* X = new double[N];
    X[N - 1] = matrix[N - 1][N] / matrix[N - 1][N - 1];
    for (int i = N - 2; i >= 0; i--) {
        double right_side = matrix[i][N];
        for (int k = N - 1; k > i; k--)
            right_side -= matrix[i][k] * X[k];
        X[i] = right_side / matrix[i][i];
    }

    // Showing results
    print_matrix(matrix, N);
    print_X(X, N);

    for (int i = 0; i < N; i++)
        delete[] matrix[i];
    delete[] matrix;
    delete[] X;
    return 0;
}