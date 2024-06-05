#include <iostream>
#include <iomanip>
#include <random>
#include "mpi.h"
using namespace std;

void swap_rows(double** matrix, int row1, int row2) {
    double* tmp = matrix[row1];
    matrix[row1] = matrix[row2];
    matrix[row2] = tmp;
}

void swap_rows(double* arr, int row1, int row2) {
    double tmp = arr[row1];
    arr[row1] = arr[row2];
    arr[row2] = tmp;
}

void print_matrix(double** matrix, int N_rows, int N_cols = -1, int limit_cols = 10, int limit_rows = 10) {
    if (N_cols == -1) N_cols = N_rows;
    for (int i = 0; i < N_rows; i++) {
        if (i != limit_rows) {
            cout << "|\t";
            for (int j = 0; j < N_cols + 1; j++) {
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
    cout << endl << endl << endl;
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
    cout << "X" << N << " = " << X[N - 1] << ".\n\n";
}

void print_time(double const& start, double const& stop) {
    cout << "Time taken: " << setprecision(5) << stop - start << " seconds.\n";
}

int main() {
    int N, local_N, mod, rank, size, root = 0;
    double min = -100, max = 100;
    double start, stop;
    double** matrix, ** local_matrix, * matrix_data, * local_matrix_data;
    MPI_Comm comm = MPI_COMM_WORLD;

    // MPI Initialization
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    // Matrix initialization
    if (rank == root) {
        cout << "Select number of equations: ";
        cin >> N;
    }

    MPI_Bcast(&N, 1, MPI_INT, root, comm);
    local_N = N / size, mod = N % size;

    matrix_data = new double[N * (N + 1)];
    matrix = new double* [N];
    for (int i = 0; i < N; i++)
        matrix[i] = &(matrix_data[i * (N + 1)]);

    local_matrix_data = new double[local_N * (N + 1)];
    local_matrix = new double* [local_N];
    for (int i = 0; i < local_N; i++)
        local_matrix[i] = &(local_matrix_data[i * (N + 1)]);

    if (rank == root) {
        // Generating random double in range
        random_device rd; // obtain a random number from hardware
        mt19937 gen(rd()); // seed the generator
        uniform_real_distribution<> distr(min, max); // define the range

        // Filling the matrix with random float numbers
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N + 1; j++)
                matrix[i][j] = distr(gen);
        
        // Timer start
        start = MPI_Wtime();
    }

    MPI_Scatter(matrix_data, local_N * (N + 1), MPI_DOUBLE, 
      &(local_matrix[0][0]), local_N * (N + 1), MPI_DOUBLE, root, comm);

    // Forward elimination
    int shift = 0;  // Will be used for keeping track of already changed rows
    for (int k = 0; k < N; k++) {
        // Searching for the maximum nonzero first multiplier
        double local_max_elem = 0;
        int local_index = 0;
        for (int i = shift; i < local_N; i++) {
            double first_elem = abs(local_matrix[i][k]);
            if (local_max_elem < first_elem) {
                local_max_elem = first_elem;
                local_index = i;
            }
        }

        struct Package { double max_elem; int rank; };

        Package p_global { 0, 0 };
        Package p_local { local_max_elem, rank };

        MPI_Allreduce(&p_local, &p_global, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);

        if (rank == root && p_global.max_elem == 0)
            return -1;

        // Pivot row that will be send to all processes and used for subtraction
        double* pivot_row = new double[N + 1 - k];
        if (rank == p_global.rank)
            for (int i = 0; i < N + 1 - k; i++)
                pivot_row[i] = local_matrix[local_index][i + k];

        MPI_Bcast(pivot_row, N + 1 - k, MPI_DOUBLE, p_global.rank, comm);

        // Subtracting the pivot row from the every other
        if (rank == p_global.rank) {
            for (int i = shift; i < local_N; i++)
                if (i != local_index)
                    for (int j = N; j >= k; j--) {
                        local_matrix[i][j] /= local_matrix[i][k];
                        local_matrix[i][j] -= pivot_row[j - k] / pivot_row[0];
                    }
                else
                    for (int j = N; j >= k; j--)
                        local_matrix[i][j] /= pivot_row[0];
            swap_rows(local_matrix, local_index, shift);
            shift++;
        }
        else {
            for (int i = shift; i < local_N; i++)
                for (int j = N; j >= k; j--) {
                    local_matrix[i][j] /= local_matrix[i][k];
                    local_matrix[i][j] -= pivot_row[j - k] / pivot_row[0];
                }
        }

        delete[] pivot_row;
    }

    MPI_Gather(local_matrix_data, local_N * (N + 1), MPI_DOUBLE,
        &(matrix[0][0]), local_N * (N + 1), MPI_DOUBLE, root, comm);

    // Back substitution
    double* X = new double[N];
    double* local_X = new double[local_N];

    double x;
    bool found_x = false;
    shift = local_N - 1;

    for (int i = shift; i >= 0; i--) {
        local_X[i] = local_matrix[i][N];
        if (local_matrix[i][N - 1] == 1.0)
            found_x = true;
    }

    if (found_x) {
        x = local_X[shift--];
        for (int receiver = 0; receiver < size; receiver++)
            if (rank != receiver)
                MPI_Send(&x, 1, MPI_DOUBLE, receiver, 0, comm);
        found_x = false;
    }
    else
        MPI_Recv(&x, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, comm, MPI_STATUS_IGNORE);

    X[N - 1] = x;

    for (int k = N - 1; k > 0; k--) {
        for (int i = shift; i >= 0; i--) {
            local_X[i] -= local_matrix[i][k] * X[k];
            if (local_matrix[i][k - 1] == 1.0)
                found_x = true;
        }

        if (found_x) {
            x = local_X[shift--];
            //cout << k - 1 << ": \t" << x << endl;
            for (int receiver = 0; receiver < size; receiver++)
                if (rank != receiver)
                    MPI_Send(&x, 1, MPI_DOUBLE, receiver, 0, comm);
            found_x = false;
        }
        else
            MPI_Recv(&x, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, comm, MPI_STATUS_IGNORE);

        MPI_Barrier(comm);
        X[k - 1] = x;
    }

    // Showing results
    if (rank == root) {
        // Timer end
        stop = MPI_Wtime();

        // Sorting the matrix (just to make it look better)
        for (int i = 0; i < N - 1; i++)
            for (int j = i; j < N; j++)
                if (matrix[j][i] == 1.0) {
                    swap_rows(matrix, j, i);
                    break;
                }

        print_matrix(matrix, N);
        print_X(X, N);
        print_time(start, stop);
    }

    delete[] X;
    delete[] local_X;
    delete[] local_matrix_data;
    delete[] local_matrix;
    delete[] matrix_data;
    delete[] matrix;
    MPI_Finalize();
    return 0;
}