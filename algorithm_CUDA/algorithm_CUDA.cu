#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <random>

using namespace std;

__device__ void swap_rows(double* device_matrix_data, int N, int index, int k, int size) {
    int t_index = threadIdx.x;
    double tmp;
    for (int i = t_index; i < N + 1; i += size) {
        tmp = device_matrix_data[index + i];
        device_matrix_data[index + i] = device_matrix_data[k * (N + 1) + i];
        device_matrix_data[k * (N + 1) + i] = tmp;
    }
}

__global__ void find_pivot_row(double* device_matrix_data, int N, int k, int N_rows_per_thread, int size) {
    // Searching for the row with maximum nonzero first element and moving this row to the top
    extern __shared__ double s[];
    double* max_elems = s;
    int* indexes = (int*)&s[size];

    int t_index = threadIdx.x;
    int row_start = (k + N_rows_per_thread * t_index) * (N + 1);
    int row_stop = (k + N_rows_per_thread * (t_index + 1)) * (N + 1);
    if (row_start > (N - 1) * (N + 1))
        return;
    if (row_stop > N * (N + 1))
        row_stop = N * (N + 1);

    double max_elem = 0;
    int index = 0;
    for (int i = row_start; i < row_stop; i += N + 1) {
        double first_elem = abs(device_matrix_data[i + k]);
        if (max_elem < first_elem) {
            max_elem = first_elem;
            index = i;
        }
    }

    max_elems[t_index] = max_elem;
    indexes[t_index] = index;
    __syncthreads();
    
    for (int i = size / 2; i > 0; i /= 2) {
        if (t_index < i) {
            if (max_elems[i] < max_elems[t_index + i]) {
                max_elems[i] = max_elems[t_index + i];
                indexes[i] = indexes[t_index + i];
            }
        }
        __syncthreads();
    }

    double a = device_matrix_data[indexes[0] + k];
    for (int i = t_index; i < N + 1; i += size)
        device_matrix_data[indexes[0] + k + i] /= a;
    swap_rows(device_matrix_data, N, indexes[0], k, size);
}

__global__ void subtract_rows(double* device_matrix_data, int N, int k, int size) {
    // Each thread calculates one row at a time
    int t_index = threadIdx.x + blockDim.x * blockIdx.x;
    if (t_index > N) return;

    for (int i = k + 1 + t_index; i < N; i += size) {
        //normalize_row(device_matrix_data, N, i * (N + 1), k, 1);
        double a = device_matrix_data[i * (N + 1) + k];
        for (int j = k; j < N + 1; j++) {
            device_matrix_data[i * (N + 1) + j] /= a;
            device_matrix_data[i * (N + 1) + j] -= device_matrix_data[k * (N + 1) + j];
        }
    }
}

__global__ void init_X(double* device_matrix_data, double* device_X, int N, int size) {
    int t_index = threadIdx.x + blockDim.x * blockIdx.x;
    if (t_index > N) return;

    for (int i = t_index; i < N; i += size)
        device_X[i] = device_matrix_data[i * (N + 1) + N];
}

__global__ void calculate_X(double* device_matrix_data, double* device_X, int N, int k, int size) {
    int t_index = threadIdx.x + blockDim.x * blockIdx.x;
    if (t_index > N) return;
    
    for (int i = t_index; i < k; i += size)
        device_X[i] -= device_matrix_data[i * (N + 1) + k] * device_X[k];
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
    cout << "X" << N << " = " << X[N - 1] << ".\n";
}

void print_time(double const& time_taken) {
    cout << "Time taken: " << setprecision(5) << 1e-3 * time_taken << " seconds.\n";
}

int main()
{
    int N;
    double min = -100, max = 100;
    double** matrix, * matrix_data, ** device_matrix, * device_matrix_data;
    float time_taken;
    int size;

    cout << "Select number of equations: ";
    cin >> N;

    // Generating random double in range
    random_device rd; // obtain a random number from hardware
    mt19937 gen(0);//rd()); // seed the generator
    uniform_real_distribution<> distr(min, max); // define the range

    // Filling the matrix with random float numbers
    matrix_data = new double[N * (N + 1)];
    matrix = new double* [N];
    for (int i = 0; i < N; i++)
        matrix[i] = &(matrix_data[i * (N + 1)]);

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N + 1; j++)
            matrix[i][j] = distr(gen);

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    // Timer start
    cudaEventRecord(start, 0);

    cudaMalloc((void**)&device_matrix_data, N * (N + 1) * sizeof(double));
    cudaMemcpy(device_matrix_data, matrix_data, N * (N + 1) * sizeof(double), cudaMemcpyHostToDevice);

    dim3 N_threads(128);
    int N_blocks;

    // Forward elimination
    for (int k = 0; k < N; k++) {
        const int N_rows_per_thread = 16;
        int N_threads_per_row = (N - k) / N_rows_per_thread;
        if ((N - k) % N_rows_per_thread != 0) N_threads_per_row++;
        int smem_size = N_threads_per_row * sizeof(double) + N_threads_per_row * sizeof(int);
        N_blocks = 1;
        int size = N_threads_per_row;


        find_pivot_row<<< N_blocks, N_threads_per_row, smem_size >>>(device_matrix_data, N, k, N_rows_per_thread, size);

        N_blocks = (N - k) / N_threads.x;
        if ((N - k) % N_threads.x != 0) N_blocks++;
        size = N_blocks * N_threads.x;

        subtract_rows<<< N_blocks, N_threads.x >>>(device_matrix_data, N, k, size);

    }

    // Back substitution
    double* X = new double[N];
    double* device_X;
    cudaMalloc((void**)&device_X, N * sizeof(double));

    N_blocks = N / N_threads.x;
    if (N % N_threads.x != 0) N_blocks++;
    size = N_blocks * N_threads.x;
    init_X <<< N_blocks, N_threads.x >>> (device_matrix_data, device_X, N, size);

    for (int k = N - 1; k >= 0; k--) {
        N_blocks = k / N_threads.x;
        if (k % N_threads.x != 0) N_blocks++;
        size = N_blocks * N_threads.x;
        calculate_X<<< N_blocks, N_threads.x >>>(device_matrix_data, device_X, N, k, size);
    }

    cudaMemcpy(X, device_X, N * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(matrix_data, device_matrix_data, N * (N + 1) * sizeof(double), cudaMemcpyDeviceToHost);

    // Timer end
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&time_taken, start, stop);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);


    // Showing results
    print_matrix(matrix, N);
    print_X(X, N);
    print_time(time_taken);

    delete[] X;
    cudaFree(device_X);
    cudaFree(device_matrix_data);
    delete[] matrix_data;
    delete[] matrix;
    return 0;
}
