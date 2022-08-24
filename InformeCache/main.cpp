#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <stdlib.h>

using namespace std;
#define MAX_SIZE 1000

int main()
{
    vector<vector<double>> A(MAX_SIZE, vector<double>(MAX_SIZE, 0));
    vector<double> x(MAX_SIZE, 0), y(MAX_SIZE, 0);

    auto t_start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < MAX_SIZE; ++i) {
        for (int j = 0; j < MAX_SIZE; ++j) {
            y[i] += A[i][j] * x[j];
        }
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    double elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();
    std::cout << "Primer Caso: " << elapsed_time_ms << " ms" << "\n";

    t_start = std::chrono::high_resolution_clock::now();

    for (int j = 0; j < MAX_SIZE; ++j) {
        for (int i = 0; i < MAX_SIZE; ++i) {
            y[i] += A[i][j] * x[j];
        }
    }
    t_end = std::chrono::high_resolution_clock::now();
    elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();
    std::cout << "Segundo Caso: " << elapsed_time_ms << " ms"<<"\n";

    return 0;
}