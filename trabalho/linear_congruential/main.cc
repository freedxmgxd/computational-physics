#include <iostream>
#include <vector>
#include <cstdlib> // For system()

using namespace std;

int main() {
    const int N = 100;
    const unsigned long long a = 1664525;
    const unsigned long long c = 1013904223;
    const unsigned long long m = 4294967296;
    unsigned long long x = 1;

    vector<double> results;

    for (int i = 0; i < N; ++i) {
        x = (a * x + c) % m;
        results.push_back(static_cast<double>(x) / m);
    }

    FILE* gnuplotPipe = popen("gnuplot -persistent", "w");
    if (!gnuplotPipe) {
        cerr << "Error: Could not open Gnuplot." << endl;
        return 1;
    }

    fprintf(gnuplotPipe, "set title 'Linear Congruential Generator Results'\n");
    fprintf(gnuplotPipe, "set xlabel 'Iteration'\n");
    fprintf(gnuplotPipe, "set ylabel 'Random Value'\n");
    fprintf(gnuplotPipe, "plot '-' using 1:2 with points title 'Random Numbers'\n");

    for (auto i = 0; i < results.size(); ++i) {
        fprintf(gnuplotPipe, "%d %f\n", i, results[i]);
    }
    fprintf(gnuplotPipe, "e\n");

    pclose(gnuplotPipe);

    return 0;
}
