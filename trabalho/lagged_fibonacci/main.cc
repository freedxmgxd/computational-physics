#include <iostream>
#include <vector>
#include <cstdlib> // For system()

using namespace std;

class LaggedFibonacciGenerator {
private:
    vector<unsigned long> laggedFib;
    int j, k;
    unsigned long m;
    int currentIndex;

public:
    LaggedFibonacciGenerator(int j, int k, unsigned long m, unsigned long seed, int initSize)
        : j(j), k(k), m(m), currentIndex(0) {
        laggedFib.resize(k, seed);
        // Initialize the sequence with the seed value
        for (int i = 1; i < k; ++i) {
            laggedFib[i] = (laggedFib[i - 1] * 69069 + 1) % m; // Arbitrary initializer
        }
    }

    unsigned long next() {
        int indexJ = (currentIndex + k - j) % k;
        int indexK = currentIndex;
        unsigned long newValue = (laggedFib[indexJ] + laggedFib[indexK]) % m;
        laggedFib[currentIndex] = newValue;
        currentIndex = (currentIndex + 1) % k;
        return newValue;
    }
};

int main() {
    // Parameters for Lagged Fibonacci Generator
    const int N = 100;
    const int j = 24; // J parameter
    const int k = 55; // K parameter
    const unsigned long m = 4294967296; // Modulus (2^32 for 32-bit numbers)
    unsigned long seed = 123456789; // Initial seed

    // Instantiate the Lagged Fibonacci Generator
    LaggedFibonacciGenerator lfg(j, k, m, seed, k);

    // Vector to store results
    vector<double> results;

    // Generate random numbers
    for (int i = 0; i < N; ++i) {
        unsigned long randomValue = lfg.next();
        results.push_back(static_cast<double>(randomValue) / m);
    }

    // Open Gnuplot process
    FILE* gnuplotPipe = popen("gnuplot -persistent", "w");
    if (!gnuplotPipe) {
        cerr << "Error: Could not open Gnuplot." << endl;
        return 1;
    }

    // Send commands directly to Gnuplot
    fprintf(gnuplotPipe, "set title 'Lagged Fibonacci Generator Results'\n");
    fprintf(gnuplotPipe, "set xlabel 'Iteration'\n");
    fprintf(gnuplotPipe, "set ylabel 'Random Value'\n");
    fprintf(gnuplotPipe, "plot '-' using 1:2 with points title 'Random Numbers'\n");

    // Send data to Gnuplot
    for (int i = 0; i < results.size(); ++i) {
        fprintf(gnuplotPipe, "%d %f\n", i, results[i]);
    }
    fprintf(gnuplotPipe, "e\n");

    // Close Gnuplot pipe
    pclose(gnuplotPipe);

    return 0;
}
