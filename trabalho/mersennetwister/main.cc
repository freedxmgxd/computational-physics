#include <iostream>
#include <vector>
#include <cstdlib> // For system()

using namespace std;

// Mersenne Twister function implementation
class MersenneTwister {
private:
    static const int n = 624;
    static const int m = 397;
    static const unsigned long a = 0x9908B0DF;
    static const unsigned long upperMask = 0x80000000;
    static const unsigned long lowerMask = 0x7FFFFFFF;
    unsigned long mt[n];
    int index;

    unsigned long tempering(unsigned long y) {
        y ^= (y >> 11);
        y ^= (y << 7) & 0x9D2C5680;
        y ^= (y << 15) & 0xEFC60000;
        y ^= (y >> 18);
        return y;
    }

public:
    MersenneTwister(unsigned long seed) {
        index = n;
        mt[0] = seed;
        for (int i = 1; i < n; ++i) {
            mt[i] = (1812433253 * (mt[i - 1] ^ (mt[i - 1] >> 30)) + i);
            mt[i] &= 0xFFFFFFFF; // Ensure 32-bit word
        }
    }

    void generateNumbers() {
        for (int i = 0; i < n; ++i) {
            unsigned long y = (mt[i] & upperMask) | (mt[(i + 1) % n] & lowerMask);
            mt[i] = mt[(i + m) % n] ^ (y >> 1);
            if (y % 2 != 0) {
                mt[i] ^= a;
            }
        }
    }

    unsigned long next() {
        if (index >= n) {
            generateNumbers();
            index = 0;
        }

        unsigned long y = mt[index++];
        return tempering(y);
    }
};

int main() {
    // Parameters
    const int N = 100; // Number of random numbers to generate

    // Initialize custom Mersenne Twister generator
    MersenneTwister mt(5489); // Seed value

    // Vector to store results
    vector<double> results;

    // Generate random numbers
    for (int i = 0; i < N; ++i) {
        results.push_back(static_cast<double>(mt.next()) / 0xFFFFFFFF);
    }

    // Open Gnuplot process
    FILE* gnuplotPipe = popen("gnuplot -persistent", "w");
    if (!gnuplotPipe) {
        cerr << "Error: Could not open Gnuplot." << endl;
        return 1;
    }

    // Send commands directly to Gnuplot
    fprintf(gnuplotPipe, "set title 'Custom Mersenne Twister Generator Results'\n");
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
