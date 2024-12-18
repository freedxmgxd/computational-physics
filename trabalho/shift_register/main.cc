#include <iostream>
#include <vector>
#include <cstdlib> // For system()

using namespace std;

unsigned long shiftRegisterGenerator(unsigned long& state, int taps[], int numTaps, int bitWidth) {
    unsigned long newBit = 0;
    for (int i = 0; i < numTaps; ++i) {
        newBit ^= (state >> taps[i]) & 1;
    }
    state = ((state << 1) | newBit) & ((1UL << bitWidth) - 1);
    return state;
}

int main() {
    // Parameters
    const int N = 100;
    const int bitWidth = 16; // 16-bit shift register
    unsigned long state = 0xACE1; // Initial state (seed)
    int taps[] = {15, 13, 12, 10}; // Example taps for 16-bit LFSR
    const int numTaps = sizeof(taps) / sizeof(taps[0]);

    // Vector to store results
    vector<double> results;

    // Generate random numbers using shift-register generator
    for (int i = 0; i < N; ++i) {
        unsigned long randomValue = shiftRegisterGenerator(state, taps, numTaps, bitWidth);
        results.push_back(static_cast<double>(randomValue) / (1UL << bitWidth));
    }

    // Open Gnuplot process
    FILE* gnuplotPipe = popen("gnuplot -persistent", "w");
    if (!gnuplotPipe) {
        cerr << "Error: Could not open Gnuplot." << endl;
        return 1;
    }

    // Send commands directly to Gnuplot
    fprintf(gnuplotPipe, "set title 'Shift Register Generator Results'\n");
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
