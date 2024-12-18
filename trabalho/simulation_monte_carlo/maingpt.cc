#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <cstdlib>
#include <fstream>

const double PI = 3.141592653589793;

int main() {
    const double T = 10.0;
    const int N = 1000;
    const int steps = 250000;

    // Create a 2D array to store the quantum numbers
    std::vector<std::vector<int>> n(N, std::vector<int>(3, 1));

    // Main loop variables
    std::vector<double> eplot;
    double E = 3.0 * N * PI * PI / 2.0;

    // Random number generators
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    std::uniform_int_distribution<> particle_dis(0, N - 1);
    std::uniform_int_distribution<> quantum_dis(0, 2);

    for (int k = 0; k < steps; ++k) {
        // Choose the particle and the move
        int i = particle_dis(gen);
        int j = quantum_dis(gen);

        double dE;
        int dn;
        if (dis(gen) < 0.5) {
            dn = 1;
            dE = (2 * n[i][j] + 1) * PI * PI / 2.0;
        } else {
            dn = -1;
            dE = (-2 * n[i][j] + 1) * PI * PI / 2.0;
        }

        // Decide whether to accept the move
        if (n[i][j] > 1 || dn == 1) {
            if (dis(gen) < exp(-dE / T)) {
                n[i][j] += dn;
                E += dE;
            }
        }

        eplot.push_back(E);
    }

    // Write data to a file for plotting with gnuplot
    std::ofstream dataFile("eplot.dat");
    for (size_t i = 0; i < eplot.size(); ++i) {
        dataFile << i << " " << eplot[i] << "\n";
    }
    dataFile.close();

    // Open gnuplot directly
    FILE* gnuplot = popen("gnuplot", "w");
    if (gnuplot) {
        fprintf(gnuplot, "set xlabel 'Step'\n");
        fprintf(gnuplot, "set ylabel 'Energy'\n");
        fprintf(gnuplot, "plot '-' with lines title 'Energy'\n");
        for (size_t i = 0; i < eplot.size(); ++i) {
            fprintf(gnuplot, "%lu %f\n", i, eplot[i]);
        }
        fprintf(gnuplot, "e\n");
        fflush(gnuplot);
        std::cout << "Press Enter to exit gnuplot...";
        std::cin.ignore();
        pclose(gnuplot);
    } else {
        std::cerr << "Failed to open gnuplot." << std::endl;
    }

    return 0;
}