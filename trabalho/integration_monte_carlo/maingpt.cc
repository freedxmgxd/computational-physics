#include <functional>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <unistd.h> // For sleep function

using namespace std;

double integrationMonteCarlo(double (*f)(double), double a, double b, int N, double h, double t, vector<pair<double, double>>& points, int& count) {
    // count = 0;
    double A = abs((a - b) * (h - t));

    for (int i = 0; i < N; ++i) {
        double x = a + (b - a) * rand() / RAND_MAX;
        double y = h + (t - h) * rand() / RAND_MAX;
        points.emplace_back(x, y);
        if (y <= f(x)) {
            count++;
        }
    }

    return A * count / N;
}

int main() {
    // Define the function (quarter-circle)
    auto quarterCircle = [](double x) { return sqrt(1 - x * x); };

    // Parameters for Monte Carlo integration
    double a = 0;
    double b = 1;
    double h = 0;
    double t = 1;
    int totalIterations = 1000;
    int step = 10; // Update every 10 points

    // Vector to store the points and a counter for points inside the function
    vector<pair<double, double>> points;
    int insideCount = 0;

    // Open Gnuplot process
    FILE* gnuplotPipe = popen("gnuplot -persistent", "w");
    if (!gnuplotPipe) {
        cerr << "Error: Could not open Gnuplot." << endl;
        return 1;
    }

    // Plot settings
    fprintf(gnuplotPipe, "set title 'Monte Carlo Integration - Animated Quarter Circle'\n");
    fprintf(gnuplotPipe, "set xlabel 'x'\n");
    fprintf(gnuplotPipe, "set ylabel 'y'\n");
    fprintf(gnuplotPipe, "set xrange [%f:%f]\n", a, b);
    fprintf(gnuplotPipe, "set yrange [%f:%f]\n", h, t);
    fprintf(gnuplotPipe, "set key off\n");

    // Perform Monte Carlo integration iteratively
    for (int i = 1; i <= totalIterations; ++i) {
        double x = a + (b - a) * rand() / RAND_MAX;
        double y = h + (t - h) * rand() / RAND_MAX;
        points.emplace_back(x, y);
        if (y <= quarterCircle(x)) {
            insideCount++;
        }

        if (i % step == 0 || i == totalIterations) {
            // Clear plot
            fprintf(gnuplotPipe, "clear\n");

            // Plot the function
            fprintf(gnuplotPipe, "plot '-' with lines lw 2 title 'f(x) = sqrt(1 - x^2)', '-' with points pt 7 lc rgb 'blue' title 'Points Inside', '-' with points pt 7 lc rgb 'red' title 'Points Outside'\n");

            // Send function data to Gnuplot
            for (double x = a; x <= b; x += 0.01) {
                fprintf(gnuplotPipe, "%f %f\n", x, quarterCircle(x));
            }
            fprintf(gnuplotPipe, "e\n");

            // Send inside points to Gnuplot
            for (const auto& point : points) {
                if (point.second <= quarterCircle(point.first)) {
                    fprintf(gnuplotPipe, "%f %f\n", point.first, point.second);
                }
            }
            fprintf(gnuplotPipe, "e\n");

            // Send outside points to Gnuplot
            for (const auto& point : points) {
                if (point.second > quarterCircle(point.first)) {
                    fprintf(gnuplotPipe, "%f %f\n", point.first, point.second);
                }
            }
            fprintf(gnuplotPipe, "e\n");

            // Update label with current estimate
            double currentEstimate = abs((a - b) * (h - t)) * insideCount / i;
            fprintf(gnuplotPipe, "set label 1 sprintf('Iteration: %d, Estimate: %.5f', %d, %f) at graph 0.02, graph 0.95\n", i, currentEstimate, i, currentEstimate);

            fflush(gnuplotPipe);
            usleep(500000); // Pause for 50ms to create animation effect
        }
    }

    cout << "Final estimate: " << abs((a - b) * (h - t)) * insideCount / totalIterations << endl;
    cout << "Actual value: " << M_PI / 4 << endl;
    cout << "Error: " << abs(M_PI / 4 - abs((a - b) * (h - t)) * insideCount / totalIterations) << endl;

    // Close Gnuplot pipe
    pclose(gnuplotPipe);

    return 0;
}
