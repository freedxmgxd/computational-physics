#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>

// Function 1: f(x) = x^2 - cos(4 * pi * x)
double function1(double x) {
    return x * x - cos(4 * M_PI * x);
}

// Function 2: f(x) = cos(x) + cos(sqrt(2)*x) + cos(sqrt(3)*x)
double function2(double x) {
    return cos(x) + cos(sqrt(2) * x) + cos(sqrt(3) * x);
}

// Gaussian random number generator (Box-Muller transform)
double gaussian_random() {
    static bool hasSpare = false;
    static double spare;

    if (hasSpare) {
        hasSpare = false;
        return spare;
    }

    hasSpare = true;
    double u, v, s;
    do {
        u = 2.0 * ((double)rand() / RAND_MAX) - 1.0;
        v = 2.0 * ((double)rand() / RAND_MAX) - 1.0;
        s = u * u + v * v;
    } while (s >= 1.0 || s == 0.0);

    s = sqrt(-2.0 * log(s) / s);
    spare = v * s;
    return u * s;
}

// Metropolis criterion
bool metropolis(double f_current, double f_new, double T) {
    if (f_new < f_current) {
        return true;
    } else {
        double probability = exp((f_current - f_new) / (0.2*T));
        return ((double)rand() / RAND_MAX) < probability;
    }
}

// Simulated Annealing
void simulated_annealing(double (*func)(double), double x_start, double T_start, double T_end, double cooling_rate, double x_min, double x_max) {
    double x = x_start;
    double T = T_start;
    double best_x = x;
    double best_f = func(x);

    // Create a file to store data for plotting
    std::ofstream dataFile("annealing_data.dat");
    if (!dataFile.is_open()) {
        std::cerr << "Error: Could not open file for writing." << std::endl;
        return;
    }

    // Create a file to store the function values for plotting
    std::ofstream funcFile("function_data.dat");
    if (!funcFile.is_open()) {
        std::cerr << "Error: Could not open function file for writing." << std::endl;
        return;
    }

    // Write the function values for plotting
    for (double x_plot = x_min; x_plot <= x_max; x_plot += 0.1) {
        funcFile << x_plot << " " << func(x_plot) << "\n";
    }
    funcFile.close();

    // Main loop
    while (T > T_end) {
        double delta = gaussian_random();
        double x_new = x + delta;

        // Ensure x_new stays within bounds
        if (x_new < x_min) x_new = x_min;
        if (x_new > x_max) x_new = x_max;

        double f_new = func(x_new);

        // Use Metropolis criterion
        if (metropolis(func(x), f_new, T)) {
            x = x_new;
            best_x = x_new;
            best_f = f_new;
        }

        // Update temperature
        T *= cooling_rate;

        // Write data to file
        dataFile << x << " " << best_f << "\n";
    }

    dataFile.close();

    // Output final result
    std::cout << "Final minimum found at x = " << best_x << ", f(x) = " << best_f << std::endl;

    // Call gnuplot to plot the data
    std::ofstream gpFile("plot_commands.gp");
    if (gpFile.is_open()) {
        gpFile << "set title 'Simulated Annealing Trajectory with Function'\n";
        gpFile << "set xlabel 'x'\n";
        gpFile << "set ylabel 'f(x)'\n";
        gpFile << "plot 'function_data.dat' with lines title 'Function', \\\n";
        gpFile << "     'annealing_data.dat' with points pointtype 7 title 'Annealing'\n";
        gpFile.close();
        system("gnuplot -persist plot_commands.gp");
    } else {
        std::cerr << "Error: Could not create gnuplot commands file." << std::endl;
    }

    // close the gnu plot window

}

int main() {
    srand(time(0));

    // Solve for function 1
    std::cout << "Minimizing f(x) = x^2 - cos(4*pi*x)..." << std::endl;
    simulated_annealing(function1, 2.0, 10.0, 1e-3, 0.9, -10.0, 10.0);

    // Solve for function 2
    std::cout << "Minimizing f(x) = cos(x) + cos(sqrt(2)*x) + cos(sqrt(3)*x)..." << std::endl;
    simulated_annealing(function2, 25.0, 20000.0, 1e-3, 0.9, 0.0, 50.0);

    std::cout << "Minimizing f(x) = cos(x) + cos(sqrt(2)*x) + cos(sqrt(3)*x)..." << std::endl;
    simulated_annealing(function2, 25.0, 100000.0, 1e-3, 0.9, 0.0, 50.0);


    return 0;
}
