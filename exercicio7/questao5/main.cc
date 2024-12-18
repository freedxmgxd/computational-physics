#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>
#include <functional>
#include <matplot/matplot.h>

using namespace std;
using namespace matplot;

// Função de Bessel Jn usando std::cyl_bessel_j
double bessel_j(int n, double x) {
    return std::cyl_bessel_j(n, x);
}

// Método de Bisseção
double bisection_method(std::function<double(double)> func, double a, double b, double tol, int max_iter) {
    int count = 0;
    while ((b - a) / 2.0 > tol && count < max_iter) {
        ++count;
        double midpoint = (a + b) / 2.0;
        if (func(midpoint) == 0) {
            return midpoint;
        } else if (func(a) * func(midpoint) < 0) {
            b = midpoint;
        } else {
            a = midpoint;
        }
    }
    return (a + b) / 2.0;
}

// Função para calcular e exibir os zeros das funções de Bessel
void calculate_bessel_zeros() {
    vector<double> tolerances = {0.1, 0.001, 0.0001, 0.00001, 0.000001};
    vector<function<double(double)>> bessel_funcs = {
        [](double x) { return bessel_j(0, x); },
        [](double x) { return bessel_j(1, x); },
        [](double x) { return bessel_j(2, x); },
        [](double x) { return bessel_j(3, x); }
    };

    cout << "\nTabela de tempos e zeros para cada função de Bessel e precisão:\n";
    cout << "Função de Bessel | Tolerância | Zero Encontrado | Tempo (s)\n";
    cout << "-------------------------------------------------------------\n";

    for (size_t i = 0; i < bessel_funcs.size(); ++i) {
        cout << "Calculando para J_" << i << "(x):\n";
        for (double tol : tolerances) {
            auto start_time = chrono::high_resolution_clock::now();
            double zero = bisection_method(bessel_funcs[i], 0, 10, tol, 100);
            auto end_time = chrono::high_resolution_clock::now();
            chrono::duration<double> elapsed_time = end_time - start_time;

            cout << "J_" << i << "(x)           | " << tol << "        | " 
                 << fixed << setprecision(5) << zero << "         | " 
                 << elapsed_time.count() << " segundos\n";
        }
    }
}

// Função para plotar as 4 primeiras funções de Bessel
void plot_bessel_functions() {
    vector<double> x(1000);
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] = i * 20.0 / (x.size() - 1);
    }

    vector<double> bessel_0, bessel_1, bessel_2, bessel_3;
    for (double xi : x) {
        bessel_0.push_back(bessel_j(0, xi));
        bessel_1.push_back(bessel_j(1, xi));
        bessel_2.push_back(bessel_j(2, xi));
        bessel_3.push_back(bessel_j(3, xi));
    }

    // Criar gráfico
    plot(x, bessel_0, "-b")->line_width(2).display_name("J_0(x)");
    hold(on);
    plot(x, bessel_1, "-r")->line_width(2).display_name("J_1(x)");
    hold(on);
    plot(x, bessel_2, "-g")->line_width(2).display_name("J_2(x)");
    hold(on);
    plot(x, bessel_3, "-m")->line_width(2).display_name("J_3(x)");

    xlabel("x");
    ylabel("J_n(x)");
    title("As primeiras 4 Funções de Bessel de primeira espécie");
    grid(true);
    show();
}

int main() {
    // Plotar as funções de Bessel
    plot_bessel_functions();
    
    // Calcular e exibir os zeros das funções de Bessel
    calculate_bessel_zeros();

    return 0;
}
