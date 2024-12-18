#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <utility>
#include <vector>

using namespace std;

// // Derivative
// double derivative(function<double(double)> f, double x, double h) {
//   auto aux = f(x + h);
//   auto aux2 = f(x - h);
//   return (aux - aux2) / (2 * h);
// }

// double derivative_2(function<double(double)> f, double x, double h) {
//   return (derivative(f, x + h, h) - derivative(f, x, h)) / h;
// }

// double derivative_3(function<double(double)> f, double x, double h) {
//   return (derivative_2(f, x + h, h) - derivative_2(f, x, h)) / h;
// }

// Regra trapezoidal
double trapezoidal(function<double(double)> f, double a, double b, int N) {
  double h = (b - a) / N;
  double sum = 0.5 * (f(a) + f(b));
  for (int i = 1; i < N; ++i) {
    sum += f(a + i * h);
  }

  //   double e = 1 / 12 * powf(h, 2) * (derivative(f, a, h) - derivative(f, b,
  //   h));

  return sum * h;
}

// Regra de Simpson
double simpson(function<double(double)> f, double a, double b, int N) {
  if (N % 2 != 0)
    N++; // Garantir que N seja par
  double h = (b - a) / N;
  double sum = f(a) + f(b);
  for (int i = 1; i < N / 2; i++) {
    sum += 4 * f(a + (2 * i - 1) * h);
  }
  for (int i = 2; i < N / 2; i++) {
    sum += 2 * f(a + 2 * i * h);
  }

  //   double e =
  //       1 / 180 * powf(h, 4) * (derivative_3(f, a, h) - derivative_3(f, b,
  //       h));

  return sum * h / 3.0;
}

// Método adaptativo de Simpson
double simpson_adaptativo(function<double(double)> f, double a, double b, int N,
                          double epsilon) {
  if (N % 2 != 0)
    N++; // Garantir que N seja par
  double h = (b - a) / N;
  double sum = f(a) + f(b);
  for (int i = 1; i < N / 2; i++) {
    sum += 4 * f(a + (2 * i - 1) * h);
  }
  for (int i = 2; i < N / 2; i++) {
    sum += 2 * f(a + 2 * i * h);
  }
  auto Ii = sum * h / 3.0;

  do {
    Ii = sum * h / 3.0;
    N *= 2;
    h = (b - a) / N;
    sum = f(a) + f(b);
    for (int i = 1; i < N / 2; i++) {
      sum += 4 * f(a + (2 * i - 1) * h);
    }
    for (int i = 2; i < N / 2; i++) {
      sum += 2 * f(a + 2 * i * h);
    }
    cout << "Ii: " << Ii << " N: " << N << endl;

  } while (abs(Ii - sum * h / 3.0) / 15.0 > epsilon);

  //   double e =
  //       1 / 180 * powf(h, 4) * (derivative_3(f, a, h) - derivative_3(f, b,
  //       h));

  return Ii;
}

// Método de Romberg
double romberg(function<double(double)> f, double a, double b, double epsilon) {
  vector<vector<double>> R(1, vector<double>(1));
  int k = 0;
  R[k][0] = trapezoidal(f, a, b, 1);
  while (true) {
    k++;
    R.emplace_back(k + 1);
    R[k][0] = trapezoidal(f, a, b, 1 << k); // 2^k subintervalos
    cout << "R[" << k << "][0]: " << R[k][0];
    for (int j = 1; j <= k; j++) {
      R[k][j] = (pow(4, j) * R[k][j - 1] - R[k - 1][j - 1]) / (pow(4, j) - 1);
      cout << " R[" << k << "][" << j << "]: " << R[k][j];
      cout << "e = " << fabs((R[k][j - 1] - R[k - 1][j - 1]) / (pow(4, j) - 1));
      if (k > 1 &&
          fabs((R[k][j - 1] - R[k - 1][j - 1]) / (pow(4, j) - 1)) < epsilon) {
        return R[k][k];
      }
    }
    cout << endl;
  }
}

int main() {
  // Definição da função como um lambda
  auto func = [](double x) { return pow(sin(sqrt(200 * fabs(x))), 2); };

  auto integral_by_wolfram =
      0.5028512714128629086973240134980234735337341059602528410035843691643110437548531357209388898179912121;

  double a = 0;
  double b = 1;
  int N = 10;
  auto epsilon = 1e-6;

  std::ofstream dataFile("results.dat");
  dataFile << "# Method Slices Result Error\n";

  double result, error;

  // Wolfram
  result = integral_by_wolfram;
  error = 0;
  dataFile << "'Wolfram 0' " << result << " " << error << "\n";

  cout << "Wolfram " << result << " " << error << "\n";

  // Trapezoidal (10 slices)
  result = trapezoidal(func, a, b, N);
  error = abs(integral_by_wolfram - result);
  dataFile << "\"Trapezoidal 10\" " << 10 << " " << result << " " << error
           << "\n";

  cout << "Trapezoidal 10 " << result << " " << error << "\n";

  // Simpson (10 slices)
  result = simpson(func, a, b, N);
  error = abs(integral_by_wolfram - result);
  dataFile << "\"Simpson 10\" " << 10 << " " << result << " " << error << "\n";

  cout << "Simpson 10 " << result << " " << error << "\n";

  // Trapezoidal (100 slices)
  result = trapezoidal(func, a, b, 100);
  error = abs(integral_by_wolfram - result);
  dataFile << "\"Trapezoidal 100\" " << 100 << " " << result << " " << error
           << "\n";

  cout << "Trapezoidal 100 " << result << " " << error << "\n";

  // Simpson (100 slices)
  result = simpson(func, a, b, 100);
  error = abs(integral_by_wolfram - result);
  dataFile << "\"Simpson 100\" " << 100 << " " << result << " " << error
           << "\n";

  cout << "Simpson 100 " << result << " " << error << "\n";

  // Trapezoidal (1000 slices)
  result = trapezoidal(func, a, b, 1000);
  error = abs(integral_by_wolfram - result);
  dataFile << "\"Trapezoidal 1000\" " << 1000 << " " << result << " " << error
           << "\n";

  cout << "Trapezoidal 1000 " << result << " " << error << "\n";

  // Simpson (1000 slices)
  result = simpson(func, a, b, 1000);
  error = abs(integral_by_wolfram - result);
  dataFile << "\"Simpson 1000\" " << 1000 << " " << result << " " << error
           << "\n";

  cout << "Simpson 1000 " << result << " " << error << "\n";

  dataFile.close();

  // Create GNU Plot script
  std::ofstream gnuplotScript("plot_results.gp");
  gnuplotScript
      << "set terminal qt\n"; // Use qt terminal for interactive window
  gnuplotScript << "set title 'Integration Results Comparison'\n";
  gnuplotScript << "set xlabel 'Slices'\n";
  gnuplotScript << "set ylabel 'Result'\n";
  gnuplotScript << "set style data histograms\n";
  gnuplotScript << "set style fill solid 1.0 border -1\n";
  gnuplotScript << "set boxwidth 0.5\n";
  gnuplotScript << "set key outside\n";
  gnuplotScript << "plot 'results.dat' using 3:xtic(1) title 'Result', \\\n";
  gnuplotScript << "     '' using 4:xtic(1) title 'Error'\n";
  gnuplotScript.close();

  // Execute GNU Plot script
  system("gnuplot -p plot_results.gp");

  // // d) Método adaptativo de Simpson
  cout << "Adaptativo: " << simpson_adaptativo(func, a, b, 2, epsilon) << endl;

  // // e) Integração de Romberg
  cout << "Romberg: " << romberg(func, a, b, epsilon) << endl;

  return 0;
}
