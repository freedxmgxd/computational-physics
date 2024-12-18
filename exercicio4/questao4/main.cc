#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

double S_1(double N) {
  double S = 0.0;
  for (double n = 1; n <= 2 * N; n++) {
    S += (powf(-1, n) * (n / (n + 1)));
  }
  return S;
}
double S_2(int N) {
  double S_1 = 0.0;
  double S_2 = 0.0;

  for (double n = 1; n <= N; n++) {
    S_1 += (2 * n - 1) / (2 * n);
  }

  for (double n = 1; n <= N; n++) {
    S_2 += (2 * n) / (2 * n + 1);
  }

  return (-S_1 + S_2);
}

double S_3(int N) {
  double S = 0.0;
  for (double n = 1; n <= N; n++) {
    S += 1.0 / (2 * n * (2 * n + 1));
  }
  return S;
}

int main() {

  std::vector<double> S_1_values;
  std::vector<double> S_2_values;
  std::vector<double> S_3_values;

  std::vector<double> S_1_diff_values;
  std::vector<double> S_2_diff_values;

  std::vector<u_int64_t> N_values;
  for (u_int64_t i = 1; i <= 1000000; i++) {
    N_values.push_back(i);
  }

  for (int N : N_values) {
    auto aux_S_1 = S_1(N);
    auto aux_S_2 = S_2(N);
    auto aux_S_3 = S_3(N);

    S_1_values.push_back(aux_S_1);
    S_2_values.push_back(aux_S_2);
    S_3_values.push_back(aux_S_3);

    S_1_diff_values.push_back(fabs((aux_S_1 - aux_S_3) / aux_S_3));
    S_2_diff_values.push_back(fabs((aux_S_2 - aux_S_3) / aux_S_3));
  }

    // Save data to a file for gnuplot
    std::ofstream data_file("output/bessel_data.dat");
    for (size_t i = 0; i < N_values.size(); i++) {
        data_file << N_values[i] << " " << S_1_diff_values[i] << " " << S_2_diff_values[i] << "\n";
    }
    data_file.close();


  // Generate the gnuplot script
  std::ofstream gnuplot_script("output/plot_bessel.gp");
  gnuplot_script << "set logscale xy\n";
  gnuplot_script << "set xlabel 'N'\n";
  gnuplot_script << "set ylabel 'D_S(N)'\n";
  gnuplot_script << "set title 'D_S_1 and D_S_2 vs N'\n";
  gnuplot_script << "set grid\n";
  gnuplot_script << "set key top left\n";
  gnuplot_script
      << "plot 'output/bessel_data.dat' using 1:2 with points pointtype 7 title "
         "'D_S_1', "
      << "'output/bessel_data.dat' using 1:3 with points pointtype 5 title 'D_S_2'\n";
  gnuplot_script.close();

  // Run gnuplot with the generated script
  std::system("gnuplot -p output/plot_bessel.gp");

  // Export CSV with the results
  std::ofstream file("results.csv");
  file << "N,S_1,S_2,S_3, S_1_diff, S_2_diff\n";
  for (size_t i = 0; i < N_values.size(); i++) {
    file << N_values[i] << "," << S_1_values[i] << "," << S_2_values[i] << ","
         << S_3_values[i] << "," << S_1_diff_values[i] << ","
         << S_2_diff_values[i] << "\n";
  }
  file.close();

  return 0;
}