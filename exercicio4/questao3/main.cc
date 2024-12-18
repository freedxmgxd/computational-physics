#include <cmath>
#include <iostream>
#include <matplot/matplot.h>
#include <vector>

double calculate_S_up(int N) {
  double S_up = 0.0;
  for (int n = 1; n <= N; n++) {
    S_up += 1.0 / n;
  }
  return S_up;
}

double calculate_S_down(int N) {
  double S_down = 0.0;
  for (int n = N; n >= 1; n--) {
    S_down += 1.0 / n;
  }
  return S_down;
}

int main() {
  std::vector<double> S_up_values;
  std::vector<double> S_down_values;

  std::vector<double> S_diff_values;

  std::vector<u_int64_t> N_values;
  for (u_int64_t i = 1; i <= 10e9; i++) {
    N_values.push_back(i);
  }

  for (int N : N_values) {
    auto aux_S_up = calculate_S_up(N);
    auto aux_S_down = calculate_S_down(N);

    S_up_values.push_back(aux_S_up);
    S_down_values.push_back(aux_S_down);

    S_diff_values.push_back((aux_S_up - aux_S_down)/ (std::fabs(aux_S_up) + std::fabs(aux_S_down)));
  }

   // Save data to a file for gnuplot
    std::ofstream data_file("output/bessel_data_q3.dat");
    for (size_t i = 0; i < N_values.size(); i++) {
        data_file << N_values[i] << " " << S_diff_values[i] << "\n";
    }
    data_file.close();


  // Generate the gnuplot script
  std::ofstream gnuplot_script("output/plot_bessel_q3.gp");
  gnuplot_script << "set logscale xy\n";
  gnuplot_script << "set xlabel 'N'\n";
  gnuplot_script << "set ylabel 'D_S(N)'\n";
  gnuplot_script << "set title 'D_S vs N'\n";
  gnuplot_script << "set grid\n";
  gnuplot_script << "set key top left\n";
  gnuplot_script
      << "plot 'output/bessel_data_q3.dat' using 1:2 with points pointtype 7 title "
         "'D_S_1'\n";
  gnuplot_script.close();

  // Run gnuplot with the generated script
  std::system("gnuplot -p output/plot_bessel_q3.gp");

  // export as csv the S_up_values, S_down_values, S_diff_values and N_values
  std::ofstream file("data.csv");
  file << "N,S_up,S_down,S_diff\n";
  for (std::size_t i = 0; i < N_values.size(); i++) {
    file << N_values[i] << "," << S_up_values[i] << "," << S_down_values[i] << "," << S_diff_values[i] << "\n";
  }
  file.close();

  return 0;
}
