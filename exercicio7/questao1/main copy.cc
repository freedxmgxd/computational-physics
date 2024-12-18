#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

// // nomear Matrix como Matrix
using Matrix = vector<vector<double>>;

// // Função para multiplicar duas matrizes
Matrix multiplyMatrices(const Matrix &A, const Matrix &B) {
  int n = A.size();
  Matrix result(n, vector<double>(n, 0.0));
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < n; ++k) {
        result[i][j] += A[i][k] * B[k][j];
      }
    }
  }
  return result;
}

// Function to compute the norm of a vector
double vectorNorm(const vector<double> &vec) {
  double sum = 0.0;
  for (double val : vec) {
    sum += val * val;
  }
  return sqrt(sum);
}

// Function to perform the QR decomposition using Householder reflections
void qrDecomposition(const vector<vector<double>> &A, vector<vector<double>> &Q,
                     vector<vector<double>> &R) {
  int m = A.size();    // Number of rows
  int n = A[0].size(); // Number of columns

  // Initialize Q as an identity matrix and R as a copy of A
  Q = vector<vector<double>>(m, vector<double>(m, 0.0));
  R = A;

  for (int i = 0; i < m; i++) {
    Q[i][i] = 1.0;
  }

  for (int j = 0; j < n; j++) {
    // -- Find H = I - tau * w * w^T to zero out elements below R[j][j]
    vector<double> w(m - j, 0.0); // Vector to store the Householder vector
    double normx = 0.0;

    // Compute the norm of the column starting at R[j][j]
    for (int i = j; i < m; i++) {
      normx += R[i][j] * R[i][j];
    }
    normx = sqrt(normx);

    double s = -copysign(1.0, R[j][j]);
    double u1 = R[j][j] - s * normx;
    w[0] = 1.0;

    for (int i = j + 1; i < m; i++) {
      w[i - j] = R[i][j] / u1;
    }

    double tau = -s * u1 / normx;

    // -- Update R = HR
    for (int col = j; col < n; col++) {
      double dot = 0.0;
      for (int row = 0; row < m - j; row++) {
        dot += w[row] * R[j + row][col];
      }
      dot *= tau;
      for (int row = 0; row < m - j; row++) {
        R[j + row][col] -= dot * w[row];
      }
    }

    // -- Update Q = QH
    for (int col = 0; col < m; col++) {
      double dot = 0.0;
      for (int row = 0; row < m - j; row++) {
        dot += Q[col][j + row] * w[row];
      }
      dot *= tau;
      for (int row = 0; row < m - j; row++) {
        Q[col][j + row] -= dot * w[row];
      }
    }
  }
}

// Function to print a matrix
void printMatrix(const vector<vector<double>> &matrix) {
  for (const auto &row : matrix) {
    for (double val : row) {
      cout << fixed << setprecision(6) << val << "\t";
    }
    cout << endl;
  }
}

Matrix huckelMatrix(double x) {
  // auto a_0 = 0.0;
  auto b_0 = -2.5; // eV

  // auto h_n = 1.5;
  // auto h_c = 0.0;

  auto k_cc = 1.0;
  auto k_cc2 = 1.1;
  auto k_nc1 = 0.8;
  auto k_nc = 1.0;

  // auto a_n = a_0 + h_n * b_0;
  // auto a_c = a_0 + h_c * b_0;

  auto b_nc1 = k_nc1 * b_0;
  auto b_nc = k_nc * b_0;
  auto b_cc = k_cc * b_0;
  auto b_cc2 = k_cc2 * b_0;

  // a_n - E, b_n-c, 0, 0, 0, 0, b_n-c
  // b_n-c, a_c - E, b_c=c, 0, 0, 0, 0
  // 0, b_c=c, a_c - E, b_c-n, 0, 0, 0
  // 0, 0, b_nc, a_n - E, b_nc, 0 , 0
  // 0, 0, 0, b_nc, a_c - E, b_cc, 0
  // 0, 0, 0, 0, b_cc, a_c - E, b_cc
  // b_nc, 0, 0, 0, 0, b_cc, a_c - E

  Matrix huckel_matrix = {
      {0.0, b_nc1, 0, 0, 0, 0, b_nc1}, {b_nc1, 0.0, b_cc2, 0, 0, 0, 0},
      {0, b_cc2, 0.0, b_nc1, 0, 0, 0}, {0, 0, b_nc, 0.0, b_nc, 0, 0},
      {0, 0, 0, b_nc, 0.0, b_cc, 0},   {0, 0, 0, 0, b_cc, 0.0, b_cc},
      {b_nc, 0, 0, 0, 0, b_cc, 0.0}};

  for (int i = 0; i < huckel_matrix.size(); i++) {
    for (int j = 0; j < huckel_matrix[i].size(); j++) {
      huckel_matrix[i][j] = huckel_matrix[i][j] / b_0;
    }
  }

  // \alpha_x = alpha_0 + h_x \beta_0
  // \beta_xx = k_xx \beta_0

  // x_c = (a_c - E) / b_0
  // x_c = (a_0 + h_c \beta_0 - E) / b_0
  // x_c = (a_0 + 0   \beta_0 - E) / b_0 = (a_0 - E) / b_0 => E = a_0 - b_0 x_c
  // x_n = (a_n - E) / b_0
  // x_n = (a_0 + h_n \beta_0 - E) / b_0
  // x_n = (a_0 + 1.5 \beta_0 - E) / b_0 = (a_0 - E) / b_0 + 1.5 = x_c + 1.5

  huckel_matrix[0][0] = x + 1.5;
  huckel_matrix[1][1] = x;
  huckel_matrix[2][2] = x;
  huckel_matrix[3][3] = x + 1.5;
  huckel_matrix[4][4] = x;
  huckel_matrix[5][5] = x;
  huckel_matrix[6][6] = x;

  // cout << "Matriz de Huckel para x = " << x << endl;
  // printMatrix(huckel_matrix);

  return huckel_matrix;
}

double triangular_matrix_determinant(const Matrix &matrix) {
  double det = 1.0;
  for (int i = 0; i < matrix.size(); i++) {
    det *= matrix[i][i];
  }
  return det;
}

double find_root(double x0, double x1, double tol = 1e-6) {
  Matrix A0 = huckelMatrix(x0);
  Matrix R0;
  Matrix Q0;
  Matrix A1 = huckelMatrix(x1);
  Matrix R1;
  Matrix Q1;
  qrDecomposition(A0, Q0, R0);
  qrDecomposition(A1, Q1, R1);
  double f0 = triangular_matrix_determinant(R0);
  double f1 = triangular_matrix_determinant(R1);

  double x = x1;
  double f = f1;

  if (f0 * f1 > 0) {
    cerr << "Erro: Não há uma raiz no intervalo fornecido!" << endl;
    return NAN;
  }

  while ((x1 - x0) / 2.0 > tol) {
    x = (x1 + x0) / 2.0;

    Matrix A = huckelMatrix(x);
    Matrix R;
    Matrix Q;
    qrDecomposition(A, Q, R);

    f = triangular_matrix_determinant(R);

    if (fabs(f) < tol) {
      return x;
    }

    if (f1 * f < 0) {
      x1 = x;
      f1 = f;
    } else {
      x0 = x;
      f0 = f;
    }
  }
  return x;
}

// Função para buscar raízes em um intervalo dado usando divisão adaptativa
void refine_intervals(vector<pair<double, double>> &intervals,
                      double tol = 1e-3) {
  vector<pair<double, double>> refined_intervals;

  for (const auto &interval : intervals) {
    double x0 = interval.first;
    double x1 = interval.second;

    Matrix A0 = huckelMatrix(x0);
    Matrix R0, Q0;
    qrDecomposition(A0, Q0, R0);
    double f0 = triangular_matrix_determinant(R0);

    Matrix A1 = huckelMatrix(x1);
    Matrix R1, Q1;
    qrDecomposition(A1, Q1, R1);
    double f1 = triangular_matrix_determinant(R1);

    // Verifica se há uma mudança de sinal no intervalo
    if (f0 * f1 <= 0) {
      refined_intervals.push_back({x0, x1});
    } else {
      // Dividir em subintervalos
      double mid = (x0 + x1) / 2.0;

      Matrix Am = huckelMatrix(mid);
      Matrix Rm, Qm;
      qrDecomposition(Am, Qm, Rm);
      double fm = triangular_matrix_determinant(Rm);

      // Verifica os subintervalos
      if (f0 * fm <= 0) {
        refined_intervals.push_back({x0, mid});
      }
      if (fm * f1 <= 0) {
        refined_intervals.push_back({mid, x1});
      }
    }
  }

  // Substituir os intervalos refinados
  intervals = refined_intervals;
}

int main() {
  cout << "Resolvendo o determinante ..." << endl;

  // Intervalos iniciais para buscar raízes (autovalores)
  vector<pair<double, double>> intervals = {{-3.0, -2.0}, {-2.0, -1.0},
                                            {-1.0, 0.0},  {0.0, 1.0},
                                            {1.0, 2.0},   {2.0, 3.0}};

  // Refinar os intervalos com base na mudança de sinal
  cout << "Refinando os intervalos..." << endl;
  refine_intervals(intervals, 1e-3);

  cout << "Intervalos refinados: " << endl;
  for (const auto &interval : intervals) {
    cout << "[" << interval.first << ", " << interval.second << "]" << endl;
  }

  auto roots = vector<double>();

  for (const auto &interval : intervals) {
    double root = find_root(interval.first, interval.second);
    if (!isnan(root)) {
      roots.push_back(root);
    }
  }

  // Imprimir os resultados
  cout << "Níveis de energia):" << endl;
  for (double x : roots) {
    // E = a_0 - b_0 x_c
    cout << setw(10) << setprecision(6) << 0 - 2.5 * x << endl;
  }

  return 0;
}