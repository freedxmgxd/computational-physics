#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include <matplot/matplot.h>


using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

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

int main() {
  auto a_0 = 0.0;
  auto b_0 = -2.5; // eV

  auto h_n = 1.5;
  auto h_c = 0.0;

  auto k_cc = 1.0;
  auto k_cc2 = 1.1;
  auto k_nc1 = 0.8;
  auto k_nc = 1.0;

  auto a_n = a_0 + h_n * b_0;
  auto a_c = a_0 + h_c * b_0;

  auto b_nc1 = k_nc1 * b_0;
  auto b_nc = k_nc * b_0;
  auto b_cc = k_cc * b_0;
  auto b_cc2 = k_cc2 * b_0;

  Matrix huckel_matrix = {
      {a_n, b_nc1, 0, 0, 0, 0, b_nc1}, {b_nc1, a_c, b_cc2, 0, 0, 0, 0},
      {0, b_cc2, a_c, b_nc1, 0, 0, 0}, {0, 0, b_nc, a_n, b_nc, 0, 0},
      {0, 0, 0, b_nc, a_c, b_cc, 0},   {0, 0, 0, 0, b_cc, a_c, b_cc},
      {b_nc, 0, 0, 0, 0, b_cc, a_c}};

  MatrixXd A(7, 7);
  for (int i = 0; i < 7; i++) {
    for (int j = 0; j < 7; j++) {
      A(i, j) = huckel_matrix[i][j];
    }
  }

  cout << "Matriz de Huckel (Eigen):" << endl;
  cout << A << endl;

  // Resolver autovalores e autovetores
  Eigen::SelfAdjointEigenSolver<MatrixXd> solver_1(A);

  VectorXd omega_squared = solver_1.eigenvalues(); // ω²
  MatrixXd displacements = solver_1.eigenvectors();

  cout << "Autovalores (ω²):" << endl;
  cout << omega_squared << endl;

  cout << "Autovetores (deslocamentos):" << endl;
  cout << displacements << endl;

  // QR decomposition
  vector<vector<double>> Q, R;
  qrDecomposition(huckel_matrix, Q, R);
  for (int i = 0; i < 7000; i++) {
    qrDecomposition(multiplyMatrices(Q, R), Q, R);
  }

  cout << "Matriz de Huckel (QR):" << endl;
  printMatrix(Q);
  cout << endl;
  printMatrix(R);

  Matrix result = multiplyMatrices(Q, R);

  cout << "Matriz de Huckel (QR):" << endl;
  printMatrix(result);

  return 0;
}