#include <vector>
#include <Eigen/Dense>
#include <unordered_map>

using namespace Eigen;
using namespace std;

typedef vector<vector<int>> AdjList;
typedef MatrixXd Matrix;
typedef VectorXd Vector;

void update_parameters(const AdjList &adjacency_list, const vector<int> &Z, const Matrix &A, const Matrix &X, const Vector &Y,
                       Matrix &p, Matrix &mu, Vector &nu, int K) {
    int N = Z.size();

    // Initialize
    Matrix Z_T_Z = Matrix::Zero(K, K);
    Matrix Z_T_A_Z = Matrix::Zero(K, K);
    Matrix Z_T_X_Z = Matrix::Zero(K, K);
    Vector Z_T_Y = Vector::Zero(K);

    // Compute Z^T Z, Z^T A Z, Z^T X Z, and Z^T Y
    for (int i = 0; i < N; i++) {
        int k = Z[i];
        Z_T_Y(k) += Y(i);

        for (int j = 0; j < N; j++) {
            int l = Z[j];
            Z_T_Z(k, l) += 1;

            if (i != j) {
                Z_T_A_Z(k, l) += A(i, j);
                Z_T_X_Z(k, l) += X(i, j);
            }
        }
    }

    // Compute inverses
    Matrix Z_T_Z_inv = Z_T_Z.inverse();
    Matrix Z_T_A_Z_inv = Z_T_A_Z.inverse();

    // Update p, mu, and nu
    p = Z_T_Z_inv * Z_T_A_Z * Z_T_Z_inv;
    mu = Z_T_A_Z_inv * Z_T_X_Z;
    nu = Z_T_Z_inv * Z_T_Y;
}

int main() {
    // Example usage of update_parameters
    int N = 5;  // Number of nodes
    int K = 3;  // Number of clusters

    // Adjacency list representation of the graph
    AdjList adjacency_list(N);
    adjacency_list[0] = {1, 2};
    adjacency_list[1] = {0, 2, 3};
    adjacency_list[2] = {0, 1};
    adjacency_list[3] = {1, 4};
    adjacency_list[4] = {3};

    // Example input matrices and vectors
    vector<int> Z = {0, 1, 2, 0, 1};  // Cluster assignments
    Matrix A = Matrix::Random(N, N);
    Matrix X = Matrix::Random(N, N);
    Vector Y = Vector::Random(N);

    // Initialize output parameters
    Matrix p(K, K);
    Matrix mu(K, K);
    Vector nu(K);

    // Update parameters
    update_parameters(adjacency_list, Z, A, X, Y, p, mu, nu, K);

    // Output the updated parameters
    cout << "Updated p:\n" << p << endl;
    cout << "Updated mu:\n" << mu << endl;
    cout << "Updated nu:\n" << nu << endl;

    return 0;
}
