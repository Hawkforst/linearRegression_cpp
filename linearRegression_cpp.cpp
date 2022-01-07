#include <cmath>
#include <cstdio>
#include <vector>
#include <iostream>
#include <algorithm>
using namespace std;



class Matrix
{
public:
    Matrix(size_t rows, size_t cols);
    double& operator()(size_t i, size_t j);
    double operator()(size_t i, size_t j) const;
    Matrix transpose();
    bool printMatrix();
    vector<Matrix> decompLU();
    int determinant(size_t ignoreCol, size_t ignoreRow);
    Matrix matrixMul(Matrix m);
    Matrix invertLU(Matrix L, Matrix U);
    Matrix invertLU();
    size_t get_mRows();
    size_t get_mCols();

private:
    size_t mRows;
    size_t mCols;
    vector<double> mData;
};

Matrix::Matrix(size_t rows, size_t cols)
    : mRows(rows), mCols(cols), mData(rows* cols)
{
}

size_t Matrix::get_mRows() {
    return mRows;
}
size_t Matrix::get_mCols() {
    return mCols;
}

double& Matrix::operator()(size_t i, size_t j) {
    return mData.at(i * mCols + j);
}

double Matrix::operator()(size_t i, size_t j) const {
    return mData[i * mCols + j];
}

Matrix Matrix::transpose() {
    Matrix tM(mCols, mRows);
    for (size_t i = 0; i < mRows; i++) {
        for (size_t j = 0; j < mCols; j++) {
            tM(j, i) = (*this)(i, j);
        }
    }
    return tM;
}

int determinant(size_t ignoreCol, size_t ignoreRow) {
    return 0;
}

vector<Matrix> Matrix::decompLU() {
    vector<Matrix> matrixes;
    Matrix l(mRows, mCols);
    Matrix u(mRows, mCols);
    // Check that matrix is square
    if (mRows != mCols) {
        cout << "Matrix is not square. Can't do LU decompozition" << endl;
        return matrixes;
    }
    size_t i = 0, j = 0, k = 0;
    for (i = 0; i < mRows; i++) {
        // Calculate L (the lower matrix)
        for (j = 0; j < mCols; j++) {
            if (j < i)
                l(j, i) = 0;
            else {
                l(j, i) = (*this)(j, i);
                for (k = 0; k < i; k++) {
                    l(j, i) = l(j, i) - l(j, k) * u(k, i);
                }
            }
        }
        // Calculate U (the upper matrix)
        for (j = 0; j < mRows; j++) {
            if (j < i)
                u(i, j) = 0;
            else if (j == i)
                u(i, j) = 1;
            else {
                u(i, j) = (*this)(i, j) / l(i, i);
                for (k = 0; k < i; k++) {
                    u(i, j) = u(i, j) - ((l(i, k) * u(k, j)) / l(i, i));
                }
            }
        }
    }
    // Pack up matrix L and U
    matrixes.push_back(l);
    matrixes.push_back(u);

    return matrixes;
}

Matrix Matrix::matrixMul(Matrix m) {
    double sum = 0;
    Matrix mulRes(mRows, m.get_mCols());
    for (size_t i = 0; i < mRows; i++) {
        for (size_t j = 0; j < m.get_mCols(); j++) {
            for (size_t k = 0; k < mCols; k++) {
                sum += (*this)(i, k) * m(k, j);
            }
            mulRes(i, j) = sum;
            sum = 0;
        }
    }
    return mulRes;
}

bool Matrix::printMatrix() {
    for (size_t i = 0; i < mRows; i++) {
        for (size_t j = 0; j < mCols; j++) {
            std::cout << mData.at(i * mCols + j) << " ";
        }
        std::cout << std::endl;
    }
    return true;
}

struct LU_matrix
{
    Matrix L;
    Matrix U;
};
// This method calculates the inverse of a Matrix using its LU decompozition
// This is currently linked to an instance of Matrix but matrices L and U can be of any matrix
Matrix Matrix::invertLU(Matrix L, Matrix U) {
    Matrix eigenMatrix(L.get_mRows(), L.get_mCols());
    for (size_t i = 0; i < eigenMatrix.get_mRows(); i++) {
        for (size_t j = 0; j < eigenMatrix.get_mCols(); j++) {
            if (i == j) {
                eigenMatrix(i, j) = 1;
            }
            else {
                eigenMatrix(i, j) = 0;
            }
        }
    }
    // calculate [L][Z] = [I] -> get three vectors
    Matrix Z(L.get_mRows(), L.get_mCols());
    double s;
    for (size_t Z_col = 0; Z_col < L.get_mCols(); Z_col++) {
        // Iterate for each element (row) in Z
        for (size_t row = 0; row < L.get_mRows(); row++) {
            s = 0;
            // On each row, iterate through existing Z elements
            for (size_t j = 0; j < row; j++) {
                s += Z(j, Z_col) * L(row, j);
            }
            Z(row, Z_col) = (eigenMatrix(row, Z_col) - s) / L(row, row);
        }
    }
    // calculate [U][X] = [Z] -> get three vectors X, X = inv(A)
    Matrix X(U.get_mRows(), U.get_mCols());
    for (size_t X_col = 0; X_col < U.get_mCols(); X_col++) {
        // Iterate for each element (row) in Z, go from bottom up
        size_t row = 0;
        for (size_t temp_row = 0; temp_row < U.get_mRows(); temp_row++) {
            row = U.get_mRows() - temp_row - 1;
            //for (size_t row = U.get_mRows() - 1; row >= 0; row--) {
            s = 0;
            // On each row, iterate through existing Z elements
            for (size_t j = row; j < U.get_mRows(); j++) {
                s += X(j, X_col) * U(row, j);
            }
            X(row, X_col) = (Z(row, X_col) - s) / U(row, row);
        }
    }

    return X;
}
// This method calculates the inverse of a Matrix using its LU decompozition
// for the current Matrix instance.
Matrix Matrix::invertLU() {
    vector<Matrix> LU = (*this).decompLU();
    auto L = LU.at(0);
    auto U = LU.at(1);
    Matrix eigenMatrix(L.get_mRows(), L.get_mCols());
    for (size_t i = 0; i < eigenMatrix.get_mRows(); i++) {
        for (size_t j = 0; j < eigenMatrix.get_mCols(); j++) {
            if (i == j) {
                eigenMatrix(i, j) = 1;
            }
            else {
                eigenMatrix(i, j) = 0;
            }
        }
    }

    // calculate [L][Z] = [I] -> get three vectors
    Matrix Z(L.get_mRows(), L.get_mCols());
    double s;
    for (size_t Z_col = 0; Z_col < L.get_mCols(); Z_col++) {
        // Iterate for each element (row) in Z
        for (size_t row = 0; row < L.get_mRows(); row++) {
            s = 0;
            // On each row, iterate through existing Z elements
            for (size_t j = 0; j < row; j++) {
                s += Z(j, Z_col) * L(row, j);
            }
            Z(row, Z_col) = (eigenMatrix(row, Z_col) - s) / L(row, row);
        }
    }
    // calculate [U][X] = [Z] -> get three vectors X, X = inv(A)
    Matrix X(U.get_mRows(), U.get_mCols());
    for (size_t X_col = 0; X_col < U.get_mCols(); X_col++) {
        // Iterate for each element (row) in Z, go from bottom up
        size_t row = 0;
        for (size_t temp_row = 0; temp_row < U.get_mRows(); temp_row++) {
            row = U.get_mRows() - temp_row - 1;
            //for (size_t row = U.get_mRows() - 1; row >= 0; row--) {
            s = 0;
            // On each row, iterate through existing Z elements
            for (size_t j = row; j < U.get_mRows(); j++) {
                s += X(j, X_col) * U(row, j);
            }
            X(row, X_col) = (Z(row, X_col) - s) / U(row, row);
        }
    }
    return X;
}

int main() {
    bool p_intercept = true;
    size_t m, n;
    cin >> n >> m;
    double tmp;

    if (p_intercept) n++;

    Matrix X(m, n);
    Matrix Y(m, 1);
    // load data matrix
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            if (p_intercept && j == 0) {
                X(i, j) = 1.0;
            }
            else {
                cin >> tmp;
                X(i, j) = tmp;
            }
        }
        cin >> tmp;
        Y(i, 0) = tmp;
    }
    // load test matrix
    cin >> m;
    Matrix X_test(m, n);
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            if (p_intercept && j == 0) {
                X_test(i, j) = 1.0;
            }
            else {
                cin >> tmp;
                X_test(i, j) = tmp;
            }
        }
    }
    // calculate regression weights using pseudoinverse
    auto regressionWeights = ((X.transpose()).matrixMul(X).invertLU()).matrixMul(X.transpose()).matrixMul(Y);
    cout << "Testing regression weights" << endl;
    regressionWeights.printMatrix();
    // predict values using calculated weights
    auto regressionPredict = X_test.matrixMul(regressionWeights);

    cout << "Testing fitting" << endl;
    regressionPredict.printMatrix();

}
