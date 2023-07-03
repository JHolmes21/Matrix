#include <iostream>
#include <iomanip>
#include <utility>

class Matrix {
public:
    // Parameterized constructor
    Matrix(int rows, int cols) : rows_(rows), cols_(cols), data_(new double* [rows]) {
        for (int i = 0; i < rows; ++i) {
            data_[i] = new double[cols];
            for (int j = 0; j < cols; ++j) {
                data_[i][j] = 0.0;
            }
        }
    }

    // Destructor
    ~Matrix() {
        for (int i = 0; i < rows_; ++i) {
            delete[] data_[i];
        }
        delete[] data_;
    }

    // Copy constructor
    Matrix(const Matrix& other) : Matrix(other.rows_, other.cols_) {
        for (int i = 0; i < rows_; ++i) {
            for (int j = 0; j < cols_; ++j) {
                data_[i][j] = other.data_[i][j];
            }
        }
    }

    // Copy assignment operator
    Matrix& operator=(const Matrix& other) {
        if (this != &other) {
            // Deallocate existing data
            for (int i = 0; i < rows_; ++i) {
                delete[] data_[i];
            }
            delete[] data_;

            // Allocate new data
            rows_ = other.rows_;
            cols_ = other.cols_;
            data_ = new double* [rows_];
            for (int i = 0; i < rows_; ++i) {
                data_[i] = new double[cols_];
                for (int j = 0; j < cols_; ++j) {
                    data_[i][j] = other.data_[i][j];
                }
            }
        }
        return *this;
    }
    void set(int i, int j, double val) {
        data_[i][j] = val;
    }
    // Move constructor
    Matrix(Matrix&& other) noexcept : rows_(other.rows_), cols_(other.cols_), data_(other.data_) {
        other.rows_ = 0;
        other.cols_ = 0;
        other.data_ = nullptr;
    }

    // Move assignment operator
    Matrix& operator=(Matrix&& other) noexcept {
        if (this != &other) {
            // Deallocate existing data
            for (int i = 0; i < rows_; ++i) {
                delete[] data_[i];
            }
            delete[] data_;

            // Move data from other to this
            rows_ = other.rows_;
            cols_ = other.cols_;
            data_ = other.data_;
            other.rows_ = 0;
            other.cols_ = 0;
            other.data_ = nullptr;
        }
        return *this;
    }

    // Friend function to overload <<
    friend std::ostream& operator<<(std::ostream& os, const Matrix& matrix) {
        for (int i = 0; i < matrix.rows_; ++i) {
            for (int j = 0; j < matrix.cols_; ++j) {
                os << std::setw(8) << std::setprecision(4) << matrix.data_[i][j];
            }
            os << std::endl;
        }
        return os;
    }

    // Friend function to overload >>
    friend std::istream& operator>>(std::istream& is, Matrix& matrix) {
        for (int i = 0; i < matrix.rows_; ++i) {
            for (int j = 0; j <
                matrix.cols_; ++j) {
                is >> matrix.data_[i][j];
            }
        }
        return is;
    }

    // Overload + for matrix addition
    Matrix operator+(const Matrix& rhs) const {
        if (rows_ != rhs.rows_ || cols_ != rhs.cols_) {
            throw std::invalid_argument("Matrix dimensions must match for addition.");
        }
        Matrix result(rows_, cols_);
        for (int i = 0; i < rows_; ++i) {
            for (int j = 0; j < cols_; ++j) {
                result.data_[i][j] = data_[i][j] + rhs.data_[i][j];
            }
        }

        return result;

    }

    // Overload - for matrix subtraction
    Matrix operator-(const Matrix& rhs) const {
        if (rows_ != rhs.rows_ || cols_ != rhs.cols_) {
            throw std::invalid_argument("Matrix dimensions must match for subtraction.");
        }
        Matrix result(rows_, cols_);
        for (int i = 0; i < rows_; ++i) {
            for (int j = 0; j < cols_; ++j) {
                result.data_[i][j] = data_[i][j] - rhs.data_[i][j];
            }
        }

        return result;

    }

    // Overload * for matrix multiplication
    Matrix operator*(const Matrix& rhs) const {
        if (cols_ != rhs.rows_) {
            throw std::invalid_argument("Matrix dimensions are not compatible for multiplication.");
        }
        Matrix result(rows_, rhs.cols_);
        for (int i = 0; i < rows_; ++i) {
            for (int j = 0; j < rhs.cols_; ++j) {
                double sum = 0.0;
                for (int k = 0; k < cols_; ++k) {
                    sum += data_[i][k] * rhs.data_[k][j];
                }
                result.data_[i][j] = sum;
            }
        }

        return result;

    }
    //Function to delete a row/column
    Matrix delete_row_col(int i, int j) const {
        Matrix result(rows_ - 1, cols_ - 1);
        int r = 0;
        for (int k = 0; k < rows_; ++k) {
            if (k != i) {
                int c = 0;
                for (int l = 0; l < cols_; ++l) {
                    if (l != j) {
                        result.data_[r][c] = data_[k][l];
                        ++c;
                    }
                }
                ++r;
            }
        }
        return result;
    }
    //Function to calculate the determinant of a matrix
    double determinant() const {
        if (rows_ != cols_) {
            throw std::invalid_argument("Matrix must be square to calculate determinant.");
        }

        if (rows_ == 1) {
            return data_[0][0];
        }

        if (rows_ == 2) {
            return data_[0][0] * data_[1][1] - data_[0][1] * data_[1][0];
        }

        double det = 0.0;
        int sign = 1;

        for (int j = 0; j < cols_; ++j) {
            Matrix submatrix = delete_row_col(0, j);
            det += sign * data_[0][j] * submatrix.determinant();
            sign *= -1;
        }

        return det;
    }


private:
    int rows_;
    int cols_;
    double** data_;
};

int main() {
    // Create the matrices
    Matrix A(3, 3);
    Matrix B(3, 3);
    Matrix C(2, 3);

    // Set the values of matrix A
    A.set(0, 0, 1.0);
    A.set(0, 1, 2.0);
    A.set(0, 2, 3.0);
    A.set(1, 0, 9.0);
    A.set(1, 1, 8.0);
    A.set(1, 2, 7.0);
    A.set(2, 0, 4.0);
    A.set(2, 1, 2.0);
    A.set(2, 2, 6.0);

    // Set the values of matrix B
    B.set(0, 0, 5.0);
    B.set(0, 1, 5.0);
    B.set(0, 2, 4.0);
    B.set(1, 0, 1.0);
    B.set(1, 1, 2.0);
    B.set(1, 2, 3.0);
    B.set(2, 0, 6.0);
    B.set(2, 1, 9.0);
    B.set(2, 2, 8.0);

    // Set the values of matrix C
    C.set(0, 0, 3.0);
    C.set(0, 1, 4.0);
    C.set(0, 2, 1.0);
    C.set(1, 0, 2.0);
    C.set(1, 1, 5.0);
    C.set(1, 2, 6.0);

    // Print the matrices
    std::cout << "Matrix A:\n" << A << std::endl;
    std::cout << "Matrix B:\n" << B << std::endl;
    std::cout << "Matrix C:\n" << C << std::endl;

    // Add the matrices
    Matrix D = A + B;
    std::cout << "Matrix A + B:\n" << D << std::endl;

    // Subtract the matrices
    Matrix E = A - B;
    std::cout << "Matrix A - B:\n" << E << std::endl;

    // Multiply the matrices
    Matrix F = A * B;
    std::cout << "Matrix A * B:\n" << F << std::endl;
    Matrix G = C * B;
    std::cout << "Matrix C * B:\n" << G << std::endl;
    try {
        Matrix H = B * C;
        std::cout << "Matrix C * B:\n" << H << std::endl;
    }
    catch (const std::invalid_argument& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
    // Calculate the determinants
    double det_A = A.determinant();
    double det_B = B.determinant();
    std::cout << "Determinant of A: " << det_A << std::endl;
    std::cout << "Determinant of B: " << det_B << std::endl;
    // Use the copy constructor on A and output the original and the copy
    Matrix A_copy(A);
    std::cout << "Original Matrix A:\n" << A << std::endl;
    std::cout << "Copy of Matrix A:\n" << A_copy << std::endl;
    // Change an element of the original matrix A
    A.set(0, 0, 42.0);
    std::cout << "Modified Matrix A:\n" << A << std::endl;
    std::cout << "Copy of Matrix A (unchanged):\n" << A_copy << std::endl;
    // Demonstrate the move constructor using matrix A
    Matrix A_moved(std::move(A));
    std::cout << "New matrix after move constructor with A:\n" << A_moved << std::endl;
    return 0;
}