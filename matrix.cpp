#include <iostream>
#include <random>

using namespace std;

template <typename T>

class matrix {
private:
	size_t rows,cols;
	T** array;
	static constexpr double EPSILON = 1e-5;
public:
	matrix(const size_t row, const size_t col,const T& num = 0) {
		if (row == 0 || col == 0) {
			throw std::invalid_argument("The matrix dimensions must be greater than zero.");
		}
		array = new T*[row];
		for (int i = 0; i < row; ++i) {
			array[i]= new T[col]
				for (int j = 0; j < col; ++j) {
					array[i][j] = num;
				}
		}
	}
	matrix(const size_t row, const size_t col, const T& min, const T& max) {
		if (row == 0 || col == 0) {
			throw std::invalid_argument("The matrix dimensions must be greater than zero.");
		}
		random_device rd;
		mt19937 gen(rd());
		uniform_real_distribution<> dis(min, max);
		array = new T * [row];
		for (int i = 0; i < row; ++i) {
			array[i] = new T[col]
				for (int j = 0; j < col; ++j) {
					array[i][j] = dis(gen);
				}
		}
	}
	matrix(const matrix& m) {
		rows = m.rows;
		cols = m.cols;
		array = new T * [rows];
		for (int i = 0; i < rows; ++i) {
			array[i] = new T[cols]
				for (int j = 0; j < cols; ++j) {
					array[i][j] = m.array[i][j];
				}
		}
	}
	~matrix() {
		for (int i = 0; i < rows; ++i) {
			delete[] array[i];
		}
		delete[] array;
	}
	T& operator () (size_t row, size_t col) {
		if (i >= rows || j >= cols) {
			throw std::out_of_range("Индекс вне границ матрицы");
		}
		return array[row][col];
	}
	matrix operator + (const matrix& m) {
		if (rows != m.rows || cols != m.cols) {
			throw std::invalid_argument("The matrix sizes must match for the addition");
		}
		matrix res(rows, cols);
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				res.array[i][j] = array[i][j] + m.array[i][j];
			}
		}
		return res;
	}
	matrix operator - (const matrix& m) {
		if (rows != m.rows || cols != m.cols) {
			throw std::invalid_argument("The matrix sizes must match for the subtraction");
		}
		matrix res(rows, cols);
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				res.array[i][j] = array[i][j] - m.array[i][j];
			}
		}
		return res;
	}
	matrix operator * (const matrix& m) {
		if (cols != m.rows) {
			throw std::invalid_argument("The number of columns of the first matrix"
										"must match the number of rows of the second matrix");
		}
		matrix res(rows, m.cols);
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < m.cols; ++j) {
				for (int k = 0; k < cols; ++k) {
					res.array[i][j] += array[i][k] * m.array[k][j];
				}
			}
		}
		return res;
	}
	matrix operator * (const T& num) {
		matrix res(rows, cols);
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				res.array[i][j] = array[i][j] * num;
			}
		}
		return res;
	}
	matrix operator * (const T& num, const matrix& m) {
		matrix res(rows, cols);
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				res.array[i][j] = m.array[i][j] * num;
			}
		}
		return res;
	}
	matrix operator / (const T& num) {
		if (abs(num) < EPSILON) {
			throw std::invalid_argument("Division by zero");
		}
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				res.array[i][j] = m.array[i][j] / num;
			}
		}
	}
	matrix operator = (const matrix& m) {
		return matrix(m);
	}
	T trace() {
		if (rows != cols) {
			throw std::invalid_argument("The trace can only be calculated for a square matrix");
		}
		T sum = 0;
		for (int i = 0; i < rows; ++i) {
			sum += array[i][i]
		}
		retirn sum;
	}
	ostream& operator << (std::ostream& os, const matrix& m)
	{
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				os << array[i][j] << ' ';
			}
			os << '\n';
		}
		return os
	}
	bool coplanar(const matrix& m1, const matrix& m2) {
		if (rows != 3 || m1.rows != 3 || m2.rows != 3 || cols != 1 || m1.cols != 1 || m2.cols != 1) {
			throw std::invalid_argument("Matrices of size 3*1 are required to calculate their coplanarity");
		}
		double det = *this(0, 0) * (m1(1, 0) * m2(2, 0) - m1(2, 0) * m2(1, 0)) -
			*this(1, 0) * (m1(2, 0) * m2(0, 0) - m1(0, 0) * m2(2, 0)) +
			*this(2, 0) * (m1(0, 0) * m2(1, 0) - m1(1, 0) * m2(0, 0));
		return abs(det) < EPSILON;
	}
};