#include <iostream>
#include <random>

using namespace std;

template <typename T>

class matrix {
private:
	size_t _rows,_cols;
	T** _array;
	static constexpr double EPSILON = 1e-5;
public:
	matrix(const size_t row, const size_t col,const T& num = 0) {
		if (row == 0 || col == 0) {
			throw std::invalid_argument("The matrix dimensions must be greater than zero.");
		}
		_rows = row;
		_cols = col;
		_array = new T*[row];
		for (int i = 0; i < row; ++i) {
			_array[i] = new T[col];
			for (int j = 0; j < col; ++j) {
				_array[i][j] = num;
			}
		}
		cout << "matrix created" << '\n';
	}
	matrix(const size_t row, const size_t col, const T& min, const T& max) {
		if (row == 0 || col == 0) {
			throw std::invalid_argument("The matrix dimensions must be greater than zero.");
		}
		_rows = row;
		_cols = col;
		random_device rd;
		mt19937 gen(rd());
		uniform_real_distribution<> dis(min, max);
		_array = new T * [row];
		for (int i = 0; i < row; ++i) {
			_array[i] = new T[col];
			for (int j = 0; j < col; ++j) {
				_array[i][j] = dis(gen);
			}
		}
		cout << "matrix created" << '\n';
	}
	matrix(const matrix& m) {
		_rows = m._rows;
		_cols = m._cols;
		_array = new T * [_rows];
		for (int i = 0; i < _rows; ++i) {
			_array[i] = new T[_cols];
			for (int j = 0; j < _cols; ++j) {
				_array[i][j] = m._array[i][j];
			}
		}
		cout << "matrix created" << '\n';
	}
	~matrix() {
		for (int i = 0; i < _rows; ++i) {
			delete[] _array[i];
		}
		delete[] _array;
		cout << "matrix deleted" << '\n';
	}
	T& operator () (size_t row, size_t col) {
		if (row >= _rows || col >= _cols) {
			throw std::out_of_range("Индекс вне границ матрицы");
		}
		return _array[row][col];
	}
	const T operator () (size_t row, size_t col) const {
		if (row >= _rows || col >= _cols) {
			throw std::out_of_range("Индекс вне границ матрицы");
		}
		return _array[row][col];
	}
	matrix operator + (const matrix& m) {
		if (_rows != m._rows || _cols != m._cols) {
			throw std::invalid_argument("The matrix sizes must match for the addition");
		}
		matrix res(_rows, _cols);
		for (int i = 0; i < _rows; ++i) {
			for (int j = 0; j < _cols; ++j) {
				res(i, j) = _array[i][j] + m(i,j);
			}
		}
		return res;
	}
	matrix operator - (const matrix& m) {
		if (_rows != m._rows || _cols != m._cols) {
			throw std::invalid_argument("The matrix sizes must match for the subtraction");
		}
		matrix res(_rows, _cols);
		for (int i = 0; i < _rows; ++i) {
			for (int j = 0; j < _cols; ++j) {
				res(i,j) = _array[i][j] - m(i, j);
			}
		}
		return res;
	}
	matrix operator * (const matrix& m) {
		if (_cols != m._rows) {
			throw std::invalid_argument("The number of columns of the first matrix"
										"must match the number of rows of the second matrix");
		}
		matrix res(_rows, m._cols);
		for (int i = 0; i < _rows; ++i) {
			for (int j = 0; j < m._cols; ++j) {
				for (int k = 0; k < _cols; ++k) {
					res(i,j) += array[i][k] * m(k,j);
				}
			}
		}
		return res;
	}
	matrix operator * (const T& num) {
		matrix res(_rows, _cols);
		for (int i = 0; i < _rows; ++i) {
			for (int j = 0; j < _cols; ++j) {
				res(i,j) = array[i][j] * num;
			}
		}
		return res;
	}

	matrix operator / (const T& num) {
		if (abs(num) < EPSILON) {
			throw std::invalid_argument("Division by zero");
		}
		matrix res(_rows, _cols);
		for (int i = 0; i < _rows; ++i) {
			for (int j = 0; j < _cols; ++j) {
				res(i, j) = *this(i,j) / num;
			}
		}
		return res;
	}
	matrix operator = (const matrix& m) {
		return matrix(m);
	}
	T trace() {
		if (_rows != _cols) {
			throw std::invalid_argument("The trace can only be calculated for a square matrix");
		}
		T sum = 0;
		for (int i = 0; i < _rows; ++i) {
			sum += array[i][i];
		}
		return sum;
	}

	bool coplanar(const matrix& m1, const matrix& m2) {
		if (_rows != 3 || m1._rows != 3 || m2._rows != 3 || _cols != 1 || m1._cols != 1 || m2._cols != 1) {
			throw std::invalid_argument("Matrices of size 3*1 are required to calculate their coplanarity");
		}
		double det = *this(0, 0) * (m1(1, 0) * m2(2, 0) - m1(2, 0) * m2(1, 0)) -
			*this(1, 0) * (m1(0, 0) * m2(2, 0) - m1(2, 0) * m2(0, 0)) +
			*this(2, 0) * (m1(0, 0) * m2(1, 0) - m1(1, 0) * m2(0, 0));
		return abs(det) < EPSILON;
	}

	size_t rows() const{ return _rows; }
	size_t cols() const { return _cols; }
};

template <typename T>
ostream& operator << (std::ostream& os, const matrix<T>& m)
{
	for (int i = 0; i < m.rows(); ++i) {
		for (int j = 0; j < m.cols(); ++j) {
			os << m(i,j) << ' ';
		}
		os << '\n';
	}
	return os;
}

template <typename T>
matrix<T> operator * (const T& num, const matrix<T>& m) {
	matrix res(m.rows(), m.cols());
	for (int i = 0; i < m.rows(); ++i) {
		for (int j = 0; j < m.cols(); ++j) {
			res(i, j) = m(i,j) * num;
		}
	}
	return res;
}

int main() {
	matrix<int> m(3, 5, 10, 20);
	matrix<int> m1(3, 5, 1);
	cout << m;
	cout << m1;
	m1 = m;
	cout << m1;
}