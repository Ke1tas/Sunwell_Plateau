#include <iostream>
#include <random>

using namespace std;

template <typename T>

class matrix {
private:
	size_t rows,cols;
	T** array;
public:
	matrix(const size_t row, const size_t col,const T& num) {
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
	size_t rows() {
		return rows;
 	}
	size_t cols() {
		return cols;
	}
	matrix operator + (const matrix& m) {
		if (rows != m.rows || cols != m.cols) {
			throw std::invalid_argument("Размеры матриц должны совпадать для сложения");
		}
		matrix
	}
};