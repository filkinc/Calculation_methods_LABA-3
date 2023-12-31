#pragma once

#include <iostream>
#include <fstream>
#include "Matrix.h"

template<class T>
class QuadMatrix : public Matrix<T>
{
public:
	QuadMatrix(size_t n) : Matrix<T>(n, n) {}

	QuadMatrix(const QuadMatrix<T>& copy) : Matrix<T>(0, 0) {
		this->elems.resize(copy.elems.size());

		for (size_t i = 0; i < this->elems.size(); ++i) {
			this->elems[i] = copy.elems[i];
		}
	}

	QuadMatrix(vector<vector<T>> data) : Matrix<T>(data[0].size(), data[0].size()) {
		if (data[0].size() != data.size()) {
			cerr << "QuadMatrix(vector<vector<T>> data)!";
		}

		this->elems = vector<vector<T>>(data[0].size());

		for (int i = 0; i < this->elems.size(); ++i) {
			this->elems[i] = data[i];
		}
	}

	//friend vector<T> operator* (const QuadMatrix<T>& A, vector<T> b);
	template<typename L>
	friend QuadMatrix<L> operator* (const QuadMatrix<L>& A, const QuadMatrix<L>& B);
	template<typename L>
	friend QuadMatrix<L> operator* (L coef, const QuadMatrix<L>& A);
	template<typename L>
	friend QuadMatrix<L> operator- (const QuadMatrix<L>& A, const QuadMatrix<L>& B);
	template<typename L>
	friend QuadMatrix<L> operator+ (const QuadMatrix<L>& A, const QuadMatrix<L>& B);
	template<typename L>
	friend vector<L> operator* (const QuadMatrix<L>& A, const vector<L>& v);

	void resize(size_t newN);
	size_t order() const;
	QuadMatrix<T> transposed() const; 
	void print() const;
	QuadMatrix<T> inv() const;
};

template<class T>
vector<T> mul(const QuadMatrix<T>& A, const vector<T>& b) {
	size_t n = A.order();
	vector<T> c(n, 0);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			c[i] += A(i, j) * b[j];
		}
	}
	return c;
}

template<typename L>
vector<L> operator* (const QuadMatrix<L>& A, const vector<L>& v) {
	return mul(A, v);
}

template<typename L>
QuadMatrix<L> operator* (const QuadMatrix<L>& A, const QuadMatrix<L>& B) {
	size_t n = A.order();

	QuadMatrix<L> res(n);

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			res(i, j) = 0;
			for (int k = 0; k < n; ++k) {
				res(i, j) += A(i, k) * B(k, j);
			}
		}
	}

	return res;
}

template<typename L>
QuadMatrix<L> operator* (L coef, const QuadMatrix<L>& A) {
	size_t n = A.order();
	QuadMatrix<L> res(n);

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			res(i, j) = coef * A(i, j);
		}
	}

	return res;
}

template<typename L>
QuadMatrix<L> operator- (const QuadMatrix<L>& A, const QuadMatrix<L>& B) {
	size_t n = A.order();

	QuadMatrix<L> res(n);

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			res(i, j) = A(i, j) - B(i, j);
		}
	}

	return res;
}

template<typename L>
QuadMatrix<L> operator+ (const QuadMatrix<L>& A, const QuadMatrix<L>& B) {
	size_t n = A.order();

	QuadMatrix<L> res(n);

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			res(i, j) = A(i, j) + B(i, j);
		}
	}

	return res;
}

template<class T>
size_t QuadMatrix<T>::order() const {
	return this->elems.size();
}

template<class T>
void QuadMatrix<T>::resize(size_t newN) {
	this->elems.resize(newN);
	for (size_t i = 0; i < newN; ++i) {
		this->elems[i].resize(newN);
	}
}

template<class T>
QuadMatrix<T> QuadMatrix<T>::transposed() const {
	QuadMatrix A(*this);
	for (int i = 0; i < A.order(); ++i) {
		for (int j = 0; j < A.order(); ++j) {
			A(i, j) = A(j, i);
		}
	}
	return A;
}

template<class T>
void QuadMatrix<T>::print() const {
	QuadMatrix A(*this);
	size_t n = A.order();
	for (int i = 0; i < n; ++i) {
		std::cout << std::endl;
		for (int j = 0; j < n; ++j) {
			std::cout << A(i, j) << ' ';
		}
	}
	return;
}

template<class T>
QuadMatrix<T> QuadMatrix<T>::inv() const{
	size_t size = (*this).order();
	QuadMatrix A(*this), E(size);
	for (size_t i = 0; i < size; ++i) {
		E(i, i) = 1;
	}
	for (size_t i = 0; i < size; ++i) {
		int m = i;
		for (int k = i; k < size; ++k) {
			if (fabs(A(k, i)) > fabs(A(m, i))) {
				m = k;
			};
		}
		for (int k = 0; k < size; ++k) {
			swap(A(m, k), A(i, k));
			swap(E(m, k), E(i, k));
		}

		if (A(i, i) == 0) {
			cout << "������� ������� �� ������� ���������!";
			system("pause");
			exit(1);
		}

		T a_ii = A(i, i);
		for (size_t j = 0; j < size; ++j) {
			E(i, j) /= a_ii;
			A(i, j) /= a_ii;
		}
		for (int k = i + 1; k < size; ++k) {
			T a_ki = A(k, i);
			for (size_t j = 0; j < size; ++j) {
				A(k, j) += A(i, j) * (-a_ki);
				E(k, j) += E(i, j) * (-a_ki);
			}
		}
	}
	for (int i = size - 1; i >= 0; i--) {
		for (int k = i - 1; k >= 0; k--) {
			T a_ki = A(k, i);
			for (int j = size - 1; j >= 0; j--) {
				A(k, j) += A(i, j) * (-a_ki);
				E(k, j) += E(i, j) * (-a_ki);
			}
		}
	}
	return E;
}