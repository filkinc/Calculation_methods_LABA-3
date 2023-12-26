#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <cmath>
#include "QuadMatrix.h"
#include "LinSolveAlgs.h"

template<typename T>
using vectorNorm = T(*)(const vector<T>&);

template<typename T>
using matrixNorm = T(*)(const QuadMatrix<T>&);

template <class T>
T norm_2(const vector<T>& x){
	size_t n = x.size();
	T xSum = 0;

	for (int i = 0; i < n; ++i) {
		xSum += x[i] * x[i];
	}
	return sqrt(xSum);
}

template <class T>
T norm_inf(const vector<T>& x) {
	size_t n = x.size();
	T xMax = 0;

	for (int i = 0; i < n; ++i) {
		if (fabs(x[i]) > xMax) {
			xMax = fabs(x[i]);
		}
	}
	return xMax;
}
template <class T>
T norm_inf(const QuadMatrix<T>& A) {
	size_t n = A.order();
	T  aSum, aMax = 0;

	for (int i = 0; i < n; ++i) {
		aSum = 0;
		for (int j = 0; j < n; j++) {
			aSum += fabs(A(i, j));
		}
		if (aSum > aMax) {
			aMax = aSum;
		}
	}
	return aMax;
}

template<class T>
T norm_1(const vector<T>& x) {
	size_t n = x.size();
	T xSum = 0;
	for (int i = 0; i < n; ++i) {
		xSum += fabs(x[i]);
	}
	return xSum;
}
template<class T>
T norm_1(const QuadMatrix<T>& A) {
	size_t n = A.order();
	T aSum, aMax = 0;
	for (int j = 0; j < n; ++j) {
		aSum = 0;
		for (int i = 0; i < n; ++i) {
			aSum += fabs(A(i, j));
		}
		if (aSum > aMax) {
			aMax = aSum;
		}
	}
	return aMax;
}