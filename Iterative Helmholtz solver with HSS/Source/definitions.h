#pragma once

/*****************************
Preprocessor definitions and
declaration of used structures
*****************************/

// C
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>

// C++
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <complex>

typedef std::complex<double> dtype;
#define MKL_Complex16 dtype

#include "C:\Program Files (x86)\IntelSWTools\compilers_and_libraries\windows\mkl\include\mkl.h"

//#define DEBUG

#define EPS 0.00000001

#define min(a,b) ((a) < (b)) ? (a) : (b)

struct size_m {
	int l;
	int n;
	double h;
};

struct BinaryMatrixTreeNode {

	int p = 0;
	double *U = NULL;
	double *VT = NULL;
	double *A = NULL;
	struct BinaryMatrixTreeNode *left;
	struct BinaryMatrixTreeNode *right;
};

typedef struct BinaryMatrixTreeNode mnode;

struct ComplexBinaryMatrixTreeNode 
{
	int n2 = 0;
	int p = 0;
	int n1 = 0;
	dtype *U = NULL;
	dtype *VT = NULL;
	dtype *A = NULL;
	struct ComplexBinaryMatrixTreeNode *left;
	struct ComplexBinaryMatrixTreeNode *right;
};

typedef struct ComplexBinaryMatrixTreeNode cmnode;

struct ComplexBinaryUnsymmetricMatrixTreeNode 
{
	cmnode *A21;
	cmnode *A12;
	struct ComplexBinaryUnsymmetricMatrixTreeNode *left = NULL;
	struct ComplexBinaryUnsymmetricMatrixTreeNode *right = NULL;
};

typedef struct ComplexBinaryUnsymmetricMatrixTreeNode cumnode;

struct MatrixCSR {

	int *ia = NULL;
	int *ja = NULL;
	dtype *values = NULL;
};

typedef struct MatrixCSR dcsr;

struct list {
	cmnode* node;
	struct list* next;
};

typedef struct list qlist;

struct list2 {
	cumnode* node;
	struct list2* next;
};

typedef struct list2 qlist2;

struct my_queue {
	struct list *first, *last;
};

struct my_queue2 {
	struct list2 *first, *last;
};



#define STRUCT_CSR

#ifdef STRUCT_CSR
#define ONLINE
#endif

//#define FULL_SVD

//#define DIM_3D

//#define COL_UPDATE
//#define COL_ADD

#define max(a, b) ((a > b) ? (a) : (b)) 

#define PI 3.141592653589793238462643

// parameters of Helmholtz equation

//#define HELMHOLTZ

#ifdef HELMHOLTZ
#define omega 4
#define ky 1.8
#define pml 15 //max(15, c0(1, 1) / omega)
#define beta_eq 2.3
#define PML
#else
#define omega 0
#define ky 0
#define pml 0
#define beta_eq 1
#endif

#define ind(x,y)        ((x) * (y))

// Функция выделения памяти под массив

template<typename T>
T* alloc_arr(int n)
{
	T *f = (T*)malloc(n * sizeof(T));

#pragma omp parallel for simd schedule(simd:static)
	for (int i = 0; i < n; i++)
		f[i] = 0.0;

	return f;
}

template<typename T>
T* alloc_arr2(int n)
{
	T *f = (T*)malloc(n * sizeof(T));
#if 0
#pragma omp parallel for simd schedule(simd:static)
	for (int i = 0; i < n; i++)
		f[i] = 0.0;
#endif

	return f;
}

template<typename T>
void free_arr(T* &arr)
{
	free(arr);
}

template <typename MatrixType>
double RelError(double(*LANGE)(const char *, const int*, const int*, const MatrixType*, const int*, double *),
	int m, int n, const MatrixType *Hrec, int ldh1, const MatrixType *Hinit, int ldh2, double eps)
{
	double norm = 0;
	MatrixType *Hdiff = alloc_arr<MatrixType>(m * n);
	int ldh = m;

	// Norm of residual
#pragma omp parallel for schedule(static)
	for (int j = 0; j < n; j++)
#pragma omp simd
		for (int i = 0; i < m; i++)
			Hdiff[i + ldh * j] = Hrec[i + ldh1 * j] - Hinit[i + ldh2 * j];

	norm = LANGE("Frob", &m, &n, Hdiff, &ldh, NULL);
	norm = norm / LANGE("Frob", &m, &n, Hinit, &ldh2, NULL);

	free_arr(Hdiff);

#if 0
	if (norm < eps) printf("Norm %12.10e < eps %12.10lf: PASSED\n", norm, eps);
	else printf("Norm %12.10lf > eps %12.10lf : FAILED\n", norm, eps);
#endif

	return norm;
}

template <typename MatrixType>
double RelErrorPart(double(*LANGE)(const char *, const int*, const int*, const MatrixType*, const int*, double *),
	char part, int m, int n, const MatrixType *Hrec, int ldh1, const MatrixType *Hinit, int ldh2, double eps)
{
	double norm = 0;
	MatrixType *Hdiff = alloc_arr<MatrixType>(m * n);
	int ldh = m;

	if (part == 'L')
	{
		// Norm of residual
#pragma omp parallel for schedule(static)
		for (int j = 0; j < n; j++)
#pragma omp simd
			for (int i = j; i < m; i++)
				Hdiff[i + ldh * j] = Hrec[i + ldh1 * j] - Hinit[i + ldh2 * j];
	}
	else if (part == 'U')
	{
#pragma omp parallel for schedule(static)
		for (int i = 0; i < m; i++)
#pragma omp simd
			for (int j = i; j < n; j++)
				Hdiff[i + ldh * j] = Hrec[i + ldh1 * j] - Hinit[i + ldh2 * j];
	}
	else
	{
		// Norm of residual
#pragma omp parallel for schedule(static)
		for (int j = 0; j < n; j++)
#pragma omp simd
			for (int i = 0; i < m; i++)
				Hdiff[i + ldh * j] = Hrec[i + ldh1 * j] - Hinit[i + ldh2 * j];
	}

//#define PRINT
#ifdef PRINT
	printf("diff\n");
	PrintMat(m, n, Hdiff, ldh);
#endif

	norm = LANGE("Frob", &m, &n, Hdiff, &ldh, NULL);
	norm = norm / LANGE("Frob", &m, &n, Hinit, &ldh2, NULL);

	free_arr(Hdiff);

#if 0
	if (norm < eps) printf("Norm %12.10e < eps %12.10lf: PASSED\n", norm, eps);
	else printf("Norm %12.10lf > eps %12.10lf : FAILED\n", norm, eps);
#endif

	return norm;
}

template <typename MatrixType>
double AbsError(double(*LANGE)(const char *, const int*, const int*, const MatrixType*, const int*, double *),
	int m, int n, const MatrixType *Hrec, int ldh1, const MatrixType *Hinit, int ldh2, double eps)
{
	double norm = 0;
	MatrixType *Hdiff = alloc_arr<MatrixType>(m * n);
	int ldh = m;

	// Norm of residual
#pragma omp parallel for schedule(static)
	for (int j = 0; j < n; j++)
#pragma omp simd
		for (int i = 0; i < m; i++)
			Hdiff[i + ldh * j] = Hrec[i + ldh1 * j] - Hinit[i + ldh2 * j];

	norm = LANGE("Frob", &m, &n, Hdiff, &ldh, NULL);

	free_arr(Hdiff);

#if 0
	if (norm < eps) printf("Norm %12.10e < eps %12.10lf: PASSED\n", norm, eps);
	else printf("Norm %12.10lf > eps %12.10lf : FAILED\n", norm, eps);
#endif

	return norm;
}

template<typename T>
void MultVectorConst(int n, T* v1, T alpha, T* v2)
{
#pragma omp parallel for simd schedule(static)
	for (int i = 0; i < n; i++)
		v2[i] = v1[i] * alpha;
}






