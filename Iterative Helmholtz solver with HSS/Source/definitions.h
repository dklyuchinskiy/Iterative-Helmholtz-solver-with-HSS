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

#include "C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2018\windows\mkl\include\mkl.h"

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
#pragma omp parallel for simd schedule(simd:static)
	for (int i = 0; i < n; i++)
		f[i] = 0.0;

	return f;
}

template<typename T>
void free_arr(T* &arr)
{
	free(arr);
}






