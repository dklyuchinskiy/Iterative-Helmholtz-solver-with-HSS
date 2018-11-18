#include "definitions.h"
#include "templates.h"
#include "TestSuite.h"

/***************************************************
2D Helmholtz solver with HSS

The known solution is generated by function u_ex().
Then we check how correct we has constructed
coefficient matrix as : ||A * u_ex - F_ex|| < eps.

If test is passed, we run the solver for matrix A
and right hand side F.

In the output we compare the relative norm of
the solutuon as:
||u_ex - u_sol|| / ||u_ex|| < eps.

The two version of solver is enabled:
a) with storing the result of factorization in
the array G of doulbe
b) with storing the result of factorization
in the set of structures Gstr, which is defined 
in the definitions.h

The second variant, also, is supported by
storing the initial coefficient matrix A
in sparse CSR format to save memory.

*************************************************/

#ifdef DIM_3D
int main()
{
	TestAll();
	system("pause");
#if 1
	int n1 = 40;		    // number of point across the directions
	int n2 = 40;
	int n3 = 40;
	int n = n1 * n2;		// size of blocks
	int NB = n3;			// number of blocks
	int size = n * NB;		// size of vector x and f: n1 * n2 * n3
	int smallsize = 400;
	double thresh = 1e-6;	// stop level of algorithm by relative error
	int ItRef = 200;		// Maximal number of iterations in refirement
	char bench[255] = "display"; // parameter into solver to show internal results
	int sparse_size = n + 2 * (n - 1) + 2 * (n - n1);
	int non_zeros_in_3diag = n + (n - 1) * 2 + (n - n1) * 2 - (n1 - 1) * 2;

	size_m x, y, z;

	x.n = n1;
	y.n = n2;
	z.n = n3;

	x.l = y.l = z.l = n1 + 1;
	x.h = x.l / (double)(x.n + 1);
	y.h = y.l / (double)(y.n + 1);
	z.h = z.l / (double)(z.n + 1);

	dtype *D;
	dtype *B_mat;

	// Memory allocation for coefficient matrix A
	// the size of matrix A: n^3 * n^3 = n^6
#ifndef ONLINE
	D = alloc_arr(size * n); // it's a matrix with size n^3 * n^2 = size * n
	B_mat = alloc_arr((size - n) * n); 
	int ldd = size;
	int ldb = size - n;
#else
	D = alloc_arr<dtype>(n * n); // it's a matrix with size n^3 * n^2 = size * n
	B_mat = alloc_arr<dtype>(n * n);
	int ldd = n;
	int ldb = n;
#endif

	// Factorization matrix
#ifndef STRUCT_CSR
	double *G = alloc_arr(size * n);
	int ldg = size;
#else
	cmnode **Gstr;
#endif

	// Solution, right hand side and block B
	dtype *B = alloc_arr<dtype>(size - n); // vector of diagonal elementes
	dtype *x_orig = alloc_arr<dtype>(size);
	dtype *x_sol = alloc_arr<dtype>(size);
	dtype *f = alloc_arr<dtype>(size);

#ifdef STRUCT_CSR
	// Memory for CSR matrix
	dcsr *Dcsr;
	int non_zeros_in_block3diag = (n + (n - 1) * 2 + (n - x.n) * 2 - (x.n - 1) * 2) * z.n + 2 * (size - n);
	Dcsr = (dcsr*)malloc(sizeof(dcsr));
	Dcsr->values = alloc_arr<dtype>(non_zeros_in_block3diag);
	Dcsr->ia = alloc_arr<int>(size + 1);
	Dcsr->ja = alloc_arr<int>(non_zeros_in_block3diag);
	Dcsr->ia[size] = non_zeros_in_block3diag + 1;
#endif

	int success = 0;
	int itcount = 0;
	double RelRes = 0;
	double norm = 0;
	int nthr = omp_get_max_threads();
	
	printf("Run in parallel on %d threads\n", nthr);
		
	printf("Grid steps: hx = %lf hy = %lf hz = %lf\n", x.h, y.h, z.h);

#ifndef STRUCT_CSR
	// Generation matrix of coefficients, vector of solution (to compare with obtained) and vector of RHS
	GenMatrixandRHSandSolution(n1, n2, n3, D, ldd, B, x_orig, f);
#else

	// Generation of vector of solution (to compare with obtained), vector of RHS and block B
	GenRHSandSolution(x, y, z, B, x_orig, f);

	// Generation of sparse coefficient matrix
#ifndef ONLINE
	GenSparseMatrix(x, y, z, B_mat, ldb, D, ldd, B_mat, ldb, Dcsr);
#else
	GenSparseMatrixOnline(x, y, z, B_mat, n, D, n, B_mat, n, Dcsr);
	free_arr(D);
#endif
	free_arr(B_mat);

	printf("Non_zeros in block-tridiagonal matrix: %d\n", non_zeros_in_block3diag);

	//	Test_CompareColumnsOfMatrix(n1, n2, n3, D, ldd, B, Dcsr, thresh);
	Test_TransferBlock3Diag_to_CSR(n1, n2, n3, Dcsr, x_orig, f, thresh);
#endif

	printf("Solving %d x %d x %d Laplace equation\n", n1, n2, n3);
	printf("The system has %d diagonal blocks of size %d x %d\n", n3, n1*n2, n1*n2);
	printf("Compressed blocks method\n");
	printf("Parameters: thresh = %g, smallsize = %d \n", thresh, smallsize);

	// Calling the solver
	
#ifndef STRUCT_CSR
	Block3DSPDSolveFast(n1, n2, n3, D, ldd, B, f, thresh, smallsize, ItRef, bench, G, ldg, x_sol, success, RelRes, itcount);
#else

#ifndef ONLINE
	Block3DSPDSolveFastStruct(x, y, z, D, ldd, B, f, Dcsr, thresh, smallsize, ItRef, bench, Gstr, x_sol, success, RelRes, itcount);
#else
	Block3DSPDSolveFastStruct(x, y, z, NULL, ldd, B, f, Dcsr, thresh, smallsize, ItRef, bench, Gstr, x_sol, success, RelRes, itcount);
#endif

#endif
	printf("success = %d, itcount = %d\n", success, itcount);
	printf("-----------------------------------\n");

	printf("Computing error ||x_{exact}-x_{comp}||/||x_{exact}||\n");
	norm = rel_error_complex(n, 1, x_sol, x_orig, size, thresh);

	if (norm < thresh) printf("Norm %12.10e < eps %12.10lf: PASSED\n", norm, thresh);
	else printf("Norm %12.10lf > eps %12.10lf : FAILED\n", norm, thresh);


#ifdef STRUCT_CSR
	Test_DirFactFastDiagStructOnline(x, y, z, Gstr, B, thresh, smallsize);
	//Test_DirSolveFactDiagStructConvergence(x, y, z, Gstr, thresh, smallsize);
	//Test_DirSolveFactDiagStructBlockRanks(x, y, z, Gstr);

	for (int i = z.n - 1; i >= 0; i--)
		FreeNodes(n, Gstr[i], smallsize);

	free(Gstr);
#endif


#ifndef ONLINE
	free_arr(D);
	free_arr(B);
#endif
	free_arr(x_orig);
	free_arr(x_sol);
	free_arr(f);

	system("pause");

	return 0;
#endif
}

#else

int main()
{
	TestAll();

//	Test_PermutLowRankApprox(5, 7, 1e-3, "SVD");
//	Test_PermutLowRankApprox(20, 20, 1e-3, "SVD");
	printf("-----------\n");
//	Test_UnsymmLUfact(5, 1e-4, "SVD", 3);
	printf("-----------\n");
//	Test_UnsymmLUfact(6, 1e-4, "SVD", 3);
	printf("-----------\n");
//	Test_UnsymmLUfact(7, 1e-4, "SVD", 3);
//	Test_UnsymmLUfact(12, 1e-8, "SVD", 5);
//	Test_UnsymmLUfact(12, 1e-8, "SVD", 5);
#if 0
	Test_UnsymmLUfact(5, 1e-8, "SVD", 5);
	Test_UnsymmLUfact(6, 1e-8, "SVD", 5);
	Test_UnsymmLUfact(10, 1e-8, "SVD", 5);
	Test_UnsymmLUfact(12, 1e-8, "SVD", 5);
	Test_UnsymmLUfact(20, 1e-8, "SVD", 5);
	Test_UnsymmLUfact(24, 1e-8, "SVD", 5);
	Test_UnsymmLUfact(40, 1e-8, "SVD", 5);
#endif

//	Test_UnsymmLUfact2(7, 1e-4, "SVD", 3);
//	Test_MyLURecFact(7, 1e-4, "SVD", 3);
//	Test_MyLURecFact(20, 1e-8, "SVD", 5);
//	TestRowInterchange(5, 7, 1e-8);


	// performance test
//	Test_UnsymmLUfact(5000, 1e-8, "SVD", 200);

	system("pause");

#if 1
	int n1 = 400;		    // number of point across the directions
	int n2 = 400;


	int n = n1 + 2 * pml;				// size of blocks
	int NB = n2 + 2 * pml;			// number of blocks
	int size = n * NB;		// size of vector x and f: n1 * n2
	int smallsize = 100;
	double thresh = 1e-6;	// stop level of algorithm by relative error

	int ItRef = 200;		// Maximal number of iterations in refinement
	char bench[255] = "print_time"; // parameter into solver to show internal results
	int sparse_size = n + 2 * (n - 1) + 2 * (n - n1);
	int non_zeros_in_3diag = n + (n - 1) * 2 + (n - n1) * 2 - (n1 - 1) * 2;

	size_m x, y, z;

	x.n = n;
	y.n = NB;

	x.l = y.l = n1 + 1.0;
	x.h = x.l / (double)(n1 + 1);
	y.h = y.l / (double)(n2 + 1);

	dtype *D;
	dtype *B_mat;

	// Memory allocation for coefficient matrix A
	// the size of matrix A: n^3 * n^3 = n^6
#ifndef ONLINE
	D = alloc_arr(size * n); // it's a matrix with size n^3 * n^2 = size * n
	B_mat = alloc_arr((size - n) * n);
	int ldd = size;
	int ldb = size - n;
#else
	D = alloc_arr<dtype>(n * n); // it's a matrix with size n^3 * n^2 = size * n
	B_mat = alloc_arr<dtype>(n * n);
	int ldd = n;
	int ldb = n;
#endif

	// Factorization matrix
#ifndef STRUCT_CSR
	double *G = alloc_arr(size * n);
	int ldg = size;
#else
	cmnode **Gstr;
#endif

	// Solution, right hand side and block B
	dtype *B = alloc_arr<dtype>(size - n); // vector of diagonal elementes
	dtype *x_orig = alloc_arr<dtype>(size);
	dtype *x_sol = alloc_arr<dtype>(size);
	dtype *f = alloc_arr<dtype>(size);

#ifdef STRUCT_CSR
	// Memory for CSR matrix
	dcsr *Dcsr;
	//int non_zeros_in_block3diag = (n + (n - 1) * 2 + (n - x.n) * 2 - (x.n - 1) * 2) * NB + 2 * (size - n); 3D
	int non_zeros_in_block3diag = (n + (n - 1) * 2) * NB + 2 * (size - n);
	Dcsr = (dcsr*)malloc(sizeof(dcsr));
	Dcsr->values = alloc_arr<dtype>(non_zeros_in_block3diag);
	Dcsr->ia = alloc_arr<int>(size + 1);
	Dcsr->ja = alloc_arr<int>(non_zeros_in_block3diag);
	Dcsr->ia[size] = non_zeros_in_block3diag + 1;  // one based indexing
#endif

	int success = 0;
	int itcount = 0;
	double RelRes = 0;
	double norm = 0;
	double timer = 0;
	int nthr = 1;

#ifdef _OPENMP
	nthr = omp_get_max_threads();
#endif

	printf("Run in parallel on %d threads\n", nthr);

	printf("Grid steps: hx = %lf hy = %lf\n", x.h, y.h);

#ifndef STRUCT_CSR
	// Generation matrix of coefficients, vector of solution (to compare with obtained) and vector of RHS
	GenMatrixandRHSandSolution(n1, n2, n3, D, ldd, B, x_orig, f);
#else

	// Generation of sparse coefficient matrix
#ifndef ONLINE
	GenSparseMatrix(x, y, z, B_mat, ldb, D, ldd, B_mat, ldb, Dcsr);
#else
	GenSparseMatrixOnline2D(x, y, B, B_mat, n, D, n, B_mat, n, Dcsr);

	// Generation of vector of solution (to compare with obtained) and vector of RHS
#ifdef HELMHOLTZ
	GenRHSandSolution2D_Syntetic(x, y, Dcsr, x_orig, f);
#else
	GenRHSandSolution2D(x, y, x_orig, f);
#endif

	free_arr(D);
#endif
	free_arr(B_mat);

	printf("Non_zeros in block-tridiagonal matrix: %d\n", non_zeros_in_block3diag);

	//	Test_CompareColumnsOfMatrix(n1, n2, n3, D, ldd, B, Dcsr, thresh);
	Test_TransferBlock3Diag_to_CSR(n1, n2, Dcsr, x_orig, f, thresh);
#endif

	printf("Solving %d x %d Helmholtz equation\n", n1, n2);
	printf("The system has %d diagonal blocks of size %d x %d\n", n, n, n);
	printf("Compressed blocks method\n");
	printf("Parameters: thresh = %g, smallsize = %d \n", thresh, smallsize);

	//system("pause");

	// Calling the solver

#ifndef STRUCT_CSR
	Block3DSPDSolveFast(n1, n2, n3, D, ldd, B, f, thresh, smallsize, ItRef, bench, G, ldg, x_sol, success, RelRes, itcount);
#else

	timer = omp_get_wtime();
#ifndef ONLINE
	Block3DSPDSolveFastStruct(x, y, D, ldd, B, f, Dcsr, thresh, smallsize, ItRef, bench, Gstr, x_sol, success, RelRes, itcount);
#else
	Block3DSPDSolveFastStruct(x, y, NULL, ldd, B, f, Dcsr, thresh, smallsize, ItRef, bench, Gstr, x_sol, success, RelRes, itcount);
#endif
	timer = omp_get_wtime() - timer;
	printf("Time HSS solver: %lf\n", timer);

#endif
	printf("success = %d, itcount = %d\n", success, itcount);
	printf("-----------------------------------\n");

	printf("Computing error ||x_{exact}-x_{HSS}||/||x_{exact}||\n");
	norm = rel_error_complex(size, 1, x_sol, x_orig, size, thresh);

	if (norm < thresh) printf("Norm %12.10e < eps %12.10lf: PASSED\n", norm, thresh);
	else printf("Norm %12.10lf > eps %12.10lf : FAILED\n", norm, thresh);


#ifdef STRUCT_CSR
	Test_DirFactFastDiagStructOnline(x, y, Gstr, B, thresh, smallsize);
	//Test_DirSolveFactDiagStructConvergence(x, y, z, Gstr, thresh, smallsize);
	//Test_DirSolveFactDiagStructBlockRanks(x, y, Gstr);
	//Test_NonZeroElementsInFactors(x, y, Gstr, B, thresh, smallsize);

	for (int i = NB - 1; i >= 0; i--)
		FreeNodes(n, Gstr[i], smallsize);

	free(Gstr);

	//system("pause");
#endif

	// Check Pardiso
	int mtype = 13;
	int *iparm = alloc_arr<int>(64);
	int *perm = alloc_arr<int>(size);
	dtype *x_sol_prd = alloc_arr<dtype>(size);
	size_t *pt = alloc_arr<size_t>(64);

	int maxfct = 1;
	int mnum = 1;
	int phase = 0;
	int rhs = 1;
	int msglvl = 1;
	int error = 0;

	iparm[0] = 0;
	pardisoinit(pt, &mtype, iparm);
	iparm[26] = 1;
	printf("Calling Pardiso...\n");
#if 0
	phase = 11;
	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &size, Dcsr->values, Dcsr->ia, Dcsr->ja, perm, &rhs, iparm, &msglvl, f, x_sol_prd, &error);

	printf("Error Pardiso: %d\n", error);
	phase = 22;
	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &size, Dcsr->values, Dcsr->ia, Dcsr->ja, perm, &rhs, iparm, &msglvl, f, x_sol_prd, &error);

	phase = 33;
	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &size, Dcsr->values, Dcsr->ia, Dcsr->ja, perm, &rhs, iparm, &msglvl, f, x_sol_prd, &error);
#else
	phase = 13;
	timer = omp_get_wtime();
	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &size, Dcsr->values, Dcsr->ia, Dcsr->ja, perm, &rhs, iparm, &msglvl, f, x_sol_prd, &error);
	timer = omp_get_wtime() - timer;
#endif
	printf("Error Pardiso: %d\n", error); fflush(0);
	printf("Time PARDISO: %lf\n", timer); fflush(0);

	printf("Computing error ||x_{exact}-x_{PRD}||/||x_{exact}||\n");

	//compare_vec(size, x_sol_prd, x_orig);
	norm = rel_error_complex(size, 1, x_sol_prd, x_orig, size, thresh);

	if (norm < thresh) printf("Norm %12.10e < eps %12.10lf: PASSED\n", norm, thresh);
	else printf("Norm %12.10lf > eps %12.10lf : FAILED\n", norm, thresh);


#ifndef ONLINE
	free_arr(D);
	free_arr(B);
#endif
	free_arr(x_orig);
	free_arr(x_sol);
	free_arr(f);

	system("pause");

	return 0;
#endif
}
#endif