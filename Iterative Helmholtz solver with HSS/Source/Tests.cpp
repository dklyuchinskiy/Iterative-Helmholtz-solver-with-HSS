#include "definitions.h"
#include "templates.h"
#include "TestSuite.h"
#include "TestFramework.h"

/***************************************
Source file contains tests for testing
functionalites, implemented in
functions.cpp and BinaryTrees.cpp.

The interface is declared in TestSuite.h
****************************************/

void TestAll()
{
	TestRunner runner;
	//void(*pt_func)(int&) = NULL;
	//pt_func = &Shell_SymRecCompress;

	printf("***** TEST LIBRARY FUNCTIONS *******\n");
	printf("****Complex precision****\n");
	runner.RunTest(Shell_LowRankApprox, Test_LowRankApproxStruct, "Test_LowRankApprox");
	runner.RunTest(Shell_SymRecCompress, Test_SymRecCompressStruct, "Test_SymRecCompress");
	runner.RunTest(Shell_DiagMult, Test_DiagMultStruct, "Test_DiagMult");
	runner.RunTest(Shell_RecMultL, Test_RecMultLStruct, "Test_RecMultL");
	runner.RunTest(Shell_Add, Test_AddStruct, "Test_Add");
	runner.RunTest(Shell_SymCompUpdate2, Test_SymCompUpdate2Struct, "Test_SymCompUpdate2");
	runner.RunTest(Shell_SymCompRecInv, Test_SymCompRecInvStruct, "Test_SymCompRecInv");
	runner.RunTest(Shell_CopyStruct, Test_CopyStruct,  "Test_CopyStruct");


	printf("********************\n");
	printf("ALL TESTS: %d\nPASSED: %d \nFAILED: %d\n", runner.GetAll(), runner.GetPassed(), runner.GetFailed());

	printf("***** THE END OF TESTING*******\n\n");

}

void Shell_LowRankApprox(ptr_test_low_rank func, const string& test_name, int &numb, int &fail_count)
{
	char method[255] = "SVD";

	for (double eps = 1e-2; eps > 1e-8; eps /= 10)
		for (int m = 3; m <= 10; m++)
			for (int n = 1; n <= 10; n++)
			{
				try
				{
					numb++;
					func(m, n, eps, method);
				}
				catch (runtime_error& e)
				{
					++fail_count;
					cerr << test_name << " fail: " << e.what() << endl;
				}
				catch (...) {
					++fail_count;
					cerr << "Unknown exception caught" << endl;
				}
			}
}

void Shell_SymRecCompress(ptr_test_sym_rec_compress func, const string& test_name, int &numb, int &fail_count)
{
	char method[255] = "SVD";
	int smallsize = 3;

	for (int n = 3; n <= 10; n++)
		for (double eps = 1e-2; eps > 1e-8; eps /= 10)
		{
			try
			{
				numb++;
				func(n, eps, method, smallsize);
			}
			catch (runtime_error& e)
			{
				++fail_count;
				cerr << test_name << " fail: " << e.what() << endl;
			}
			catch (...) {
				++fail_count;
				cerr << "Unknown exception caught" << endl;
			}
		}
}

void Shell_DiagMult(ptr_test_sym_rec_compress func, const string& test_name, int &numb, int &fail_count)
{
	char method[255] = "SVD";
	int smallsize = 3;

	for (int n = 3; n <= 10; n++)
		for (double eps = 1e-2; eps > 1e-8; eps /= 10)
		{
			try
			{
				numb++;
				func(n, eps, method, smallsize);
			}
			catch (runtime_error& e)
			{
				++fail_count;
				cerr << test_name << " fail: " << e.what() << endl;
			}
			catch (...) {
				++fail_count;
				cerr << "Unknown exception caught" << endl;
			}
		}

}

void Shell_RecMultL(ptr_test_mult_diag func, const string& test_name, int &numb, int &fail_count)
{
	char method[255] = "SVD";
	int smallsize = 3;

	for (double eps = 1e-2; eps > 1e-8; eps /= 10)
		for (int n = 3; n <= 10; n++)
			for (int k = 1; k <= 10; k++)
			{
				try
				{
					numb++;
					func(n, k, eps, method, smallsize);
				}
				catch (runtime_error& e)
				{
					++fail_count;
					cerr << test_name << " fail: " << e.what() << endl;
				}
				catch (...) {
					++fail_count;
					cerr << "Unknown exception caught" << endl;
				}
			}

}

void Shell_Add(ptr_test_add func, const string& test_name, int &numb, int &fail_count)
{
	char method[255] = "SVD";
	int smallsize = 3;

	for (double eps = 1e-4; eps > 1e-9; eps /= 10)
		for (int n = 3; n <= 10; n++)
			for (int alpha = -10; alpha < 10; alpha += 2)
				for (int beta = -10; beta < 10; beta += 2)
				{
					if (alpha != 0 && beta != 0)
					{
						try
						{
							numb++;
							func(n, alpha, beta, eps, method, smallsize);
						}
						catch (runtime_error& e)
						{
							++fail_count;
							cerr << test_name << " fail: " << e.what() << endl;
						}
						catch (...) {
							++fail_count;
							cerr << "Unknown exception caught" << endl;
						}
					}
				}
}

void Shell_SymCompUpdate2(ptr_test_update func, const string& test_name, int &numb, int &fail_count)
{
	char method[255] = "SVD";
	int smallsize = 3;

	for (double eps = 1e-3; eps > 1e-9; eps /= 10)
		for (double alpha = -10; alpha < 10; alpha += 2)
			for (int n = 3; n <= 10; n++)
				for (int k = 1; k <= 10; k++)
				{
					try
					{
						numb++;
						func(n, k, alpha, eps, method, smallsize);
					}
					catch (runtime_error& e)
					{
						++fail_count;
						cerr << test_name << " fail: " << e.what() << endl;
					}
					catch (...) {
						++fail_count;
						cerr << "Unknown exception caught" << endl;
					}
				}
}


void Shell_SymCompRecInv(ptr_test_sym_rec_compress func, const string& test_name, int &numb, int &fail_count)
{
	char method[255] = "SVD";
	int smallsize = 3;

	for (double eps = 1e-2; eps > 1e-8; eps /= 10)
		for (int n = 3; n <= 10; n++)
		{
			try
			{
				numb++;
				func(n, eps, method, smallsize);
			}
			catch (runtime_error& e)
			{
				++fail_count;
				cerr << test_name << " fail: " << e.what() << endl;
			}
			catch (...) {
				++fail_count;
				cerr << "Unknown exception caught" << endl;
			}
		}

}

void Shell_CopyStruct(ptr_test_sym_rec_compress func, const string& test_name, int &numb, int &fail_count)
{
	char method[255] = "SVD";
	int smallsize = 3;

	for (double eps = 1e-2; eps > 1e-8; eps /= 10)
		for (int n = 3; n <= 10; n++)
		{
			try
			{
				numb++;
				func(n, eps, method, smallsize);
			}
			catch (runtime_error& e)
			{
				++fail_count;
				cerr << test_name << " fail: " << e.what() << endl;
			}
			catch (...) {
				++fail_count;
				cerr << "Unknown exception caught" << endl;
			}
		}
}


void Test_LowRankApproxStruct(int m, int n, double eps, char *method)
{
	// A - matrix in dense order
	dtype *A = alloc_arr2<dtype>(m * n);
	dtype *A_init = alloc_arr2<dtype>(m * n);
	dtype *A_rec = alloc_arr2<dtype>(m * n);
	char str[255];

	int lda = m;

	double norm = 0;
	dtype alpha = 1.0;
	dtype beta = 0.0;

	for (int j = 0; j < n; j++)
		for (int i = 0; i < m; i++)
		{
			A[i + lda * j] = dtype(1.0 / (n + i + j + 1), 1.0 / (i + j + 1));
			A_init[i + lda * j] = dtype(1.0 / (n + i + j + 1), 1.0 / (i + j + 1));
		}

#if 1
	cmnode *Astr = (cmnode*)malloc(sizeof(cmnode));
	//printf("Test for LowRankApproximationStruct m = %d n = %d ", m, n);
	LowRankApproxStruct(m, n, A, lda, Astr, eps, "SVD"); // memory allocation for Astr inside function
#else
	mnode *Astr;
	//printf("Test for LowRankApproximationStruct2 using return m = %d n = %d ", m, n);
	Astr = LowRankApproxStruct2(m, n, A, lda, eps, "SVD"); // memory allocation for Astr inside function
#endif

	zgemm("no", "no", &m, &n, &Astr->p, &alpha, Astr->U, &m, Astr->VT, &Astr->p, &beta, A_rec, &lda);

	norm = rel_error_complex(m, n, A_rec, A_init, lda, eps);
	sprintf(str, "Struct: n = %d m = %d ", n, m);
	AssertLess(norm, eps, str);


	free_arr<dtype>(A);
	free_arr<dtype>(A_init);
	free_arr<dtype>(A_rec);
	free_arr<dtype>(Astr->U);
	free_arr<dtype>(Astr->VT);
	free_arr<cmnode>(Astr);
}

void Test_SymRecCompressStruct(int n, double eps, char *method, int smallsize)
{
	//printf("*****Test for SymRecCompressStruct  n = %d eps = %e ******* ", n, eps);
	char frob = 'F';
	double norm = 0;

	dtype *H = alloc_arr2<dtype>(n * n); // init
	dtype *H1 = alloc_arr2<dtype>(n * n); // compressed
	dtype *H2 = alloc_arr<dtype>(n * n); // recovered init

	int ldh = n;

	Hilbert(n, n, H, ldh);
	Hilbert(n, n, H1, ldh);

#ifdef DEBUG
	print(n, n, H1, ldh, "H1");
#endif

	cmnode *H1str; // pointer to the tree head
	SymRecCompressStruct(n, H1, ldh, H1str, smallsize, eps, "SVD"); // recursive function means recursive allocation of memory for structure fields
	SymResRestoreStruct(n, H1str, H2, ldh, smallsize);

#ifdef DEBUG
	//print(n, n, H1, ldh, "H1 compressed");
	print(n, n, H2, ldh, "H recovered");
#endif

	// Norm of residual || A - L * U ||
	norm = rel_error_complex(n, n, H2, H, ldh, eps);

#ifdef DEBUG
	print(n, n, H, ldh, "H init");
	print(n, n, H2, ldh, "diff");
#endif

	char str[255];
	sprintf(str, "Struct: n = %d ", n);
	AssertLess(norm, eps, str);

	FreeNodes(n, H1str, smallsize);
	free_arr(H);
	free_arr(H2);
	free_arr(H1);
}

void Test_DiagMultStruct(int n, double eps, char *method, int smallsize)
{
	//printf("*****Test for DiagMultStruct  n = %d ******* ", n);
	dtype *Hd = alloc_arr2<dtype>(n * n); // diagonal Hd = D * H * D
	dtype *H1 = alloc_arr2<dtype>(n * n); // compressed H
	dtype *H2 = alloc_arr2<dtype>(n * n); // recovered H after D * H1 * D
	dtype *d = alloc_arr2<dtype>(n);
	char str[255];

	double norm = 0;
	int ldh = n;

	for (int j = 0; j < n; j++)
	{
		d[j] = j + 1;
	}

	Hilbert(n, n, H1, ldh);
	Hilbert(n, n, Hd, ldh);

	for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
		{
			Hd[i + ldh * j] *= d[j];
			Hd[i + ldh * j] *= d[i];
		}
#ifdef DEBUG
	print(n, n, H1, ldh, "Initial H");
#endif

	cmnode *HCstr;
	// Compress H1 to structured form
	SymRecCompressStruct(n, H1, ldh, HCstr, smallsize, eps, method);

	// Compressed H1 = D * H * D
	DiagMultStruct(n, HCstr, d, smallsize);

	// Recove H1 to uncompressed form
	SymResRestoreStruct(n, HCstr, H2, ldh, smallsize);

#ifdef DEBUG
	print(n, n, Hd, ldh, "Initial Hd = D * H * D");
	print(n, n, H2, ldh, "Recovered H2 = (D * H * D)comp");
#endif

	// Compare Hd and H2
	norm = rel_error_complex(n, n, H2, Hd, ldh, eps);

	sprintf(str, "Struct: n = %d ", n);
	AssertLess(norm, eps, str);

	FreeNodes(n, HCstr, smallsize);
	free_arr(Hd); // diagonal Hd = D * H * D
	free_arr(H1); // compressed H
	free_arr(H2); // recovered H after D * H1 * D
	free_arr(d);
}

/* Тест на сравнение результатов умножения Y = H * X сжимаемой матрицы H на произвольную X.
Сравниваются результаты со сжатием и без */
void Test_RecMultLStruct(int n, int k, double eps, char *method, int smallsize)
{
	//printf("*****Test for RecMultLStruct  n = %d k = %d ******* ", n, k);
	dtype *H = alloc_arr2<dtype>(n * n); // init and compressed
	dtype *X = alloc_arr2<dtype>(n * k);
	dtype *Y = alloc_arr2<dtype>(n * k); // init Y
	dtype *Y1 = alloc_arr2<dtype>(n * k); // after multiplication woth compressed
	char str[255];

	double norm = 0;
	dtype alpha = 1.0;
	dtype beta = 0.0;

	int ldh = n;
	int ldy = n;
	int ldx = n;

	Hilbert(n, n, H, ldh);
	Hilbert(n, k, X, ldx);

	zgemm("No", "No", &n, &k, &n, &alpha, H, &ldh, X, &ldx, &beta, Y, &ldy);

#ifdef DEBUG
	print(n, n, H, ldy, "H init");
	print(n, k, X, ldy, "X init");
	print(n, k, Y, ldy, "Y init");
#endif

	cmnode *Hstr;
	// Compress H
	SymRecCompressStruct(n, H, ldh, Hstr, smallsize, eps, method);

	// RecMult Y1 = comp(H) * X
	RecMultLStruct(n, k, Hstr, X, ldx, Y1, ldy, smallsize);

	norm = rel_error_complex(n, k, Y1, Y, ldy, eps);
	sprintf(str, "Struct: n = %d k = %d ", n, k);
	AssertLess(norm, eps, str);

#ifdef DEBUG
	print(n, n, H, ldy, "H comp");
	print(n, k, Y1, ldy, "Y1 rec");
#endif

	FreeNodes(n, Hstr, smallsize);
	free_arr(H);
	free_arr(X);
	free_arr(Y);
	free_arr(Y1);
}

void Test_AddStruct(int n, dtype alpha, dtype beta, double eps, char *method, int smallsize)
{
	//printf("*****Test for Add n = %d ******* ", n);
	dtype *H1 = alloc_arr2<dtype>(n * n);
	dtype *H2 = alloc_arr2<dtype>(n * n);
	dtype *G = alloc_arr2<dtype>(n * n);
	dtype *H1c = alloc_arr2<dtype>(n * n);
	dtype *H2c = alloc_arr2<dtype>(n * n);
	dtype *GcR = alloc_arr2<dtype>(n * n);
	char str[255];

	int ldh = n;
	int ldg = n;
	double norm = 0;

	Hilbert(n, n, H1, ldh);
	Hilbert(n, n, H1c, ldh);

#pragma omp parallel for simd schedule(simd:static)
	for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
		{
			H2[i + ldh * j] = 1.0 / (i*i + j*j + 1);
			H2c[i + ldh * j] = 1.0 / (i*i + j*j + 1);
		}

#ifdef DEBUG
	print(n, n, H1, ldh, "H1");
	print(n, n, H2, ldh, "H2");
#endif

	cmnode *H1str, *H2str;
	SymRecCompressStruct(n, H1c, ldh, H1str, smallsize, eps, method);
	SymRecCompressStruct(n, H2c, ldh, H2str, smallsize, eps, method);

#ifdef DEBUG
	print(n, n, H1c, ldh, "H1c");
	print(n, n, H2c, ldh, "H2c");
#endif

	cmnode *Gstr;
	Add_dense(n, n, alpha, H1, ldh, beta, H2, ldh, G, ldg);
	AddStruct(n, alpha, H1str, beta, H2str, Gstr, smallsize, eps, method);

#ifdef DEBUG
	print(n, n, G, ldg, "res_dense");
	print(n, n, Gc, ldg, "res_comp");
#endif

	SymResRestoreStruct(n, Gstr, GcR, ldg, smallsize);

#ifdef DEBUG
	print(n, n, GcR, ldg, "res_comp_restore");
#endif
	// |GcR - G| / |G|
	norm = rel_error_complex(n, n, GcR, G, ldg, eps);
	sprintf(str, "Struct: n = %d n = %d alpha = %lf", n, n, alpha, beta);
	AssertLess(norm, eps, str);

	FreeNodes(n, H1str, smallsize);
	FreeNodes(n, H2str, smallsize);
	FreeNodes(n, Gstr, smallsize);
	free_arr(H1);
	free_arr(H2);
	free_arr(G);
	free_arr(H1c);
	free_arr(H2c);
	free_arr(GcR);
}

// B = H - V * Y * VT
void Test_SymCompUpdate2Struct(int n, int k, dtype alpha, double eps, char* method, int smallsize)
{
	//printf("*****Test for SymCompUpdate2Struct  n = %d k = %d ***** ", n, k);
	dtype *B_rec = alloc_arr<dtype>(n * n);
	dtype *Y = alloc_arr2<dtype>(k * k); int ldy = k;
	dtype *V = alloc_arr2<dtype>(n * k); int ldv = n; int ldvtr = k;
	dtype *HC = alloc_arr2<dtype>(n * n); int ldh = n;
	dtype *H = alloc_arr2<dtype>(n * n);
	dtype *C = alloc_arr2<dtype>(n * k); int ldc = n;
	char str[255];

	dtype alpha_one = 1.0;
	dtype beta_zero = 0.0;
	dtype beta_one = 1.0;
	double norm = 0;


	Hilbert(n, n, HC, ldh);
	Hilbert(n, n, H, ldh);

	Clear(k, k, Y, ldy);
	Clear(n, k, V, ldv);

#pragma omp parallel for simd schedule(simd:static)
	for (int i = 0; i < k; i++)
		Y[i + ldy * i] = i + 1;

#pragma omp parallel for simd schedule(simd:static)
	for (int j = 0; j < k; j++)
		for (int i = 0; i < n; i++)
			V[i + ldv * j] = (i + j + 1);

	// C = V * Y
	zsymm("Right", "Up", &n, &k, &alpha_one, Y, &ldy, V, &ldv, &beta_zero, C, &ldc);

	// H = H + alpha * C * VT
	zgemm("No", "Trans", &n, &n, &k, &alpha, C, &ldc, V, &ldv, &beta_one, H, &ldh);

	cmnode *HCstr;
	// Compressed update
	SymRecCompressStruct(n, HC, ldh, HCstr, smallsize, eps, method);

	cmnode *Bstr;
	SymCompUpdate2Struct(n, k, HCstr, alpha, Y, ldy, V, ldv, Bstr, smallsize, eps, method);
	SymResRestoreStruct(n, Bstr, B_rec, ldh, smallsize);

#ifdef DEBUG
	print(n, n, B_rec, ldb, "B_rec");
	print(n, n, H, ldh, "H");
#endif

	// || B_rec - H || / || H ||
	norm = rel_error_complex(n, n, B_rec, H, ldh, eps);
	sprintf(str, "Struct: n = %d k = %d alpha = %lf", n, k, alpha);
	AssertLess(norm, eps, str);

	FreeNodes(n, Bstr, smallsize);
	FreeNodes(n, HCstr, smallsize);
	free_arr(B_rec);
	free_arr(H);
	free_arr(HC);
	free_arr(Y);
	free_arr(C);
	free_arr(V);
}

void Test_SymCompRecInvStruct(int n, double eps, char *method, int smallsize)
{
	//printf("***** Test_SymCompRecInvStruct n = %d eps = %lf ****", n, eps);
	dtype *H = alloc_arr2<dtype>(n * n);
	dtype *Hc = alloc_arr2<dtype>(n * n);
	dtype *Brec = alloc_arr2<dtype>(n * n);
	dtype *Y = alloc_arr2<dtype>(n * n);
	char str[255];

	int ldh = n;
	int ldb = n;
	int ldy = n;

	dtype alpha_mone = -1.0;
	dtype beta_one = 1.0;
	double norm = 0;

	Hilbert(n, n, H, ldh);
	Hilbert(n, n, Hc, ldh);

	// for stability
	for (int i = 0; i < n; i++)
	{
		H[i + ldh * i] += 1.0;
		Hc[i + ldh * i] += 1.0;
	}

	cmnode *HCstr, *BCstr;
	SymRecCompressStruct(n, Hc, ldh, HCstr, smallsize, eps, method);
	SymCompRecInvStruct(n, HCstr, BCstr, smallsize, eps, method);
	SymResRestoreStruct(n, BCstr, Brec, ldb, smallsize);

	Eye(n, Y, ldy);

	// Y = Y - H * Brec
	zgemm("No", "No", &n, &n, &n, &alpha_mone, H, &ldh, Brec, &ldb, &beta_one, Y, &ldy);

	norm = zlange("Frob", &n, &n, Y, &ldy, NULL);
	sprintf(str, "Struct: n = %d", n);
	AssertLess(norm, eps, str);

	//if (norm < eps) printf("Norm %10.8e < eps %10.8lf: PASSED\n", norm, eps);
	//else printf("Norm %10.8lf > eps %10.8e : FAILED\n", norm, eps);

	FreeNodes(n, HCstr, smallsize);
	FreeNodes(n, BCstr, smallsize);
	free_arr(H);
	free_arr(Hc);
	free_arr(Brec);
	free_arr(Y);
}

void Test_CopyStruct(int n, double eps, char *method, int smallsize)
{
	dtype *H = alloc_arr2<dtype>(n * n);
	dtype *H1 = alloc_arr2<dtype>(n * n);
	dtype *H2 = alloc_arr2<dtype>(n * n);
	char str[255];

	double norm = 0;
	int ldh = n;

	//printf("***Test CopyStruct n = %d ", n);

	Hilbert(n, n, H, ldh);
	Hilbert(n, n, H1, ldh);

	cmnode* Hstr, *Hcopy_str;
	SymRecCompressStruct(n, H, ldh, Hstr, smallsize, eps, method);
	CopyStruct(n, Hstr, Hcopy_str, smallsize);
	SymResRestoreStruct(n, Hcopy_str, H2, ldh, smallsize);

	norm = rel_error_complex(n, n, H2, H1, ldh, eps);
	sprintf(str, "Struct: n = %d", n);
	AssertLess(norm, eps, str);

	FreeNodes(n, Hstr, smallsize);
	FreeNodes(n, Hcopy_str, smallsize);
	free_arr(H2);
	free_arr(H1);
	free_arr(H);
}

void Test_DirFactFastDiagStructOnline(size_m x, size_m y, cmnode** Gstr, dtype *B, double eps, int smallsize)
{
	printf("Testing factorization...\n");
	int n = x.n;
	int size = n * y.n;
	int nbr = y.n;
	char bench[255] = "No";
	dtype *DD = alloc_arr<dtype>(n * n); int lddd = n;
	dtype *DR = alloc_arr<dtype>(n * n); int lddr = n;
	dtype *alpX = alloc_arr<dtype>(n + 2);
	dtype *alpY = alloc_arr<dtype>(n + 2);
	double norm = 0;

	double timer = 0;
	timer = omp_get_wtime();

	SetPml(0, x, y, n, alpX, alpY);
	GenerateDiagonal1DBlock(0, x, y, DD, lddd, alpX, alpY);

	cmnode *DCstr;
	SymCompRecInvStruct(n, Gstr[0], DCstr, smallsize, eps, "SVD");
	SymResRestoreStruct(n, DCstr, DR, lddr, smallsize);

	printf("Block %d. ", 0);
	norm = rel_error_complex(n, n, DR, DD, lddd, eps);

	if (norm < eps) printf("Norm %12.10e < eps %12.10lf: PASSED\n", norm, eps);
	else printf("Norm %12.10lf > eps %12.10lf : FAILED\n", norm, eps);

	free_arr(DR);
	free_arr(DD);
	FreeNodes(n, DCstr, smallsize);

	for (int k = 1; k < nbr; k++)
	{
		dtype *DR = alloc_arr<dtype>(n * n); int lddr = n;
		dtype *HR = alloc_arr<dtype>(n * n); int ldhr = n;
		dtype *DD = alloc_arr<dtype>(n * n); int lddd = n;
		cmnode *DCstr, *Hstr;

		printf("Block %d. ", k);

		SymCompRecInvStruct(n, Gstr[k], DCstr, smallsize, eps, "SVD");
		SymResRestoreStruct(n, DCstr, DR, lddr, smallsize);

		CopyStruct(n, Gstr[k - 1], Hstr, smallsize);
		DiagMultStruct(n, Hstr, &B[ind(k - 1, n)], smallsize);
		SymResRestoreStruct(n, Hstr, HR, ldhr, smallsize);

		// Norm of residual
#pragma omp parallel for schedule(static)
		for (int j = 0; j < n; j++)
#pragma omp simd
			for (int i = 0; i < n; i++)
				HR[i + ldhr * j] = HR[i + ldhr * j] + DR[i + lddr * j];

		SetPml(k, x, y, n, alpX, alpY);
		GenerateDiagonal1DBlock(k, x, y, DD, lddd, alpX, alpY);

		norm = rel_error_complex(n, n, HR, DD, lddd, eps);

		if (norm < eps) printf("Norm %12.10e < eps %12.10lf: PASSED\n", norm, eps);
		else printf("Norm %12.10lf > eps %12.10lf : FAILED\n", norm, eps);

		FreeNodes(n, DCstr, smallsize);
		FreeNodes(n, Hstr, smallsize);
		free_arr(DR);
		free_arr(HR);
		free_arr(DD);
	}
	timer = omp_get_wtime() - timer;
	printf("Time: %lf\n", timer);

}

void Test_TransferBlock3Diag_to_CSR(int n1, int n2, dcsr* Dcsr, dtype* x_orig, dtype *f, double eps)
{
	int n = n1;
	int size = n * n2;
	double RelRes = 0;
	dtype *g = alloc_arr2<dtype>(size);
	ResidCSR(n1, n2, Dcsr, x_orig, f, g, RelRes);

	if (RelRes < eps) printf("Norm %10.8e < eps %10.8lf: PASSED\n", RelRes, eps);
	else printf("Norm %10.8lf > eps %10.8e : FAILED\n", RelRes, eps);

	free_arr(g);
}


void Test_DirSolveFactDiagStructBlockRanks(size_m x, size_m y, cmnode** Gstr)
{
	printf("----------Trees information-----------\n");
	int *size = alloc_arr<int>(y.n);
	int *depth = alloc_arr<int>(y.n);

	double time = omp_get_wtime();
#pragma omp parallel
	{
#pragma omp single
		for (int i = 0; i < y.n; i++)
		{
			size[i] = TreeSize(Gstr[i]);
			depth[i] = MaxDepth(Gstr[i]);
		}
	}
	double result = omp_get_wtime() - time;
	printf("Computational time of TreeSize and MaxDepth for all %d trees: %lf\n", y.n, result);

	for (int i = 0; i < y.n; i++)
	{
		printf("For block %2d. Size: %d, MaxDepth: %d, Ranks: ", i, size[i], depth[i]);
		PrintRanksInWidthList(Gstr[i]);
		printf("\n");
	}

	free(size);
	free(depth);

}

void Test_NonZeroElementsInFactors(size_m x, size_m y, cmnode **Gstr, dtype* B, double thresh, int smallsize)
{
	long long non_zeros_exact = 0;
	long long non_zeros_HSS = 0;

	non_zeros_exact = (x.n * x.n) * y.n + 2 * x.n * (y.n - 1);

	for(int k = 0; k < y.n; k++)
		non_zeros_HSS += CountElementsInMatrixTree(x.n, Gstr[k]);

	non_zeros_HSS += 2 * (y.n - 1) * x.n;

	int compr_size = x.n / smallsize;
	int compr_level = log(compr_size) / log(2);
	int loc_size = x.n;
	long long zeros_ideal = 0;

	for (int j = 0; j < compr_level; j++)
	{
		loc_size = ceil(loc_size / 2.0);
		compr_size = ceil(x.n / loc_size);
		zeros_ideal += (loc_size * loc_size * compr_size) * y.n;
	}

	printf("Compression level: %d Compression size: %d\n", compr_level, compr_size);
	printf("non_zeros_exact: %ld\nnon_zeros_HSS: %ld\n", non_zeros_exact, non_zeros_HSS);
	printf("coefficient of compression: %lf (ideal: %lf)\n",  (double)non_zeros_exact / non_zeros_HSS, (double)non_zeros_exact/(non_zeros_exact - zeros_ideal));
	
}




#if 0




void Test_DirSolveFactDiagStructConvergence(size_m x, size_m y, size_m z, mnode** Gstr, double thresh, int smallsize)
{
	printf("------------Test convergence-----------\n");
	int n = x.n * y.n;
	double norm = 0;
	for (int i = 1; i < z.n; i++)
	{
		double *GR1 = alloc_arr(n * n); int ldg1 = n;
		double *GR2 = alloc_arr(n * n); int ldg2 = n;
		double *GRL = alloc_arr(n * n); int ldgl = n;
		printf("For block %2d. \n", i);

		SymResRestoreStruct(n, Gstr[i], GR2, ldg2, smallsize);
		SymResRestoreStruct(n, Gstr[i - 1], GR1, ldg1, smallsize);
		norm = rel_error(n, n, GR2, GR1, ldg1, thresh);
		printf("norm |G[%2d] - G[%2d]|/|G[%2d]| = %12.10lf\n", i - 1, i, i - 1, norm);

		SymResRestoreStruct(n, Gstr[z.n - 1], GRL, ldgl, smallsize);
		norm = rel_error(n, n, GR1, GRL, ldgl, thresh);
		printf("norm |G[%2d] - G[%2d]|/|G[%2d]| = %12.10lf\n\n", i - 1, z.n - 1, z.n - 1, norm);

		free_arr(&GR1);
		free_arr(&GR2);
		free_arr(&GRL);
	}
}

void Test_RankEqual(mnode *Astr, mnode *Bstr)
{
	try
	{
		char str[255] = "Rank(A01) = Rank(B01) ";
		AssertEqual(Astr->p, Bstr->p, str);
	}
	catch (exception &ex)
	{
		cout << ex.what();
	}

	if (Astr->left != NULL || Bstr->left != NULL)
	{
		Test_RankEqual(Astr->left, Bstr->left);
	}

	if (Astr->right != NULL || Bstr->right != NULL)
	{
		Test_RankEqual(Astr->right, Bstr->right);
	}
}

void Test_RankAdd(mnode *Astr, mnode *Bstr, mnode *Cstr)
{
	try
	{
		char str[255] = "Rank(C01) <= Rank(A01) + Rank(B01) ";
		AssertLess(Astr->p + Bstr->p, Cstr->p, str);
	}
	catch (exception &ex)
	{
		cout << ex.what();
	}

}

void Test_Dense_to_CSR(size_m x, size_m y, size_m z, int non_zeros_in_3diag, double *D, int ldd)
{
	int n = x.n * y.n;
	printf("non_zeros: %d\n", non_zeros_in_3diag);
	double *values = alloc_arr(non_zeros_in_3diag);
	int *ia = (int*)malloc(non_zeros_in_3diag * sizeof(int));
	int *ja = (int*)malloc(non_zeros_in_3diag * sizeof(int));
	map<vector<int>, double> CSR;
	CSR = dense_to_CSR(n, n, &D[0], ldd, ia, ja, values);
	print(n, n, &D[0], ldd, "D[0]");
	print_map(CSR);
	free(ia);
	free(ja);
	free(values);

	/*
	print(n, n, &B_mat[0], ldb, "B[0]");
	print(size, n, D, ldd, "D");

	printf("all non_zero elements: %d\n", non_zeros_in_block3diag);
	for (int i = 0; i < size + 1; i ++)
	printf("%d: ia = %d value(ia) = %lf diff = %d\n", i, Dcsr->ia[i], Dcsr->values[Dcsr->ia[i] - 1], Dcsr->ia[i+1]- Dcsr->ia[i]);

	print_vec(non_zeros_in_block3diag, Dcsr->ja, Dcsr->values, "ja and values");
	print_map(CSR);*/

}

void Test_QueueList(int n, double eps, char* method, int smallsize)
{
	printf("----Test for queue implementation with LIST---\n n = %d \n", n);
	// A - matrix in dense order
	double *A = alloc_arr(n * n);
	double *A_init = alloc_arr(n * n);
	int lda = n;

	for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
		{
			A[i + lda * j] = 1.0 / (i + j + 1);
			A_init[i + lda * j] = 1.0 / (i + j + 1);
		}

	mnode *ACstr;
	SymRecCompressStruct(n, A, lda, ACstr, smallsize, eps, method);

	//printf("first level. p = %d\n", ACstr->p);
	//printf("second level. left: %d right: %d\n", Astr->left->p, Astr->right->p);

	printf("Size: %d, MaxDepth: %d, Ranks: ", TreeSize(ACstr), MaxDepth(ACstr));
	//	PrintRanksInWidth(ACstr);
	printf("List:\n");
	PrintRanksInWidthList(ACstr);

	FreeNodes(n, ACstr, smallsize);
	free_arr(&A);
	free_arr(&A_init);
}

void Test_CompareColumnsOfMatrix(int n1, int n2, int n3, double* D, int ldd, double* B, dcsr* Dcsr, double thresh)
{
	int n = n1 * n2;
	int size = n * n3;
	double RelRes = 0;
	double *g = alloc_arr(size);
	double *f1 = alloc_arr(size);
	double *f2 = alloc_arr(size);
	double Res1 = 0;
	double Res2 = 0;

	for (int i = 0; i < size; i++)
	{
		double *x_test = alloc_arr(size);
		x_test[i] = 1;
		Mult_Au(n1, n2, n3, D, ldd, B, x_test, f1);
		mkl_dcsrgemv("No", &size, Dcsr->values, Dcsr->ia, Dcsr->ja, x_test, f2);
		//	print_vec(size, f1, f2, "f1 and f2");
		rel_error(size, 1, f1, f2, size, thresh);
		free(x_test);
	}

	free_arr(&f1);
	free_arr(&f2);
	free_arr(&g);
}
#endif


