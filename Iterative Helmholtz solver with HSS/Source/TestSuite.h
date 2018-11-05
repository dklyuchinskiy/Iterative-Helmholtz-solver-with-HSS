#pragma once

/**********************************
Prototypes for all test functions.
**********************************/

typedef void(*ptr_test_low_rank)(int, int, double, char *);
typedef void(*ptr_test_sym_rec_compress)(int, double, char *, int);
typedef void(*ptr_test_mult_diag)(int, int, double, char *, int);
typedef void(*ptr_test_add)(int, dtype, dtype, double, char *, int);
typedef void(*ptr_test_update)(int, int, dtype, double, char*, int);
typedef void(*ptr_test_update2)(int, int, int, dtype, double, char*, int);
typedef void(*ptr_test)(int, double, double, int, double, char *);
typedef void(*ptr_test_shell)(ptr_test, const string&, int &, int &);

// Tests
void TestAll();

// Tests - BinaryTrees
void Test_LowRankApproxStruct(int m, int n, double eps, char *method);
void Test_SymRecCompressStruct(int n, double eps, char *method, int smallsize);
void Test_DiagMultStruct(int n, double eps, char *method, int smallsize);
void Test_RecMultLStruct(int n, int k, double eps, char *method, int smallsize);
void Test_UnsymmRecMultLStruct(int n, int k, double eps, char *method, int smallsize);
void Test_UnsymmRecMultRStruct(int n, int k, double eps, char *method, int smallsize);
void Test_AddStruct(int n, dtype alpha, dtype beta, double eps, char *method, int smallsize);
void Test_UnsymmAddStruct(int n, dtype alpha, dtype beta, double eps, char *method, int smallsize);
void Test_SymCompUpdate2Struct(int n, int k, dtype alpha, double eps, char* method, int smallsize);
void Test_UnsymmCompUpdate3Struct(int n, int k1, int k2, dtype alpha, double eps, char* method, int smallsize);
void Test_UnsymmCompUpdate2Struct(int n, int k, dtype alpha, double eps, char* method, int smallsize);
void Test_SymCompRecInvStruct(int n, double eps, char *method, int smallsize);
void Test_UnsymmCompRecInvStruct(int n, double eps, char *method, int smallsize);
void Test_CopyStruct(int n, double eps, char *method, int smallsize);
void Test_CopyUnsymmStruct(int n, double eps, char *method, int smallsize);
void Test_UnsymmLUfact(int n, double eps, char* method, int smallsize);
void Test_UnsymmRecCompressStruct(int n, double eps, char *method, int smallsize);


// Tests Shells
void Shell_LowRankApprox(ptr_test_low_rank func, const string& test_name, int &numb, int &fail_count);
void Shell_SymRecCompress(ptr_test_sym_rec_compress func, const string& test_name, int &numb, int &fail_count);
void Shell_DiagMult(ptr_test_sym_rec_compress func, const string& test_name, int &numb, int &fail_count);
void Shell_RecMultL(ptr_test_mult_diag func, const string& test_name, int &numb, int &fail_count);
void Shell_Add(ptr_test_add func, const string& test_name, int &numb, int &fail_count);
void Shell_SymCompUpdate2(ptr_test_update func, const string& test_name, int &numb, int &fail_count);
void Shell_UnsymmCompUpdate3(ptr_test_update2 func, const string& test_name, int &numb, int &fail_count);
void Shell_SymCompRecInv(ptr_test_sym_rec_compress func, const string& test_name, int &numb, int &fail_count);
void Shell_CopyStruct(ptr_test_sym_rec_compress func, const string& test_name, int &numb, int &fail_count);


// Solver
void Test_DirFactFastDiagStructOnline(size_m x, size_m y, cmnode** Gstr, dtype *B, double thresh, int smallsize);
void Test_TransferBlock3Diag_to_CSR(int n1, int n2, dcsr* Dcsr, dtype* x_orig, dtype *f, double eps);
void Test_DirSolveFactDiagStructBlockRanks(size_m x, size_m y, cmnode** Gstr);

#if 0

void Test_QueueList(int n, double eps, char* method, int smallsize);
void Test_Dense_to_CSR(size_m x, size_m y, size_m z, int non_zeros_in_3diag, double *D, int ldd);
void Test_CompareColumnsOfMatrix(int n1, int n2, int n3, double* D, int ldd, double* B, dcsr* Dcsr, double thresh);

void Test_DirSolveFactDiagStructConvergence(size_m x, size_m y, size_m z, mnode** Gstr, double thresh, int smallsize);



// TestOnline

void Test_RankEqual(mnode *Astr, mnode *AIstr);
void Test_RankAdd(mnode *Astr, mnode *Bstr, mnode* Cstr);
#endif