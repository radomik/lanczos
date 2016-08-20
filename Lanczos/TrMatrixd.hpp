#ifndef TRMATRIXD_HPP
#define TRMATRIXD_HPP

#include "Vecd.hpp"

/**
 * Operations on linearized reduced symmetric square matrices of form:
 * [(0,0)...(0,n-1);(1,1)...(1,n-1);(2,2)...(2,n-1);...;(n-2,n-2);(n-2,n-1);(n-1,n-1)]
 * where each row starts with diagonal element.
 */
namespace TrMatrixd {
	/**
	 * Print linearized reduced symmetric square matrix to file using Matlab syntax.
	 *
	 * @param A    Linearized reduced symmetric square matrix
	 * @param n    Size of matrix
	 * @param file File to be printed to
	 * @param name Matlab variable name
	 */
	void print(const double* A, size_t n, FILE* file, const char* name);
	
	/**
	 * Print linearized reduced symmetric square matrix to CSV file.
	 *
	 * @param A    Linearized reduced symmetric square matrix
	 * @param n    Size of matrix
	 * @param filename File name to be printed to
	 */
	int saveCsv(const char* filename, const double* A, size_t n);

	/**
	 * Multiply linearized reduced symmetric square matrix A by vector b.
	 *
	 * @param res Product vector of size n
	 * @param A   Multiplicand - linearized reduced symmetric square matrix
	 * @param b   Multiplier vector of size n
	 * @param n   Size of operands
	 */
	void mulByVector(double* res, const double* A, const double* b, size_t n);
	inline void mulByVector(Vecd& res, const double* A, const Vecd& b) {
		mulByVector(res.data(), A, b.dataConst(), b.size());
	}

	double* readCsv(const char* filename, const char* separator, size_t* out_n = NULL, size_t* out_nelem = NULL);
	
	int     saveBin(const char* filename, const double* A, size_t n);
	double* readBin(const char* filename, size_t* out_n = NULL, size_t* out_nelem = NULL);
	
	double* readCsvOrBin(const std::string& work_dir, const std::string& filename_without_ext, size_t* out_n = NULL, size_t* out_nelem = NULL);
	
	inline size_t numElemFromN(size_t n) { return (n*(n+1))>>1; }
	bool equals(const double* X, size_t Xn, const double* Y, size_t Yn, float eps, size_t* first_diff_index = NULL, float* first_diff = NULL);
};

#endif
