#ifndef MATRIXD_HPP
#define MATRIXD_HPP

#include <cstdio>
#include <cstdlib>
#include <string>

/**
 * Operations on linearized square matrices.
 */
namespace Matrixd {
	/**
	 * Print linearized square matrix to file using Matlab syntax.
	 *
	 * @param A    Linearized square matrix
	 * @param n    Size of matrix
	 * @param file File to be printed to
	 * @param name Matlab variable name
	 */
	void print(const double* A, size_t n, FILE* file, const char* name);
	
	/**
	 * Print linearized symmetric square matrix to CSV file.
	 *
	 * @param A    Linearized symmetric square matrix
	 * @param n    Size of matrix
	 * @param filename File name to be printed to
	 */
	int saveCsv(const char* filename, const double* A, size_t n);

	double* readCsv(const char* filename, const char* separator, size_t* out_n, size_t* out_nelem);

	int saveBin(const char* filename, const double* A, size_t n);

	double* readBin(const char* filename, size_t* out_n, size_t* out_nelem);

	double* readCsvOrBin(const std::string& work_dir, const std::string& filename_without_ext, size_t* out_n, size_t* out_nelem);
	
	inline size_t numElemFromN(size_t n) { return n*n; }
};

#endif
