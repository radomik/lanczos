#include "TrMatrixd.hpp"
#include "string_utils.hpp"
#include <algorithm>
#include <fstream>
#include <cstring>
#include <climits>
#include <errno.h>
#include <file_utils.hpp>

/**
 * Print linearized reduced symmetric square matrix to file using Matlab syntax.
 *
 * @param A    Linearized reduced symmetric square matrix
 * @param n    Size of matrix
 * @param file File to be printed to
 * @param name Matlab variable name
 */
void TrMatrixd::print(const double* A, size_t n, FILE* file, const char* name) {
	const size_t t = (n<<1) + 1;
	const size_t u = n - 1;

	fprintf(file, "%s = [\n", name);
	for (size_t i = 0; i < n; i++) {
		size_t ai = ((t - i) * i)>>1;

		if (i != 0) { // elements left from diagonal
			size_t step = u;
			size_t sta_curr = i;

			fprintf(file, ";\n");

			for ( ; sta_curr < ai; ) {
				fprintf(file, "%15.10f ", A[sta_curr]);
				sta_curr += step;
				step--;
			}
		}

		// elements right of and including diagonal
		const size_t ai_end = ai + (n-i);
		for ( ; ai < ai_end; ) {
			fprintf(file, "%15.10f ", A[ai++]);
		}
	}
	fprintf(file, "]\n");
}

/**
 * Print linearized reduced symmetric square matrix to CSV file.
 *
 * @param A    Linearized reduced symmetric square matrix
 * @param n    Size of matrix
 * @param filename File name to be printed to
 */
int TrMatrixd::saveCsv(const char* filename, const double* A, size_t n) {
	FILE* f = fopen(filename, "wb");
	if (! f) {
		fprintf(stderr, "%s: Error opening file '%s' [%s]\n", __FUNCTION__, filename, strerror(errno));
		return -1;
	}
	
	const size_t t = (n<<1) + 1;
	const size_t u = n - 1;

	for (size_t i = 0; i < n; i++) {
		size_t ai = ((t - i) * i)>>1;

		if (i != 0) { // elements left from diagonal
			size_t step = u;
			size_t sta_curr = i;

			fprintf(f, "\n");

			for ( ; sta_curr < ai; ) {
				fprintf(f, "%lf,", A[sta_curr]);
				sta_curr += step;
				step--;
			}
		}

		// elements right of and including diagonal
		const size_t ai_end = ai + (n-i);
		for ( ; ; ) {
			fprintf(f, "%lf", A[ai]);
			if (++ai == ai_end) break;
			fprintf(f, ",");
		}
	}
	
	fclose(f);
	fprintf(stderr, "%s: Saved file '%s'\n", __FUNCTION__, filename);
	return 0;
}

/**
 * Multiply linearized reduced symmetric square matrix A by vector b.
 *
 * @param res Product vector of size n
 * @param A   Multiplicand - linearized reduced symmetric square matrix
 * @param b   Multiplier vector of size n
 * @param n   Size of operands
 */
void TrMatrixd::mulByVector(double* res, const double* A, const double* b, size_t n) {
	const size_t t = (n<<1) + 1;
	const size_t u = n - 1;

	// loop is optimized so each iteration is independent and may be parallelized
	for (size_t i = 0; i < n; i++) {
		size_t ai = ((t - i) * i)>>1, bi = 0;
		double res_i = 0.0;

		if (i != 0) { // elements left from diagonal
			size_t step = u;
			size_t sta_curr = i;
			for ( ; sta_curr < ai; ) {
				res_i += b[bi++] * A[sta_curr];
				sta_curr += step;
				step--;
			}
		}

		// elements right of and including diagonal
		const size_t ai_end = ai + (n-i);
		for ( ; ai < ai_end; ) {
			res_i += b[bi++] * A[ai++];
		}

		res[i] = res_i;
	}
}

double* TrMatrixd::readCsv(const char* filename, const char* separator, size_t* out_n, size_t* out_nelem) {
	std::ifstream file(filename);

	if (! file) {
		fprintf(stderr, "%s: Failed to open '%s' [%s]\n", __FUNCTION__, filename, strerror(errno));
		return NULL;
	}

	std::string delim = separator;
	std::string line;
	double* data = NULL;
	size_t  i = 0;
	size_t  di = 0;
	size_t  n;
	bool    error = false;
	n = 0;

	while (std::getline(file, line)) {
		if (i == 0) {
			n = std::count(line.begin(), line.end(), delim[0]);
			if (n == 0) {
				break;
			}
			n++;
			if (out_n) {
				*out_n = n;
			}
			size_t nelem = TrMatrixd::numElemFromN(n);
			if (out_nelem) {
				*out_nelem = nelem;
			}
			data = new double[nelem];
		}

		size_t j = 0;
		size_t pos;
		while ((pos = line.find(delim)) != std::string::npos) {
			if (j >= i) {
				if (j >= n) {
					fprintf(stderr, "%s: [1] Invalid column index %zu at row %zu\n", __FUNCTION__, j, i);
					error = true;
					break;
				}

				std::string token = line.substr(0, pos);
				string_utils::trim(token);

				if (sscanf(token.c_str(), "%lf", &data[di++]) != 1) {
					fprintf(stderr, "%s: [1] Cannot read value at row index: %zu and column index: %zu\n", __FUNCTION__, i, j);
					error = true;
					break;
				}
			}

			line.erase(0, pos + delim.length());
			j++;
		}

		if (error) {
			break;
		}

		string_utils::trim(line);
		if (! line.empty()) {
			if (j >= n) {
				fprintf(stderr, "%s: [2] Invalid column index %zu at row %zu\n", __FUNCTION__, j, i);
				error = true;
				break;
			}

			if (sscanf(line.c_str(), "%lf", &data[di++]) != 1) {
				fprintf(stderr, "%s: [2] Cannot read value at row index: %zu and column index: %zu\n", __FUNCTION__, i, j);
				error = true;
				break;
			}
		}

		i++;
	}

	if (n == 0) {
		fprintf(stderr, "%s: No data in file '%s'\n", __FUNCTION__, filename);
	}
	else {
		if (error) {
			delete[] data;
			data = NULL;
		}
	}

	file.close();
	return data;
}

int TrMatrixd::saveBin(const char* filename, const double* A, size_t n) {
	if (n > UINT_MAX) {
		fprintf(stderr, "%s: Error. Matrix size n=%zu is above UINT_MAX=%u\n", __FUNCTION__, n, UINT_MAX);
		return -1;
	}

	uint32_t nu = (uint32_t)n;
	size_t nelem = TrMatrixd::numElemFromN(n);

	fprintf(stderr, "%s: Saving matrix of size n=%u, nelem=%zu, sizeof(A[0])=%zu\n", __FUNCTION__, nu, nelem, sizeof(A[0]));

	FILE* bin = fopen(filename, "wb");
	if (! bin) {
		fprintf(stderr, "%s: Error opening output file '%s' [%s]\n", __FUNCTION__, filename, strerror(errno));
		return -1;
	}
	
	int ret = (fwrite(&nu, sizeof(nu), 1, bin) != 1) || (fwrite(A, sizeof(A[0]), nelem, bin) != nelem);
	
	if (ret) {
		fprintf(stderr, "%s: Error writing file '%s' [%s]\n", __FUNCTION__, filename, strerror(errno));
		ret = -1;
	}
	else {
		fprintf(stderr, "%s: Saved file '%s'\n", __FUNCTION__, filename);
	}
	
	fclose(bin);
	return ret;
}

double* TrMatrixd::readBin(const char* filename, size_t* out_n, size_t* out_nelem) {
	FILE* bin = fopen(filename, "rb");
	if (! bin) {
		fprintf(stderr, "%s: Error opening input file '%s' [%s]\n", __FUNCTION__, filename, strerror(errno));
		return NULL;
	}
	
	uint32_t nu;
	size_t ret;
	if (((ret=fread(&nu, sizeof(nu), 1, bin)) != 1) || (nu == 0)) {
		fprintf(stderr, "%s: Error reading matrix size (n=%u, ret=%zu) [%s]\n", __FUNCTION__, nu, ret, strerror(errno));
		fclose(bin);
		return NULL;
	}
	
	double* A = NULL;
	size_t nelem = TrMatrixd::numElemFromN(nu);
	
	if (out_n) {
		*out_n = nu;
	}
	if (out_nelem) {
		*out_nelem = nelem;
	}
	
	fprintf(stderr, "%s: Reading matrix of size n=%u, nelem=%zu, sizeof(A[0])=%zu\n", __FUNCTION__, nu, nelem, sizeof(A[0]));
	
	A = new double[nelem];
	
	if ((ret=fread(A, sizeof(A[0]), nelem, bin)) != nelem) {
		fprintf(stderr, "%s: Error reading matrix content (nelem=%zu, ret=%zu) [%s]\n", __FUNCTION__, nelem, ret, strerror(errno));
		fclose(bin);
		delete[] A;
		return NULL;
	}
	
	fclose(bin);
	return A;
}

double* TrMatrixd::readCsvOrBin(const std::string& work_dir, const std::string& filename_without_ext, size_t* out_n, size_t* out_nelem) {
	CsvBinPaths p(work_dir, filename_without_ext);
	double* A;
	size_t n;
	if (! (A=TrMatrixd::readBin(p.bin(), &n, out_nelem))) {
		if (! (A=TrMatrixd::readCsv(p.csv(), ",", &n, out_nelem))) {
			fprintf(stderr, "%s: Error reading matrix [csv=%s, bin=%s]\n", __FUNCTION__, p.csv(), p.bin());
			return NULL;
		}
		TrMatrixd::saveBin(p.bin(), A, n);
	}
	if (out_n) {
		*out_n = n;
	}
	return A;
}

bool TrMatrixd::equals(const double* X, size_t Xn, const double* Y, size_t Yn, float eps, size_t* first_diff_index, float* first_diff) {
	if (!(X && Y)) {
		fprintf(stderr, "%s: Warning at least one matrix is null (X: %p, Y: %p)\n", __FUNCTION__, X, Y);
		return false;
	}
	if (Xn != Yn) {
		return false;
	}
	
	size_t nelem = TrMatrixd::numElemFromN(Xn);
	
	for (size_t i = 0; i < nelem; i++) {
		const float diff = fabsf(X[i]-Y[i]);
		if (diff > eps) {
			if (first_diff_index) {
				*first_diff_index = i;
			}
			if (first_diff) {
				*first_diff = diff;
			}
			return false;
		}
	}
	
	if (first_diff_index) {
		*first_diff_index = 0;
	}
	if (first_diff) {
		*first_diff = 0.0f;
	}
	return true;
}
