#include "Matrixd.hpp"
#include <file_utils.hpp>
#include <string_utils.hpp>
#include <algorithm>
#include <fstream>
#include <cstring>
#include <climits>
#include <errno.h>
#include <stdint.h>

/**
 * Print linearized square matrix to file using Matlab syntax.
 *
 * @param A    Linearized square matrix
 * @param n    Size of matrix
 * @param file File to be printed to
 * @param name Matlab variable name
 */
void Matrixd::print(const double* A, size_t n, FILE* file, const char* name) {
	fprintf(file, "%s = [\n", name);
	for (size_t i = 0; i < n; i++) {
		size_t i_off = i*n;
		if (i != 0) {
			fprintf(file, ";\n");
		}
		for (size_t j = 0; j < n; j++) {
			fprintf(file, "%15.10f ", A[i_off+j]);
		}
	}
	fprintf(file, "]\n");
}

/**
 * Print linearized symmetric square matrix to CSV file.
 *
 * @param A    Linearized symmetric square matrix
 * @param n    Size of matrix
 * @param filename File name to be printed to
 */
int Matrixd::saveCsv(const char* filename, const double* A, size_t n) {
	FILE* f = fopen(filename, "wb");
	if (! f) {
		fprintf(stderr, "%s: Error opening file '%s' [%s]\n", __FUNCTION__, filename, strerror(errno));
		return -1;
	}
	
	for (size_t i = 0; i < n; i++) {
		size_t i_off = i*n;
		if (i != 0) {
			fprintf(f, "\n");
		}
		for (size_t j = 0; j < n; j++) {
			fprintf(f, "%s%15.10f ", j?",":"", A[i_off+j]);
		}
	}
	
	fclose(f);
	return 0;
}

double* Matrixd::readCsv(const char* filename, const char* separator, size_t* out_n, size_t* out_nelem) {
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
			size_t nelem = Matrixd::numElemFromN(n);
			if (out_nelem) {
				*out_nelem = nelem;
			}
			data = new double[nelem];
		}

		size_t j = 0;
		size_t pos;
		while ((pos = line.find(delim)) != std::string::npos) {
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

int Matrixd::saveBin(const char* filename, const double* A, size_t n) {
	if (n > UINT_MAX) {
		fprintf(stderr, "%s: Error. Matrix size n=%zu is above UINT_MAX=%u\n", __FUNCTION__, n, UINT_MAX);
		return -1;
	}

	uint32_t nu = (uint32_t)n;
	size_t nelem = Matrixd::numElemFromN(n);

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
		printf("File '%s' saved successfully\n", filename);
	}
	
	fclose(bin);
	return ret;
}

double* Matrixd::readBin(const char* filename, size_t* out_n, size_t* out_nelem) {
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
	size_t nelem = Matrixd::numElemFromN(nu);
	
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

double* Matrixd::readCsvOrBin(const std::string& work_dir, const std::string& filename_without_ext, size_t* out_n, size_t* out_nelem) {
	CsvBinPaths p(work_dir, filename_without_ext);
	double* A;
	size_t n;
	if (! (A=Matrixd::readBin(p.bin(), &n, out_nelem))) {
		if (! (A=Matrixd::readCsv(p.csv(), ",", &n, out_nelem))) {
			fprintf(stderr, "%s: Error reading matrix [csv=%s, bin=%s]\n", __FUNCTION__, p.csv(), p.bin());
			return NULL;
		}
		Matrixd::saveBin(p.bin(), A, n);
		if (out_n) *out_n = n;
	}
	return A;
}
