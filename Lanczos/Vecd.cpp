#include "Vecd.hpp"
#include <cstring>
#include <cassert>
#include <ctime>
#include <cstring>
#include <climits>
#include <errno.h>
#include <file_utils.hpp>

#define EACH_I(v, expr) \
for (size_t i = 0; i < (v).m_size; i++) { \
	expr; \
}

#define EACH_SI(v, expr) \
const double* end = (v).m_data + (v).m_size; \
for (double* si = (v).m_data; si < end; si++) { \
	expr; \
}

#define EACH_CONST_SI(v, expr) \
const double* end = (v).m_data + (v).m_size; \
for (const double* si = (v).m_data; si < end; si++) { \
	expr; \
}

void Vecd::set(const Vecd& src, size_t count, size_t startIndex) {
	if (count == ((size_t)-1)) {
		count = m_size;
	}
	memcpy(m_data+startIndex, src.m_data+startIndex, count*sizeof(m_data[0]));
}

void Vecd::random() {
	assert(sizeof(m_data[0]) % sizeof(int) == 0);
	srand(time(NULL));

	const int* vi_end = (int*)(m_data+m_size);
	for (int* vi = (int*)m_data; vi != vi_end; vi++) {
		*vi = rand();
	}
}

double Vecd::getNormPow() const {
	double s = 0.0;
	EACH_SI(*this, s += (*si) * (*si));
	return s;
}

void Vecd::normalize() {
	const double norm = getNorm();
	EACH_SI(*this, *si /= norm);
}

void Vecd::zero() {
	EACH_SI(*this, *si = 0.0);
}

void Vecd::print(const double* data, size_t size, FILE* file, const char* name) {
	fprintf(file, "%s = [\n", name);
	for (size_t i = 0; i < size; i++) {
		fprintf(file, "%s%.4f ", i ? " ;\n" : "", data[i]);
	}
	fprintf(file, " ]\n\n");
}

void Vecd::mulByScalar(Vecd& res, const Vecd& v, double a) {
	EACH_I(res, res[i] = v[i] * a);
}

void Vecd::divByScalar(Vecd& res, const Vecd& v, double a) {
	EACH_I(res, res[i] = v[i] / a);
}

void Vecd::mulInvByScalar(Vecd& res, const Vecd& v, double a) {
	EACH_I(res, res[i] = 1.0 / (v[i] * a));
}

void Vecd::sub(Vecd& res, const Vecd& v, const Vecd& u) {
	EACH_I(res, res[i] = v[i] - u[i]);
}

void Vecd::add(Vecd& res, const Vecd& v, const Vecd& u) {
	EACH_I(res, res[i] = v[i] + u[i]);
}

double Vecd::dotProduct(const Vecd& v, const Vecd& u) {
	double s = 0.0;
	EACH_I(v, s += v[i] * u[i]);
	return s;
}

int Vecd::readCsv(const char* filename, bool countLines) {
	FILE* f = fopen(filename, "rb");
	if (! f) {
		fprintf(stderr, "%s: Error opening file '%s' [%s]\n", __FUNCTION__, filename, strerror(errno));
		return -1;
	}
	
	char line[32];
	if (countLines) {
		size_t n = 0;
		while (fgets(line, sizeof(line), f) != NULL) {
			n++;
		}
		
		fprintf(stderr, "%s: File '%s' line count: %zu\n", __FUNCTION__, filename, n);
		
		if (_resize(n) != 0) {
			fprintf(stderr, "%s: Failed to resize vector to %zu\n", __FUNCTION__, n);
			fclose(f);
			return -1;
		}
		
		rewind(f);
	}
	
	size_t i = 0;
	while ((i < m_size) && (fgets(line, sizeof(line), f) != NULL)) {
		if (sscanf(line, "%lf", &m_data[i]) != 1) {
			fprintf(stderr, "%s: Error reading value at row index %zu\n", __FUNCTION__, i);
			fclose(f);
			return -1;
		}
		i++;
	}

	fclose(f);

	if (i != m_size) {
		fprintf(stderr, "%s: Premature end of file (read: %zu, expected: %zu)\n", __FUNCTION__, i, m_size);
		return -1;
	}

	return 0;
}

int Vecd::saveCsv(const char* filename) {
	FILE* f = fopen(filename, "wb");
	if (! f) {
		fprintf(stderr, "%s: Error opening output file '%s' [%s]\n", __FUNCTION__, filename, strerror(errno));
		return -1;
	}
	EACH_CONST_SI(*this, fprintf(f, "%lf\n", *si));
	fclose(f);
	fprintf(stderr, "%s: Saved file '%s'\n", __FUNCTION__, filename);
	return 0;
}

int Vecd::saveBin(const char* filename) const {
	if ((!m_data) || (m_size==0)) {
		fprintf(stderr, "%s: Cannot save empty vector\n", __FUNCTION__);
		return -1;
	}
	
	if (m_size > UINT_MAX) {
		fprintf(stderr, "%s: Error. Vector size n=%zu is above UINT_MAX=%u\n", __FUNCTION__, m_size, UINT_MAX);
		return -1;
	}
	
	uint32_t nu = (uint32_t)m_size;
	
	FILE* bin = fopen(filename, "wb");
	if (! bin) {
		fprintf(stderr, "%s: Error opening output file '%s' [%s]\n", __FUNCTION__, filename, strerror(errno));
		return -1;
	}
	
	size_t ret;
	if ((ret=fwrite(&nu, sizeof(nu), 1, bin)) != 1) {
		fprintf(stderr, "%s: Error writing vector size (ret=%zu) [%s]\n", __FUNCTION__, ret, strerror(errno));
		fclose(bin);
		return -1;
	}
	
	if ((ret=fwrite(m_data, sizeof(m_data[0]), nu, bin)) != nu) {
		fprintf(stderr, "%s: Error writing vector content (n=%u, ret=%zu) [%s]\n", __FUNCTION__, nu, ret, strerror(errno));
		fclose(bin);
		return -1;
	}
	
	fclose(bin);
	fprintf(stderr, "%s: Saved file '%s'\n", __FUNCTION__, filename);
	return 0;
}

int Vecd::readBin(const char* filename) {
	FILE* bin = fopen(filename, "rb");
	if (! bin) {
		fprintf(stderr, "%s: Error opening input file '%s' [%s]\n", __FUNCTION__, filename, strerror(errno));
		return -1;
	}
	
	uint32_t nu;
	size_t ret;
	if (((ret=fread(&nu, sizeof(nu), 1, bin)) != 1) || (nu == 0)) {
		fprintf(stderr, "%s: Error reading vector size (n=%u, ret=%zu) [%s]\n", __FUNCTION__, nu, ret, strerror(errno));
		fclose(bin);
		return -1;
	}
	
	if (_resize(nu) != 0) {
		fprintf(stderr, "%s: Failed to resize vector to %u\n", __FUNCTION__, nu);
		fclose(bin);
		return -1;
	}
	
	if ((ret=fread(m_data, sizeof(m_data[0]), nu, bin)) != nu) {
		fprintf(stderr, "%s: Error reading vector content (n=%u, ret=%zu) [%s]\n", __FUNCTION__, nu, ret, strerror(errno));
		fclose(bin);
		return -1;
	}
	
	fclose(bin);
	return 0;
}

int Vecd::readCsvOrBin(const std::string& work_dir, const std::string& filename_without_ext) {
	CsvBinPaths p(work_dir, filename_without_ext);
	
	if (readBin(p.bin()) != 0) {
		if (readCsv(p.csv(), true) != 0) {
			fprintf(stderr, "%s: Error reading vector [csv=%s, bin=%s]\n", __FUNCTION__, p.csv(), p.bin());
			return -1;
		}
		saveBin(p.bin());
	}
	return 0;
}

int Vecd::_resize(size_t size) {
	bool alloc = true;
	
	if (m_data) {
		if ((alloc = (size != m_size))) {
			if (! m_free) {
				fprintf(stderr, "%s: Cannot dispose vector data with m_free=false\n", __FUNCTION__);
				return -1;
			}
			delete[] m_data;
		}
	}
	
	if (alloc) {
		m_data = new double[size];
		m_size = size;
		m_free = true; 
	}
	
	return 0;
}

bool Vecd::equals(const Vecd& other, float eps, size_t* first_diff_index, float* first_diff) const {
	if ((! m_data) || (! other.m_data)) {
		fprintf(stderr, "%s: Warning at least one vector is null (m_data: %p, other.m_data: %p)\n", __FUNCTION__, m_data, other.m_data);
		return false;
	}
	
	if (m_size != other.m_size) {
		return false;
	}
	
	for (size_t i = 0; i < m_size; i++) {
		const float diff = fabsf(m_data[i]-other.m_data[i]);
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
