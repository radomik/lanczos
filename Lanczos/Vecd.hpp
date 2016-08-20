#ifndef VECD_HPP
#define VECD_HPP

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <stdint.h>
#include <string>

class Vecd {
public:
	Vecd() : m_data(NULL), m_size(0), m_free(true) { }

	Vecd(size_t size) {
		m_data = new double[size];
		m_size = size;
		m_free = true;
	}

	Vecd(double* data, size_t size)
		: m_data(data), m_size(size), m_free(false) { }

	~Vecd() {
		if (m_data && m_free) {
			delete[] m_data;
			m_data = NULL;
		}
	}

	void random();
	void zero();
	double getNormPow() const;
	double getNorm() const { return sqrt(getNormPow()); }
	void   normalize();
	void print(FILE* file, const char* name) const {
		print(m_data, m_size, file, name);
	}
	static void print(const double* data, size_t size, FILE* file, const char* name);

	void set(const Vecd& src, size_t count = ((size_t)-1), size_t startIndex = 0);

	const double& operator[](size_t index) const { return m_data[index]; }
	double& operator[](size_t index) { return m_data[index]; }

	double* data() { return m_data; }
	const double* dataConst() const { return m_data; }
	size_t size() const { return m_size; }

	static void mulByScalar(Vecd& res, const Vecd& v, double a);
	static void divByScalar(Vecd& res, const Vecd& v, double a);

	/**
	 * @return EACH res[i] = 1 / (a * v[i])
	 */
	static void mulInvByScalar(Vecd& res, const Vecd& v, double a);
	static void sub(Vecd& res, const Vecd& v, const Vecd& u);
	static void add(Vecd& res, const Vecd& v, const Vecd& u);

	/**
	 *
	 * @return SUM v[i]*u[i]
	 */
	static double dotProduct(const Vecd& v, const Vecd& u);
	
	int readCsv(const char* filename, bool countLines = false);
	int saveCsv(const char* filename);
	
	int saveBin(const char* filename) const;
	int readBin(const char* filename);
	
	int readCsvOrBin(const std::string& work_dir, const std::string& filename_without_ext);
	
	bool equals(const Vecd& other, float eps, size_t* first_diff_index = NULL, float* first_diff = NULL) const;
private:
	double* m_data;
	size_t  m_size;
	bool    m_free;
	
	int _resize(size_t size);
};

#endif
