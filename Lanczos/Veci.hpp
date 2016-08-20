#ifndef VECI_HPP
#define VECI_HPP

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <stdint.h>

class Veci {
public:
	Veci() : m_data(NULL), m_size(0), m_free(true) { }

	Veci(size_t size) {
		m_data = new int[size];
		m_size = size;
		m_free = true;
	}

	Veci(int* data, size_t size)
		: m_data(data), m_size(size), m_free(false) { }

	~Veci() {
		if (m_data && m_free) {
			delete[] m_data;
			m_data = NULL;
		}
	}

	void random();
	void zero();
	void print(FILE* file, const char* name) const {
		print(m_data, m_size, file, name);
	}
	static void print(const int* data, size_t size, FILE* file, const char* name);

	void set(const Veci& src, size_t count = ((size_t)-1), size_t startIndex = 0);

	const int& operator[](size_t index) const { return m_data[index]; }
	int& operator[](size_t index) { return m_data[index]; }

	int* data() { return m_data; }
	const int* dataConst() const { return m_data; }
	size_t size() const { return m_size; }

private:
	int*    m_data;
	size_t  m_size;
	bool    m_free;
};

#endif
