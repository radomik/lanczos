#ifndef VECB_HPP
#define VECB_HPP

#include "utils.hpp"

class Vecb {
public:
	Vecb() : m_data(NULL), m_size(0), m_free(true) { }

	Vecb(size_t size) {
		m_data = new uint8_t[size];
		m_size = size;
		m_free = true;
	}

	Vecb(uint8_t* data, size_t size) : m_data(data), m_size(size), m_free(false) { }

	~Vecb() {
		if (m_data && m_free) {
			delete[] m_data;
			m_data = NULL;
		}
	}

	void zero();
	void print(FILE* file, const char* name) const {
		print(m_data, m_size, file, name);
	}
	static void print(const uint8_t* data, size_t size, FILE* file, const char* name);

	void set(const Vecb& src, size_t count = ((size_t)-1), size_t startIndex = 0);

	const uint8_t& operator[](size_t index) const { return m_data[index]; }
	uint8_t& operator[](size_t index) { return m_data[index]; }

	uint8_t* data() { return m_data; }
	const uint8_t* dataConst() const { return m_data; }
	size_t size() const { return m_size; }

private:
	uint8_t*    m_data;
	size_t  m_size;
	bool    m_free;
};

#endif
