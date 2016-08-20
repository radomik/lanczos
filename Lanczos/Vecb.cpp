#include "Vecb.hpp"

void Vecb::set(const Vecb& src, size_t count, size_t startIndex) {
	if (count == ((size_t)-1)) {
		count = m_size;
	}
	memcpy(m_data+startIndex, src.m_data+startIndex, count*sizeof(m_data[0]));
}

void Vecb::zero() {
	memset(m_data, 0, m_size*sizeof(m_data[0]));
}

void Vecb::print(const uint8_t* data, size_t size, FILE* file, const char* name) {
	fprintf(file, "%s = [\n", name);
	for (size_t i = 0; i < size; i++) {
		fprintf(file, "%s%u ", i ? " ;\n" : "", data[i]);
	}
	fprintf(file, " ]\n\n");
}
