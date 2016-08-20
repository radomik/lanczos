#include "Veci.hpp"
#include <cstring>
#include <cassert>
#include <ctime>

void Veci::set(const Veci& src, size_t count, size_t startIndex) {
	if (count == ((size_t)-1)) {
		count = m_size;
	}
	memcpy(m_data+startIndex, src.m_data+startIndex, count*sizeof(m_data[0]));
}

void Veci::random() {
	assert(sizeof(m_data[0]) % sizeof(int) == 0);
	srand(time(NULL));

	const int* vi_end = (int*)(m_data+m_size);
	for (int* vi = (int*)m_data; vi != vi_end; vi++) {
		*vi = rand();
	}
}

void Veci::zero() {
	memset(m_data, 0, m_size*sizeof(m_data[0]));
}

void Veci::print(const int* data, size_t size, FILE* file, const char* name) {
	fprintf(file, "%s = [\n", name);
	for (size_t i = 0; i < size; i++) {
		fprintf(file, "%s%d ", i ? " ;\n" : "", data[i]);
	}
	fprintf(file, " ]\n\n");
}
