#include <cstdio>

#define ARR_LEN(a) (sizeof(a)/sizeof((a)[0]))

void printReducedArray(const char* name, const double* A, size_t n) {
	const size_t t = (n<<1) + 1;
	const size_t u = n - 1;

	printf("%s = [\n", name);
	for (size_t i = 0; i < n; i++) {
		size_t ai = ((t - i) * i)>>1;

		if (i != 0) { // elements left from diagonal
			size_t step = u;
			size_t sta_curr = i;

			printf(";\n");

			for ( ; sta_curr < ai; ) {
				printf("%8.2f ", A[sta_curr]);
				sta_curr += step;
				step--;
			}
		}

		// elements right of and including diagonal
		const size_t ai_end = ai + (n-i);
		for ( ; ai < ai_end; ) {
			printf("%8.2f ", A[ai++]);
		}
	}
	printf("]\n");
}

void printVector(const char* name, const double* b, size_t n) {
	printf("%s = [ ", name);
	for (size_t i = 0; i < n; i++) {
		printf("%s%8.2f ", i ? "; " : "", b[i]);
	}
	printf("]\n");
}

/**
 * Multiply linearized reduced symmetric square matrix A by vector b.
 *
 * @param res Product vector of size n
 * @param A   Multiplicand matrix linearized by rows and each row starts with element at diagonal:
 *            [(0,0)...(0,n-1);(1,1)...(1,n-1);(2,2)...(2,n-1);...;(n-2,n-2);(n-2,n-1);(n-1,n-1)]
 * @param b   Multiplier vector of size n
 * @param n   Size of operands
 */
void matmulvec(double* res, const double* A, const double* b, size_t n) {
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

int main(int argc, char** argv) {
	const double b[] = {
		1.0,  2.0,  3.0,  4.0,  5.0
	};

	const double A[((ARR_LEN(b))*(ARR_LEN(b)+1)) >> 1] = {
		1.0,  2.0,  3.0,  4.0,  5.0,
			  6.0,  7.0,  8.0,  9.0,
				   10.0, 11.0, 12.0,
						 13.0, 14.0,
							   15.0
	};

	double res[ARR_LEN(b)];

	printReducedArray("A", A, ARR_LEN(b));
	printVector("b", b, ARR_LEN(b));

	matmulvec(res, A, b, ARR_LEN(b));

	printVector("res", res, ARR_LEN(b));

	return 0;
}
