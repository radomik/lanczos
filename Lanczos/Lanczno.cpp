#include "Lanczno.hpp"
#include "Vecd.hpp"
#include "Veci.hpp"
#include "TrMatrixd.hpp"
#include "Matrixd.hpp"
#include "mrrr.hpp"
#include "timing.h"
#include <algorithm>

/**
 * Lanczos with no reorthogonalization.
 *
 * @param ritz Preallocated array of size m
 * @param S Preallocated array of size m*m
 * @param A symmetric, square matrix (linearized reduced symmetric square matrix)
 * @param n size of matrix A
 * @param m Lanczos algorithm iteration count
 *
 * @return 0 on success, non-zero value on error
 */
int lanczno(double* ritz, double* S, const double* A, size_t n, uint32_t m, const Vecd* initStartVec) {
	assert(m > 1);

	Vecd startvec(n); //TODO: create constructor accepting double* data, size_t n which creates vector from preallocated buffer and does not free it on its own

	if (initStartVec) {
		if (initStartVec->size() != n) {
			fprintf(stderr, "%s: Invalid initStartVec size: %zu (expected: %zu)\n", __FUNCTION__, initStartVec->size(), n);
			return -1;
		}

		startvec.set(*initStartVec);
	}
	else {
		startvec.random();     // z = randn(n,1);
		startvec.normalize();  // startvec = z/norm(z);
	}

	Vecd v(n);
	Vecd v2(n);
	Vecd vt(n);
	Vecd r(n);
	Vecd rt(n);
	Vecd a(m);
	Vecd b(m);
	double anorm;

	v.set(startvec);
	a.zero();
	b.zero();

	v.print(stderr, "v");
	a.print(stderr, "a");
	b.print(stderr, "b");

	reset_and_start_timer();

	// initial Lanczos iteration, m steps
	for (uint32_t k = 0; ; ) {
		fprintf(stderr, "for: k == %u\n", k);
		if (k == 0) {
			TrMatrixd::mulByVector(r, A, v);  // r{n} = a{n}{n} * v{n}
		} else {
			TrMatrixd::mulByVector(rt, A, v);     // rt{n} = a{n}{n} * v{n}
			Vecd::mulByScalar(vt, v2, b[k-1]); // vt{n} = b[k-1]{1} * v2{n}
			Vecd::sub(r, rt, vt);              // r{n}  = rt{n} - vt{n}
		}

		r.print(stderr, "r");

		a[k] = Vecd::dotProduct(v, r);  // a[k]{1} = SUM v[i]*r[i]
		a.print(stderr, "a");

		Vecd::mulByScalar(vt, v, a[k]); // vt{n}   = a[k]{1} * v{n}
		Vecd::sub(r, r, vt);            // r{n}    = r{n} - vt{n}
		r.print(stderr, "r");
		b[k] = r.getNorm();             // b[k]{1} = |r{n}|
		b.print(stderr, "b");

		// estimate |A|_2 by |T|_1
		if (k == 0) {
			anorm = fabs(a[0] + b[0]);
		} else {
			anorm = std::max(anorm, b[k-1]+fabs(a[k])+b[k]);
		}

		fprintf(stderr, "anorm = %f\n\n", anorm);

		if (++k == m) break;

		// prepare next step, k = {1 ... m-1}
		v2.set(v);                          // v2{n} = v{n}
		v2.print(stderr, "v2");
		Vecd::divByScalar(v, r, b[k-1]); // EACH_i v[i] = r[i]/b[k-1]
		v.print(stderr, "v");
	}

	double ts = get_elapsed_mcycles();
	double ts_total = ts;

	fprintf(stderr, "%s: Initial Lanczos iteration completed in %f mcycles [total: %f]\n", __FUNCTION__, ts, ts_total);

	// tridiagonal matrix
	reset_and_start_timer();
	int info = mrrr(ritz, S, a, b, true);

	ts = get_elapsed_mcycles();
	ts_total += ts;
	fprintf(stderr, "%s: mrrr() completed in %f mcycles [total: %f]\n", __FUNCTION__, ts, ts_total);

	if (info != 0) {
		fprintf(stderr, "%s: mrrr() failed with status %d\n", __FUNCTION__, info);
		return info;
	}

	Vecd::print(ritz, m, stderr, "ritz");
	Matrixd::print(S, m, stderr, "S");

	//~ // residual estimation
	//~ const double* evec_first = evec;
	//~ const double* evec_last = evec + (m-1)*m; // offset is a desired row index multiplied by column count
	//~ const double  b_last = b[m-1];
	//~ const double tol = DOUBLE_EPS * anorm;
	//~ Vecd lres(m);
	//~ Vecd cul(m);
	//~ Veci idx(m);

	//~ for (size_t i = 0; i < m; i++) {
		//~ double lres_i = fabs(b_last * evec_last[i]);
		//~ double cul_i = fabs(evec_first[i]);

		//~ // remove non-converged and spurious Ritz values
		//~ idx[i] = (lres_i < tol) && (cul_i > tol);

		//~ lres[i] = lres_i;
		//~ cul[i] = cul_i;
	//~ }

	//~ ts = get_elapsed_mcycles();
	//~ ts_total += ts;

	//~ fprintf(stderr, "%s: mrrr() and residual estimation steps completed in %f mcycles\n", __FUNCTION__, ts);


	//~ int idx = (lres < tol) && (cul > tol);
	//~ double ritz = eval[idx];

	return 0;
}
