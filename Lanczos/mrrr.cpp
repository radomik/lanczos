#include "mrrr.hpp"
#include "Veci.hpp"
#include <clapack/clapack.h>

const double DOUBLE_EPS = 1e-13;

/**
 * Multiple Relatively Robust Representations for Tridiagonals.
 * 
 * @param ritz Preallocated output array of size a.size() == m
 * @param S    Preallocated output array of size m*m
 * @param a    Main diagonal vector of size m
 * @param b    Super-/Subdiagonal vector of size m
 * @param dbg  Print debug info
 *
 * @return 0 on success, non-zero value on error
 */
int mrrr(double* ritz, double* S, const Vecd& a, const Vecd& b, bool dbg) {
	// matrix size
	int m = a.size(); // m

	char jobz[] = "V";
	char range[] = "A";

	// main diagonal
	Vecd d(m); d.set(a);

	// sub-/superdiagonal
	Vecd e(m); e.set(b);

	double vl = 0, vu = 0; // unused
	int il = 0, iu = 0;    // unused
	double abstol = DOUBLE_EPS; // unused?
	int numEval;                // number of eigenvalues found

	// misc arguments
	double v_work1 = 0.0;
	int v_iwork1 = 0;
	int ldz = m;
	int lwork = -1;
	int liwork = -1;
	int info;

	Veci isuppz(m<<1);
	double* work1 = &v_work1;
	int* iwork1 = &v_iwork1;

	isuppz.zero();

	// query for workspace size
	dstegr_( jobz, range, &m, d.data(), e.data(), &vl, &vu, &il, &iu, &abstol,
		&numEval, ritz, S, &ldz, isuppz.data(), work1, &lwork, iwork1, &liwork,
		&info );

	if (info < 0) {
		fprintf(stderr, "%s: [ERROR 1] Input to dstegr had an illegal value [info: %d]\n", __FUNCTION__, info);
		return info;
	}

	if (dbg) {
		fprintf(stderr, "%s: [1] lwork: %d (%f), liwork: %d (%d)\n", __FUNCTION__, lwork, work1[0], liwork, iwork1[0]);
	}

	lwork  = (int)work1[0];
	liwork = iwork1[0];

	if (dbg) {
		fprintf(stderr, "%s: [2] lwork: %d (%f), liwork: %d (%d)\n", __FUNCTION__, lwork, work1[0], liwork, iwork1[0]);
	}
	
	Vecd work(lwork);
	Veci iwork(liwork);
	work.zero();
	iwork.zero();

	// call LAPACK routine
	dstegr_( jobz, range, &m, d.data(), e.data(), &vl, &vu, &il, &iu, &abstol,
		&numEval, ritz, S, &ldz, isuppz.data(), work.data(), &lwork, iwork.data(), &liwork,
		&info );

	if (info != 0) {
		fprintf(stderr, "%s: [ERROR 2] Input to dstegr had an illegal value [info: %d]\n", __FUNCTION__, info);
	}
	return info;
}
