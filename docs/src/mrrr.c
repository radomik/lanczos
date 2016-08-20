#include "mex.h"
#include <clapack/clapack.h>
#include <stdio.h>

void mexFunction(
	int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[]
) {
	if ((nrhs != 2) || (nlhs != 2)) {
		mexErrMsgTxt("Expects two input and output arguments");
	}
	
	char s[128];
	
	/* matrix size */
	int n1 = mxGetM(prhs[0]);
	int n2 = mxGetN(prhs[0]);
	int n = (n1 > n2 ? n1 : n2);
	
	char jobz[] = "V";
	char range[] = "A";
	
	/* main diagonal */
	double* d = mxCalloc(n, sizeof(double));
	double* temp = mxGetPr(prhs[0]);
	int i;
	for (i = 0; i < n; ++i) {
		d[i] = temp[i];
	}
	
	/* sub-/superdiagonal */
	double* e = mxCalloc(n, sizeof(double));
	temp = mxGetPr(prhs[1]);
	for (i = 0; i < n; ++i) {
		e[i] = temp[i];
	}
	
	double vl = 0, vu = 0; /* unused */
	int il = 0, iu = 0;    /* unused */
	double abstol = 1e-13; /* unused? */
	int m; /* number of eigenvalues found */
	
	/* computed eigenvalues */
	plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
	double* w = mxGetPr(plhs[0]);
	
	/* computed eigenvectors */
	plhs[1] = mxCreateDoubleMatrix(n, n, mxREAL);
	double* z = mxGetPr(plhs[1]);
	
	/* misc arguments */
	int ldz = n;
	int* isuppz = (int*) mxCalloc(2*n, sizeof(int));
	double* work1 = (double*) mxCalloc(1, sizeof(double));
	int lwork = -1;
	int* iwork1 = (int*) mxCalloc(1, sizeof(int));
	int liwork = -1;
	int info;
	
	snprintf(s, sizeof(s), "z = %p\n", z);
	mexPrintf(s);
	
	/* query for workspace size */
	dstegr_(jobz, range, &n, d, e, &vl, &vu, &il, &iu, &abstol,
		&m, w, z, &ldz, isuppz, work1, &lwork, iwork1, &liwork,
		&info);
	if (info < 0) {
		snprintf(s, sizeof(s), "[1] Input to dstegr had an illegal value (info: %d)", info);
		mexErrMsgTxt(s);
	}
	
	lwork = (int) work1[0];
	double* work = (double*) mxCalloc(lwork, sizeof(double));
	liwork = (int) iwork1[0];
	int* iwork = (int*) mxCalloc(liwork, sizeof(int));
	
	/* call LAPACK routine */
	dstegr_(jobz, range, &n, d, e, &vl, &vu, &il, &iu, &abstol,
		&m, w, z, &ldz, isuppz, work, &lwork, iwork, &liwork,
		&info);
	if (info != 0) {
		snprintf(s, sizeof(s), "[2] Input to dstegr had an illegal value (info: %d)", info);
		mexErrMsgTxt(s);
	}
	
	mxFree(d);
	mxFree(e);
	mxFree(isuppz);
	mxFree(work1);
	mxFree(work);
	mxFree(iwork1);
	mxFree(iwork);
}
