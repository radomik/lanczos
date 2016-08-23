
typedef uint32_t u32;
typedef uint64_t u64;
typedef double flt;

///TODO: Fully rewrite this in C++, it is only draft for now

// distinguish the clusters of eigenvalues
vector<pair<u32,u32>> ci;
for (u32 k = 1, i = 0; k < mm; k++) {
	d = fabs(ritz[k] - ritz[i]);
	if (d > tol) {
		ci.push_back(pair(i, k-1));
		i = k;
		k++;
	}
}
ci.push_back(pair(i, mm));

// reflect using Householder
vector<flt> e(ci.size()); //TODO: use doubles everywhere, not floats, use u32 instead of u32 everywhere but size_t for u32*u32
vector<flt> c(ci.size());
vector<flt> x;
vector<flt> u;

flt* S2 = new flt[m * ci.size()]; // array of m rows and ci.size() columns
for (u32 i = 0; i < ci.size(); i++) {
	pair<u32,u32> cii = ci[i];
	max_cul_idx = 0;               // j
	max_cul_value = cul[cii.left]; // y
	for (u32 j = cii.left+1, ix = 1; j <= cii.right; j++, ix++) {
		if (cul[j] > max_cul_value) {
			max_cul_value = cul[j];
			max_cul_idx = ix;
		}
	}
	ji = cii.left + max_cul_idx;
	// x and u vectors are of size no more than mm
	// x = {S[0][cii.left], S[0][cii.left+1], ..., S[0][cii.right]}
	x.resize(cii.right - cii.left + 1);
	memcpy(x.data(), S+cii.left, x.size());
	u = x;
	u[max_cul_idx] += sign(x[max_cul_idx]) * x.norm();
	
	// s = S(:,ji)-2/(u'*u)*(S(:,cii)*u)*u(j);
	// S2(:,i) = s;
	flt un = u.normPow();
	flt uj = u[max_cul_idx];
	for (u32 p = 0; p < m; p++) { //TODO: paralellize this loop with ISPC/OpenMP
		flt st = 0.0;
		u64 pm = p*m; 
		for (u32 q = cii.left, r = 0; q <= cii.right; q++, r++) {
			st += S[pm + q] * u[r];
		}
		S2[p*ci.size()+i] = S[pm + ji] - 2.0 / un*st*u[j];
	}
	
	for (u32 q = cii.left; q <= cii.right; q++) {
		idx[q] = false; // idx is a 0/1 vector obtained from Rescon step
	}
	idx[ji] = true;
	e[i] = ritz[ji];
	c[i] = lres[ji];
}

// e is a vector of computed eigen values

// compute eigen vectors
vector<flt> v = startvec; // startvec from initial Lanczos iteration input
vector<flt> v2;
flt* X = new flt[n * ci.size()];
flt a, b, a2, b2 = 0.0;

for (u32 k = 0; ; k++) {
	flt* S2k = &S2[k * ci.size()];
	for (u32 i = 0; i < n; i++) {
		flt* Xi = &X[i*ci.size()];
		for (u32 j = 0; j < ci.size(); j++) {
			Xi[j] += v[i] * S2k[j];
		}
	}
	
	if (k == m-1) {
		break;
	}
	
	if (k == 0) {
		TrMatrixd::mulByVector(r, A, v);  // r{n} = A{n}{n} * v{n}
	} else {
		TrMatrixd::mulByVector(rt, A, v);  // rt{n} = A{n}{n} * v{n}
		Vecd::mulByScalar(vt, v2, b2);     // vt{n} = b2{1} * v2{n}
		Vecd::sub(r, rt, vt);              // r{n}  = rt{n} - vt{n}
	}
	
	a = Vecd::dotProduct(v, r);  // a{1} = SUM v[i]*r[i]

	Vecd::mulByScalar(vt, v, a); // vt{n}   = a{1} * v{n}
	Vecd::sub(r, r, vt);         // r{n}    = r{n} - vt{n}
	b = r.getNorm();
	v2 = v;
	a2 = a;
	b2 = b;
	Vecd::divByScalar(v, r, b); // EACH_i v[i] = r[i]/b
}

vector<flt> xv(ci.size(), 0.0);

// xv is a vector of sums of squared elements in each column of X
for (u32 i = 0; i < n; i++) {
	flt* Xi = &X[i*ci.size()];
	for (u32 j = 0; j < ci.size(); j++) {
		xv[j] += Xi[j]*Xi[j];
	}
}

for (u32 i = 0; i < ci.size(); i++) {
	// xv[i] = 1.0 / diag(X)[i] / sqrt(xv[i])
	xv[i] = 1.0 / X[i*ci.size() + i] / sqrt(xv[i]);
}

Matrixd::mulByVector(X, X, xv);
