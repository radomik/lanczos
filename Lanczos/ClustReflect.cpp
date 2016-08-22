
typedef uint32_t u32;

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
vector<double> e; //TODO: use doubles everywhere, not floats, use u32 instead of u32 everywhere but size_t for u32*u32
vector<double> c;
double* S2 = new double[m * ci.size()]; // array of m rows and ci.size() columns
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
	x = {S[0][cii.left], S[0][cii.left+1], ..., S[0][cii.right]}
	u = x;
	u[max_cul_idx] += sign(x[max_cul_idx]) * x.norm();
	
	// s = S(:,ji)-2/(u'*u)*(S(:,cii)*u)*u(j);
	// S2(:,i) = s;
	double un = u.normPow();
	double uj = u[max_cul_idx];
	for (u32 p = 0; p < m; p++) { //TODO: paralellize this loop with ISPC/OpenMP
		double st = 0.0;
		size_t pm = p*m; 
		for (u32 q = cii.left, r = 0; q <= cii.right; q++, r++) {
			st += S[pm + q] * u[r];
		}
		S2[p*ci.size()+i] = S[pm + ji] - 2.0 / un*st*u[j];
	}
	
	for (u32 q = cii.left; q <= cii.right; q++) {
		idx[q] = false; // idx is a 0/1 vector obtained from Rescon step
	}
	idx[ji] = true;
	e.push_back(ritz[ji]);
	c.push_back(lres[ji]);
}
