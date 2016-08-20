
///TODO: Fully rewrite this in C++, it is only draft for now

// distinguish the clusters of eigenvalues
vector<pair<size_t,size_t>> ci;
for (size_t k = 1, i = 0; k < mm; k++) {
	d = fabs(ritz[k] - ritz[i]);
	if (d > tol) {
		ci.push_back(pair(i, k-1));
		i = k;
		k++;
	}
}
ci.push_back(pair(i, mm));

// reflect using Householder
for (size_t i = 0; i < ci.size(); i++) {
	pair<int,int> cii = ci[i];
	max_cul_idx = 0;           // j
	max_cul_value = cul[cii.left]; // y
	for (size_t j = cii.left+1, idx = 1; j <= cii.right; j++, idx++) {
		if (cul[j] > max_cul_value) {
			max_cul_value = cul[j];
			max_cul_idx = idx;
		}
	}
	ji = cii.left + max_cul_idx;
	// x and u vextors are of size no more than mm
	x = {S[0][cii.left], S[0][cii.left+1], ..., S[0][cii.right]}
	u = x
	u[max_cul_idx] += sign(x[j]) * x.norm();
	
}
