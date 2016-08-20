#include "Rescon.hpp"

/**
 * Residual estimation and removing non-converged and spurious Ritz values.
 * 
 * @param s     Input non-symmetric matrix of size m*m obtained in MRRR step
 * @param ritz  Input ritz vector of size m obtained in MRRR step
 * @param b     Input b vector of size m obtained in LANCNO_INIT step
 * @param anorm Input value of anorm obtained in LANCNO_INIT step
 * @param eps   Input value of double epsilon
 * @param dbg   Print debug info
 * 
 * @return 0 on success, 1 when all RITZ values has been removed
 * 
 */
int Rescon::rescon(
	const double* s,
	const Vecd&   ritz,
	const Vecd&   b,
	const double  anorm,
	const double  eps,
	bool          dbg
) {
	size_t m = b.size();
	
	const double* s_last_row = s + m*(m-1);
	const double  b_last = b[m-1];
	const double  tol = eps * anorm;
	
	if (dbg) {
		fprintf(stderr, "%s: eps=%g, anorm=%g, tol=%g\n", __FUNCTION__, eps, anorm, tol);
	}
	
	_free();
	m_idx = new Vecb(m);
	size_t count_one = 0;
	Vecd* temp_lres = new Vecd(m);
	Vecd* temp_cul = new Vecd(m);
	
	for (size_t i = 0; i < m; i++) {
		const double lres_i = fabs(b_last * s_last_row[i]); // lres = abs(b(m)*S(m,:))'; % residual estimation
		const double cul_i = fabs(s[i]);                    // cul = abs(S(1,:))';
		
		// remove non-converged and spurious Ritz values
		(*m_idx)[i] = (lres_i < tol) && (cul_i > tol);
		
		if ((*m_idx)[i]) {
			count_one++;
		}
		
		(*temp_lres)[i] = lres_i;
		(*temp_cul)[i] = cul_i;
	}
	
	if (dbg) {
		fprintf(stderr, "%s: m=%zu, count_one=%zu\n", __FUNCTION__, m, count_one);
	}
	
	if (count_one == 0) {
		fprintf(stderr, "%s: [ERROR] All Ritz values are marked to be removed\n", __FUNCTION__);
		return 1;
	}
	
	m_lres = new Vecd(count_one);
	m_cul  = new Vecd(count_one);
	m_ritz = new Vecd(count_one);
	m_s    = new double[m * count_one];
	m_k    = count_one;
	for (size_t i = 0, k = 0; i < m; i++) {
		if ((*m_idx)[i]) {
			m_lres[k] = (*temp_lres)[i];
			m_cul[k]  = (*temp_cul)[i];
			m_ritz[k] = ritz[i];
			k++;
		}
	}
	
	// rewrite only columns s[i][k] which has idx[k] = 1
	for (size_t i = 0; i < m; i++) {
		const size_t sd = i * count_one;
		const size_t ss = i * m;
		for (size_t j = 0, k = 0; j < m; j++) {
			if ((*m_idx)[j]) {
				m_s[sd + k] = s[ss + j];
				k++;
			}
		}
	}
	
	delete(temp_lres);
	delete(temp_cul);
	return 0;
}

void Rescon::_free() {
	if (m_lres) {
		delete(m_lres); m_lres = NULL;
	}
	if (m_cul) {
		delete(m_cul); m_cul = NULL;
	}
	if (m_idx) {
		delete(m_idx); m_idx = NULL;
	}
	if (m_s) {
		delete[](m_s); m_s = NULL;
	}
}
