#ifndef RESCON_HPP
#define RESCON_HPP

#include "utils.hpp"
#include "Vecd.hpp"
#include "Vecb.hpp"

class Rescon {
public:
	Rescon() : m_lres(NULL), m_cul(NULL), m_ritz(NULL), m_idx(NULL), m_s(NULL) { }
	~Rescon() { _free(); }
	
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
	int rescon(
		const double* s,
		const Vecd&   ritz,
		const Vecd&   b,
		const double  anorm,
		const double  eps,
		bool          dbg
	);
	
	const Vecd& lres() const {
		assert(m_lres != NULL);
		return *m_lres;
	}
	
	const Vecd& cul() const {
		assert(m_cul != NULL);
		return *m_cul;
	}
	
	const Vecd& ritz() const {
		assert(m_ritz != NULL);
		return *m_ritz;
	}
	
	const Vecb& idx() const {
		assert(m_idx != NULL);
		return *m_idx;
	}
	
	const double* s() const {
		assert(m_s != NULL);
		return m_s;
	}
	
	size_t k() const {
		return m_k;
	}

private:
	Vecd*   m_lres; ///< output vector lres of size k, k <= m
	Vecd*   m_cul;  ///< output vector cul of size k, k <= m
	Vecd*   m_ritz; ///< output vector ritz of size k, k <= m
	Vecb*   m_idx;  ///< output 0/1 vector idx of size m (maybe use local variable)
	double* m_s;    ///< output s matrix of size m*k, k <= m
	double  m_k;    ///< output k value
	
	void _free();
};

#endif
