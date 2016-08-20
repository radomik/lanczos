#ifndef MRRR_HPP
#define MRRR_HPP

#include "utils.hpp"
#include "Vecd.hpp"

extern const double DOUBLE_EPS;

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
int mrrr(double* ritz, double* S, const Vecd& a, const Vecd& b, bool dbg);

#endif
