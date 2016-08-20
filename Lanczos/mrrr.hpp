#ifndef MRRR_HPP
#define MRRR_HPP

#include "Vecd.hpp"

extern const double DOUBLE_EPS;

/**
 * @param ritz Preallocated output array of size a.size() == m
 * @param S    Preallocated output array of size m*m
 * @param a    Main diagonal vector of size m
 * @param b    Super-/Subdiagonal vector of size m
 *
 * @return 0 on success, non-zero value on error
 */
int mrrr(double* ritz, double* S, const Vecd& a, const Vecd& b);

#endif
