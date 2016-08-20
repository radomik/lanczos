#ifndef LANCZNO_HPP
#define LANCZNO_HPP

#include "utils.hpp"
#include "Vecd.hpp"

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
int lanczno(double* ritz, double* S, const double* A, size_t n, uint32_t m, const Vecd* initStartVec = NULL);

#endif
