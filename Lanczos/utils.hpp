#ifndef UTILS_HPP
#define UTILS_HPP

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cassert>
#include <ctime>
#include <climits>
#include <stdint.h>
#include <string>
#include <errno.h>

struct bool_index {
	const char* const str;
	const bool value;
	bool_index(const char* const _str, bool _value) : str(_str), value(_value) { }
};

extern const bool_index BOOLS[];

#define sbool(x) ((x) ? (BOOLS[0].str) : (BOOLS[1].str))
#define syesno(x) ((x) ? (BOOLS[2].str) : (BOOLS[3].str))

#define array_len(arr) (sizeof(arr)/sizeof((arr)[0]))

int arg2bool(bool& retval, int argc, char** argv, int argIndex);

#endif
